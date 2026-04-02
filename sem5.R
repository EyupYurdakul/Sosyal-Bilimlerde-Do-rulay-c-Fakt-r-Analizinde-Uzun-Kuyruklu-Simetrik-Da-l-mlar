library(MASS)
library(lavaan)
library(dplyr)
# ==========================================================
# 1. VERİ ÜRETİMİ (DÜZELTİLDİ)
# ==========================================================
n <- 10000
gamma_true <- 1
var_F1 <- 1
var_F2 <- 1
set.seed(456)
# Latent değişkenler
F1 <- rnorm(n, mean = 0, sd = sqrt(var_F1))
F2 <- gamma_true*F1 + sqrt(1 - gamma_true^2)*rnorm(n) # DÜZELTME
# Gözlenen değişkenler
theta_epsilon <- diag(c(1, 2, 3, 4, 5))
epsilon <- mvrnorm(n, mu = rep(0,5), Sigma = theta_epsilon)
X1 <- 0.6*F1 + epsilon[,1]
X2 <- 0.7*F1 + epsilon[,2]
Y1 <- 0.8*F2 + epsilon[,3]
Y2 <- 0.9*F2 + epsilon[,4]
Y3 <- 1.0*F2 + epsilon[,5]
df <- data.frame(X1, X2, Y1, Y2, Y3)
# ==========================================================
# 2. SEM MODELİ (LAVAAN)
# ==========================================================
sem_model <- '
  F1 =~ 1*X1 + X2
  F2 =~ 1*Y1 + Y2 + Y3
  F1 ~~ NA*F1 # Varyans sabitleme
  F2 ~~ NA*F2 # Varyans sabitleme
  F2 ~ gamma*F1 # Yapısal ilişki
'
fit_sem <- sem(model = sem_model, 
               data = df,
               estimator = "WLS", std.lv = FALSE)
summary(fit_sem)
# ==========================================================
# 3. MANUEL WLS TAHMİN (DÜZELTİLDİ)
# ==========================================================
# Kovaryans matrisi ve vektörleştirme fonksiyonu
S <- cov(df)
vech <- function(M) M[lower.tri(M, diag = TRUE)]
s_obs <- vech(S)
# Teorik kovaryans matrisi fonksiyonu
sigma_theta <- function(params) {
  lambda <- params[1:5]
  gamma <- params[6]
  theta <- params[7:11]
  # Yapısal parametreler
  cov_F1F2 <- gamma
  # Beklenen kovaryans matrisi
  Sigma <- matrix(0, nrow = 5, ncol = 5)
  # F1 yüklemeleri (X1, X2)
  Sigma[1:2, 1:2] <- tcrossprod(lambda[1:2]) + diag(theta[1:2])
  # F2 yüklemeleri (Y1, Y2, Y3)
  Sigma[3:5, 3:5] <- tcrossprod(lambda[3:5]) + diag(theta[3:5])
  # Çapraz kovaryanslar (F1 ve F2 arasındaki ilişkiler)
  Sigma[1:2, 3:5] <- tcrossprod(lambda[1:2], lambda[3:5]) * cov_F1F2
  Sigma[3:5, 1:2] <- t(Sigma[1:2, 3:5])
  return(vech(Sigma))
}
# Amaç fonksiyonu (WLS)
W <- diag(length(s_obs)) # Diyagonal ağırlık matrisi
#W=solve(diag(s_obs))
WLS_objective <- function(params) {
  sigma_model <- sigma_theta(params)
  #
  diff <- s_obs - sigma_model
  t(diff) %*% W %*% diff
}
W
s_obs
# Optimizasyon başlangıç parametreleri ve sınırlar
start_params <- c(0.6,0.7,0.8,0.9,1.0, # Lambda değerleri (yüklemeler)
                  gamma_true, # Gamma değeri (yapısal regresyon katsayısı)
                  c(1,2,3,4,5)) # Theta değerleri (hata varyansları)
lower_bounds <- c(rep(-Inf,6), rep(0.01,5)) # Theta için alt sınır pozitif olmalı
upper_bounds <- c(rep(Inf,11))
result_manual <- optim(
  par = start_params,
  fn = WLS_objective,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(maxit = 1000, reltol = 1e-12)
)
param_names <- c("lambda1", "lambda2", "lambda3", "lambda4", "lambda5", 
                 "gamma", "theta1", "theta2", "theta3", "theta4", "theta5")
names(result_manual$par) <- param_names
# ==========================================================
# 4. PARAMETRELERİ SINIFLANDIR VE RAPORLA (VEKTÖR FORMATINDA)
# ==========================================================
Lambda_manual <- c(
  result_manual$par["lambda1"],
  result_manual$par["lambda2"],
  result_manual$par["lambda3"],
  result_manual$par["lambda4"],
  result_manual$par["lambda5"]
)
Theta_manual <- c(
  result_manual$par["theta1"],
  result_manual$par["theta2"],
  result_manual$par["theta3"],
  result_manual$par["theta4"],
  result_manual$par["theta5"]
)
Gamma_manual <- result_manual$par["gamma"]
# Lavaan'dan parametreleri çekme
lambda_sem_est <- inspect(fit_sem, what = "est")$lambda
lambda_F1 <- lambda_sem_est[1:2, "F1"]
lambda_F2 <- lambda_sem_est[3:5, "F2"]
lambda_sem_vector <- c(lambda_F1, lambda_F2)
theta_sem_est <- diag(inspect(fit_sem, what = "est")$theta)
psi_sem_est <- inspect(fit_sem, what = "est")$psi
psi_sem_vector <- c(psi_sem_est[1,1], psi_sem_est[2,2], psi_sem_est[1,2])
Gamma_sem <- tryCatch({
  g <- lavaan::inspect(fit_sem, what = "est")$regression["F2", "F1"]
  if (is.null(g) || is.na(g)) NA else as.numeric(g)
}, error = function(e) NA)
# ==========================================================
# 5. KARŞILAŞTIRMA RAPORU (Lambda, Theta, Gamma)
# ==========================================================
lambda_comp <- data.frame(
  YUKLEME = c("X1=~F1", "X2=~F1", "Y1=~F2", "Y2=~F2", "Y3=~F2"),
  TRUE_LAMBDA = c(0.6, 0.7, 0.8, 0.9, 1.0),
  SEM_EST = round(lambda_sem_vector, 4),
  MANUAL_WLS = round(Lambda_manual, 4)
)
theta_comp <- data.frame(
  Degisken = c("X1", "X2", "Y1", "Y2", "Y3"),
  TRUE_Theta = c(1, 2, 3, 4, 5),
  SEM_EST = round(theta_sem_est, 4),
  MANUAL_WLS = round(Theta_manual, 4)
)
gamma_comp <- data.frame(
  Katsayi = "F2 ~ F1",
  TRUE_GAMMA = gamma_true,
  SEM_EST = round(Gamma_sem, 4),
  MANUAL_WLS = round(Gamma_manual, 4)
)
# ==========================================================
# 6. SONUÇLARI YAZDIR
# ==========================================================
cat("\n==== LAMBDA (YÜKLEME) KARŞILAŞTIRMASI ====\n")
print(lambda_comp)
cat("\n==== THETA (RESIDUAL VARIANCES) KARŞILAŞTIRMASI ====\n")
print(theta_comp)
cat("\n==== GAMMA (STRUCTURAL REGRESSION) KARŞILAŞTIRMASI ====\n")
print(gamma_comp)
# Optimizasyon sonucu
cat("\nOptimizasyon sonucu:", result_manual$message, 
    "\nToplam iterasyon:", result_manual$counts[["function"]])
# LİB GEN Lİ sem kaynaklı kitaplar vd kaglee araştır
#w ve lavaan değerini aynı yapmak 
#amos ile uyarlılık 
# modeli genellemek ve matris formatına çevirmek
