library(MASS)
library(lavaan)
library(dplyr)
# ==========================================================
# 1. VERİ ÜRETİMİ (DÜZELTİLDİ)
# ==========================================================
n <- 10000
gamma_true <- 0.3
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
  F1 =~ NA*X1 + X2
  F2 =~ NA*Y1 + Y2 + Y3
  F1 ~~ 1*F1    # Varyans sabitleme
  F2 ~~ 1*F2    # Varyans sabitleme
  F2 ~ gamma*F1 # Yapısal ilişki
'
fit_sem <- sem(model = sem_model, 
               data = df, 
               estimator = "WLS",
               std.lv = FALSE)
summary(fit_sem)
# ==========================================================
# 4. MANUEL WLS TAHMİN (DÜZELTİLDİ)
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
WLS_objective <- function(params) {
  sigma_model <- sigma_theta(params)
  
  diff <- s_obs - sigma_model
  t(diff) %*% W %*% diff
}
# Optimizasyon başlangıç parametreleri ve sınırlar
start_params <- c(0.6,0.7,0.8,0.9,1.0,   # Lambda değerleri (yüklemeler)
                  gamma_true,            # Gamma değeri (yapısal regresyon katsayısı)
                  c(1,2,3,4,5))          # Theta değerleri (hata varyansları)

lower_bounds <- c(rep(-Inf,6), rep(0.01,5))   # Theta için alt sınır pozitif olmalı
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
                 "gamma",
                 "theta1", "theta2", "theta3", "theta4", "theta5")
names(result_manual$par) <- param_names
# Sonuçları yazdırma
cat("\n==== MANUEL WLS OPTİMİZASYON SONUÇLARI ====\n")
print(result_manual$par)
# ==========================================================
# 4. SONUÇ KARŞILAŞTIRMA (GÜNCELLENDİ)
# ==========================================================
# Lavaan sonuçları
sem_est <- c(
  parameterEstimates(fit_sem) %>% 
    filter(op == "=~") %>% pull(est),
  parameterEstimates(fit_sem) %>% 
    filter(op == "~") %>% pull(est),
  diag(inspect(fit_sem, "est")$theta)
)
# Manuel WLS sonuçları
manual_est <- result$par
# Karşılaştırma tablosu
comparison <- data.frame(
  Parametre = c(paste0("λ",1:5), "γ", paste0("θ",1:5)),
  Gerçek = c(0.6,0.7,0.8,0.9,1.0,0.3,1:5),
  Lavaan = round(sem_est, 4),
  Manuel_WLS = round(manual_est, 4)
)
# Sonuçları göster
cat("================ PARAMETRE KARŞILAŞTIRMASI ================\n")
print(comparison)
# Yakınsama kontrolü
cat("\nOptimizasyon sonucu:", result$message, 
    "\nToplam iterasyon:", result$counts[1])
