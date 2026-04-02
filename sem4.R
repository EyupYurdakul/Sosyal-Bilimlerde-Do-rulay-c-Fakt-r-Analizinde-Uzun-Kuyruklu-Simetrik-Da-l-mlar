rm(list=ls())
library(MASS)
library(lavaan)
library(dplyr)
library(numDeriv)

# ==========================================================
# 1. VERİ ÜRETİMİ
# ==========================================================
n <- 10000
gamma_true <- 1
var_F1 <- 1
var_F2 <- 1

set.seed(456)

# Latent değişkenler
F1 <- rnorm(n, mean = 0, sd = sqrt(var_F1))
F2 <- gamma_true*F1 + sqrt(1 - gamma_true^2)*rnorm(n)

# Gözlenen değişkenler
theta_epsilon <- diag(c(1, 2, 3, 4, 5))
epsilon <- mvrnorm(n, mu = rep(0,5), Sigma = theta_epsilon)

X1 <- 0.6*F1 + epsilon[,1]
X2 <- 0.7*F1 + epsilon[,2]
Y1 <- 0.8*F2 + epsilon[,3]
Y2 <- 0.9*F2 + epsilon[,4]
Y3 <- 1.0*F2 + epsilon[,5]

df <- data.frame(X1, X2, Y1, Y2, Y3)
X1[1] <- 100  # Adding an outlier
X1[1]
head(df,10)
# ==========================================================
# 2. SEM MODELİ (LAVAAN)
# ==========================================================
sem_model <- '
  F1 =~ NA*X1 + X2
  F2 =~ NA*Y1 + Y2 + Y3
  F1 ~~ 1*F1
  F2 ~~ 1*F2
  F2 ~ gamma*F1
'

fit_sem <- sem(model = sem_model, data = df, estimator = "WLS", std.lv = FALSE)
fit_sem_robust <- sem(model = sem_model, data = df, estimator = "WLSMV")

# ==========================================================
# 3. MANUEL ROBUST WLS TAHMİN
# ==========================================================
# Kovaryans matrisi ve vektörleştirme
S <- cov(df)
vech <- function(M) M[lower.tri(M, diag = TRUE)]
s_obs <- vech(S)

# Asimptotik kovaryans matrisi fonksiyonu
calculate_Gamma <- function(S, n) {
  p <- nrow(S)
  Gamma <- matrix(0, nrow = p*(p+1)/2, ncol = p*(p+1)/2)
  idx <- 1
  
  for (i in 1:p) {
    for (j in i:p) {
      idy <- 1
      for (k in 1:p) {
        for (l in k:p) {
          Gamma[idx, idy] <- (n-1)/n * (S[i,k]*S[j,l] + S[i,l]*S[j,k])
          idy <- idy + 1
        }
      }
      idx <- idx + 1
    }
  }
  return(Gamma)
}

# Teorik kovaryans matrisi fonksiyonu
sigma_theta <- function(params) {
  lambda <- params[1:5]
  gamma <- params[6]
  theta <- params[7:11]
  
  Sigma <- matrix(0, nrow = 5, ncol = 5)
  Sigma[1:2, 1:2] <- tcrossprod(lambda[1:2]) + diag(theta[1:2])
  Sigma[3:5, 3:5] <- tcrossprod(lambda[3:5]) + diag(theta[3:5])
  Sigma[1:2, 3:5] <- tcrossprod(lambda[1:2], lambda[3:5]) * gamma
  Sigma[3:5, 1:2] <- t(Sigma[1:2, 3:5])
  
  return(vech(Sigma))
}

# Robust ağırlık matrisi
Gamma <- calculate_Gamma(S, n)
W_robust <- solve(Gamma)

# Robust amaç fonksiyonu
WLS_robust_objective <- function(params) {
  sigma_model <- sigma_theta(params)
  diff <- s_obs - sigma_model
  t(diff) %*% W_robust %*% diff
}

# Optimizasyon parametreleri
start_params <- c(0.6,0.7,0.8,0.9,1.0, gamma_true, c(1,2,3,4,5))
param_names <- c("lambda1", "lambda2", "lambda3", "lambda4", "lambda5", 
                 "gamma", "theta1", "theta2", "theta3", "theta4", "theta5")

lower_bounds <- c(rep(-Inf,6), rep(0.01,5))
upper_bounds <- rep(Inf,11)

# Robust optimizasyon
result_robust <- optim(
  par = start_params,
  fn = WLS_robust_objective,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(maxit = 1000, factr = 1e-12)
)
names(result_robust$par) <- param_names

# Robust standart hatalar
jacobian <- jacobian(sigma_theta, result_robust$par)
cov_robust <- solve(t(jacobian) %*% W_robust %*% jacobian)
robust_se <- sqrt(diag(cov_robust))
names(robust_se) <- param_names

# ==========================================================
# 4. SONUÇLARI KARŞILAŞTIRMA
# ==========================================================
# Parametreleri çek
get_params <- function(result) {
  list(
    lambda = result$par[1:5],
    gamma = result$par["gamma"],
    theta = result$par[7:11]
  )
}

params_robust <- get_params(result_robust)

# Lavaan'dan parametreler
lambda_sem <- c(inspect(fit_sem, "est")$lambda[1:2,1], inspect(fit_sem, "est")$lambda[3:5,2])
theta_sem <- diag(inspect(fit_sem, "est")$theta)
gamma_sem <- inspect(fit_sem, "est")$beta[2,1]

lambda_sem_robust <- c(inspect(fit_sem_robust, "est")$lambda[1:2,1], inspect(fit_sem_robust, "est")$lambda[3:5,2])
theta_sem_robust <- diag(inspect(fit_sem_robust, "est")$theta)
gamma_sem_robust <- inspect(fit_sem_robust, "est")$beta[2,1]

# Karşılaştırma tabloları
lambda_comp <- data.frame(
  Parametre = c("X1=~F1", "X2=~F1", "Y1=~F2", "Y2=~F2", "Y3=~F2"),
  Gerçek = c(0.6, 0.7, 0.8, 0.9, 1.0),
  SEM_WLS = round(lambda_sem, 4),
  SEM_WLSMV = round(lambda_sem_robust, 4),
  Manuel_Robust = round(params_robust$lambda, 4),
  SE_Robust = round(robust_se[1:5], 4)
)

theta_comp <- data.frame(
  Parametre = paste0("theta", 1:5),
  Gerçek = 1:5,
  SEM_WLS = round(theta_sem, 4),
  SEM_WLSMV = round(theta_sem_robust, 4),
  Manuel_Robust = round(params_robust$theta, 4),
  SE_Robust = round(robust_se[7:11], 4)
)

gamma_comp <- data.frame(
  Parametre = "F2~F1",
  Gerçek = gamma_true,
  SEM_WLS = round(gamma_sem, 4),
  SEM_WLSMV = round(gamma_sem_robust, 4),
  Manuel_Robust = round(params_robust$gamma, 4),
  SE_Robust = round(robust_se["gamma"], 4)
)

# ==========================================================
# 5. SONUÇLARI GÖSTER
# ==========================================================
cat("=== LAMBDA YÜKLEME KATSAYILARI ===\n")
print(lambda_comp)

cat("\n=== THETA HATA VARYANSLARI ===\n")
print(theta_comp)

cat("\n=== GAMMA YAPISAL KATSAYI ===\n")
print(gamma_comp)

cat("\nOptimizasyon sonucu:", result_robust$message, 
    "\nToplam iterasyon:", result_robust$counts[["function"]])