library(MASS)
library(lavaan)
library(haven)

data <- read_sav("Downloads/tez4.sav")

# Veri setini kontrol et

if (is.null(data)) {
  stop("Veri seti yüklenmedi. Lütfen 'data' değişkenini uygun bir veri seti ile doldurun.")
}

# Lavaan model tanımı (F1 ve F2 gizli değişken olarak ele alınıyor)
model <- '
  # Latent de??i??ken tan??mlar??
  F1 =~ NA*X1 + X2
  F2 =~ NA*Y1 + Y2 + Y3
  F1 ~~ 2*F1
  F2 ~~ 1*F2
 F2 ~ gamma*F1 
'

# Modeli lavaan ile tahmin etme
fit <- cfa(model, data = data, estimator = "ML")
# Sonuçların özeti
summary(fit, standardized = TRUE, fit.measures = TRUE)







# Parametre tahminlerini alma
param_estimates <- parameterEstimates(fit)
# Tahmin edilen parametreleri görüntüleme
cat("Lavaan'dan tahmin edilen gamma:", 
    param_estimates[param_estimates$lhs == "F2" & param_estimates$op == "~" & param_estimates$rhs == "F1", "est"], "\n")

cat("Lavaan'dan tahmin edilen beta_F1 katsayıları:", 
    param_estimates[param_estimates$lhs == "F1" & param_estimates$op == "=~", "est"], "\n")

cat("Lavaan'dan tahmin edilen beta_F2 katsayıları:", 
    param_estimates[param_estimates$lhs == "F2" & param_estimates$op == "=~", "est"], "\n")

summary(fit, standardized = TRUE)