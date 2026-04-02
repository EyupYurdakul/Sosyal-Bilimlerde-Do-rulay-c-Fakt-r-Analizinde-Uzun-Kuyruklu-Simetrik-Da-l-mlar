rm(list=ls())
library(MASS)
library(lavaan)
library(haven)
# Model tanımı
model <- '
# Latent değişken tanımları
  F1 =~ NA*X1 + NA*X2
  F2 =~ NA*Y1 + NA*Y2 + NA*Y3
  F2 ~~ F1
  F2 ~~ F1
  
  # Latent değişkenler arası ilişki

'
data <- read_sav("Downloads/tez4.sav")
View(tez4)
# Veri seti kontrolü
if (is.null(data)) {
  stop("Veri seti yüklenmedi. Lütfen 'data' değişkenini uygun bir veri seti ile doldurun.")
}
fit <- sem(model, data = data)
fit
# Sonuçları özetle
summary(fit, standardized = TRUE)

