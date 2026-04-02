rm(list=ls())
library(MASS)
library(lavaan)
library(haven)
# Model tanımı
 model <- '
  # Latent değişken tanımları
  F1 =~ 1*X1 + X2
  F2 =~ 1*Y1 + Y2 + Y3
  F1 ~~ F1
  F2 ~~ F2

  F2 ~ gamma*F1
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

