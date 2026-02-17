# scripts/01_acquisition.R
# Objectif : Charger les données réelles du package BGLR (blé) et simuler des spectres NIR

library(BGLR)
library(tidyverse)

# Créer les dossiers nécessaires
dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# 1. Charger les données wheat
data(wheat)
# wheat.X : matrice génotypique (599 lignées × 1279 SNPs) codée en -1, 0, 1
# wheat.Y : matrice phénotypique (599 × 4) avec NA pour validation croisée
#   colonnes : Env1, Env2, Env4, Env5
# wheat.A : matrice de parenté pedigree (non utilisée ici)

genotypes <- wheat.X
phenotypes <- wheat.Y
colnames(phenotypes) <- c("Env1", "Env2", "Env4", "Env5")

cat("Dimensions génotypes :", dim(genotypes), "\n")
cat("Dimensions phénotypes :", dim(phenotypes), "\n")
cat("Taux de NA dans phénotypes :", sum(is.na(phenotypes)) / length(phenotypes), "\n")

# 2. Simulation de données SPIR (phénomique) pour les 599 lignées
set.seed(123)
n_lines <- nrow(genotypes)
n_wavelengths <- 200

# On simule des spectres corrélés aux phénotypes réels
# On utilise la première composante principale des phénotypes comme signal commun
pca_pheno <- prcomp(phenotypes, center = TRUE, scale. = TRUE, na.action = na.omit)
signal <- predict(pca_pheno, phenotypes)[, 1]  # première PC
signal <- scale(signal) * 2  # amplitude

# Création de la matrice spectrale avec structure par blocs
nir_spectra <- matrix(NA, nrow = n_lines, ncol = n_wavelengths)

for (i in 1:n_lines) {
  # Bruit de fond
  baseline <- rnorm(n_wavelengths, 0, 0.5)
  # Trois pics informatifs
  peak1 <- dnorm(1:n_wavelengths, mean = 50, sd = 15) * signal[i] * 5
  peak2 <- dnorm(1:n_wavelengths, mean = 120, sd = 20) * signal[i] * 3
  peak3 <- dnorm(1:n_wavelengths, mean = 180, sd = 10) * signal[i] * 4
  nir_spectra[i, ] <- baseline + peak1 + peak2 + peak3 + rnorm(n_wavelengths, 0, 0.2)
}

colnames(nir_spectra) <- paste0("wave_", 1:n_wavelengths)

# Ajouter un peu de bruit pour plus de réalisme
nir_spectra <- nir_spectra + matrix(rnorm(n_lines * n_wavelengths, 0, 0.1), nrow = n_lines)

# Sauvegarde
save(genotypes, phenotypes, nir_spectra, 
     file = "D:/Treatment/R_WorkSpace/projet-ble-selection/data/wheat_raw.RData")
cat("Données brutes sauvegardées dans data/wheat_raw.RData\n")
