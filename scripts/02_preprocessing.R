# scripts/02_preprocessing.R
# Objectif : Contrôle qualité, calcul de la matrice G, prétraitement des spectres NIR

load("D:/Treatment/R_WorkSpace/projet-ble-selection/data/wheat_raw.RData")
library(tidyverse)
library(prospectr)
library(Matrix)

# 1. Prétraitement génotypique
# Les données BGLR sont déjà nettoyées, mais on peut vérifier MAF et taux de missing
# Calcul de la matrice de parenté génomique G (VanRaden)
X <- genotypes
# Centrage-réduction des SNPs (recommandé pour GBLUP)
X_scaled <- scale(X, center = TRUE, scale = TRUE)
# G = (X_scaled %*% t(X_scaled)) / ncol(X_scaled)
G_matrix <- tcrossprod(X_scaled) / ncol(X_scaled)

# Vérification que G est définie positive (petite correction si nécessaire)
if (min(eigen(G_matrix, symmetric = TRUE, only.values = TRUE)$values) < 0) {
  G_matrix <- as.matrix(nearPD(G_matrix)$mat)
}
rownames(G_matrix) <- colnames(G_matrix) <- rownames(phenotypes)

# 2. Prétraitement spectral
# Vérification des dimensions initiales
cat("Dimensions de nir_spectra :", dim(nir_spectra), "\n")

# Standard Normal Variate (SNV)
nir_snv <- standardNormalVariate(X = nir_spectra)
cat("Dimensions après SNV :", dim(nir_snv), "\n")

# Dérivée première (Savitzky-Golay, fenêtre 11, ordre 2)
# Note : cette fonction conserve le même nombre de colonnes
nir_deriv1 <- savitzkyGolay(X = nir_snv, p = 2, w = 11, m = 1)
cat("Dimensions après dérivée :", dim(nir_deriv1), "\n")

# Note: Les fonctions comme savitzkyGolay ou gapDer appliquent une fenêtre glissante 
# sur le spectre. Selon la méthode de gestion des extrémités (padding, troncature, etc.), 
# il peut y avoir une perte de points au début et à la fin du spectre, réduisant ainsi
# le nombre total de colonnes. Par exemple, avec une fenêtre de taille w, si aucune 
# extrapolation n'est faite, les (w-1)/2 premières et dernières longueurs d'onde 
# peuvent être supprimées.

# On garde la version traitée
nir_processed <- nir_deriv1

# S'assurer que le nombre de colonnes n'a pas changé
if (ncol(nir_processed) != ncol(nir_spectra)) {
  warning("Le nombre de colonnes a changé après prétraitement. Adaptation des noms.
          Voir la note précédente en commentaire dans le script pour mieux comprendre")
  # On génère des noms génériques si nécessaire
  colnames(nir_processed) <- paste0("proc_wave_", 1:ncol(nir_processed))
} else {
  colnames(nir_processed) <- paste0("proc_", colnames(nir_spectra))
}

cat("Dimensions finales de nir_processed :", dim(nir_processed), "\n")

# Calcul de la matrice de parenté phénomique P (corrélation entre individus)
# On utilise la matrice de corrélation des spectres prétraités
P_matrix <- cor(t(nir_processed))
# S'assurer qu'elle est définie positive
if (min(eigen(P_matrix, symmetric = TRUE, only.values = TRUE)$values) < 0) {
  P_matrix <- as.matrix(nearPD(P_matrix)$mat)
}
rownames(P_matrix) <- colnames(P_matrix) <- rownames(phenotypes)

# Sauvegarde
save(genotypes, phenotypes, X_scaled, G_matrix, nir_processed, P_matrix,
     file = "D:/Treatment/R_WorkSpace/projet-ble-selection/data/wheat_processed.RData")
cat("Données prétraitées sauvegardées dans data/wheat_processed.RData\n")
