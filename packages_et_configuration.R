# Installation des packages nécessaires
required_packages <- c("BGLR", "rrBLUP", "tidyverse", "tidymodels", "bonsai",
                       "prospectr", "viridis", "corrplot", "ggcorrplot", "vcfR",
                       "data.table", "sommer", "here", "targets", "tarchetypes")
install_if_missing <- function(p) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

# Création des défférents sous répertoires du projet

# projet_ble_selection/
# ├── data/                     # Données brutes et transformées
# ├── scripts/
# │ ├── 01_acquisition.R
# │ ├── 02_preprocessing.R
# │ ├── 03_gwas.R
# │ ├── 04_prediction_benchmark.R
# │ ├── 05_visualisation.R
# │ └── 06_shiny_app/
# │     ├── app.R
# │     └── Dockerfile
# ├── outputs/                   # Graphiques et résultats
# ├── reports/                   # Rapport Rmd
# └── renv.lock                  # Gestion des versions (optionnel)
