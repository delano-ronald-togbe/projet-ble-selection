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

# Créer une structure de projet
# dir.create("projet_ble_selection")
# setwd("projet_ble_selection")
# dir.create("data")
# dir.create("scripts")
# dir.create("outputs")
# dir.create("reports")