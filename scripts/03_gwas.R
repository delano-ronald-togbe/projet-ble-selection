# scripts/03_gwas.R
# Objectif : Réaliser une GWAS pour chaque environnement

load("D:/Treatment/R_WorkSpace/projet-ble-selection/data/wheat_processed.RData")
library(sommer)
library(tidyverse)
library(ggrepel)

# 1. Préparation des données
pheno_df <- as.data.frame(phenotypes)
pheno_df$Line <- rownames(phenotypes)

# Ajout de composantes principales pour contrôler la structure de population
pca_G <- prcomp(G_matrix, center = FALSE, scale. = FALSE)
pc_scores <- as.data.frame(pca_G$x[, 1:5])
colnames(pc_scores) <- paste0("PC", 1:5)
pc_scores$Line <- rownames(phenotypes)

pheno_df <- left_join(pheno_df, pc_scores, by = "Line")
pheno_df$Line <- as.character(pheno_df$Line)  # s'assurer que c'est un vecteur de caractères

# 2. Alignement des marqueurs avec pheno_df
# Vérifier si genotypes a des noms de lignes
if (!is.null(rownames(genotypes))) {
  # Aligner par les noms
  idx <- match(pheno_df$Line, rownames(genotypes))
  if (any(is.na(idx))) {
    stop("Certaines lignées de pheno_df ne sont pas dans rownames(genotypes)")
  }
  markers <- genotypes[idx, , drop = FALSE]
} else {
  # Pas de noms de lignes, on suppose le même ordre
  if (nrow(genotypes) != nrow(pheno_df)) {
    stop(
      "Le nombre de lignées dans genotypes et pheno_df diffère et genotypes n'a pas de noms de lignes."
    )
  }
  markers <- genotypes
}

#rownames(markers) <- colnames(G_matrix)

cat("Dimensions de markers après alignement :", dim(markers), "\n")
cat(
  "Vérification : les premières lignées de pheno_df et les premières lignes de markers correspondent ?\n"
)
print(head(pheno_df$Line))
if (!is.null(rownames(markers))) {
  print(head(rownames(markers)))
} else {
  print("Pas de noms de lignes dans markers")
}

# 3. Boucle sur les traits
traits <- c("Env1", "Env2", "Env4", "Env5")
gwas_results_list <- list()

for (trait in traits) {
  cat("\n========== Analyse GWAS pour", trait, "==========\n")
  
  # 2.1 Exécution de GWAS en utilisant le modèle ajusté
  # Important : spécifier gTerm = "u:Line" pour indiquer quel effet aléatoire correspond aux génotypes
  gwas_out <- GWAS(
    fixed = as.formula(paste(trait, "~ 1 + PC1 + PC2 + PC3 + PC4 + PC5")),
    random = ~ vsr(Line, Gu = G_matrix),
    rcov = ~ units,
    M = markers,
    gTerm = "u:Line",
    # nom de l'effet aléatoire dans le modèle
    n.PC = 0,
    # les PC sont déjà dans le modèle fixe
    P3D = TRUE,
    # réutilisation des variances (plus rapide)
    verbose = TRUE,
    min.MAF = 0.00,
    # on garde tous les marqueurs
    data = pheno_df   # Ajout explicite des données
  )
  
  # 2.2 Extraction des p-values
  # Selon la version, les p-values peuvent être dans gwas_out$p.value ou gwas_out$Pval
  if (!is.null(gwas_out$p.value)) {
    pvals <- gwas_out$p.value
  } else if (!is.null(gwas_out$pvals)) {
    pvals <- gwas_out$pvals
  } else {
    stop("Impossible de trouver les p-values dans l'objet retourné par GWAS")
  }
  
  # 2.3 Création du tableau de résultats
  gwas_tab <- data.frame(
    Trait = trait,
    SNP = colnames(markers),
    Chromosome = paste0("chr_", rep(1:7, length.out = ncol(markers))),
    # simulation (7 chromosomes)
    Position = 1:ncol(markers),
    P_value = pvals,
    stringsAsFactors = FALSE
  )
  
  gwas_results_list[[trait]] <- gwas_tab
  
  # 2.4 Manhattan plot
  
  # La list des SNPs to mettre en evidence sont dans l'4=''objet snpsOfInterest
  
  # Préparer le jeu de données
  don <- gwas_tab %>%
    
    # Calculer la taille de chaque chromosome
    group_by(Chromosome) %>%
    summarise(chr_len = max(Position)) %>%
    
    # Calculer la position cumulée de chaque chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Ajouter cette information au jeu de données initial
    left_join(gwas_tab, ., by = c("Chromosome" = "Chromosome")) %>%
    
    # Ajouter une position cumulée pour chaque SNP
    arrange(Chromosome, Position) %>%
    mutate(Position_cum = Position + tot) %>%
    
    # Ajouter des informations de surbrillance et d'annotation
    mutate(is_highlight = ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate(is_annotate = ifelse(-log10(P_value) > -log10(0.05 / ncol(markers)), "yes", "no"))
  
  # Préparer l'axe X
  axisdf <- don %>%
    group_by(Chromosome) %>%
    summarize(center = (max(Position_cum) + min(Position_cum)) / 2)
  
  # Créer le graphique
  manhattan <- ggplot(don, aes(x = Position_cum, y = -log10(P_value))) +
    
    # Afficher tous les points
    geom_point(aes(color = as.factor(Chromosome)), alpha = 0.8, size = 1) +
    scale_color_manual(values = rep(c("steelblue", "tomato"), 22)) +
    
    # Personnaliser l'axe X :
    scale_x_continuous(label = axisdf$Chromosome, breaks = axisdf$center) +
    # scale_y_continuous(expand = c(0, 0)) +     # supprimer l'espace entre la zone de tracé et l'axe X
    
    # Ajouter les points en surbrillance
    geom_point(
      data = subset(don, is_highlight == "yes"),
      color = "green",
      size = 2
    ) +
    
    # Ajouter des étiquettes avec ggrepel pour éviter les chevauchements
    geom_label_repel(data = subset(don, is_annotate == "yes"),
                     aes(label = SNP),
                     size = 2.5) +
    
    geom_hline(
      yintercept = -log10(0.05 / ncol(markers)),
      linetype = "dashed",
      color = "red"
    ) +
    
    # Personnaliser le thème :
    theme_bw() +
    theme(
      legend.position = "none",
      # panel.border = element_blank(),
      panel.grid.major.x = element_blank()
      # panel.grid.minor.x = element_blank()
    ) +
    labs(title = paste("GWAS -", trait),
         x = "Position cumulée",
         y = "-log10(p)")
  
  ggsave(
    paste0(
      "D:/Treatment/R_WorkSpace/projet-ble-selection/outputs/gwas_manhattan_",
      trait,
      ".png"
    ),
    manhattan,
    width = 8,
    height = 4
  )
  
  # 2.5 QQ-plot
  gwas_qq <- gwas_tab %>%
    arrange(P_value) %>%
    mutate(expected = -log10(ppoints(n())),
           observed = -log10(P_value))
  
  qqplot <- ggplot(gwas_qq, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.5) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed"
    ) +
    labs(title = paste("QQ-plot -", trait),
         x = "-log10(p) attendue",
         y = "-log10(p) observée") +
    theme_minimal()
  
  ggsave(
    paste0(
      "D:/Treatment/R_WorkSpace/projet-ble-selection/outputs/gwas_qqplot_",
      trait,
      ".png"
    ),
    qqplot,
    width = 6,
    height = 6
  )
  
  cat("Terminé pour", trait, "\n")
}

# 3. Sauvegarde des résultats combinés
save(gwas_results_list, file = "D:/Treatment/R_WorkSpace/projet-ble-selection/data/gwas_results_all.RData")

# 4. Tableau des SNPs significatifs (correction de Bonferroni)
threshold <- 0.05 / ncol(markers)
significant_snps <- bind_rows(gwas_results_list) %>%
  filter(P_value < threshold) %>%
  arrange(Trait, P_value)

write.csv(
  significant_snps,
  "D:/Treatment/R_WorkSpace/projet-ble-selection/outputs/significant_snps.csv",
  row.names = FALSE
)

cat(
  "\nToutes les analyses GWAS sont terminées. Résultats sauvegardés dans data/gwas_results_all.RData et outputs/\n"
)
