# Load required libraries
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(RColorBrewer)
library(seqinr)
library(ggrepel)

# 1. PHYLOGENETIC TREE ANALYSIS # Read bootstrap support values
bootstrap_tree <- read.tree("T14.raxml.supportTBE")
tree_with_support <- ladderize(midpoint(bootstrap_tree))
tree_with_support$node.label <- round(100*as.numeric(tree_with_support$node.label),0)

# Define genetic groups based on sample names - corrected based on tree structure
group_A_samples <- c("KC951854.1", "MN072621.1", "MH381810.1",
  "MW020570.1", "MN072620.1", "PV167794.1")  # Group A (red cluster)
group_B_samples <- c("AY077836.1", "AY077835.1", "KX576657.1",
  "NC_004003.1", "MN072622.1")  # Group B (blue cluster)
group_C_samples <- c("MN072623.1", "MN072625.1", "MN072624.1")  # Group C (green cluster)

# Create grouping function
assign_groups <- function(tip_labels) {
  groups <- rep("Other", length(tip_labels))
  groups[tip_labels %in% group_A_samples] <- "I"
  groups[tip_labels %in% group_B_samples] <- "II" 
  groups[tip_labels %in% group_C_samples] <- "III"
  return(groups) }

# Assign groups to tips
tip_groups <- assign_groups(tree_with_support$tip.label)
names(tip_groups) <- tree_with_support$tip.label
group_colors <- c("I" = "grey", "II" = "purple", "III" = "orange")

# Create the phylogenetic tree plot
p1 <- ggtree(tree_with_support, layout="rectangular") +
  geom_tiplab(aes(color = tip_groups[label]), 
              size = 3, hjust = -0.1) +
  scale_color_manual(values = group_colors, name = "") +
  theme_tree2() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 5)) +
  xlim(0, max(node.depth.edgelength(tree_with_support)) * 1.3)

# Add scale bar (0.05 substitutions per site)
p1 <- p1 + geom_treescale(x = 0, y = 1, width = 0.05, 
                 label = "0.05 substitutions per site", fontsize = 1.5)

# Add bootstrap support for nodes ≥90%
p1 <- p1 + geom_text2(aes(subset = !isTip & as.numeric(label) >= 90, 
                         label = paste0(label, "%")), 
                     hjust = 1.2, vjust = -0.3, size = 1.9)
ggsave("goatpox_phylogeny.png", p1, width = 5, height = 5, dpi = 300)
ggsave("goatpox_phylogeny.pdf", p1, width = 5, height = 5)

# 2. PCA ANALYSIS
if (file.exists("genomes.aln")) {
  # Read alignment using seqinr
  alignment <- read.alignment("genomes.aln", format = "fasta")
  
  # Convert alignment to matrix
  seq_matrix <- as.matrix(alignment)
  
  # Convert DNA characters to numeric (A=1, T=2, G=3, C=4, gap/N=0)
  convert_to_numeric <- function(x) {
    x_clean <- toupper(as.character(x))
    x_clean[x_clean == "A"] <- 1
    x_clean[x_clean == "T"] <- 2
    x_clean[x_clean == "G"] <- 3
    x_clean[x_clean == "C"] <- 4
    x_clean[x_clean == "-" | x_clean == "N" | x_clean == "?" | x_clean == "X"] <- 0
    return(suppressWarnings(as.numeric(x_clean)))  }
  
  # Apply conversion to each sequence
  numeric_matrix <- t(apply(seq_matrix, 1, convert_to_numeric))
  
  # Remove columns with all NAs or constant values
  valid_cols <- apply(numeric_matrix, 2, function(x) {
    x_clean <- x[!is.na(x) & x != 0]  # Remove gaps and NAs
    length(unique(x_clean)) > 1 && length(x_clean) > 0  })
  
  numeric_matrix <- numeric_matrix[, valid_cols]
  
  # Replace remaining NAs with 0
  numeric_matrix[is.na(numeric_matrix)] <- 0
  
  # Perform PCA with scaling only if there's variance
  col_vars <- apply(numeric_matrix, 2, var, na.rm = TRUE)
  constant_cols <- col_vars == 0 | is.na(col_vars)
  
  # Remove constant columns
  numeric_matrix <- numeric_matrix[, !constant_cols]
  
  # Perform PCA without scaling if some columns still have zero variance
  pca_result <- prcomp(numeric_matrix, center = TRUE, scale. = FALSE)
  
  # Extract PC1 and PC2
  pca_data <- data.frame(    Sample = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],    PC2 = pca_result$x[, 2]  )
  
  # Assign groups to PCA data
  pca_data$Group <- assign_groups(pca_data$Sample)
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  # Create PCA plot
  p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.4) +
    geom_text_repel(aes(x = PC1,   y = PC2, label = Sample),
       nudge_x = 0.06, nudge_y = 0.06, size=3, max.overlaps = 999) +
    scale_color_manual(values = group_colors) +
    labs(      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = ""    ) +    theme_minimal() +
    theme(      legend.position = "bottom",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)    )
  ggsave("goatpox_pca.pdf", p2, width =var_explained[1]/10,
     height = var_explained[2]/10) } 

# 3. COMBINED VISUALIZATION (optional)
if (exists("p2")) {
  library(gridExtra)
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  ggsave("goatpox_combined_analysis.png", combined_plot, 
         width = 10, height = 4, dpi = 300)
  ggsave("goatpox_combined_analysis.pdf", combined_plot, 
         width = 10, height = 4) }
