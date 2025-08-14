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
region_names <- c("LSDV009", "LSDV013", "LSDV026", "LSDV136", "LSDV132")

# 2. Loop through each region name
for (region in region_names) {

# 3. Construct the filenames for the current region
  tree_file <- paste0(region, ".raxml.supportTBE")
  genomes_aln <- paste0("region_", region, ".aln")
  bootstrap_tree <- read.tree(tree_file)

tree_with_support <- ladderize(midpoint(bootstrap_tree))
tree_with_support$node.label <- round(100*as.numeric(tree_with_support$node.label),0)

# Define genetic groups based on sample names - corrected based on tree structure
group_A_samples <- c("AY077833.1",
"AY077834.1",
"MN072630.1",
"MN072631.1",
"MT137384.1",
"MW020571.1",
"MW167070.1",
"ON961655.1",
"ON961656.1",
"ON961657.1",
"OQ434235.1",
"OQ434236.1",
"OQ434237.1",
"OQ434238.1",
"OQ434239.1",
"PQ014465.1",
"PV167793.1",
"PV434148.1"
)
group_B_samples <- c("MN072629.1",
"MW167071.1",
"OR239060.1",
"AY077832.1",
"NC_004002.1",
"PP886236.1",
"PP886237.1",
"PP886238.1",
"PP886239.1"
)
group_C_samples <- c("MN072626.1", "MN072627.1")

# Create grouping function
assign_groups <- function(tip_labels) {
  groups <- rep("Other", length(tip_labels))
  groups[tip_labels %in% group_A_samples] <- "3.3"
  groups[tip_labels %in% group_B_samples] <- "3.2" 
  groups[tip_labels %in% group_C_samples] <- "3.1"
  return(groups) }

# Assign groups to tips
tip_groups <- assign_groups(tree_with_support$tip.label)
names(tip_groups) <- tree_with_support$tip.label
group_colors <- c("3.3" = "red", "3.2" = "blue", "3.1" = "darkgreen")

# Create the phylogenetic tree plot
p1 <- ggtree(tree_with_support, layout="rectangular") +
  geom_tiplab(aes(color = tip_groups[label]), 
              size = 3, hjust = -0.1) +
  scale_color_manual(values = group_colors, name = "") +
  theme_tree2() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 5)) +
  xlim(0, max(node.depth.edgelength(tree_with_support)) * 1.3)

# Add bootstrap support for nodes ≥90%
p1 <- p1 + geom_text2(aes(subset = !isTip & as.numeric(label) >= 90, 
                         label = paste0(label, "%")), 
                     hjust = 1.2, vjust = -0.3, size = 1.9)
ggsave(paste0(region,"_GTPV_phylogeny.pdf"), p1, width = 5, height = 5)

# 2. PCA ANALYSIS
if (file.exists(genomes_aln)) {
  # Read alignment using seqinr
  alignment <- read.alignment(genomes_aln, format = "fasta")
  
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
  ggsave(paste0(region,"_GTPV_pca.pdf"), p2, width =var_explained[1]/10,
     height = var_explained[2]/10) } 

# 3. COMBINED VISUALIZATION (optional)
if (exists("p2")) {
  library(gridExtra)
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  ggsave(paste0(region, "_GTPV_combined_.pdf"), combined_plot, 
         width = 10, height = 4) }

}