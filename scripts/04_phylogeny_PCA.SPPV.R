library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(RColorBrewer)
library(seqinr)
library(ggrepel)
library(gridExtra)

# 1. PHYLOGENETIC TREE ANALYSIS
region_names <- paste0("LSDV", sprintf("%03d", 1:156))

trim_names <- function(names) {
  sapply(names, function(name) {
    pos <- regexpr("\\.1", name)
    if (pos > 0) {
      substr(name, 1, pos + 1)
    } else {
      name   }
  }, USE.NAMES = FALSE)
}

# 2. Loop through each region name
for (region in region_names) {
  # 3. Construct the filenames for the current region
  tree_file <- paste0("GENE_TREE/", region, ".raxml.bestTree")
  genomes_aln <- paste0( "GENE_FASTA/",  region, ".fasta")
  
  # Check if files exist
  if (!file.exists(tree_file)) {
    cat("Tree file not found for region:", region, "\n")
    next  }
  
  bootstrap_tree <- read.tree(tree_file)
  bootstrap_tree$tip.label <- trim_names(bootstrap_tree$tip.label)
  tree_with_support <- ladderize(midpoint(bootstrap_tree))
  
  # Handle node labels safely
  if (!is.null(tree_with_support$node.label)) {
    tree_with_support$node.label <- round(100 * as.numeric(tree_with_support$node.label), 0)
    tree_with_support$node.label[is.na(tree_with_support$node.label)] <- 0  }
  
  # Define genetic groups based on sample names
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
    return(groups)  }
  
  # Assign groups to tips
  tip_groups <- assign_groups(tree_with_support$tip.label)
  names(tip_groups) <- tree_with_support$tip.label
  group_colors <- c("3.3" = "red", "3.2" = "blue", "3.1" = "darkgreen")
  
  # Create the phylogenetic tree plot
  p1 <- ggtree(tree_with_support, layout = "rectangular") +
    geom_tiplab(aes(color = tip_groups[label]),
                size = 3, hjust = -0.1) +
    scale_color_manual(values = group_colors, name = "") +
    theme_tree2() + ggtitle(paste("GTPV gene", region)) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 5)) +
    xlim(0, max(node.depth.edgelength(tree_with_support)) * 1.3)
  
  # Add bootstrap support for nodes ≥90%
  if (!is.null(tree_with_support$node.label)) {
    p1 <- p1 + geom_text2(aes(subset = !isTip & as.numeric(label) >= 90,
                             label = paste0(label, "%")),
                         hjust = 1.2, vjust = -0.3, size = 1.9)  }
  
  # 2. PCA ANALYSIS
  if (file.exists(genomes_aln)) {
    cat("Processing alignment for region:", region, "\n")
    
    # Read alignment using seqinr
    tryCatch({
      alignment <- read.alignment(genomes_aln, format = "fasta")
      alignment$nam <- trim_names(alignment$nam)
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
        return(suppressWarnings(as.numeric(x_clean)))
      }
      
      # Apply conversion to each sequence
      numeric_matrix <- t(apply(seq_matrix, 1, convert_to_numeric))
      
      # Remove columns with all NAs or constant values
      valid_cols <- apply(numeric_matrix, 2, function(x) {
        x_clean <- x[!is.na(x) & x != 0]  # Remove gaps and NAs
        length(unique(x_clean)) > 1 && length(x_clean) > 0
      })
      
      numeric_matrix <- numeric_matrix[, valid_cols]
      
      # Replace remaining NAs with 0
      numeric_matrix[is.na(numeric_matrix)] <- 0
      
      # Check if we have enough variable columns
      if (ncol(numeric_matrix) < 2) {
        cat("Not enough variable positions for PCA in region:", region, "\n")
        next
      }
      
      # Perform PCA with scaling only if there's variance
      col_vars <- apply(numeric_matrix, 2, var, na.rm = TRUE)
      constant_cols <- col_vars == 0 | is.na(col_vars)
      
      # Remove constant columns
      if (any(!constant_cols)) {
        numeric_matrix <- numeric_matrix[, !constant_cols]
      } else {
        cat("No variable columns found for region:", region, "\n")
        next
      }
      
      # Perform PCA
      pca_result <- prcomp(numeric_matrix, center = TRUE, scale. = FALSE)
      
      # Extract PC1 and PC2
      pca_data <- data.frame(
        Sample = rownames(pca_result$x),
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2]
      )
      
      # Assign groups to PCA data
      pca_data$Group <- assign_groups(pca_data$Sample)
      var_explained <- summary(pca_result)$importance[2, ] * 100
      
      # Create PCA plot
      p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_text_repel(aes(label = Sample),
                       nudge_x = 0.06, nudge_y = 0.06, size = 3, max.overlaps = 999) +
        scale_color_manual(values = group_colors) +
        labs( title=paste("GTPV gene", region), 
          x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
          y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
          color = ""
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      # Save PCA plot with fixed dimensions (minimum 4x4 inches)
      pca_width <- max(4, min(12, var_explained[1] / 5))
      pca_height <- max(4, min(12, var_explained[2] / 5))
      
      # 3. COMBINED VISUALIZATION
      combined_plot <- grid.arrange(p1, p2, ncol = 2)
      ggsave(paste0("GENOME/", region, "/", region, "_SPPV_combined.pdf"), combined_plot,
             width = 10, height = 4)
      
    }, error = function(e) {
      cat("Error processing alignment for region:", region, ":", e$message, "\n")
    })
  } else {
    cat("Alignment file not found for region:", region, "\n")
  }
  
  cat("Completed region:", region, "\n\n")
}