# R script to compare SNVs/Kb between SPPV and GTPV samples
# Reads VCF files and GenBank files, creates a comparative plot with CDS annotations

library(ggplot2)
library(dplyr)
library(readr)
library(genbankr)
library(patchwork)

# Function to read and process VCF file
read_vcf_snvs <- function(vcf_path, sample_name) {
  # Read VCF file, skipping header lines
  vcf_lines <- readLines(vcf_path)

  # Find the header line (starts with #CHROM)
  header_line <- grep("^#CHROM", vcf_lines)

  vcf_data <- read.table(text = vcf_lines[(header_line+1):length(vcf_lines)],
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

  # Filter for SNVs only
  snv_data <- vcf_data[grepl("TYPE=snv", vcf_data$INFO), ]

  # Add sample name
  snv_data$Sample <- sample_name

  return(snv_data[, c("CHROM", "POS", "Sample")])
}

# Function to calculate SNVs per kb using sliding windows
calculate_snv_density <- function(snv_data, window_size_kb = 2, step_size_bp = 50) {
  window_size_bp <- window_size_kb * 1000

  # Get the range of positions for this sample
  pos_range <- range(snv_data$POS)
  min_pos <- pos_range[1]
  max_pos <- pos_range[2]

  # Create sliding windows
  window_starts <- seq(from = min_pos, to = max_pos - window_size_bp + 1, by = step_size_bp)

  # Calculate SNV counts for each window
  window_results <- data.frame()

  for(i in 1:length(window_starts)) {
    window_start <- window_starts[i]
    window_end <- window_start + window_size_bp - 1
    window_center <- window_start + (window_size_bp / 2)

    # Count SNVs in this window
    snvs_in_window <- sum(snv_data$POS >= window_start & snv_data$POS <= window_end)
    snvs_per_kb <- snvs_in_window / window_size_kb

    # Store results
    window_results <- rbind(window_results, data.frame(
      CHROM = snv_data$CHROM[1],
      window_center = window_center,
      snv_count = snvs_in_window,
      snvs_per_kb = snvs_per_kb,
      Sample = snv_data$Sample[1]
    ))
  }
  return(window_results)
}

# Function to calculate SNVs per kb for specific genomic regions
calculate_region_snv_density <- function(snv_data, regions_df, sample_name) {
  results <- data.frame()
  
  for(i in 1:nrow(regions_df)) {
    start_pos <- regions_df$start[i]
    end_pos <- regions_df$end[i]
    gene_label <- regions_df$label[i]
    
    # Count SNVs within this region
    snvs_in_region <- sum(snv_data$POS >= start_pos & snv_data$POS <= end_pos)
    
    # Calculate region length in kb
    region_length_bp <- end_pos - start_pos + 1
    region_length_kb <- region_length_bp / 1000
    
    # Calculate SNVs per kb
    snvs_per_kb <- snvs_in_region / region_length_kb
    
    # Store results
    results <- rbind(results, data.frame(
      Gene = paste0("LSDV", sprintf("%03d", as.numeric(gene_label))),
      Start = start_pos,
      End = end_pos,
      Length_bp = region_length_bp,
      Length_kb = round(region_length_kb, 3),
      SNV_count = snvs_in_region,
      SNVs_per_kb = round(snvs_per_kb, 2),
      Sample = sample_name,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# Function to perform random resampling of genomic regions
random_region_sampling <- function(snv_data, genome_length, n_samples = 1000, seed = 123) {
  set.seed(seed)
  
  # Calculate average gene length from our highlight regions
  avg_gene_length <- mean(all_highlight_regions$end - all_highlight_regions$start + 1)
  
  # Store results
  random_snvs_per_kb <- numeric(n_samples)
  
  for(i in 1:n_samples) {
    # Generate random start position (ensuring we don't go beyond genome)
    max_start <- genome_length - avg_gene_length + 1
    random_start <- sample(1:max_start, 1)
    random_end <- random_start + avg_gene_length - 1
    
    # Count SNVs in this random region
    snvs_in_region <- sum(snv_data$POS >= random_start & snv_data$POS <= random_end)
    
    # Calculate SNVs per kb
    region_length_kb <- avg_gene_length / 1000
    random_snvs_per_kb[i] <- snvs_in_region / region_length_kb
  }
  
  return(random_snvs_per_kb)
}

# Main analysis
tryCatch({
  # File paths
  sppv_file <- "results_SPPV/vcf/gfavariants.vcf"
  gtpv_file <- "results_GTPV/vcf/gfavariants.vcf"

  # Read VCF files
  cat("Reading SPPV VCF file...\n")
  sppv_snvs <- read_vcf_snvs(sppv_file, "SPPV")

  cat("Reading GTPV VCF file...\n")
  gtpv_snvs <- read_vcf_snvs(gtpv_file, "GTPV")

  # Combine data
  all_snvs <- rbind(sppv_snvs, gtpv_snvs)

  # Print summary statistics
  cat("\nSummary statistics:\n")
  cat("SPPV SNVs:", nrow(sppv_snvs), "\n")
  cat("GTPV SNVs:", nrow(gtpv_snvs), "\n")

  # Get position ranges for each sample
  sppv_range <- range(sppv_snvs$POS)
  gtpv_range <- range(gtpv_snvs$POS)
  cat("SPPV position range:", sppv_range[1], "-", sppv_range[2], "\n")
  cat("GTPV position range:", gtpv_range[1], "-", gtpv_range[2], "\n")

  # Calculate SNV density using sliding windows (1kb windows, 100bp steps)
  cat("\nCalculating SNV density using sliding windows (1kb windows, 100bp steps)...\n")

  # Process each sample separately to handle different position ranges
  sppv_density <- calculate_snv_density(sppv_snvs, window_size_kb = 2, step_size_bp = 50)
  gtpv_density <- calculate_snv_density(gtpv_snvs, window_size_kb = 2, step_size_bp = 50)

  # Combine results
  snv_density <- rbind(sppv_density, gtpv_density)

  # Read GenBank files for CDS annotations
  cat("Reading GenBank files...\n")
  sppv_gb <- genbankr::readGenBank("nc_004002.gb")
  gtpv_gb <- genbankr::readGenBank("nc_004003.gb")

  # Define regions of interest to highlight
  highlight_regions <- data.frame(
  start = c(390, 1768, 5427, 7422, 17199, 118888, 128515, 147380, 148325),
  end = c(1231, 2342, 6182, 8557, 18470, 119899, 129451, 147889, 149466),
  label = c("002", "004", "009", "013", "026", "132", "136", "153","155"))

  highlight_regions2 <- data.frame(
  start = c(7759, 11558, 56791, 121407, 132499, 134652, 136343, 139663, 141206, 144149),
  end = c(8594, 13246, 57381, 127430, 133176, 136295, 138247, 141159, 142549, 145807),
  label = c("012", "019", "067", "134", "141", "144", "145", "147", "148", "151")
)

  # Combine both highlight regions for complete analysis
  all_highlight_regions <- rbind(highlight_regions, highlight_regions2)

  # Create the main SNV density plot
  cat("Creating SNV density plot...\n")
  snv_plot <- ggplot(snv_density, aes(x = window_center/1000, y = snvs_per_kb, color=Sample)) +
    geom_rect(data = highlight_regions,
              aes(xmin = start/1000, xmax = end/1000, ymin = -Inf, ymax = Inf),
              fill = "darkorange", alpha = 0.5, inherit.aes = FALSE) +
    geom_rect(data = highlight_regions2,
              aes(xmin = start/1000, xmax = end/1000, ymin = -Inf, ymax = Inf),
              fill = "purple", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(linewidth = 1.4, alpha = 0.8) +
    scale_color_manual(values = c("SPPV" = "#E31A1C", "GTPV" = "#1F78B4")) +
    scale_x_continuous(breaks = seq(0, max(snv_density$window_center/1000), by = 10),
                       minor_breaks = seq(0, max(snv_density$window_center/1000), by = 5)) +
    labs(x = "Genome position (kb)", y = "SNPs/Kb") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position.inside = c(0.94, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      legend.margin = margin(8, 12, 8, 12),
      panel.grid.minor = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.6),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

  # Function to create CDS plot
  create_cds_plot <- function(gb_data, sample_name, genome_length) {
    cds_data <- gb_data@cds %>%
      as.data.frame() %>%
      mutate(start = start,
             end = end,
             strand = as.character(strand),
             y = ifelse(strand == "+", 1, -1),
             highlight = ifelse(0, "darkorange", "grey50")      )

    cds_plot <- ggplot(cds_data, aes(xmin = start, xmax = end, ymin = 0, ymax = y)) +
      geom_rect(data = highlight_regions,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "darkorange", alpha = 0.7, inherit.aes = FALSE) +
      geom_rect(data = highlight_regions2,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "purple", alpha = 0.7, inherit.aes = FALSE) +
      geom_rect(aes(fill = highlight), color = "black", linewidth = 0.15) +
      scale_fill_identity() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
      theme_minimal(base_size = 10) +
      labs(x = "Genome position (Kb)", y = sample_name) +
      scale_x_continuous(breaks = seq(0, genome_length, by = 10000),
                         minor_breaks = seq(0, genome_length, by = 5000),
                         labels = seq(0, genome_length, by = 10000) / 1000) +
      scale_y_continuous(breaks = c(-1, 1), labels = c("-", "+")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_text(size = 10, face = "bold"),
            panel.grid.minor = element_line(color = "gray85", linewidth = 0.3),
            panel.grid.major = element_line(color = "gray70", linewidth = 0.6),
            panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))
    return(cds_plot)
  }

  # Create CDS plots for both samples
  max_genome_length <- 151000
  sppv_cds_plot <- create_cds_plot(sppv_gb, "SPPV", max_genome_length)
  gtpv_cds_plot <- create_cds_plot(gtpv_gb, "GTPV", max_genome_length)

  # Combine all plots
  cat("Combining plots...\n")
  final_plot <- snv_plot / sppv_cds_plot / gtpv_cds_plot +
    plot_layout(heights = c(4, 1.1, 1.1))

  # Save the combined plot
  ggsave("snv_density_with_cds.png", plot = final_plot, width = 10, height = 6, dpi = 300)
  ggsave("snv_density_with_cds.pdf", plot = final_plot, width = 10, height = 6)

  density_summary <- snv_density %>%
    group_by(Sample) %>%
    summarise(
      mean_snvs_per_kb = mean(snvs_per_kb),
      median_snvs_per_kb = median(snvs_per_kb),
      max_snvs_per_kb = max(snvs_per_kb),
      total_bins = n(),
      .groups = 'drop'
    )
  print(density_summary)

  # ========== Gene-based SNV Analysis ==========

  # Calculate SNV density for each gene region for both samples
  cat("\nCalculating SNV density for gene regions...\n")

  # GTPV analysis
  gtpv_gene_snvs <- calculate_region_snv_density(gtpv_snvs, all_highlight_regions, "GTPV")

  # SPPV analysis  
  sppv_gene_snvs <- calculate_region_snv_density(sppv_snvs, all_highlight_regions, "SPPV")

  # Create formatted tables
  cat("\n=== Table 1: GTPV Gene SNV Analysis ===\n")
  gtpv_table <- gtpv_gene_snvs[, c("Gene", "Start", "End", "Length_kb", "SNV_count", "SNVs_per_kb")]
  gtpv_table <- gtpv_table[order(as.numeric(gsub("LSDV", "", gtpv_table$Gene))), ]
  print(gtpv_table, row.names = FALSE)

  cat("\n=== Table 2: SPPV Gene SNV Analysis ===\n")
  sppv_table <- sppv_gene_snvs[, c("Gene", "Start", "End", "Length_kb", "SNV_count", "SNVs_per_kb")]
  sppv_table <- sppv_table[order(as.numeric(gsub("LSDV", "", sppv_table$Gene))), ]
  print(sppv_table, row.names = FALSE)

  # Save tables to CSV files
  write.csv(gtpv_table, "GTPV_gene_SNV_analysis.csv", row.names = FALSE)
  write.csv(sppv_table, "SPPV_gene_SNV_analysis.csv", row.names = FALSE)

  # Create a comparison table
  cat("\n=== Comparison Table: GTPV vs SPPV SNVs per kb ===\n")
  comparison_table <- merge(
    gtpv_table[, c("Gene", "SNVs_per_kb")], 
    sppv_table[, c("Gene", "SNVs_per_kb")], 
    by = "Gene", 
    suffixes = c("_GTPV", "_SPPV")
  )
  comparison_table$Difference <- comparison_table$SNVs_per_kb_SPPV - comparison_table$SNVs_per_kb_GTPV
  comparison_table <- comparison_table[order(as.numeric(gsub("LSDV", "", comparison_table$Gene))), ]
  print(comparison_table, row.names = FALSE)

  # Save comparison table
  write.csv(comparison_table, "SNV_comparison_GTPV_vs_SPPV.csv", row.names = FALSE)

  # Summary statistics for gene regions
  cat("\n=== Summary Statistics ===\n")
  cat("GTPV genes - Mean SNVs/kb:", round(mean(gtpv_table$SNVs_per_kb), 2), "\n")
  cat("GTPV genes - Median SNVs/kb:", round(median(gtpv_table$SNVs_per_kb), 2), "\n")
  cat("SPPV genes - Mean SNVs/kb:", round(mean(sppv_table$SNVs_per_kb), 2), "\n")
  cat("SPPV genes - Median SNVs/kb:", round(median(sppv_table$SNVs_per_kb), 2), "\n")

  # Identify genes with highest SNV density
  if(nrow(gtpv_table) > 0) {
    cat("\nTop 5 genes with highest SNV density in GTPV:\n")
    top_gtpv <- gtpv_table[order(gtpv_table$SNVs_per_kb, decreasing = TRUE)[1:min(5, nrow(gtpv_table))], c("Gene", "SNVs_per_kb")]
    print(top_gtpv, row.names = FALSE)
  }

  if(nrow(sppv_table) > 0) {
    cat("\nTop 5 genes with highest SNV density in SPPV:\n")
    top_sppv <- sppv_table[order(sppv_table$SNVs_per_kb, decreasing = TRUE)[1:min(5, nrow(sppv_table))], c("Gene", "SNVs_per_kb")]
    print(top_sppv, row.names = FALSE)
  }

  # ========== Random Resampling Analysis ==========
  
  # Calculate overall SNV density (total SNVs / genome length in kb)
  genome_length <- 151000  # Approximate genome length in bp
  genome_length_kb <- genome_length / 1000
  
  gtpv_overall_density <- nrow(gtpv_snvs) / genome_length_kb
  sppv_overall_density <- nrow(sppv_snvs) / genome_length_kb
  
  cat("\n=== Overall Genome SNV Density ===\n")
  cat("GTPV overall SNVs/Kb:", round(gtpv_overall_density, 3), "\n")
  cat("SPPV overall SNVs/Kb:", round(sppv_overall_density, 3), "\n")
  
  # Perform random resampling
  cat("\n=== Random Region Resampling Analysis ===\n")
  cat("Performing random resampling (1000 iterations)...\n")
  
  gtpv_random_snvs <- random_region_sampling(gtpv_snvs, genome_length)
  sppv_random_snvs <- random_region_sampling(sppv_snvs, genome_length)
  
  # Calculate statistics for random regions
  gtpv_random_mean <- mean(gtpv_random_snvs)
  gtpv_random_sd <- sd(gtpv_random_snvs)
  gtpv_random_median <- median(gtpv_random_snvs)
  
  sppv_random_mean <- mean(sppv_random_snvs)
  sppv_random_sd <- sd(sppv_random_snvs)
  sppv_random_median <- median(sppv_random_snvs)
  
  cat("\nGTPV Random Regions (n=1000):\n")
  cat("  Mean SNVs/Kb:", round(gtpv_random_mean, 3), "\n")
  cat("  Standard Deviation:", round(gtpv_random_sd, 3), "\n")
  cat("  Median SNVs/Kb:", round(gtpv_random_median, 3), "\n")
  cat("  Range:", round(min(gtpv_random_snvs), 3), "-", round(max(gtpv_random_snvs), 3), "\n")
  
  cat("\nSPPV Random Regions (n=1000):\n")
  cat("  Mean SNVs/Kb:", round(sppv_random_mean, 3), "\n")
  cat("  Standard Deviation:", round(sppv_random_sd, 3), "\n")
  cat("  Median SNVs/Kb:", round(sppv_random_median, 3), "\n")
  cat("  Range:", round(min(sppv_random_snvs), 3), "-", round(max(sppv_random_snvs), 3), "\n")
  
  # Compare gene regions to random expectation
  cat("\n=== Comparison: Gene Regions vs Random Expectation ===\n")
  
  # Calculate how many standard deviations away from random mean
  gtpv_gene_mean <- mean(gtpv_table$SNVs_per_kb)
  sppv_gene_mean <- mean(sppv_table$SNVs_per_kb)
  
  gtpv_z_score <- (gtpv_gene_mean - gtpv_random_mean) / gtpv_random_sd
  sppv_z_score <- (sppv_gene_mean - sppv_random_mean) / sppv_random_sd
  
  cat("GTPV:\n")
  cat("  Gene regions mean:", round(gtpv_gene_mean, 3), "SNVs/Kb\n")
  cat("  Random regions mean:", round(gtpv_random_mean, 3), "SNVs/Kb\n")
  cat("  Difference:", round(gtpv_gene_mean - gtpv_random_mean, 3), "SNVs/Kb\n")
  cat("  Z-score:", round(gtpv_z_score, 2), "\n")
  cat("  Interpretation:", if(abs(gtpv_z_score) > 2) "Significantly different" else "Not significantly different", "from random\n")
  
  cat("\nSPPV:\n")
  cat("  Gene regions mean:", round(sppv_gene_mean, 3), "SNVs/Kb\n")
  cat("  Random regions mean:", round(sppv_random_mean, 3), "SNVs/Kb\n")
  cat("  Difference:", round(sppv_gene_mean - sppv_random_mean, 3), "SNVs/Kb\n")
  cat("  Z-score:", round(sppv_z_score, 2), "\n")
  cat("  Interpretation:", if(abs(sppv_z_score) > 2) "Significantly different" else "Not significantly different", "from random\n")
  
  # Save random sampling results
  random_results <- data.frame(
    Sample = rep(c("GTPV", "SPPV"), each = 1000),
    Random_SNVs_per_Kb = c(gtpv_random_snvs, sppv_random_snvs)
  )
  write.csv(random_results, "random_sampling_results.csv", row.names = FALSE)
  
  # Summary table
  summary_comparison <- data.frame(
    Sample = c("GTPV", "SPPV"),
    Overall_SNVs_per_Kb = c(gtpv_overall_density, sppv_overall_density),
    Gene_Regions_Mean = c(gtpv_gene_mean, sppv_gene_mean),
    Random_Regions_Mean = c(gtpv_random_mean, sppv_random_mean),
    Random_Regions_SD = c(gtpv_random_sd, sppv_random_sd),
    Z_Score = c(gtpv_z_score, sppv_z_score),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== Summary Comparison Table ===\n")
  print(summary_comparison, row.names = FALSE)
  write.csv(summary_comparison, "SNV_density_summary_comparison.csv", row.names = FALSE)

}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
})