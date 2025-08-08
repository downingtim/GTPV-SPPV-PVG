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
  # FIX: Use readGenBank() to parse the files into S4 objects, not just store the filenames
  sppv_gb <- genbankr::readGenBank("nc_004002.gb")
  gtpv_gb <- genbankr::readGenBank("nc_004003.gb")

  # Define regions of interest to highlight
  highlight_regions <- data.frame(
    start = c(32000, 65000, 81000, 118000),
    end = c(34000, 67000, 84000, 121000),
    label = c("32-34 Kb", "65-67 Kb", "81-84 Kb", "118-121 Kb")
  )

  # Create the main SNV density plot
  cat("Creating SNV density plot...\n")
  snv_plot <- ggplot(snv_density, aes(x = window_center/1000, y = snvs_per_kb, color=Sample)) +
    geom_rect(data = highlight_regions,
              aes(xmin = start/1000, xmax = end/1000, ymin = -Inf, ymax = Inf),
              fill = "darkorange", alpha = 0.5, inherit.aes = FALSE) +
    # FIX: Use 'linewidth' instead of the deprecated 'size' argument
    geom_line(linewidth = 1.4, alpha = 0.8) +
    scale_color_manual(values = c("SPPV" = "#E31A1C", "GTPV" = "#1F78B4")) +
    scale_x_continuous(breaks = seq(0, max(snv_density$window_center/1000), by = 10),
                       minor_breaks = seq(0, max(snv_density$window_center/1000), by = 5)) +
    labs(x = "Genome position (kb)", y = "SNVs per kb") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      # FIX: Use 'legend.position.inside' for numeric coordinates
      legend.position.inside = c(0.94, 0.98),
      legend.justification = c(1, 1),
      # FIX: Use 'linewidth' for rect borders
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      legend.margin = margin(8, 12, 8, 12),
      # FIX: Use 'linewidth' for grid lines
      panel.grid.minor = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.6),
      # FIX: Use 'linewidth' for the panel border
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
             # NOTE: This logic is static; it always colors CDS grey.
             # You may want to implement logic to check for overlaps with 'highlight_regions'.
             highlight = ifelse(0, "darkorange", "grey50")
      )

    cds_plot <- ggplot(cds_data, aes(xmin = start, xmax = end, ymin = 0, ymax = y)) +
      geom_rect(data = highlight_regions,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "darkorange", alpha = 0.7, inherit.aes = FALSE) +
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
            # FIX: Use 'linewidth' for grid lines
            panel.grid.minor = element_line(color = "gray85", linewidth = 0.3),
            panel.grid.major = element_line(color = "gray70", linewidth = 0.6),
            # FIX: Use 'linewidth' for the panel border
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

}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
})