# R script to create a focused comparison of SNV density in the
# first and last 2.5 kb of SPPV and GTPV genomes.

library(ggplot2)
library(dplyr)
library(readr)
library(genbankr)
library(patchwork)

# --- Data Loading and Processing Functions (Unchanged) ---

# Function to read and process VCF file
read_vcf_snvs <- function(vcf_path, sample_name) {
  vcf_lines <- readLines(vcf_path)
  header_line <- grep("^#CHROM", vcf_lines)
  vcf_data <- read.table(text = vcf_lines[(header_line+1):length(vcf_lines)],
                         sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
  snv_data <- vcf_data[grepl("TYPE=snv", vcf_data$INFO), ]
  snv_data$Sample <- sample_name
  return(snv_data[, c("CHROM", "POS", "Sample")])
}

# Function to calculate SNVs per kb using sliding windows
calculate_snv_density <- function(snv_data, window_size_kb = 0.5, step_size_bp = 50) {
  window_size_bp <- window_size_kb * 1000
  pos_range <- range(snv_data$POS)
  min_pos <- pos_range[1]
  max_pos <- pos_range[2]
  window_starts <- seq(from = min_pos, to = max_pos - window_size_bp + 1, by = step_size_bp)
  
  # Using purrr::map_dfr for a faster alternative to the for-loop with rbind
  results <- purrr::map_dfr(window_starts, function(start) {
    end <- start + window_size_bp - 1
    center <- start + (window_size_bp / 2)
    snvs_in_window <- sum(snv_data$POS >= start & snv_data$POS <= end)
    data.frame(
      CHROM = snv_data$CHROM[1],
      window_center = center,
      snvs_per_kb = snvs_in_window / window_size_kb,
      Sample = snv_data$Sample[1]
    )
  })
  return(results)
}

# --- Main Analysis ---

tryCatch({
  # --- 1. Load and Process All Data (as before) ---
  
  # File paths
  sppv_file <- "results_SPPV/vcf/gfavariants.vcf"
  gtpv_file <- "results_GTPV/vcf/gfavariants.vcf"
  
  # VCF data
  cat("Reading VCF files...\n")
  sppv_snvs <- read_vcf_snvs(sppv_file, "SPPV")
  gtpv_snvs <- read_vcf_snvs(gtpv_file, "GTPV")
  
  # SNV Density data
  cat("Calculating SNV density...\n")
  sppv_density <- calculate_snv_density(sppv_snvs)
  gtpv_density <- calculate_snv_density(gtpv_snvs)
  snv_density <- rbind(sppv_density, gtpv_density)
  
  # GenBank data
  cat("Reading GenBank files...\n")
  sppv_gb <- genbankr::readGenBank("nc_004002.gb")
  gtpv_gb <- genbankr::readGenBank("nc_004003.gb")
  
  
  # --- 2. Master Plotting Function for Zoomed Regions ---
  
  create_zoomed_plot <- function(region_coords, region_title) {
    
    # region_coords should be a vector like c(start, end)
    start_pos <- region_coords[1]
    end_pos <- region_coords[2]
    
    # A. SNV Density Plot for the region
    snv_plot <- ggplot(snv_density, aes(x = window_center / 1000, y = snvs_per_kb, color = Sample)) +
      geom_line(linewidth = 1.2, alpha = 0.8) +
      scale_color_manual(values = c("SPPV" = "#E31A1C", "GTPV" = "#1F78B4")) +
      labs(        x = "Genome position (Kb)",       y = "SNVs per Kb"
      ) + ylim(0,600) + 
      # ZOOM: Set coordinate limits for this specific region
      coord_cartesian(xlim = c(start_pos / 1000, end_pos / 1000)) +
      theme_minimal() +
      theme(
        legend.position = "none", # Legend will be added to the combined plot
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)      )
      
    # B. Function to create a CDS plot for the region
    create_cds_subplot <- function(gb_data, sample_name) {
      cds_data <- gb_data@cds %>%
        as.data.frame() %>%
        mutate(          y = ifelse(strand == "+", 1, -1)        )
      
      ggplot(cds_data, aes(xmin = start, xmax = end, ymin = 0, ymax = y)) +
        geom_rect(fill = "grey50", color = "black", linewidth = 0.2) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(y = sample_name) + ylim(0,600) + 
        # ZOOM: Set coordinate limits for this specific region
        coord_cartesian(xlim = region_coords) +
        scale_x_continuous(labels = scales::label_number(scale = 1/1000, suffix = "k")) +
        scale_y_continuous(breaks = c(-1, 1), labels = c("-", "+"), name = sample_name) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.title.y = element_text(size = 10, face = "bold", angle = 0, vjust = 0.5),
          panel.grid = element_blank()
        )
    }
    
    # C. Create CDS plots for SPPV and GTPV
    sppv_cds_plot <- create_cds_subplot(sppv_gb, "SPPV")
    gtpv_cds_plot <- create_cds_subplot(gtpv_gb, "GTPV")
    
    # D. Combine the plots for this region using patchwork
    combined_regional_plot <- snv_plot / sppv_cds_plot / gtpv_cds_plot +
      plot_layout(heights = c(4, 1, 1))
      
    return(combined_regional_plot)
  }
  
  # --- 3. Generate and Combine Plots for Both Regions ---
  
  max_genome_length <- 151000 # Define max length to calculate end region
  
  cat("Creating plot for the first 2.5kb...\n")
  plot_first_2500 <- create_zoomed_plot(
    region_coords = c(1, 5200),
    region_title = "First 5 Kb"  )
  
  cat("Creating plot for the last 2.5kb...\n")
  plot_last_2500 <- create_zoomed_plot(
    region_coords = c(max_genome_length - 6000, max_genome_length),
    region_title = "Last 5 Kb" )
  
  final_plot <- plot_first_2500 | plot_last_2500 &
    theme(legend.position = "bottom")

  final_plot <- final_plot +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))    )
  
  ggsave("snv_density_terminal_regions.png", plot = final_plot, width = 10, height = 5, dpi = 300)
  ggsave("snv_density_terminal_regions.pdf", plot = final_plot, width = 10, height = 5)
  
}, error = function(e) {
  cat("An error occurred:", e$message, "\n")
})