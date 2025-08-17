# VCF Analysis: Computing SNP density and diversity metrics from sites-only VCF
# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Function to read VCF file manually
read_vcf_manual <- function(file_path) {
  # Read all lines
  lines <- readLines(file_path)
  
  # Find header line (starts with #CHROM)
  header_idx <- which(grepl("^#CHROM", lines))
  
  # Extract header
  header <- strsplit(lines[header_idx], "\t")[[1]]
  header <- gsub("^#", "", header)  # Remove # from CHROM
  
  # Extract data lines (non-comment lines after header)
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]
  
  # Parse data
  vcf_data <- do.call(rbind, lapply(data_lines, function(x) strsplit(x, "\t")[[1]]))
  vcf_df <- as.data.frame(vcf_data, stringsAsFactors = FALSE)
  colnames(vcf_df) <- header
  
  # Convert POS to numeric
  vcf_df$POS <- as.numeric(vcf_df$POS)
  
  return(vcf_df)
}

# Read the VCF file
vcf_file <- "vcf/gfavariants.vcf"
cat("Reading VCF file:", vcf_file, "\n")

vcf_data <- read_vcf_manual(vcf_file)
cat("VCF data loaded successfully\n")
cat("Dimensions:", nrow(vcf_data), "variants x", ncol(vcf_data), "columns\n")
cat("Column names:", paste(colnames(vcf_data), collapse = ", "), "\n")

# Check if this is a sites-only VCF (standard columns only)
standard_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
has_samples <- any(!colnames(vcf_data) %in% standard_cols)

if (!has_samples) {
  cat("This appears to be a sites-only VCF file (no sample genotype data)\n")
  cat("Will calculate SNP density and basic diversity metrics from variant sites\n")
} else {
  cat("Sample columns found:", paste(colnames(vcf_data)[!colnames(vcf_data) %in% standard_cols], collapse = ", "), "\n")
}

# Filter for SNVs only (TYPE=snv in INFO column)
snv_indices <- grepl("TYPE=snv", vcf_data$INFO)
snv_data <- vcf_data[snv_indices, ]

cat("Filtered to", sum(snv_indices), "SNVs\n")

# Extract positions
positions <- snv_data$POS

# Define regions
max_pos <- max(positions)
cat("Position range:", min(positions), "to", max_pos, "\n")

regions <- list(
  "5_prime_accessory" = c(100, 13850),
  "core" = c(13851, 106910),
  "3_prime_accessory" = c(106911, max_pos),
  "genome_wide" = c(1, max_pos)
)

# Function to estimate theta from segregating sites (Watterson's estimator)
# For sites-only VCF, we'll estimate assuming a standard sample size
estimate_theta_watterson <- function(n_snps, region_length, n_samples = 10) {
  if (n_samples <= 1 || n_snps == 0) return(0)
  
  # Calculate harmonic number a_n
  harmonic_n <- sum(1 / (1:(n_samples - 1)))
  
  # Watterson's theta = S / (a_n * L)
  theta_w <- n_snps / (harmonic_n * region_length)
  
  return(theta_w * 1000)  # per Kb
}

# Function to estimate nucleotide diversity from biallelic sites
# This is a rough estimate based on expected heterozygosity
estimate_pi_from_sites <- function(snv_subset, region_length, assumed_maf = 0.2) {
  if (nrow(snv_subset) == 0) return(0)
  
  # For each biallelic site, estimate pi assuming Hardy-Weinberg equilibrium
  # pi = 2 * p * (1-p) where p is the minor allele frequency
  # We'll assume an average MAF for estimation
  
  n_sites <- nrow(snv_subset)
  
  # Average nucleotide diversity per site assuming MAF
  avg_pi_per_site <- 2 * assumed_maf * (1 - assumed_maf)
  
  # Total pi for the region
  total_pi <- avg_pi_per_site * n_sites / region_length
  
  return(total_pi * 1000)  # per Kb
}

# Function to analyze region
analyze_region <- function(region_name, start_pos, end_pos) {
  # Filter SNVs in this region
  region_indices <- which(positions >= start_pos & positions <= end_pos)
  
  if (length(region_indices) == 0) {
    return(data.frame(
      Region = region_name,
      Start = start_pos,
      End = end_pos,
      Length_bp = end_pos - start_pos + 1,
      N_SNPs = 0,
      SNPs_per_Kb = 0,
      Pi_per_Kb_est = 0,
      Theta_W_per_Kb_est = 0,
      Analysis_type = "sites_only"
    ))
  }
  
  region_snvs <- snv_data[region_indices, ]
  n_snps <- length(region_indices)
  region_length <- end_pos - start_pos + 1
  
  # Calculate metrics
  snps_per_kb <- (n_snps / region_length) * 1000
  
  # Estimate pi (assuming MAF = 0.2 for conservative estimate)
  pi_per_kb_est <- estimate_pi_from_sites(region_snvs, region_length, assumed_maf = 0.2)
  
  # Estimate Watterson's theta (assuming 10 samples for conservative estimate)
  theta_w_per_kb_est <- estimate_theta_watterson(n_snps, region_length, n_samples = 30)
  
  return(data.frame(
    Region = region_name,
    Start = start_pos,
    End = end_pos,
    Length_bp = region_length,
    N_SNPs = n_snps,
    SNPs_per_Kb = round(snps_per_kb, 4),
    Pi_per_Kb_est = round(pi_per_kb_est, 6),
    Theta_W_per_Kb_est = round(theta_w_per_kb_est, 6),
    Analysis_type = "sites_only"
  ))
}

# Analyze all regions
results <- data.frame()

cat("\n=== Analyzing regions ===\n")
for (region_name in names(regions)) {
  region_bounds <- regions[[region_name]]
  cat("Analyzing", region_name, "region (", region_bounds[1], "-", region_bounds[2], ")\n")
  result <- analyze_region(region_name, region_bounds[1], region_bounds[2])
  results <- rbind(results, result)
}

# Display results
cat("\n=== VCF Analysis Results (Sites-Only Analysis) ===\n")
cat("Total SNVs found:", nrow(snv_data), "\n")
cat("Note: Pi and Theta estimates are based on assumed parameters since no genotype data is available\n")
cat("Pi estimates assume MAF = 0.2; Theta estimates assume n = 14 samples\n\n")

print(results)

# Create visualization
tryCatch({
  # Reshape data for plotting (exclude genome-wide for regional comparison)
  plot_data <- results %>%
    filter(Region != "genome_wide") %>%
    select(Region, SNPs_per_Kb, Pi_per_Kb_est, Theta_W_per_Kb_est) %>%
    pivot_longer(cols = c(SNPs_per_Kb, Pi_per_Kb_est, Theta_W_per_Kb_est), 
                 names_to = "Metric", values_to = "Value")
  
  # Clean up names
  plot_data$Region <- gsub("_", " ", plot_data$Region)
  plot_data$Region <- tools::toTitleCase(plot_data$Region)
  
  plot_data$Metric <- gsub("_", " ", plot_data$Metric)
  plot_data$Metric <- gsub("per Kb est", "per Kb (est)", plot_data$Metric)
  
  # Create plot
  p1 <- ggplot(plot_data, aes(x = Region, y = Value, fill = Region)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Genomic Diversity Metrics by Region\n(Estimates from Sites-Only VCF)",
         x = "Genomic Region", y = "Value per Kb") +
    scale_fill_brewer(palette = "Set2")
  
  print(p1)
}, error = function(e) {
  cat("Could not create plot:", e$message, "\n")
})

# Summary statistics
cat("\n=== Summary Statistics ===\n")
if (nrow(results) > 0) {
  max_snp_idx <- which.max(results$SNPs_per_Kb)
  max_pi_idx <- which.max(results$Pi_per_Kb_est)
  max_theta_idx <- which.max(results$Theta_W_per_Kb_est)
  
  cat("Region with highest SNP density:", results$Region[max_snp_idx], 
      "(", results$SNPs_per_Kb[max_snp_idx], "SNPs/Kb)\n")
  cat("Region with highest estimated nucleotide diversity (π):", results$Region[max_pi_idx], 
      "(", results$Pi_per_Kb_est[max_pi_idx], "π/Kb)\n")
  cat("Region with highest estimated Watterson's θ:", results$Region[max_theta_idx], 
      "(", results$Theta_W_per_Kb_est[max_theta_idx], "θ/Kb)\n")
}

# Additional analysis: SNP distribution
cat("\n=== SNP Distribution Analysis ===\n")
regional_comparison <- results %>%
  filter(Region != "genome_wide") %>%
  arrange(desc(SNPs_per_Kb))

cat("Regional SNP density ranking:\n")
for (i in 1:nrow(regional_comparison)) {
  cat(i, ". ", regional_comparison$Region[i], ": ", 
      regional_comparison$SNPs_per_Kb[i], " SNPs/Kb\n", sep = "")
}

# Calculate SNP density ratios
core_snps <- results$SNPs_per_Kb[results$Region == "core"]
accessory_5_snps <- results$SNPs_per_Kb[results$Region == "5_prime_accessory"]
accessory_3_snps <- results$SNPs_per_Kb[results$Region == "3_prime_accessory"]

if (core_snps > 0) {
  cat("\nSNP density ratios relative to core:\n")
  cat("5' accessory / core:", round(accessory_5_snps / core_snps, 2), "\n")
  cat("3' accessory / core:", round(accessory_3_snps / core_snps, 2), "\n")
}

# Export results
write.csv(results, "vcf_analysis_results.csv", row.names = FALSE)
cat("\nResults exported to 'vcf_analysis_results.csv'\n")

# Export SNP positions by region for further analysis
cat("\nExporting SNP positions by region...\n")
for (region_name in names(regions)[1:3]) {  # Skip genome-wide
  region_bounds <- regions[[region_name]]
  region_positions <- positions[positions >= region_bounds[1] & positions <= region_bounds[2]]
  
  if (length(region_positions) > 0) {
    filename <- paste0("snp_positions_", region_name, ".txt")
    writeLines(as.character(region_positions), filename)
    cat("Exported", length(region_positions), "SNP positions for", region_name, "to", filename, "\n")
  }
}