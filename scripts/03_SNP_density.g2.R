library(vcfR)
library(dplyr)
library(tidyr)

# 1. Setup Clade Mapping (Merging 3.2.1/3.2.2 and excluding MN072628.1)
mapping_data <- "
Accession,Clade
KC951854.1,2.1; MH381810.1,2.1; MN072620.1,2.1; MN072621.1,2.1; MW020570.1,2.1; PV167794.1,2.1
AY077835.1,2.2; AY077836.1,2.2; KX576657.1,2.2; MN072622.1,2.2; NC_004003.1,2.2
MN072623.1,2.3; MN072624.1,2.3; MN072625.1,2.3
MN072626.1,3.1; MN072627.1,3.1
AY077833.1,3.3; AY077834.1,3.3; MN072630.1,3.3; MN072631.1,3.3; MT137384.1,3.3; MW020571.1,3.3; 
MW167070.1,3.3; ON961655.1,3.3; ON961656.1,3.3; ON961657.1,3.3; OQ434235.1,3.3; OQ434236.1,3.3; 
OQ434237.1,3.3; OQ434238.1,3.3; OQ434239.1,3.3; PQ014465.1,3.3; PV167793.1,3.3; PV434148.1,3.3
MN072629.1,3.2; MW167071.1,3.2; OR239060.1,3.2; AY077832.1,3.2; NC_004002.1,3.2; 
PP886236.1,3.2; PP886237.1,3.2; PP886238.1,3.2; PP886239.1,3.2"

# Clean up mapping string into a proper data frame
clade_map <- read.table(text=gsub(";", "\n", mapping_data), sep=",", 
                        strip.white=TRUE, col.names=c("Accession", "Clade"), stringsAsFactors=FALSE)

# 2. Define Genomic Regions
regions <- list(
  "Whole" = c(start = 1, end = 150000), # Adjust end to your actual max length
  "5_prime" = c(start = 100, end = 13850),
  "Core" = c(start = 13851, end = 106910),
  "3_prime" = c(start = 106911, end = 150000)
)

# 3. Calculation Helper Function
calc_stats <- function(data, group_label, region_name, L) {
  n <- n_distinct(data$Accession)
  if (n < 2) return(NULL)
  
  # Count unique segregating sites (S) in this group/region
  S <- n_distinct(data$POS)
  
  # Watterson's Theta
  a1 <- sum(1/(1:(n-1)))
  theta <- S / a1
  
  # Pi (Simplified biallelic assumption for pangenome paths)
  # Pi = Average pairwise differences
  pi <- S * (2 * (1/n) * (1 - 1/n)) * (n / (n - 1))
  
  data.frame(
    Group = group_label,
    Region = region_name,
    n = n,
    S = S,
    pi_per_kb = (pi / L) * 1000,
    theta_per_kb = (theta / L) * 1000
  )
}

# 4. Main Processing Function
process_vcf <- function(vcf_path, segment_name) {
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  fix_df <- as.data.frame(getFIX(vcf)) %>%
    mutate(POS = as.numeric(POS)) %>%
    rename(Accession = CHROM) %>%
    inner_join(clade_map, by = "Accession")
  
  results_list <- list()
  
  for (r_name in names(regions)) {
    r_range <- regions[[r_name]]
    L <- r_range["end"] - r_range["start"] + 1
    
    # Filter data for this region
    region_data <- fix_df %>% filter(POS >= r_range["start"] & POS <= r_range["end"])
    
    # A. Calculate for the Whole Segment Group (All clades combined)
    results_list[[paste(r_name, "Whole")]] <- calc_stats(region_data, paste("Whole", segment_name), r_name, L)
    
    # B. Calculate for each Clade individually
    clade_results <- region_data %>%
      group_by(Clade) %>%
      group_modify(~ calc_stats(.x, as.character(.y$Clade), r_name, L))
    
    results_list[[paste(r_name, "Clades")]] <- clade_results
  }
  
  bind_rows(results_list)
}

# 5. Run analysis for G and S
cat("Analyzing Segment G...\n")
res_g <- process_vcf("results_GTPV/vcf/gfavariants.vcf", "G")

cat("Analyzing Segment S...\n")
res_s <- process_vcf("results_SPPV/vcf/gfavariants.vcf", "S")

# 6. Merge and Display
final_table <- bind_rows(res_g, res_s)

# Format for clear viewing
print(final_table %>% arrange(Group, factor(Region, levels=names(regions))))

library(dplyr)
library(ggplot2)

# 1. Data Input
regional_data <- data.frame(
  Group = c(rep("Whole G", 3), rep("Whole S", 3)),
  Region = rep(c("5_prime", "Core", "3_prime"), 2),
  n = c(rep(14, 3), rep(29, 3)),
  S = c(2648, 14135, 8068, 2067, 7056, 6096),
  # pi_kb is actually pi per bp in your formula, but let's re-calculate from S and L
  L = c(13751, 93060, 43091, 13751, 93060, 43091)
)

# 2. Corrected Tajima's D Function (Based on S and L)
calculate_tajimas_d_raw <- function(n, S, L) {
  if (S == 0) return(0)
  
  # a1, a2, b1, b2, c1, c2, e1, e2 constants
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / (a1^2)
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  
  # Tajima's D
  # theta_w = S / a1
  # pi = pi_total / L (This requires pi_total, which we don't have, 
  # but Tajima's D can be calculated just from S, a1, a2, and a constant)
  # Actually, the standard D = (pi - theta_w) / sqrt(var_d)
  # Let's re-use the pi_kb and theta_kb to infer the actual pi/theta, 
  # but ensure they are per-site, not per-kb.
  
  # IF pi_kb IS ACTUALLY MEAN PI PER BASE (pi / site):
  # pi_total = mean_pi_per_site * L
  
  # To avoid confusion, let's stick to the standard formula using S:
  # D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))
  # To do this, we MUST know average pairwise differences (pi_total)
  # Your pi_kb seems to be "nucleotide diversity (pi) per kb".
  
  # Let's redefine based on: pi_total = (pi_kb * L) / 1000
  # And theta_w = S / a1
  
  return(NA) # Placeholder for now
}

# Revised Function based on provided pi_kb being (diversity / 1000bp)
calculate_tajimas_d_fixed <- function(n, S, pi_kb, L) {
  # 1. Total Pairwise differences (Pi)
  # If pi_kb is pi/1000bp, total pi = (pi_kb/1000) * L
  total_pi <- (pi_kb / 1000) * L
  
  # 2. Watterson's Theta (based on S)
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  theta_w <- S / a1
  
  # 3. Variance components
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / (a1^2)
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  
  var_d <- (e1 * S) + (e2 * S * (S - 1))
  
  # 4. Final D
  d <- total_pi - theta_w
  return(d / sqrt(var_d))
}

# 3. Calculate D and compare across regions
# Need to add pi_kb back to the data
regional_data$pi_kb <- c(27.51, 21.70, 26.75, 10.37, 5.23, 9.76)

regional_analysis <- regional_data %>%
  rowwise() %>%
  mutate(
    Tajimas_D = calculate_tajimas_d_fixed(n, S, pi_kb, L)  )

print(regional_analysis %>% select(Group, Region, Tajimas_D))

