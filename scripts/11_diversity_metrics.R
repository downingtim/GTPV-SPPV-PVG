library(tidyverse)
library(ggplot2)

# Read the CSV file
df <- read.csv('diversity_metrics.csv')

# Remove 'LSDV' prefix from gene names for plotting
df$Gene_short <- gsub("^LSDV", "", df$Gene)

# Calculate correlation
cor_snps <- cor(df$GTPV_SNPs.Kb, df$SPPV_SNPs.Kb, use = "complete.obs")
cor_test_snps <- cor.test(df$GTPV_SNPs.Kb, df$SPPV_SNPs.Kb)

# Print correlation results
cat("SNPs/Kb CORRELATION ANALYSIS (GTPV vs SPPV)\n")
cat("===========================================\n")
cat("Pearson r =", round(cor_snps, 4), "\n")
cat("p-value =", format(cor_test_snps$p.value, scientific = TRUE), "\n\n")

library(ggrepel)

# Log10 scale plot with non-transformed axis labels
p <- ggplot(df, aes(x = (GTPV_SNPs.Kb ), y = (SPPV_SNPs.Kb + 0.1))) +
  geom_point(size = 3, alpha = 0.4, color = "black") +
#  geom_text(aes(label = Gene_short), size = 2, hjust = -0.4, vjust = -0.4, alpha = 0.8) +
  geom_text_repel(aes(label = Gene_short), size = 2.5, alpha = 0.8,
                max.overlaps = Inf,  # Show all labels even if overlapping
                box.padding = 0.3,   # Space around labels
                point.padding = 0.3) + # Space from points
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue", alpha = 0.2) +
  annotate("text", x = 10^1.5, y = 10^2.1, 
           label = paste("r =", round(cor_snps, 3)),
           size = 5, fontface = "bold") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                labels = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50),
                labels = c(1, 2, 5, 10, 20, 50)) +
  labs(       x = "GTPV SNPs/Kb",
       y = "SPPV SNPs/Kb") + # ylim(0.1,50) + xlim(0.1,100) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11),
        panel.grid.major = element_line(color = "lightgray", size = 0.3),
        panel.grid.minor = element_line(color = "lightgray", size = 0.15, linetype = "dotted"))
ggsave('snps_kb_correlation.png', p, width =7, height =7, dpi = 300)
