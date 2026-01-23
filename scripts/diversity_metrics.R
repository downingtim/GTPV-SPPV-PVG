library(tidyverse)
library(ggplot2)
library(ggrepel)

# Read the CSV file
df <- read.csv('diversity_metrics.csv')
df$Gene_short <- gsub("^LSDV", "", df$Gene)

df2 <- subset(df, GTPV_SNPs.Kb >2 & SPPV_SNPs.Kb >1)
cor_snps <- cor(df2$GTPV_SNPs.Kb, df2$SPPV_SNPs.Kb, use = "complete.obs")
cor_test_snps <- cor.test(df2$GTPV_SNPs.Kb, df2$SPPV_SNPs.Kb)
cat("Pearson r =", round(cor_snps, 4), "\n")
cat("p-value =", format(cor_test_snps$p.value, scientific = TRUE), "\n\n")

highlight_genes <- c("002","004","009", "013", "026", "132", "136","153","155")
df$highlight <- df$Gene_short %in% highlight_genes

# Log10 scale plot with non-transformed axis labels
p <- ggplot(df, aes(x = (GTPV_SNPs.Kb + 1.01 ), y = (SPPV_SNPs.Kb + 1.01 ))) +
   geom_point(aes(color = highlight), size =2.5, alpha = 0.4) +
   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
   guides(color = "none") + 
  geom_text_repel(aes(label = Gene_short), size =3, alpha = 0.8,
                max.overlaps = Inf,  # Show all labels even if overlapping
                box.padding = 0.3,   # Space around labels
                point.padding = 0.3) + # Space from points
  geom_smooth(
  method = "lm", se = TRUE, color = "darkblue", fill = "lightblue", alpha = 0.2)  +
  annotate("text", x =4, y =1.5, 
           label = paste("r =", round(cor_snps, 3)),
           size =7, fontface = "bold") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                labels = c( 1, 2, 5, 10, 20, 50, 100)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50),
                labels = c( 1, 2, 5, 10, 20, 50)) +
  labs(       x = "GTPV SNPs/Kb",       y = "SPPV SNPs/Kb") + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(size = 15),
        panel.grid.major = element_line(color = "lightgray", size = 0.4),
        panel.grid.minor = element_line(color = "lightgray", size = 0.2, linetype = "dotted"))
	ggsave('snps_kb_correlation.png', p, width =6, height =6, dpi = 300)

p2 <- ggplot(df, aes(x = (GTPV_SNPs.Kb  ), y = (SPPV_SNPs.Kb ))) +
  geom_point(size = 3, alpha = 0.4, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue", alpha = 0.2) +
  annotate("text", x =65, y =25, 
           label = paste("r =", round(cor_snps, 3)),
           size =7, fontface = "bold") +
  labs(       x = "GTPV SNPs/Kb",
       y = "SPPV SNPs/Kb") +  ylim(0,32) + xlim(0,85) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11),
        panel.grid.major = element_line(color = "lightgray", size = 0.3),
        panel.grid.minor = element_line(color = "lightgray", size = 0.15, linetype = "dotted"))
ggsave('snps_kb_correlation.nolog.png', p2, width =6, height =6, dpi = 300)
