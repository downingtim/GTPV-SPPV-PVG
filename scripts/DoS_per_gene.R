# ================================
# DoS and dN/dS analysis + plots
# ================================

library(tidyverse)
library(ggrepel)

# ----------------
# Load data
# ----------------
df <- read_csv("DoS_per_gene.csv", show_col_types = FALSE)
df <- df %>% mutate(
    gene_num = as.integer(str_extract(gene, "\\d+")),
    dNdS = ifelse(dS > 0, dN / dS, NA_real_)  )

p1 <- ggplot(df, aes(x = gene_num)) +
  geom_point(    aes(y = DoS_G),    colour = "blue",
    size = 2,    na.rm = TRUE ) +
  geom_point(    aes(y = DoS_S),    colour = "red",
    size = 2,   na.rm = TRUE  ) +
  geom_hline(    yintercept = 0,    linetype = "dashed",
    colour = "grey50"  ) +
  labs(    x = "Genes",    y = "DoS rate",
    subtitle = "Blue: GTPV | Red: SPPV") +
  theme_bw() +  theme(panel.grid.minor = element_blank())
pdf("DoS_plot1.pdf",width=9,height=5)
print(p1)
dev.off()

p3 <- ggplot(df, aes(x = gene_num)) +
  geom_point(    aes(y = log10(dNdS)),
    colour = "black",    size = 2,    na.rm = TRUE ) +
  geom_hline(    yintercept = 0,    linetype = "dashed",
    colour = "grey50"  ) +
  labs(    x = "Genes",    y = "dN/dS") + 
  theme_bw() +  theme(panel.grid.minor = element_blank())

pdf("DoS_plot_dnds.pdf",width=9,height=5)
print(p3)
dev.off()

# ----------------
# Plot 2:
# DoS_G vs DoS_S
# ----------------
rho <- cor(  df$DoS_G,  df$DoS_S,
  method = "spearman",  use = "pairwise.complete.obs")
p2 <- ggplot(df, aes(x = DoS_G, y = DoS_S, label = gsub("LSDV","",gene))) +
  geom_hline(    yintercept = 0,   linetype = "dashed",    colour = "grey70"  ) +
  geom_vline(    xintercept = 0,    linetype = "dashed",    colour = "grey70"  ) +
  geom_point(    size =3,    colour = "black",    na.rm = TRUE, alpha=0.6) +
  geom_text_repel(    size =2,    max.overlaps = Inf,    na.rm = TRUE  ) +
  coord_equal() + labs(    x = "DoS GTPV",    y = "DoS SPPV",  ) +
   geom_smooth(    method = "lm", alpha=0.6,
    se = FALSE,    colour = "blue",    linewidth = 1.5,    na.rm = TRUE )  +
  theme_bw() + ylim(-1.001,1.001) + xlim(-1.001,1.001) 
  p2 <- p2 +  annotate(    "text",    x =-0.9,    y = 0.9,    hjust = 0,
    label = paste0("rho=", round(rho,2)),    size = 4  )
pdf("DoS_plot2.pdf",width=6,height=6)
print(p2)
dev.off()
# ----------------
# Numeric summaries
# ----------------
cat("\n=== MEAN VALUES ===\n")
medians <- df %>%
  summarise(
    dN_mean    = median(dN, na.rm = TRUE),
    dS_mean    = median(dS, na.rm = TRUE),
    dNdS_mean    = median(dN/dS, na.rm = TRUE),
    pN_G_mean  = median(pN_G, na.rm = TRUE),
    pS_G_mean  = median(pS_G, na.rm = TRUE),
    DoS_G_mean = median(DoS_G, na.rm = TRUE),
    pN_S_mean  = median(pN_S, na.rm = TRUE),
    pS_S_mean  = median(pS_S, na.rm = TRUE),
    DoS_S_mean = median(DoS_S, na.rm = TRUE)
  )
sds <- df %>%
  summarise(
    dN_sd    = sd(dN, na.rm = TRUE),
    dS_sd    = sd(dS, na.rm = TRUE),
    dNdS_sd  = sd(dN/dS, na.rm = TRUE),
    pN_G_sd  = sd(pN_G, na.rm = TRUE),
    pS_G_sd  = sd(pS_G, na.rm = TRUE),
    DoS_G_sd = sd(DoS_G, na.rm = TRUE),
    pN_S_sd  = sd(pN_S, na.rm = TRUE),
    pS_S_sd  = sd(pS_S, na.rm = TRUE),
    DoS_S_sd = sd(DoS_S, na.rm = TRUE)
  )

print(data.frame(medians))
print(data.frame(sds))
cat("\n=== DISTRIBUTION SUMMARIES ===\n")
print(  summary(    df %>%
      select(dN, dS,  pN_G, pS_G, DoS_G, pN_S, pS_S, DoS_S)  ))

cat("\n=== SPEARMAN CORRELATION: DoS_G vs DoS_S ===\n")
print(  cor.test(    df$DoS_G,    df$DoS_S,
    method = "spearman",    use = "pairwise.complete.obs"  ) )
