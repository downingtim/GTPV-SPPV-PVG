if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ape)
library(ggvenn)
library(ggplot2)

# 1. LOAD ALIGNMENT
dna <- read.dna("all.core2.aln", format = "fasta")
dna <- read.dna("all.genome.aln", format = "fasta")
dna_mat <- as.character(dna)

# 2. TARGET INDIVIDUAL SAMPLES
# Update these strings if your fasta headers are slightly different
samples <- list(
  LSDV = "LSDV_AF325528.1_NI_2490_Neethling_2490_Kenya_1958", 
  GTPV = "GTPV_NC_004003.1.fa#1#NC_004003.1",
  SPPV_Mod = "SPPV_MW167071.1.fa#1#MW167071.1",
  SPPV_Anc = "MYRZ180",
  SPPV_Pre = "CL299"
)

# Extract and name
rep_mat <- dna_mat[unlist(samples), ]
rownames(rep_mat) <- names(samples)

# 3. FIND POLYMORPHIC SITES
# Only sites where these 5 specific samples show variation
poly_sites <- which(apply(rep_mat, 2, function(x) {
  ux <- unique(x[x %in% c("a", "t", "c", "g")])
  length(ux) > 1
}))

# 4. MATCH LOGIC (Position + Allele)
get_sample_matches <- function(sites, s1_id, s2_id, s3_id, display_names) {
  s1_seq <- rep_mat[s1_id, sites]
  s2_seq <- rep_mat[s2_id, sites]
  s3_seq <- rep_mat[s3_id, sites]
  
  match_list <- list(
    V1 = paste0(sites, "_", s1_seq),
    V2 = paste0(sites, "_", s2_seq),
    V3 = paste0(sites, "_", s3_seq)
  )
  
  match_list <- lapply(match_list, function(x) x[!grepl("_[^atcg]$", x)])
  names(match_list) <- display_names
  return(match_list)
}

# 5. HIGH-VISIBILITY PLOTTING FUNCTION
save_readable_venn <- function(data, filename, title) {
  p <- ggvenn(data, 
              fill_alpha = 0.95,              # Solid colors
              fill_color = c("#D95F02", "#1B9E77", "#7570B3"), # High-contrast Palette
              stroke_color = "black",
              stroke_size = 1.2,
              set_name_size = 9,              # HUGE species names
              text_size =7) +                # HUGE numbers inside
       labs(title = title) +
       theme_void() +                         # Clean white background
       theme(
         plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b=20)),
         plot.background = element_rect(fill = "white", color = NA)
       )
  
  # Save with a tight bounding box
  ggsave(filename, plot = p, width =10, height =10, dpi = 300, bg = "white")
  message(paste("Saved highly readable file:", filename))
}

# 6. EXECUTE SAVING
save_readable_venn(
  get_sample_matches(poly_sites, "LSDV", "GTPV", "SPPV_Mod", c("LSDV", "GTPV", "Modern SPPV")),
  "Read_Modern.png", "Shared Alleles: Modern"
)

save_readable_venn(
  get_sample_matches(poly_sites, "LSDV", "GTPV", "SPPV_Pre", c("LSDV", "GTPV", "Pre-Modern SPPV")),
  "Read_Premodern.png", "Shared Alleles: Pre-Modern"
)

save_readable_venn(
  get_sample_matches(poly_sites, "LSDV", "GTPV", "SPPV_Anc", c("LSDV", "GTPV", "Ancient SPPV")),
  "Read_Ancient.png", "Shared Alleles: Ancient"
)

### part 2

library(ape)
library(genbankr)
library(dplyr)

# 1. LOAD ALIGNMENT AND DEFINE SAMPLES
#dna <- read.dna("all.core2.aln", format = "fasta")
dna <- read.dna("all.genome.aln", format = "fasta")
dna_mat <- as.character(dna)

samples <- list(
  LSDV = "LSDV_AF325528.1_NI_2490_Neethling_2490_Kenya_1958",
  GTPV = "GTPV_NC_004003.1.fa#1#NC_004003.1",
  SPPV_Mod = "SPPV_MW167071.1.fa#1#MW167071.1",
  SPPV_Anc = "MYRZ180",  SPPV_Pre = "CL299" )

# 2. EXTRACT SUBSET AND FIND POLYMORPHIC SITES
rep_mat <- dna_mat[unlist(samples), ]
rownames(rep_mat) <- names(samples)

# Find sites where at least one of these 5 samples differs
poly_sites <- which(apply(rep_mat, 2, function(x) {
  valid <- x[x %in% c("a", "t", "c", "g")]
  length(unique(valid)) > 1 }))

# 3. MAP ALIGNMENT POSITIONS TO AF325528.1 GENOME COORDINATES
# We count non-gap characters in the LSDV sequence to find the 'real' genomic position
lsdv_seq <- dna_mat[samples$LSDV, ]
gaps <- lsdv_seq == "-"
lsdv_coords <- cumsum(!gaps)
lsdv_coords[gaps] <- NA  # Mark gaps as NA

# 4. LOAD GENBANK ANNOTATIONS (AF325528.1.gb)
gb <- readGenBank("AF325528.1.gb")
genes <- as.data.frame(genes(gb))

# 5. BUILD THE TABLE
results <- data.frame(
  Alignment_Site = poly_sites,
  LSDV_Genome_Pos = lsdv_coords[poly_sites],
  LSDV_Allele = toupper(rep_mat["LSDV", poly_sites]),
  GTPV_Allele = toupper(rep_mat["GTPV", poly_sites]),
  Ancient_Allele = toupper(rep_mat["SPPV_Anc", poly_sites]),
  PreModern_Allele = toupper(rep_mat["SPPV_Pre", poly_sites]),
  Modern_Allele = toupper(rep_mat["SPPV_Mod", poly_sites])
)

# Function to find gene info for a given genome position
get_gene_info <- function(pos) {
  if (is.na(pos)) return(list(ID = "Intergenic", Name = "-"))
  
  match <- genes[genes$start <= pos & genes$end >= pos, ]
  if (nrow(match) > 0) {
    return(list(ID = match$gene_id[1], Name = ifelse(is.na(match$gene[1]), "-", match$gene[1])))
  } else {    return(list(ID = "Intergenic", Name = "-"))  }
}

# Apply mapping
gene_data <- t(sapply(results$LSDV_Genome_Pos, get_gene_info))
results$Gene_ID <- unlist(gene_data[,1])
results$Gene_Name <- unlist(gene_data[,2])

# 6. SAVE AS CSV
write.csv(results, "Capripox_Allele_Mapping.csv", row.names = FALSE)
message("Extraction complete. File saved: Capripox_Allele_Mapping.csv")


# 1. DEFINE THE CORE SETS (The Venn Centers)
# A site is only 'Center-Shared' if it's a match AND not missing/gap
get_center_ids <- function(df, sppv_col) {
  df$Alignment_Site[
    df$LSDV_Allele == df$GTPV_Allele & 
    df$LSDV_Allele == df[[sppv_col]] & 
    !df$LSDV_Allele %in% c("-", "N") & 
    !df[[sppv_col]] %in% c("-", "N")
  ]
}

# Get the IDs for the three circles you reported
ids_522 <- get_center_ids(results, "Ancient_Allele")
ids_158 <- get_center_ids(results, "PreModern_Allele")
ids_95  <- get_center_ids(results, "Modern_Allele")

# 2. EXTRACT THE TRANSITIONS (The 'Drop-outs')
# The 364 sites: Sites that were in the 522 set but are no longer in the 158 set
ids_loss_364 <- setdiff(ids_522, ids_158)

# The 63 sites: Sites that were in the 158 set but are no longer in the 95 set
ids_loss_63  <- setdiff(ids_158, ids_95)

# 3. CREATE THE DATA FRAMES
component_364 <- results[results$Alignment_Site %in% ids_loss_364, ]
component_63  <- results[results$Alignment_Site %in% ids_loss_63, ]

# 4. FINAL VERIFICATION
cat("\n--- Venn Intersection Counts ---\n")
cat("Ancient Core (Should be 522): ", length(ids_522), "\n")
cat("Pre-Modern Core (Should be 158): ", length(ids_158), "\n")
cat("Modern Core (Should be 95): ", length(ids_95), "\n")

cat("\n--- Component Counts ---\n")
cat("Medieval Loss (Should be 364): ", nrow(component_364), "\n")
cat("Modern Loss (Should be 63): ", nrow(component_63), "\n")

# 5. SAVE
write.csv(component_364, "Final_364_Medieval_Loss.csv", row.names = FALSE)
write.csv(component_63, "Final_63_Modern_Loss.csv", row.names = FALSE)

library(readxl)
library(Biostrings)

# ============================================================
# LOAD DATA
# ============================================================

mapping <- read_excel("Capripox_Allele_Mapping.xlsx")
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)

ref <- readDNAStringSet("AF325528.1.fasta")[[1]]
gb  <- readLines("AF325528.1.gb")

genome_len <- length(ref)

# ============================================================
# CLEAN POSITION COLUMN
# ============================================================

# Convert Site_LSDV to numeric, treating "-", "NA", "" as missing
mapping$Site_LSDV_num <- suppressWarnings(
  as.numeric(ifelse(mapping$Site_LSDV %in% c("-", "NA", ""), NA, mapping$Site_LSDV))
)

# ============================================================
# PARSE CDS FEATURES
# ============================================================

cds_idx <- grep("^     CDS", gb)

cds <- do.call(rbind, lapply(cds_idx, function(i){
  x <- gb[i]
  strand <- ifelse(grepl("complement", x), "-", "+")
  x <- gsub("complement\\(|join\\(|\\)|[<>]", "", x)
  x <- trimws(sub("^     CDS", "", x))
  parts <- strsplit(x, ",")[[1]]
  starts <- as.numeric(sub("\\.\\..*", "", parts))
  ends   <- as.numeric(sub(".*\\.\\.", "", parts))
  
  gene <- NA
  j <- i + 1
  while(j <= length(gb) && !grepl("^     [A-Za-z]", gb[j])){
    if(grepl("/gene=", gb[j])){
      gene <- sub('.*gene="([^"]+)".*', '\\1', gb[j])
    }
    j <- j + 1
  }
  
  data.frame(start = min(starts), end = max(ends), strand = strand, gene = gene)
}))

cds <- cds[!is.na(cds$start), ]

# ============================================================
# HELPERS
# ============================================================

safe_base <- function(x){
  if(is.na(x) || length(x)==0) return(NA_character_)
  x <- toupper(x)
  if(x %in% c("A","C","G","T")) return(x)
  return(NA_character_)
}


rc <- function(x) as.character(complement(DNAString(x)))

translate_codon <- function(x){
  if(is.na(x) || nchar(x)!=3) return(NA_character_)
  if(grepl("[^ACGT]", x)) return(NA_character_)
  as.character(Biostrings::translate(DNAString(x)))
}

mut_codon <- function(ref, alt, pos){
  if(is.na(alt) || !(alt %in% c("A","C","G","T"))) {
    return(ref)   # keep reference codon if allele is invalid
  }
  v <- strsplit(ref, "")[[1]]
  v[pos] <- alt
  paste0(v, collapse = "")
}

# ============================================================
# OUTPUT COLUMNS
# ============================================================

cols <- c("Mutation_Status","Gene","Gene_Nuc_Position","Gene_AA_Position",
          "LSDV_AA","GTPV_AA","Anc_AA","Premodern_AA","Modern_AA",
          "LSDV_codon","GTPV_codon","Anc_codon","Premodern_codon","Modern_codon")
mapping[cols] <- NA

# ============================================================
# MAIN LOOP
# ============================================================

for(i in seq_len(nrow(mapping))){
  
  pos <- mapping$Site_LSDV_num[i]
  
  if(is.na(pos) || pos < 1 || pos > genome_len){
    mapping$Mutation_Status[i] <- "Intergenic"
    next
  }
  
  hit_idx <- which(pos >= cds$start & pos <= cds$end)
  if(length(hit_idx) == 0){
    mapping$Mutation_Status[i] <- "Intergenic"
    next
  }
  
  hit <- cds[hit_idx[1], ]
  mapping$Gene[i] <- hit$gene
  
  # extract CDS sequence
  cds_seq <- subseq(ref, hit$start, hit$end)
  if(hit$strand == "-") cds_seq <- reverseComplement(cds_seq)
  cds_seq <- as.character(cds_seq)
  
  # position within CDS sequence
  if(hit$strand == "+"){
    gpos <- pos - hit$start + 1
  } else {
    gpos <- hit$end - pos + 1
  }
  
  aa_pos <- ((gpos - 1) %/% 3) + 1
  codon_pos <- ((gpos - 1) %% 3) + 1
  
  mapping$Gene_Nuc_Position[i] <- gpos
  mapping$Gene_AA_Position[i]  <- aa_pos
  
  start <- ((aa_pos - 1) * 3) + 1
  end   <- start + 2
  if(end > nchar(cds_seq)) next
  
  ref_codon <- substr(cds_seq, start, end)
  if(nchar(ref_codon)!=3) next
  
  A <- c(
    LSDV = safe_base(mapping$LSDV_Allele[i]),
    GTPV = safe_base(mapping$GTPV_Allele[i]),
    ANC  = safe_base(mapping$Ancient[i]),
    PRE  = safe_base(mapping$PreModern[i]),
    MOD  = safe_base(mapping$Modern[i])
  )
  
  if(hit$strand == "-"){
    A <- sapply(A, function(x) if(!is.na(x)) rc(x) else NA)
  }
  
  codons <- c(
    LSDV = mut_codon(ref_codon, A["LSDV"], codon_pos),
    GTPV = mut_codon(ref_codon, A["GTPV"], codon_pos),
    ANC  = mut_codon(ref_codon, A["ANC"], codon_pos),
    PRE  = mut_codon(ref_codon, A["PRE"], codon_pos),
    MOD  = mut_codon(ref_codon, A["MOD"], codon_pos)
  )
  
  aas <- sapply(codons, translate_codon)
  
  mapping[i, c("LSDV_codon","GTPV_codon","Anc_codon","Premodern_codon","Modern_codon")] <- codons
  mapping[i, c("LSDV_AA","GTPV_AA","Anc_AA","Premodern_AA","Modern_AA")] <- aas
  
  # classification
  refAA <- translate_codon(ref_codon)
  altAAs <- unique(na.omit(aas))
  
  if(length(altAAs)==0 || all(altAAs==refAA)){
    mapping$Mutation_Status[i] <- "Synonymous"
  } else {
    mapping$Mutation_Status[i] <- "Nonsynonymous"
  }
}

# ============================================================
# OUTPUT
# ============================================================

write.table(mapping,"mapping_annotated.tsv",sep="\t",quote=FALSE,row.names=FALSE)
table(mapping$Mutation_Status, useNA="ifany")
View(mapping)
mapping[10:12, c("Site_LSDV","Gene","Gene_Nuc_Position","Gene_AA_Position",
                 "LSDV_codon","LSDV_AA","GTPV_codon","GTPV_AA")]

str(mapping) 

# ancient subset
intersect_anc <- subset(mapping, Site_LSDV != "NA" & 
                    LSDV_Allele==GTPV_Allele & LSDV_Allele==Ancient )
str(intersect_anc) # 1041 alleles shared initially
View(intersect_anc)

# ancient to premodern changes
intersect_anc_premodern <- subset(intersect_anc, Ancient !=
           PreModern & Ancient != Modern & Modern==PreModern & 
             Ancient != "-" & PreModern != "-" & Modern != "-" )
str(intersect_anc_premodern) # 657
View(intersect_anc_premodern)

# premodern to modern changes
intersect_anc_pre <- subset(mapping, Site_LSDV != "NA" &
    LSDV_Allele==GTPV_Allele & LSDV_Allele==Ancient & LSDV_Allele==PreModern
    &   Ancient != "-" & PreModern != "-" & Modern != "-" &
       LSDV_Allele != "-" & GTPV_Allele != "-"  )
str(intersect_anc_pre) # 258 alleles shared initially
View(intersect_anc_pre)

# ============================================================
# Function: nonsyn_syn_ratio
# ============================================================
# Input: a data.frame like your mapping, with columns "Gene" and "Mutation_Status"
# Output: a data.frame with Gene, nonsyn_count, syn_count, ratio
# Only includes genes with nonsyn > 0 and syn > 0

nonsyn_syn_ratio <- function(df){
    counts <- table(df$Gene, df$Mutation_Status)
  counts <- as.data.frame.matrix(counts)
  # ensure columns exist
  if(!"Nonsynonymous" %in% colnames(counts)) counts$Nonsynonymous <- 0
  if(!"Synonymous" %in% colnames(counts)) counts$Synonymous <- 0
  counts$Gene <- rownames(counts)
  # filter genes with both nonsyn and syn > 0
#  counts <- counts[counts$Nonsynonymous > 0 & counts$Synonymous > 0, ]
  counts$Ratio <- round((counts$Nonsynonymous + 0.00001) 
             /(counts$Synonymous + 0.00001), 4)
  counts <- counts[, c("Gene","Nonsynonymous","Synonymous","Ratio")]
  rownames(counts) <- NULL
  return(counts)}

mapping_ratios <- nonsyn_syn_ratio(mapping) # all
intersect_anc_ratios <- nonsyn_syn_ratio(intersect_anc)
str(intersect_anc) # 1041 mutations
str(intersect_anc_ratios) # 145 genes
View(intersect_anc_ratios)
# ancient is like LSDV/GTPV but not pre & modern

intersect_anc_premodern_ratios <- nonsyn_syn_ratio(intersect_anc_premodern)
str(intersect_anc_premodern) # 657 mutations
str(subset(intersect_anc_premodern, Gene_ID=="Intergenic")) # 30
str(subset(intersect_anc_premodern, Mutation_Status=="Nonsynonymous")) # 316
str(subset(intersect_anc_premodern, Mutation_Status=="Synonymous")) # 311
str(intersect_anc_premodern_ratios) # 104 genes
View(intersect_anc_premodern_ratios)
# ancient to premodern changes that are retained in modern

intersect_anc_pre_ratios <- nonsyn_syn_ratio(intersect_anc_pre)
# premodern changes that are retained in modern
View(intersect_anc_pre) # 258 mutations
str(subset(intersect_anc_pre, Mutation_Status=="Nonsynonymous")) # 118
str(subset(intersect_anc_pre, Mutation_Status=="Synonymous")) # 127
str(subset(intersect_anc_pre, Gene_ID=="Intergenic")) # 13
str(intersect_anc_pre_ratios) # 102 genes
View(intersect_anc_pre_ratios)

intersect( intersect_anc_pre$Site_LSDV, intersect_anc_premodern$Site_LSDV)


# get weird matches

a1 <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele==GTPV_Allele &
                LSDV_Allele!=Ancient & LSDV_Allele==PreModern
             & LSDV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1) # 66 where premodern became like LSDV/GTPV & preserved in modern
View(a1)

a1m <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele==GTPV_Allele &
               LSDV_Allele!=Ancient & LSDV_Allele!=PreModern
             & LSDV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1m) # 6 where modern became like LSDV/GTPV 
View(a1m)

a1m <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele==GTPV_Allele &
                LSDV_Allele!=Ancient & LSDV_Allele!=PreModern
              & LSDV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
              &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1m) # 6 where modern became like LSDV/GTPV 
View(a1m)

b1 <- subset(mapping, Site_LSDV != "NA" & Ancient != PreModern 
             & Ancient!=Modern & PreModern==Modern & Ancient != "-"
             & PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(b1) # 851 where premodern became different & preserved in modern
str(subset(b1, Gene_ID=="Intergenic")$Site_LSDV) # 46 inter
str(subset(b1, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 382 nonsyn
str(subset(b1, Mutation_Status=="Synonymous")$Site_LSDV) # 423 syn

b1m <- subset(mapping, Site_LSDV != "NA" & Ancient == PreModern 
             & Ancient!=Modern & PreModern!=Modern & Ancient != "-"
             & PreModern != "-"  
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(b1m) # 851 where modern became different 
str(subset(b1m, Gene_ID=="Intergenic")$Site_LSDV) # 17 inter
str(subset(b1m, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 142 nonsyn
str(subset(b1m, Mutation_Status=="Synonymous")$Site_LSDV) # 157 syn

a1_lsdv <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele!=GTPV_Allele &
               LSDV_Allele!=Ancient & LSDV_Allele==PreModern
             & LSDV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1_lsdv) # 26 where premodern became like LSDV only & preserved in modern
View(a1_lsdv)
str(subset(a1_lsdv, Gene_ID=="Intergenic")$Site_LSDV) # 1 inter
str(subset(a1_lsdv, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 8 nonsyn
str(subset(a1_lsdv, Mutation_Status=="Synonymous")$Site_LSDV) # 17 syn

a1m_lsdv <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele!=GTPV_Allele &
                    LSDV_Allele!=Ancient & LSDV_Allele!=PreModern
                  & LSDV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
                  &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1m_lsdv) # 7 where modern became like LSDV only 
View(a1m_lsdv)
str(subset(a1m_lsdv, Gene_ID=="Intergenic")$Site_LSDV) # 1 inter
str(subset(a1m_lsdv, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 1 nonsyn
str(subset(a1m_lsdv, Mutation_Status=="Synonymous")$Site_LSDV) # 5 syn


a1_gtpv <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele!=GTPV_Allele &
                    GTPV_Allele!=Ancient & GTPV_Allele==PreModern
                  & GTPV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
                  &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-") 
str(a1_gtpv) # 50 where premodern became like GTPV only & preserved in modern
View(subset(a1_gtpv, Mutation_Status=="Nonsynonymous"))
str(subset(a1_gtpv, Gene_ID=="Intergenic")$Site_LSDV) # 5 inter
str(subset(a1_gtpv, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 18 nonsyn
str(subset(a1_gtpv, Mutation_Status=="Synonymous")$Site_LSDV) # 27 syn

a1m_gtpv <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele!=GTPV_Allele &
                    GTPV_Allele!=Ancient & GTPV_Allele!=PreModern
                  & GTPV_Allele==Modern &   Ancient != "-" &   PreModern != "-"
                  &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
str(a1m_gtpv) # 25 where premodern became like GTPV only & preserved in modern
View(a1m_gtpv)
str(subset(a1m_gtpv, Gene_ID=="Intergenic")$Site_LSDV) # 1 inter
str(subset(a1m_gtpv, Mutation_Status=="Nonsynonymous")$Site_LSDV) # 12 nonsyn
str(subset(a1m_gtpv, Mutation_Status=="Synonymous")$Site_LSDV) # 12 syn

# Nonsynonymous changes where both pre-modern and modern SPPV differed
# from ancient SPPV were altered more often to be similar to GTPV than
# LSDV (18 vs 8), as were synonymous ones (27 vs 17). This compares to
# 382 nonsyn and 423 sites were both pre-modern and modern SPPV differed
# from ancient SPPV were altered.
# 45 out of 805 (5.6%) vs 25 out of 805 (3.1%)

# This was also true for
# nonsyn mutations in modern SPPV differing from ancient and pre-modern SPPV
# (12 vs 1) and syn ones (12 vs 5). This compares to
# 142 nonsyn and 157 sites were both pre-modern and modern SPPV differed
# from ancient SPPV were altered.
# 24 out of 299 (8.0%) vs 6 out of 299 (2.0%)

mat <- matrix(c(18, 27, 382, 423), nrow = 2, byrow = T)
fisher.test(mat) # NS
mat <- matrix(c(8, 17, 382, 423), nrow = 2, byrow = T)
fisher.test(mat) # NS
mat <- matrix(c(8, 17, 18, 27), nrow = 2, byrow = T)
fisher.test(mat) # NS


mat <- matrix(c(12, 12, 142, 157), nrow = 2, byrow = T)
fisher.test(mat) # NS
mat <- matrix(c(1, 5, 142, 157), nrow = 2, byrow = T)
fisher.test(mat) # NS
mat <- matrix(c(12, 12, 1, 5), nrow = 2, byrow = T)
fisher.test(mat) # NS

df <- data.frame(
  period = c("Premodern","Premodern","Modern","Modern"),
  mutation_type = c("Nonsyn","Syn","Nonsyn","Syn"),
  GTPV_like = c(18,27,12,12),
  LSDV_like = c(8,17,1,5),
  total = c(382,423,142,157) )

expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i){
  data.frame(  period = df$period[i],
    mutation_type = df$mutation_type[i],
    likeness = c(rep(1, df$GTPV_like[i]), rep(0, df$LSDV_like[i]))) }))

fit <- glm(likeness ~ mutation_type * period,
           data = expanded,
           family = binomial)
summary(fit)


# visualise
# 
library(ggplot2)
df <- data.frame(
  Type = rep(c("Nonsyn","Syn"), each=2),
  Target = c(rep(c("Pre-modern LSDV-like","Pre-modern GTPV-like"), 2),
             rep(c("Modern LSDV-like","Modern GTPV-like"), 2) ),
  Count = c(8,18,17,27 , 1,12,5,12)  )

ggplot(df, aes(x=Type, y=Count, fill=Target)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Mutation type", y="Count") +
  scale_fill_manual(values=c("steelblue","tomato",
                             "blue", "red")) +
  theme_minimal()

######## oh wait we want amino acid not nucleotide  woops ###
# get weird matches - amino acid level !!

a1 <- subset(mapping, Site_LSDV != "NA" & LSDV_Allele==GTPV_AA &
               LSDV_AA!=Anc_AA & LSDV_AA==Premodern_AA
             & LSDV_AA==Modern &   Anc_AA != "NA" &   Premodern_AA != "NA"
             &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
             &   Ancient != "-" &   PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(a1) # 1 where premodern became like LSDV/GTPV & preserved in modern
View(a1)

a1m <- subset(mapping, Site_LSDV != "NA" & LSDV_AA==GTPV_AA &
                LSDV_AA!=Anc_AA & LSDV_AA!=Premodern_AA
              & LSDV_AA==Modern_AA &   Anc_AA != "NA" &   Premodern_AA != "NA"
              &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
              &   Ancient != "-" &   PreModern != "-"
              &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(a1m) # 3 where modern became like LSDV/GTPV 
View(a1m)

b1 <- subset(mapping, Site_LSDV != "NA" & Anc_AA != Premodern_AA
             & Anc_AA!=Modern_AA & Premodern_AA==Modern_AA
             &   Anc_AA != "NA" &   Premodern_AA != "NA"
             &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
             &   Ancient != "-" &   PreModern != "-"
             &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(b1) # 380 where premodern became different & preserved in modern 
dim(subset(b1, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
           & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
           & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155"))
#
b1m <- subset(mapping, Site_LSDV != "NA" & Anc_AA == Premodern_AA 
              & Anc_AA!=Modern_AA & Premodern_AA!=Modern_AA
              &   Anc_AA != "NA" &   Premodern_AA != "NA"
              &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
              &   Ancient != "-" &   PreModern != "-"
              &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
dim(b1m) # 139 where modern became different  
dim(subset(b1m, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
        & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
        & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155"))
# 110 exlcuding psuedogenes

a1_lsdv <- subset(mapping, Site_LSDV != "NA" & LSDV_AA!=GTPV_AA &
                    LSDV_AA!=Anc_AA & LSDV_AA==Premodern_AA
                  & LSDV_AA==Modern_AA &   Anc_AA != "NA" & Premodern_AA != "NA"
                  &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
                  &   Ancient != "-" &   PreModern != "-"
                  &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
a1_lsdv <- subset(a1_lsdv, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
           & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
           & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155")
dim(a1_lsdv) # 7 where premodern became like LSDV only & preserved in modern

a1m_lsdv <- subset(mapping, Site_LSDV != "NA" & LSDV_AA!=GTPV_AA &
                     LSDV_AA!=Anc_AA & LSDV_AA!=Premodern_AA
                   & LSDV_AA==Modern_AA &   Anc_AA != "NA" & Premodern_AA != "NA"
                   &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
                   &   Ancient != "-" &   PreModern != "-"
                   &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
View(a1m_lsdv) 
a1m_lsdv <- subset(a1m_lsdv, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
           & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
           & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155")
dim(a1m_lsdv) # 1 where modern became like LSDV only 

a1_gtpv <- subset(mapping, Site_LSDV != "NA" & LSDV_AA!=GTPV_AA &
                    GTPV_AA!=Anc_AA & GTPV_AA==Premodern_AA
                  & GTPV_AA==Modern_AA &   Anc_AA != "NA" & Premodern_AA != "NA"
                  &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
                  &   Ancient != "-" &   PreModern != "-"
                  &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
a1_gtpv <- subset(a1_gtpv, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
           & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
           & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155")
dim(a1_gtpv) # 16 where premodern became like GTPV only & preserved in modern

a1m_gtpv <- subset(mapping, Site_LSDV != "NA" & LSDV_AA!=GTPV_AA &
                     GTPV_AA!=Anc_AA & GTPV_AA!=Premodern_AA 
                   & GTPV_AA==Modern_AA &   Anc_AA != "NA" & Premodern_AA != "NA"
                   &   LSDV_AA != "NA"  & Modern_AA != "NA" &  GTPV_AA != "NA"
                   &   Ancient != "-" &   PreModern != "-"
                   &   LSDV_Allele != "-"  & Modern != "-" &  GTPV_Allele != "-")
a1m_gtpv <- subset(a1m_gtpv, Gene_ID != "LSDV002" & Gene_ID!="LSDV009" & Gene_ID != "LSDV013"
           & Gene_ID != "LSDV026" & Gene_ID != "LSDV132" & Gene_ID != "LSDV136"
           & Gene_ID != "LSDV004" & Gene_ID != "LSDV153" & Gene_ID != "LSDV155")
dim(a1m_gtpv) # 9 where modern became like GTPV only

intersect(a1_gtpv$Gene_ID, a1_lsdv$Gene_ID) # none
intersect(a1_gtpv$Gene_ID, a1m_gtpv$Gene_ID)
# 2 changes each at premodern LSDV032" "LSDV042" with 1 additional at modern
intersect(a1_gtpv$Gene_ID, a1m_lsdv$Gene_ID) # none
intersect(a1_lsdv$Gene_ID, a1m_lsdv$Gene_ID) # none
intersect(a1_lsdv$Gene_ID, a1m_gtpv$Gene_ID) # 1 nonsyn each "LSDV059"

a1_lsdv$Context="LSDV_ancient_to_premodern"
a1m_lsdv$Context="LSDV_premodern_to_modern"
a1_gtpv$Context="GTPV_ancient_to_premodern"
a1m_gtpv$Context="GTPV_premodern_to_modern"
table1 <- rbind(a1_lsdv, a1m_lsdv, a1_gtpv, a1m_gtpv)
table1$Site_LSDV <-as.numeric(table1$Site_LSDV) 
View(table1)
write.csv(table1, file="table8.csv")
View(mapping)
