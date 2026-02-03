library(seqinr)
library(dplyr)
#########################################################
## INPUTS
#########################################################

aln_file  <- "GTPV_KX894508.aln"
gene_file <- "../table.tsv"

ref_name <- "KX894508"
#ref_name <- "NC_004003.1"
trim_ends <- 0

dir.create("GENE_FASTA", showWarnings = FALSE)
dir.create("GENE_TREE", showWarnings = FALSE)

#########################################################
## READ GENE TABLE
#########################################################

genes <- read.table(gene_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- genes[, c("LSDV_gene_ID", "start.2", "end.2")]
colnames(genes) <- c("gene","start","end")

#########################################################
## LOAD ALIGNMENT
#########################################################

alignment  <- read.alignment(aln_file, format="fasta")
seq_matrix <- as.matrix(alignment)

#########################################################
## FIND LSDV REFERENCE
#########################################################

print(ref_name)
str(seq_matrix)
ref_row <- which(rownames(seq_matrix)==ref_name)
if(length(ref_row)!=1) stop("Reference not found")
ref_seq <- toupper(seq_matrix[ref_row,])

#########################################################
## ALIGNMENT → GENOME MAP
#########################################################

is_base <- ref_seq %in% c("A","C","G","T","N")
ref_pos <- cumsum(is_base)
ref_pos[!is_base] <- NA
genome_len <- max(ref_pos, na.rm=TRUE)

#########################################################
## FILTER EDGE GENES
#########################################################

genes <- genes %>%
  filter(start > trim_ends,
         end   < genome_len - trim_ends,
         !is.na(gene))

#########################################################
## RUN RAxML PER GENE
#########################################################

#for(i in 1:nrow(genes)){
for(i in 6: nrow(genes)){
#for(i in 1:5){ 

  g <- genes[i,]

  message("Processing ", g$gene, " ", g$start, "-", g$end)

  #######################################################
  ## MAP GENE TO ALIGNMENT
  #######################################################

  cols <- which(ref_pos >= g$start & ref_pos <= g$end)
  if(length(cols) < 100) next

  window_mat <- seq_matrix[, cols, drop=FALSE]

  #######################################################
  ## WRITE FASTA
  #######################################################

  seqs <- apply(window_mat, 1, paste, collapse="")
  names(seqs) <- rownames(window_mat)
  str(seqs[2:15])
  fasta <- paste0("GENE_FASTA/", g$gene, ".fasta")
  write.fasta(as.list(seqs[2:15]), names(seqs[2:15]), fasta)

  #######################################################
  ## RUN RAxML
  #######################################################

  prefix <- paste0("GENE_TREE/", g$gene)

  cmd <- paste(
    "raxml-ng --all",
    "--msa", fasta,
    "--model GTR+G",
    "--seed 1241",
    "--prefix", prefix,
    "--redo"  )
  system(cmd)
}
