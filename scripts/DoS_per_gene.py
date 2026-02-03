from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import Counter
import pandas as pd
import os

# -----------------------------
# CONFIG
# -----------------------------
GENBANK_FILE = "KX894508.gb.temp" # altered for LD019 and LD026
G_DIR = "results_GTPV/GENE_FASTA"
S_DIR = "results_SPPV/GENE_FASTA"
OUTFILE = "DoS_per_gene.csv"

STOP_CODONS = {"TAA", "TAG", "TGA"}
table = CodonTable.unambiguous_dna_by_name["Standard"]

# -----------------------------
# GENBANK PARSING
# -----------------------------
def load_gene_strands(genbank_file):
    gene_strand = {}

    record = SeqIO.read(genbank_file, "genbank")
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene = feature.qualifiers.get("gene", [None])[0]
        if gene is None:
            continue

        gene_id = gene.replace("LD", "LSDV")
        strand = "+" if feature.location.strand == 1 else "-"
        gene_strand[gene_id] = strand

    return gene_strand

# -----------------------------
# ORF HANDLING
# -----------------------------
def orient_sequence(seq, strand):
    return seq if strand == "+" else seq.reverse_complement()

def trim_to_orf(seq):
    seq = str(seq).upper()

    # find first ATG
    start = seq.find("ATG")
    if start == -1:
        return None

    seq = seq[start:]

    # truncate at first in-frame stop
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return Seq(seq[:i+3])

    # no stop codon found → allow full length if divisible by 3
    if len(seq) % 3 == 0:
        return Seq(seq)

    return None

# -----------------------------
# ALIGNMENT NORMALISATION
# -----------------------------
def normalise_alignment(aln, strand):
    norm_seqs = []

    for rec in aln:
        seq = orient_sequence(rec.seq, strand)
        orf = trim_to_orf(seq)
        if orf is None:
            continue
        norm_seqs.append(orf)

    if len(norm_seqs) < 2:
        return None

    min_len = min(len(s) for s in norm_seqs)
    min_len -= min_len % 3

    trimmed = [s[:min_len] for s in norm_seqs]
    return trimmed

def consensus_from_seqs(seqs):
    cons = []
    for i in range(len(seqs[0])):
        col = [s[i] for s in seqs]
        cons.append(Counter(col).most_common(1)[0][0])
    return Seq("".join(cons))

# -----------------------------
# CODON LOGIC
# -----------------------------
def codon_is_valid(c):
    return len(c) == 3 and "N" not in c and "-" not in c

def classify_change(c1, c2):
    if c1 == c2:
        return None
    if c1 not in table.forward_table or c2 not in table.forward_table:
        return None
    return "syn" if table.forward_table[c1] == table.forward_table[c2] else "nonsyn"

# -----------------------------
# COUNTS
# -----------------------------
def polymorphism_counts(seqs, consensus):
    pN = pS = 0
    for seq in seqs:
        for i in range(0, len(consensus), 3):
            c0 = str(consensus[i:i+3])
            c1 = str(seq[i:i+3])
            if not codon_is_valid(c0) or not codon_is_valid(c1):
                continue
            change = classify_change(c0, c1)
            if change == "syn":
                pS += 1
            elif change == "nonsyn":
                pN += 1
    return pN, pS

def divergence_counts(cons_G, cons_S):
    dN = dS = 0
    for i in range(0, len(cons_G), 3):
        cG = str(cons_G[i:i+3])
        cS = str(cons_S[i:i+3])
        if not codon_is_valid(cG) or not codon_is_valid(cS):
            continue
        change = classify_change(cG, cS)
        if change == "syn":
            dS += 1
        elif change == "nonsyn":
            dN += 1
    return dN, dS

def dos(dN, dS, pN, pS):
    if (dN + dS) == 0 or (pN + pS) == 0:
        return None
    return (dN / (dN + dS)) - (pN / (pN + pS))

# -----------------------------
# MAIN
# -----------------------------
gene_strands = load_gene_strands(GENBANK_FILE)
results = []

for gene in range(1, 157):
    gene_id = f"LSDV{gene:03d}"
    strand = gene_strands.get(gene_id)

    if strand is None:
        continue

    g_file = f"{G_DIR}/{gene_id}.fasta"
    s_file = f"{S_DIR}/{gene_id}.fasta"

    if not os.path.exists(g_file) or not os.path.exists(s_file):
        continue

    aln_G = AlignIO.read(g_file, "fasta")
    aln_S = AlignIO.read(s_file, "fasta")

    G_seqs = normalise_alignment(aln_G, strand)
    S_seqs = normalise_alignment(aln_S, strand)

    if G_seqs is None or S_seqs is None:
        continue

    cons_G = consensus_from_seqs(G_seqs)
    cons_S = consensus_from_seqs(S_seqs)

    pN_G, pS_G = polymorphism_counts(G_seqs, cons_G)
    pN_S, pS_S = polymorphism_counts(S_seqs, cons_S)
    dN, dS = divergence_counts(cons_G, cons_S)

    results.append({
        "gene": gene_id,
        "strand": strand,
        "dN": dN,
        "dS": dS,
        "pN_G": pN_G,
        "pS_G": pS_G,
        "DoS_G": dos(dN, dS, pN_G, pS_G),
        "pN_S": pN_S,
        "pS_S": pS_S,
        "DoS_S": dos(dN, dS, pN_S, pS_S)
    })

df = pd.DataFrame(results)
df.to_csv(OUTFILE, index=False)
