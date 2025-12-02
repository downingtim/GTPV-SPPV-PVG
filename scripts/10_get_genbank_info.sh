#!/bin/bash

# Input FASTA
FASTA="GTPV.fasta"

# Temp folder
mkdir -p gb_tmp

echo -e "Accession\tLength\tCDS_count"

# Loop through each FASTA header
grep ">" "$FASTA" | while read header; do
    
    # Example header:
    # >NC_004003.1.fa#1#NC_004003.1
    
    # Extract the *file-like* name before the first '#'
    full=$(echo "$header" | sed 's/>//; s/#.*//')
    # Remove .fa if present → NC_004003.1
    acc=$(echo "$full" | sed 's/\.fa$//')

    # Extract genome length from FASTA
    length=$(awk -v id="$header" '
        $0 ~ id {p=1; next}
        p && /^>/ {p=0}
        p {gsub(/[^ACGTacgt]/,""); L+=length($0)}
        END {print L}
    ' "$FASTA")

    # Download GenBank file if not yet downloaded
    gb="gb_tmp/${acc}.gb"
    if [ ! -f "$gb" ]; then
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=gb&retmode=text" \
            -o "$gb"
    fi

    # Count CDS entries
    cds_count=$(grep -c "^     CDS" "$gb")

    # Print result
    echo -e "${acc}\t${length}\t${cds_count}"

done
