#!/bin/bash

# Script to extract specific regions from multiple genomes in a FASTA file
# Usage: ./extract_lsdv_regions.sh

# Configuration
FASTA_FILE="pangrowth/communities.genomes.fasta"
OUTPUT_DIR="lsdv_regions"

# Define regions to extract (region_name:start-end)
declare -A REGIONS=(
    ["LSDV009"]="5455-6157"
    ["LSDV013"]="8558-9648"
    ["LSDV026"]="17599-18470"
    ["LSDV132"]="119288-119899"
    ["LSDV136"]="128916-129451"
)

# Check if required tools are available
command -v samtools >/dev/null 2>&1 || {
    echo "Error: samtools is required but not installed. Please install samtools."
    echo "On Ubuntu/Debian: sudo apt-get install samtools"
    echo "On macOS: brew install samtools"
    exit 1
}

# Check if input file exists
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: Input file $FASTA_FILE not found!"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Extracting LSDV regions from $FASTA_FILE..."
echo "Found $(grep -c '^>' "$FASTA_FILE") sequences in input file"

# Index the FASTA file if index doesn't exist
if [[ ! -f "$FASTA_FILE.fai" ]]; then
    echo "Creating FASTA index..."
    samtools faidx "$FASTA_FILE"
fi

# Get list of all sequence headers in the FASTA file
GENOMES=($(grep '^>' "$FASTA_FILE" | sed 's/^>//' | cut -d' ' -f1))

echo "Found ${#GENOMES[@]} genomes in the file"

# Extract each region for all genomes
for region_name in "${!REGIONS[@]}"; do
    coordinates="${REGIONS[$region_name]}"
    start=$(echo "$coordinates" | cut -d'-' -f1)
    end=$(echo "$coordinates" | cut -d'-' -f2)
    
    output_file="$OUTPUT_DIR/${region_name}.fa"
    
    echo "Extracting region $region_name ($coordinates) for all genomes..."
    
    # Clear output file
    > "$output_file"
    
    # Extract region from each genome
    for genome in "${GENOMES[@]}"; do
        # Try to extract the region and clean up the header
        if extracted_seq=$(samtools faidx "$FASTA_FILE" "${genome}:${start}-${end}" 2>/dev/null); then
            # Replace the modified header with the original genome name
            echo "$extracted_seq" | sed "1s/^.*$/>${genome}/" >> "$output_file"
            echo "  ✓ Extracted from $genome"
        else
            echo "  ✗ Failed to extract from $genome (sequence too short or not found)"
        fi
    done
    
    # Check if any sequences were extracted
    seq_count=$(grep -c '^>' "$output_file" 2>/dev/null || echo "0")
    echo "  → Saved $seq_count sequences to $output_file"
    echo
done

echo "Extraction complete! Check the $OUTPUT_DIR directory for output files."

# Summary
echo "Summary:"
for region_name in "${!REGIONS[@]}"; do
    output_file="$OUTPUT_DIR/${region_name}.fa"
    if [[ -f "$output_file" ]]; then
        seq_count=$(grep -c '^>' "$output_file" 2>/dev/null || echo "0")
        echo "  $region_name: $seq_count sequences"
    fi
done
