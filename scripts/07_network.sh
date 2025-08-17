#!/bin/bash

# Script to process FASTA files and generate community detection PDFs
# Usage: ./process_communities.sh [cpus] [seed_value]
# Default cpus: 4, default seed: 1500

# Set parameters
CPUS=${1:-4}
SEED=${2:-1500}

# Ensure we have the required tools and Python
command -v python >/dev/null 2>&1 || { echo "Error: python not found. Please install Python." >&2; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "Error: bgzip not found. Please install htslib." >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found. Please install samtools." >&2; exit 1; }
command -v wfmash >/dev/null 2>&1 || { echo "Error: wfmash not found. Please install wfmash." >&2; exit 1; }

# Function to process a single FASTA file
process_fasta() {
    local fasta_file="$1"
    local base_name=$(basename "$fasta_file" .fa)
    local output_dir="$base_name"
    
    echo "Processing $fasta_file -> $output_dir"
    
    # Create output directory if it doesn't exist
    if [ ! -d "$output_dir" ]; then
        echo "Error: Output directory $output_dir does not exist!"
        return 1
    fi
    
    # Store the absolute path to the script directory
    local script_dir=$(pwd)
    
    # Create a temporary working directory
    local temp_dir=$(mktemp -d)
    cd "$temp_dir" || exit 1
    
    # Copy input file to temp directory
    cp "$script_dir/$fasta_file" "./input.fa"
    
    # Set reference and segsize variables
    REFERENCE="input.fa"
    SEGSIZE=""
    
    # Preprocess the FASTA file
    echo "Preprocessing $fasta_file..."
    python "$script_dir/preprocess.py" "$REFERENCE" "NNNNNNN" $SEED > processed.fa
    
    # Check if trim.bed was created and adjust parameters
    if [ -s trim.bed ]; then
        REFERENCE="processed.fa"
        SEGSIZE="-s $SEED"
        echo "Using processed reference with segsize: $SEGSIZE"
    fi
    
    echo "Using reference: $REFERENCE"
    
    # Compress and index the reference
    echo "Compressing and indexing reference..."
    bgzip -@ 4 "$REFERENCE"
    samtools faidx "${REFERENCE}.gz"
    
    # Perform community detection based on pairwise alignments
    # 90% ID at mash level with 6 mappings per segment
    echo "Running wfmash for community detection..."
    if [ -s trim.bed ]; then
        wfmash "${REFERENCE}.gz" -p 90 -n 28 $SEGSIZE -t $CPUS -m > genomes.mapping.tmp
        python "$script_dir/editpaf.py" genomes.mapping.tmp trim.bed > genomes.mapping.paf
    else
        wfmash "${REFERENCE}.gz" -p 90 -n 28 $SEGSIZE -t $CPUS -m > genomes.mapping.paf
    fi
    
    # Convert PAF mappings into a network
    echo "Converting PAF to network..."
    python "$script_dir/paf2net.py" -p genomes.mapping.paf
    
    # Generate communities plot
    echo "Generating communities plot..."
    python "$script_dir/net2communities.py" -e genomes.mapping.paf.edges.list.txt \
                                           -w genomes.mapping.paf.edges.weights.txt \
                                           -n genomes.mapping.paf.vertices.id2name.txt \
                                           --plot
    
    # Move output files to the target directory
    echo "Moving outputs to $output_dir..."
    
    # Find and move the PDF file(s) - community detection usually creates a PDF
    find . -name "*.pdf" -exec mv {} "$script_dir/$output_dir/" \;
    
    # Optionally move other relevant output files
    if [ -f "genomes.mapping.paf.edges.list.txt" ]; then
        mv genomes.mapping.paf.edges.list.txt "$script_dir/$output_dir/"
    fi
    if [ -f "genomes.mapping.paf.edges.weights.txt" ]; then
        mv genomes.mapping.paf.edges.weights.txt "$script_dir/$output_dir/"
    fi
    if [ -f "genomes.mapping.paf.vertices.id2name.txt" ]; then
        mv genomes.mapping.paf.vertices.id2name.txt "$script_dir/$output_dir/"
    fi
    
    # Clean up
    cd "$script_dir"
    rm -rf "$temp_dir"
    
    echo "Finished processing $fasta_file"
    echo "----------------------------------------"
}

# Main execution
echo "Community Detection Pipeline"
echo "CPUs: $CPUS, Seed: $SEED"
echo "========================================"

# Process each FASTA file
failed_files=()
for fasta_file in LSDV*.fa; do
#for fasta_file in aligned*.fa; do
    if [ -f "$fasta_file" ]; then
        if ! process_fasta "$fasta_file"; then
            failed_files+=("$fasta_file")
            echo "WARNING: Failed to process $fasta_file"
        fi
    fi
done

# Summary
echo "========================================"
echo "Processing complete!"

if [ ${#failed_files[@]} -eq 0 ]; then
    echo "All files processed successfully."
    exit 0
else
    echo "Failed to process ${#failed_files[@]} file(s):"
    printf ' - %s\n' "${failed_files[@]}"
    exit 1
fi
