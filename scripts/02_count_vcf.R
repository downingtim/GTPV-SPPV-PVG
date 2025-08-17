#!/bin/bash

# Script to count mutation types in first and last 2500 RECORDS of a VCF file
# Usage: ./count_mutations.sh gfavariants.vcf

VCF_FILE="$1"

if [ -z "$VCF_FILE" ]; then
    echo "Usage: $0 <vcf_file>"
    echo "Example: $0 gfavariants.vcf"
    exit 1
fi

if [ ! -f "$VCF_FILE" ]; then
    echo "Error: File $VCF_FILE not found!"
    exit 1
fi

echo "Analyzing VCF file: $VCF_FILE"
echo "========================================"

# Extract variant lines (skip header lines starting with #)
VARIANTS=$(grep -v "^#" "$VCF_FILE")

if [ -z "$VARIANTS" ]; then
    echo "No variants found in the file!"
    exit 1
fi

# Count total records
TOTAL_RECORDS=$(echo "$VARIANTS" | wc -l)
echo "Total variant records: $TOTAL_RECORDS"
echo

# Function to count mutation types in given records
count_mutations() {
    local variants_subset="$1"
    local label="$2"
    
    echo "=== $label ==="
    
    if [ -z "$variants_subset" ]; then
        echo "No variants found in this subset."
        echo "SNV: 0"
        echo "DEL: 0" 
        echo "INS: 0"
        echo "MNP: 0"
        echo "COMPLEX: 0"
        echo "COMPOUND: 0"
        echo "Total variants: 0"
        echo
        return
    fi
    
    # Count each mutation type
    snv_count=$(echo "$variants_subset" | grep -c "TYPE=snv")
    del_count=$(echo "$variants_subset" | grep -c "TYPE=del")
    ins_count=$(echo "$variants_subset" | grep -c "TYPE=ins")
    mnp_count=$(echo "$variants_subset" | grep -c "TYPE=mnp")
    complex_count=$(echo "$variants_subset" | grep -c "TYPE=complex")
    
    # Count compound mutations (lines with multiple TYPE= entries)
    compound_count=$(echo "$variants_subset" | grep -c "TYPE=.*TYPE=")
    
    # Total variants in subset
    total_variants=$(echo "$variants_subset" | wc -l)
    
    echo "SNV: $snv_count"
    echo "DEL: $del_count"
    echo "INS: $ins_count" 
    echo "MNP: $mnp_count"
    echo "COMPLEX: $complex_count"
    echo "COMPOUND: $compound_count"
    echo "Total variants: $total_variants"
    echo
}

# Get first 2500 records (or all if fewer than 2500)
if [ "$TOTAL_RECORDS" -le 2500 ]; then
    FIRST_2500="$VARIANTS"
    echo "File has $TOTAL_RECORDS records (≤2500), using all records for 'first 2500'"
else
    FIRST_2500=$(echo "$VARIANTS" | head -2500)
fi

# Get last 2500 records (or all if fewer than 2500)
if [ "$TOTAL_RECORDS" -le 2500 ]; then
    LAST_2500="$VARIANTS"
    echo "File has $TOTAL_RECORDS records (≤2500), using all records for 'last 2500'"
else
    LAST_2500=$(echo "$VARIANTS" | tail -2500)
fi

# Count mutations in first 2500 records
count_mutations "$FIRST_2500" "First 2500 records"

# Count mutations in last 2500 records  
count_mutations "$LAST_2500" "Last 2500 records"

# Optional: Show detailed breakdown of compound mutations
echo "=== Compound Mutation Details ==="
echo "First 2500 records compound mutations:"
FIRST_COMPOUNDS=$(echo "$FIRST_2500" | grep "TYPE=.*TYPE=")
if [ -n "$FIRST_COMPOUNDS" ]; then
    echo "$FIRST_COMPOUNDS" | while IFS= read -r line; do
        pos=$(echo "$line" | awk '{print $2}')
        types=$(echo "$line" | grep -o "TYPE=[^;]*" | tr '\n' ' ')
        echo "Position $pos: $types"
    done
else
    echo "None found"
fi

echo
echo "Last 2500 records compound mutations:"
LAST_COMPOUNDS=$(echo "$LAST_2500" | grep "TYPE=.*TYPE=")
if [ -n "$LAST_COMPOUNDS" ]; then
    echo "$LAST_COMPOUNDS" | while IFS= read -r line; do
        pos=$(echo "$line" | awk '{print $2}')
        types=$(echo "$line" | grep -o "TYPE=[^;]*" | tr '\n' ' ')
        echo "Position $pos: $types"
    done
else
    echo "None found"
fi