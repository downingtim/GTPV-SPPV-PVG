#!/usr/bin/env bash

GB_FILE="$1"
ALN_FILE="align/genomes.aln"

# -------------------------------------------------------------
# 1. Extract aligned sequence for NC_004003.1
# -------------------------------------------------------------
awk '
  /^>/ {printing=0}
  /^>.*NC_004003\.1/ {printing=1; next}
  printing==1 {printf "%s", $0}
' "$ALN_FILE" > ref_aligned.tmp

if [[ ! -s ref_aligned.tmp ]]; then
    echo "ERROR: NC_004003.1 not found in $ALN_FILE"
    exit 1
fi

# -------------------------------------------------------------
# 2. Extract gene + CDS (simple version that worked for you)
# -------------------------------------------------------------
awk '
  /^ +\/gene=/ {
      g=$0; sub(/.*="/,"",g); sub(/".*/,"",g)
  }
  /^ +\/locus_tag=/ {
      g=$0; sub(/.*="/,"",g); sub(/".*/,"",g)
  }
  /^ +CDS/ {
      loc=$0
      sub(/.*CDS +/,"",loc)
      print g "\t" loc
  }
' "$GB_FILE" > cds_raw.tsv

# -------------------------------------------------------------
# 3. Adjust coordinates & clean formatting
# -------------------------------------------------------------
adjust_coords() {
    loc="$1"

    # remove complement( and )
    loc=$(echo "$loc" | sed 's/complement(//; s/)//')

    # now loc is like 1..480
    start=$(echo "$loc" | cut -d. -f1)
    end=$(echo "$loc"   | cut -d. -f3)

    # count alignment gaps
    gaps_before=$(cut -c1-"$start" ref_aligned.tmp | grep -o "-" | wc -l | tr -d ' ')
    gaps_inside=$(cut -c"$start"-"$end" ref_aligned.tmp | grep -o "-" | wc -l | tr -d ' ')

    new_start=$((start + gaps_before))
    new_end=$((end + gaps_before + gaps_inside))

    # OUTPUT FORMAT YOU REQUESTED:
    #   new_start <tab> new_end   (NO "..", NO complement)
    echo -e "${new_start}\t${new_end}"
}

export -f adjust_coords

# -------------------------------------------------------------
# 4. Output table
# -------------------------------------------------------------
echo -e "gene_ID\tstart\tend"

while IFS=$'\t' read -r gene loc; do
    vals=$(adjust_coords "$loc")
    set -- $vals
    echo -e "${gene}\t${1}\t${2}"
done < cds_raw.tsv

rm -f ref_aligned.tmp
