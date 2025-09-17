#!/bin/bash

### This pipeline will follow the first Pipeline to detect Differential Footprints between conditions or just Footprint analysis.
### If you are interested in the certain TFs and where they bind based on the footprints detected in a certain non-coding region or promoter region
### you can give this region here and see on a Plottrack where they are. Furthermore you can analyse if there is different binding detected for this TFs in the given region!

SPECIES=$1    # Human or Mouse
TF_INPUT=$2   # Either: "IRF1 IRF2 ..." OR path to file (tfs.txt)
OUTPUT_GTF=$3 # Output GTF filename (must end with .gtf)
LOCI=$4       # bed.file format with regions chr, start, end, name_of_peak, other additional columns, one line per Locus
SAMPLE_TABLE=$5
PATH_TO_BIGWIGS=${6:-}
NAMES_OF_BW_SAMPLE=${7:-}
TRACK_COLORS=${8:-}
PATH_TO_TF_DIRECTORIES=$9
SAMPLE_NAMES=$10



mkdir -p "$PWD/Gene_Annotation/"

# ----------------------------------------------------------------------
# Download GTF depending on species and remove Chr if written with Chr
# ----------------------------------------------------------------------

if [[ "$SPECIES" == "Human" ]]; then
    wget -nc -P "$PWD/Gene_Annotation/" \
        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
    INPUT_GTF="$PWD/Gene_Annotation/gencode.v49.basic.annotation.gtf.gz"
elif [[ "$SPECIES" == "Mouse" ]]; then
    wget -nc -P "$PWD/Gene_Annotation/" \
        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz
    INPUT_GTF="$PWD/Gene_Annotation/gencode.vM38.basic.annotation.gtf.gz"
else
    echo 'Error: Species must be "Mouse" or "Human" as argument.'
    exit 1
fi

# --- Normalize chromosome names (remove "chr" prefix if present) ---
first_chr=$(zcat "$INPUT_GTF" | awk '$1 !~ /^#/ {print $1; exit}')

if [[ "$first_chr" == chr* ]]; then
    echo "Detected 'chr' prefix in chromosomes → stripping it"
    zcat "$INPUT_GTF" | sed 's/^chr//' > "${INPUT_GTF%.gz}.nochr.gtf"
    INPUT_GTF="${INPUT_GTF%.gz}.nochr.gtf"
else
    echo "Chromosomes already without 'chr' prefix"
    gunzip -c "$INPUT_GTF" > "${INPUT_GTF%.gz}"
    INPUT_GTF="${INPUT_GTF%.gz}"
fi


# -----------------------------
# Parse TF list (file or inline)
# -----------------------------
if [[ -f "$TF_INPUT" ]]; then
    # Case: file provided → read one gene per line (or tab/space)
    TF_LIST=$(tr '\t' ' ' < "$TF_INPUT" | tr '\n' ' ')
else
    # Case: quoted string list
    TF_LIST=$TF_INPUT
fi

# --------------------------------------------------------------------------
# Merge TF.bed files and remove their duplicates if different conditions
# -------------------------------------------------------------------------
for Sample in "${SAMPLE_NAMES[@]}"; do
    OUTPUT="$PWD/Gene_Annotation/merged_TFs_${Sample}_bound.bed"
    > "$OUTPUT"   # create/empty output file

    echo "Processing sample: $Sample"

    for TF in "${TF_LIST[@]}"; do
        BEDFILE="$PATH_TO_TF_DIRECTORIES/_${TF}/beds/_${TF}_${Sample}_bound.bed"

        if [[ -f "$BEDFILE" ]]; then
            echo "Adding $TF to $OUTPUT"
            cat "$BEDFILE" >> "$OUTPUT"
        else
            echo "Warning: $BEDFILE not found, skipping $TF"
        fi
    done

    # Sort and remove duplicates safely (use temp file, then replace)
    TMP="${OUTPUT%.bed}.tmp"
    sort -k1,1 -k2,2n "$OUTPUT" | uniq > "$TMP" && mv "$TMP" "$OUTPUT"

    echo "Merged file written to ${OUTPUT}"
done


# ============================== CONFIGURE VARIABLES for PLOTTrack ==============================
REGIONS=$LOCI
GTF="$INPUT_GTF"

if [[ -f "$SAMPLE_TABLE" ]]; then
    BIGWIGS=($(awk '{print $1}' "$SAMPLE_TABLE"))
    LABELS=($(awk '{print $2}' "$SAMPLE_TABLE"))
    COLORS=($(awk '{print $3}' "$SAMPLE_TABLE"))
else
    BIGWIGS="${PATH_TO_BIGWIGS}"
    LABELS=${NAMES_OF_BW_SAMPLE}
    COLORS=${TRACK_COLORS}
fi

# ============================== GENERATE PLOTTracks for each Sample ==============================
for Sample in "${SAMPLE_NAMES[@]}"; do
    TFBS="merged_TFs_${Sample}_bound.bed"

    if [[ ! -f "$TFBS" ]]; then
        echo "?? Warning: TFBS file not found for $Sample: $TFBS"
        continue
    fi

    OUTDIR="PLOTTracks/${Sample}"
    mkdir -p "$OUTDIR"

    echo "Running TOBIAS PlotTracks for $Sample..."

    TOBIAS PlotTracks \
        --bigwigs "${BIGWIGS[@]}" \
        --regions "$REGIONS" \
        --sites "$TFBS" \
        --gtf "$GTF" \
        --labels "${LABELS[@]}" \
        --colors "${COLORS[@]}" \
        --width 15 \
        --max-transcripts 1 \
        --outdir "$OUTDIR"

    echo "PlotTracks completed for $Sample ? $OUTDIR"
done


# Create a regex pattern like which introduces OR (|} betweent the gene names.This allows awk to filter for each gene name in the list. The pattern looks like this in the end: e.g. (IRF1|IRF2|...)
# PATTERN=$(echo "$TF_LIST" | tr ' ' '|')

# -----------------------------
# Filter GTF ### we dont need to filter GTF but probably remove the Chr or check and remove if it is there!!!!! Instead we use htere the output fo the BINDetect files and merge them to lists
# -----------------------------


# zcat "$INPUT_GTF" | awk -v pat="gene_name \"($PATTERN)\"" '
#     $0 ~ pat { print }
# ' > "$OUTPUT_GTF"

# echo "Filtered GTF written to: $OUTPUT_GTF"
