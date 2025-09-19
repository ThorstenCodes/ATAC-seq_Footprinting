#!/bin/bash

### This pipeline will follow the first Pipeline to detect Differential Footprints between conditions or just Footprint analysis.
### If you are interested in the certain TFs and where they bind based on the footprints detected in a certain non-coding region or promoter region
### you can give this region here and see on a Plottrack where they are. Furthermore you can analyse if there is different binding detected for this TFs in the given region!

# ------------------------
#  Pipefail integration
# ------------------------
set -euo pipefail
trap 'echo "Error on line $LINENO. Exiting."; exit $?' ERR


# ----------------
# Input Variables
# ----------------

SPECIES=$1    # Human or Mouse
TF_INPUT=$2   # Either: "IRF1 IRF2 ..." OR path to file (tfs.txt)
LOCI=$3       # path to bed.file format with regions chr, start, end, name_of_peak, other additional columns, one line per Locus
SAMPLE_TABLE=${4:-} # tab delimited file with no header: Column1:path-to-bigwig files, Column 2: Sample name, Column 3: Track Color
PATH_TO_BIGWIGS=${5:-}  # PATH to bigwig you want to see in figure. Can be footprint score, can be a normalized bigwig peak file.
NAMES_OF_BW_SAMPLE=${6:-} # often the same as SAMPLE_Name. Identifier for the Sample such as Sample_1, PK, ...
TRACK_COLORS=${7:-} # Color of the peaks in the bigwig shown in the figure, requires same amount as Sample names and bigwig samples
PATH_TO_TF_DIRECTORIES=$8 # Where are the BINDetect_outputs saved you want as input to vizualize the footprints in the genome
SAMPLE_NAMES=($9) # not required when SAMPLE_TABLE is given but required if not



mkdir -p "$PWD/Gene_Annotation/"

# ----------------------------------------------------------------------
# Download GTF depending on species and remove Chr if written with Chr
# ----------------------------------------------------------------------

if [[ "$SPECIES" == "Human" ]]; then
    INPUT_GTF="$PWD/Gene_Annotation/gencode.v49.basic.annotation.gtf.gz"
    if [ -f "$INPUT_GTF" ]; then
        echo "Human GTF already exists, skipping download."
    else
        wget -nc -P "$PWD/Gene_Annotation/" \
            https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
    fi
elif [[ "$SPECIES" == "Mouse" ]]; then
    INPUT_GTF="$PWD/Gene_Annotation/gencode.vM38.basic.annotation.gtf.gz"
    if [ -f "$INPUT_GTF" ]; then
        echo "Mouse GTF already exists, skipping download."
    else
        wget -nc -P "$PWD/Gene_Annotation/" \
            https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz
    fi
else
    echo 'Error: Species must be "Mouse" or "Human".'
    exit 1
fi

# Now INPUT_GTF is always defined
echo "Using GTF: $INPUT_GTF"

# --- Normalize chromosome names (remove "chr" prefix if present) ---

# Since the output will be the chr1 (which after stripping is 1 only) we need to add the OR operator and true otherwise pipefail thinks it is an error code 1, although its literally the output!!!
first_chr=$(gunzip -c "$INPUT_GTF" | awk '$1 !~ /^#/ {print $1; exit}' || true)

if [[ -z "$first_chr" ]]; then
    echo "Error: Could not detect first chromosome from $INPUT_GTF"
    exit 1
fi


if [[ "$first_chr" == chr* ]]; then
    echo "Detected 'chr' prefix in chromosomes → stripping it (overwriting $INPUT_GTF)"
    tmpfile=$(mktemp)
    gunzip -c "$INPUT_GTF" | sed 's/^chr//' | gzip > "$tmpfile"
    mv "$tmpfile" "$INPUT_GTF"
else
    echo "Chromosomes already without 'chr' prefix → ensuring compressed file"
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

# -------------------------
# Initialize arrays
# -------------------------
BIGWIGS=()
LABELS=()
COLORS=()

# -------------------------
# Fill arrays
# -------------------------
if [[ -f "$SAMPLE_TABLE" ]]; then
    echo "Loading BIGWIGS, LABELS, COLORS from sample table: $SAMPLE_TABLE"

    # Skip header lines starting with '#' and take columns 1-3
    mapfile -t BIGWIGS < <(awk 'NR>1 && $1!~/^#/{print $1}' "$SAMPLE_TABLE")
    mapfile -t LABELS  < <(awk 'NR>1 && $1!~/^#/{print $2}' "$SAMPLE_TABLE")
    mapfile -t COLORS  < <(awk 'NR>1 && $1!~/^#/{print $3}' "$SAMPLE_TABLE")

elif [[ -n "$PATH_TO_BIGWIGS" ]] && [[ -n "$NAMES_OF_BW_SAMPLE" ]] && [[ -n "$TRACK_COLORS" ]]; then
    echo "Loading BIGWIGS, LABELS, COLORS from provided arguments"
    BIGWIGS=($PATH_TO_BIGWIGS)
    LABELS=($NAMES_OF_BW_SAMPLE)
    COLORS=($TRACK_COLORS)
else
    echo "ERROR: No valid sample table or arrays provided for BIGWIGS, LABELS, COLORS"
    exit 1
fi

# -------------------------
# Check consistency
# -------------------------
if [[ ${#BIGWIGS[@]} -ne ${#LABELS[@]} ]] || [[ ${#BIGWIGS[@]} -ne ${#COLORS[@]} ]]; then
    echo "ERROR: Number of BIGWIGS (${#BIGWIGS[@]}), LABELS (${#LABELS[@]}), and COLORS (${#COLORS[@]}) do not match"
    exit 1
fi

echo "Loaded ${#BIGWIGS[@]} samples successfully."


# ============================== GENERATE PLOTTracks for each Sample ==============================

run_plottracks() {
    local sample=$1
    local TFBS="$PWD/Gene_Annotation/merged_TFs_${sample}_bound.bed"
    local OUTDIR="PLOTTracks/${sample}"
    mkdir -p "$OUTDIR"

    echo "Running TOBIAS PlotTracks for ${sample}..."

    TOBIAS PlotTracks \
        --bigwigs "${BIGWIGS[@]}" \
        --regions "$REGIONS" \
        --sites "$TFBS" \
        --gtf "$INPUT_GTF" \
        --labels "${LABELS[@]}" \
        --colors "${COLORS[@]}" \
        --width 15 \
        --max-transcripts 1 \
        --outdir "$OUTDIR"

    echo "PlotTracks completed for ${sample} → results in ${OUTDIR}"
}

if [[ -f "$SAMPLE_TABLE" ]]; then
    for sample in $(awk '{print $2}' "$SAMPLE_TABLE"); do
        run_plottracks "$sample"
    done
elif [[ -n "${SAMPLE_NAMES:-}" ]]; then
    for sample in "${SAMPLE_NAMES[@]}"; do
        run_plottracks "$sample"
    done
else
    echo "Error: Provide either SAMPLE_TABLE (with sample names in column 2) or SAMPLE_NAMES array."
    exit 1
fi





# Create a regex pattern like which introduces OR (|} betweent the gene names.This allows awk to filter for each gene name in the list. The pattern looks like this in the end: e.g. (IRF1|IRF2|...)
# PATTERN=$(echo "$TF_LIST" | tr ' ' '|')

# -----------------------------
# Filter GTF ### we dont need to filter GTF but probably remove the Chr or check and remove if it is there!!!!! Instead we use htere the output fo the BINDetect files and merge them to lists
# -----------------------------


# zcat "$INPUT_GTF" | awk -v pat="gene_name \"($PATTERN)\"" '
#     $0 ~ pat { print }
# ' > "$OUTPUT_GTF"

# echo "Filtered GTF written to: $OUTPUT_GTF"
