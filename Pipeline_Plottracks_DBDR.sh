#!/bin/bash

### This pipeline will follow the first Pipeline to detect Differential Footprints between conditions or just Footprint analysis.
### If you are interested in the certain TFs and where they bind based on the footprints detected in a certain non-coding region or promoter region
### you can give this region here and see on a Plottrack where they are. Furthermore you can analyse if there is different binding detected for this TFs in the given region!

SPECIES=$1    # Human or Mouse
TF_INPUT=$2   # Either: "IRF1 IRF2 ..." OR path to file (tfs.txt)
OUTPUT_GTF=$3 # Output GTF filename (must end with .gtf)
LOCI=$4       # bed.file format with regions chr, start, end, name_of_peak, other additional columns, one line per Locus
PATH_TO_BIGWIGS=${5:-}
NAME_OF_BW_SAMPLE=${6:-}
TRACK_COLORS=${7:-}


mkdir -p "$PWD/Gene_Annotation/"

# -----------------------------
# Download GTF depending on species
# -----------------------------
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

# Create a regex pattern like which introduces OR (|} betweent the gene names.This allows awk to filter for each gene name in the list. The pattern looks like this in the end: e.g. (IRF1|IRF2|...)
PATTERN=$(echo "$TF_LIST" | tr ' ' '|')

# -----------------------------
# Filter GTF
# -----------------------------
zcat "$INPUT_GTF" | awk -v pat="gene_name \"($PATTERN)\"" '
    $0 ~ pat { print }
' > "$OUTPUT_GTF"

echo "Filtered GTF written to: $OUTPUT_GTF"


# ============================== PLOTTracks ==============================
REGIONS=$LOCI
GTF="$INPUT_GTF"

# BigWig files
BIGWIGS=$PATH_TO_BIGWIGS

LABELS=("Chip" "GFPATAC" "GFPFoot" "PKRFATAC" "PKRFFoot")

COLORS=(
    "green"
    "orange"
    "blue"
    "orange"
    "blue"
)

# Validate
if [[ ${#BIGWIGS[@]} -ne ${#LABELS[@]} ]] || [[ ${#BIGWIGS[@]} -ne ${#COLORS[@]} ]]; then
    echo "ERROR: Mismatch between number of bigwigs, labels, and colors"
    exit 1
fi

# # ============================== MAIN LOOP ==============================
# Conditions=("GFP0h" "KrFo24h")

# for Sample in "${Conditions[@]}"; do
#     TFBS="merged_selected_TFs_list_${Sample}_bound_unique.bed"

#     if [[ ! -f "$TFBS" ]]; then
#         echo "?? Warning: TFBS file not found for $Sample: $TFBS"
#         continue
#     fi

#     OUTDIR="TOBIAS_latest_samplespecific_final/${Sample}"
#     mkdir -p "$OUTDIR"

#     echo "Running TOBIAS PlotTracks for $Sample..."

#     TOBIAS PlotTracks \
#         --bigwigs "${BIGWIGS[@]}" \
#         --regions "$REGIONS" \
#         --sites "$TFBS" \
#         --gtf "$GTF" \
#         --labels "${LABELS[@]}" \
#         --colors "${COLORS[@]}" \
#         --width 15 \
#         --max-transcripts 1 \
#         --outdir "$OUTDIR"

#     echo "PlotTracks completed for $Sample ? $OUTDIR"
# done

# ### GTF Filtering

# #!/bin/bash

# # ========= CONFIGURE THIS =========
# INPUT_GTF="gencode.vM10.annotation.gtf.gz"                     # Path to compressed GTF
# GENE_LIST=("Foxp1" "Frmd4b" "Lmod3" )                                    # Genes to keep
# OUTPUT_GTF="filtered_TF_genes.gtf"                            # Output path
# # ==================================

# # Create a regex pattern like: (IRF1|IRF2|...)
# PATTERN=$(IFS=\| ; echo "${GENE_LIST[*]}")

# # Decompress and filter
# zcat "$INPUT_GTF" | awk -v pat="gene_name \"(${PATTERN})\"" '
#     $0 ~ pat { print }
# ' > "$OUTPUT_GTF"

# echo "Filtered GTF written to: $OUTPUT_GTF"
