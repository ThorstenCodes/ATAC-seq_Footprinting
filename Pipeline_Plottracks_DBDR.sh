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
SAMPLE_TABLE=${4:-} # tab delimited .txt file : Column1:path-to-bigwig files, Column 2: Sample name, Column 3: Track Color, no path to the file, only filename.
PATH_TO_BIGWIGS=${5:-}  # PATH to bigwig you want to see in figure. Can be footprint score, can be a normalized bigwig peak file.
NAMES_OF_BW_SAMPLE=${6:-} # often the same as SAMPLE_Name. Identifier for the Sample such as Sample_1, PK, ...
TRACK_COLORS=${7:-} # Color of the peaks in the bigwig shown in the figure, requires same amount as Sample names and bigwig samples
PATH_TO_TF_DIRECTORIES=$8 # Where are the BINDetect_outputs saved you want as input to vizualize the footprints in the genome e.g /media/group/project/Footprinting/...
SAMPLE_NAMES=($9) # not required when SAMPLE_TABLE is given but required if not
CONTROL=${10:-} # what is the control sample for the differenetial comparision, can also be empty ""
UBA=${11:-} # was TF found as bound, unbound or all



mkdir -p "$PWD/Gene_Annotation/"
mkdir -p "$PWD/selected_TF_Footprints"
LOCI=$(echo "$LOCI" | sed "s/r\$//")

# ----------------------------------------------------------------------
# Download GTF depending on species and remove Chr if written with Chr
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Download GTF depending on species and remove Chr if written with Chr
# ----------------------------------------------------------------------

if [[ "$SPECIES" == "Human" ]]; then
    INPUT_GTF_GZ="$PWD/Gene_Annotation/gencode.v49.basic.annotation.gtf.gz"
    INPUT_GTF="$PWD/Gene_Annotation/gencode.v49.basic.annotation.gtf"  # uncompressed

    if [[ ! -f "$INPUT_GTF" ]]; then
        if [[ ! -f "$INPUT_GTF_GZ" ]]; then
            wget -nc -P "$PWD/Gene_Annotation/" \
                https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
        fi
        echo "Unzipping Human GTF..."
        gunzip -c "$INPUT_GTF_GZ" | sed 's/^chr//' > "$INPUT_GTF"
    else
        echo "Human GTF already exists (uncompressed)."
    fi

elif [[ "$SPECIES" == "Mouse" ]]; then
    INPUT_GTF_GZ="$PWD/Gene_Annotation/gencode.vM38.basic.annotation.gtf.gz"
    INPUT_GTF="$PWD/Gene_Annotation/gencode.vM38.basic.annotation.gtf"  # uncompressed

    if [[ ! -f "$INPUT_GTF" ]]; then
        if [[ ! -f "$INPUT_GTF_GZ" ]]; then
            wget -nc -P "$PWD/Gene_Annotation/" \
                https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz
        fi
        echo "Unzipping Mouse GTF..."
        gunzip -c "$INPUT_GTF_GZ" | sed 's/^chr//' > "$INPUT_GTF"
    else
        echo "Mouse GTF already exists (uncompressed)."
    fi

else
    echo 'Error: Species must be "Mouse" or "Human".'
    exit 1
fi


# Now INPUT_GTF is always defined
echo "Using GTF: $INPUT_GTF"

# --- Normalize chromosome names (remove "chr" prefix if present) ---

# Since the output will be the chr1 (which after stripping is 1 only) we need to add the OR operator and true otherwise pipefail thinks it is an error code 1, although its literally the output!!!
first_chr=$(awk '$1 !~ /^#/ {print $1; exit}' "$INPUT_GTF" || true)

if [[ -z "$first_chr" ]]; then
    echo "Error: Could not detect first chromosome from $INPUT_GTF"
    exit 1
fi


if [[ "$first_chr" == chr* ]]; then
    echo "Detected 'chr' prefix in chromosomes → stripping it (overwriting $INPUT_GTF)"
    tmpfile=$(mktemp)
    gunzip -c "$INPUT_GTF" | sed 's/^chr//' > "$tmpfile"
    mv "$tmpfile" "$INPUT_GTF"
else
  tmpfile=$(mktemp)
    echo "Chromosomes already without 'chr' prefix → ensuring compressed file"
fi



# -----------------------------
# Parse TF list (file or inline)
# -----------------------------

if [[ -f "$TF_INPUT" ]]; then
    # File with TF names → split into array
    mapfile -t TF_LIST < <(tr -s ' \t' '\n' < "$TF_INPUT")
else
    # Inline string → split on spaces into array
    read -r -a TF_LIST <<< "$TF_INPUT"
fi


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
    mapfile -t BIGWIGS < <(awk 'NR>1 && $1!~/^#/ && NF>=3 {print $1}' "$SAMPLE_TABLE")
    mapfile -t LABELS  < <(awk 'NR>1 && $1!~/^#/ && NF>=3 {print $2}' "$SAMPLE_TABLE")
    mapfile -t COLORS < <(awk 'NR>1 && $1!~/^#/ && NF>=3 {gsub("\r","",$3); print $3}' "$SAMPLE_TABLE")



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


# --------------------------------------------------------------------------
# Merge TF.bed files and remove their duplicates if different conditions
# -------------------------------------------------------------------------

if [[ -f "$SAMPLE_TABLE" ]]; then
    echo "Reading sample names from table: $SAMPLE_TABLE"
    # read sample names from column 2, skipping header and comments
    mapfile -t SAMPLE_NAMES < <(awk 'NR>1 && $1!~/^#/ && NF>=3 {if(!seen[$2]++){print $2}}' "$SAMPLE_TABLE")
    echo "SAMPLE_NAMES: ${SAMPLE_NAMES[*]}"
else
    echo "Using SAMPLE_NAMES directly: ${SAMPLE_NAMES[*]}"
    echo "SAMPLE_NAMES: ${SAMPLE_NAMES[*]}"
fi

if [[ $UBA == "bound" || $UBA == "unbound" ]]; then
    # per-sample merging
    for Sample in "${SAMPLE_NAMES[@]}"; do
        OUTPUT="$PWD/selected_TF_Footprints/merged_TFs_${Sample}_${UBA}.bed"
        > "$OUTPUT"
        echo "Processing sample: $Sample"

        for TF in "${TF_LIST[@]}"; do
            BEDFILE="$PATH_TO_TF_DIRECTORIES/_${TF}/beds/_${TF}_${Sample}_${UBA}.bed"
            if [[ -f "$BEDFILE" ]]; then
                echo "Adding $TF to $OUTPUT"
                cat "$BEDFILE" >> "$OUTPUT"
            else
                echo "Warning: $BEDFILE not found, skipping $TF"
            fi
        done

        # Sort and remove duplicates per sample
        TMP="${OUTPUT%.bed}.tmp"
        sort -k1,1 -k2,2n "$OUTPUT" | uniq > "$TMP" && mv "$TMP" "$OUTPUT"
        echo "Merged file written to ${OUTPUT}"
    done

elif [[ $UBA == "all" ]]; then
    # single merged file for all TFs, no per-sample dependency
    OUTPUT="$PWD/selected_TF_Footprints/merged_TFs_all.bed"
    > "$OUTPUT"
    echo "Merging Footprints for TFs found in all samples"

    for TF in "${TF_LIST[@]}"; do
        BEDFILE="$PATH_TO_TF_DIRECTORIES/_${TF}/beds/_${TF}_all.bed"
        if [[ -f "$BEDFILE" ]]; then
            echo "Adding $TF to $OUTPUT"
            cat "$BEDFILE" >> "$OUTPUT"
        else
            echo "Warning: $BEDFILE not found, skipping $TF"
        fi
    done

    # Sort and remove duplicates for the single all-file
    TMP="${OUTPUT%.bed}.tmp"
    sort -k1,1 -k2,2n "$OUTPUT" | uniq > "$TMP" && mv "$TMP" "$OUTPUT"
    echo "Merged file written to ${OUTPUT}"

else
    echo "UBA variable needs to be bound, unbound or all !!!"
    exit 1
fi




# --------------------------------------------------------------------------
# GENERATE PLOTTracks for each Sample and Loci
# -------------------------------------------------------------------------

run_plottracks() {
    local sample=$1
    local TFBS="$PWD/selected_TF_Footprints/merged_TFs_${sample}_${UBA}.bed"
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

    echo "PlotTracks completed for ${sample}. Results in ${OUTDIR}"
}

if [[ $UBA == "bound" || $UBA == "unbound" ]]; then
    # Run PlotTracks for each sample
    if [[ -f "$SAMPLE_TABLE" ]]; then
        for sample in $(awk 'NR>1 && $1!~/^#/ {if(!seen[$2]++){print $2}}' "$SAMPLE_TABLE"); do
            run_plottracks "$sample"
        done
    elif [[ -n "${SAMPLE_NAMES:-}" ]]; then
        for sample in "${SAMPLE_NAMES[@]}"; do
            run_plottracks "$sample"
        done
    else
        echo "Error: Provide either SAMPLE_TABLE or SAMPLE_NAMES array."
        exit 1
    fi

elif [[ $UBA == "all" ]]; then
    # Run PlotTracks once using the first sample
    if [[ -f "$SAMPLE_TABLE" ]]; then
        first_sample=$(awk 'NR>1 && $1!~/^#/ {print $2; exit}' "$SAMPLE_TABLE")
    elif [[ -n "${SAMPLE_NAMES:-}" ]]; then
        first_sample="${SAMPLE_NAMES[0]}"
    else
        echo "Error: Provide either SAMPLE_TABLE or SAMPLE_NAMES array."
        exit 1
    fi

    echo "UBA=all → running PlotTracks once using sample: $first_sample"
    run_plottracks "$first_sample"

else
    echo "UBA variable must be bound, unbound or all !!!"
    exit 1
fi




# --------------------------------------------------------------------------
#  Differential Motifs of given TFs as Table output
# -------------------------------------------------------------------------

run_DiMo() {
    local sample=$1
    local control=$2
    local PATH_TO_FILES=$PWD/selected_TF_Footprints/
    local OUTDIR="Differential_Motifs/${sample}_vs_${control}"

    mkdir -p "$OUTDIR"

    # Define merged BED files for sample and control
    local SAMPLE_BED="$PATH_TO_FILES/merged_TFs_${sample}_bound.bed"
    local CONTROL_BED="$PATH_TO_FILES/merged_TFs_${control}_bound.bed"

    # Check if BED files exist
    if [[ ! -f "$SAMPLE_BED" || ! -f "$CONTROL_BED" ]]; then
        echo "Error: BED file(s) missing!"
        echo "SAMPLE_BED=$SAMPLE_BED"
        echo "CONTROL_BED=$CONTROL_BED"
        return 1
    fi

    # Define output files
    local OVERLAP_FILE="${OUTDIR}/${sample}_overlap_with_${control}.bed"
    local NON_OVERLAP_FILE="${OUTDIR}/${sample}_non_overlap_with_${control}.bed"

    # Run bedtools
    bedtools intersect -a "$CONTROL_BED" -b "$SAMPLE_BED" -wa -u > "$OVERLAP_FILE"
    echo "Overlapping regions saved to $OVERLAP_FILE"

    bedtools intersect -a "$SAMPLE_BED" -b "$CONTROL_BED" -wa -v > "$NON_OVERLAP_FILE"
    echo "Non-overlapping regions saved to $NON_OVERLAP_FILE"

    echo "Comparison completed. Results saved in: $OUTDIR"
}

# Only run differential motif analysis if UBA="bound" AND CONTROL is defined
if [[ $UBA == "bound" ]]; then
    if [[ -n $CONTROL ]]; then
        for sample in $(awk 'NR>1 && $1!~/^#/ {if(!seen[$2]++){print $2}}' "$SAMPLE_TABLE"); do
            if [[ $sample != $CONTROL ]]; then
                run_DiMo "$sample" "$CONTROL"
            fi
        done
    else
        echo "Differential Motif Analysis is not run, since CONTROL is not defined."
        exit 1
    fi
else
    echo "Differential Motif Analysis is only run for UBA='bound'. Current UBA='$UBA'."
fi


#### we also need to add the motifsearch for no footrprints but only motifs with p-value
