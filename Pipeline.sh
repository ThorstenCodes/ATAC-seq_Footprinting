#!/usr/bin/env bash
set -euo pipefail

# Arguments
PATH_BAM_FILES=$1   # Folder with BAM files
PEAKS=$2
GENOME=$3
BLACKLIST=$4
MOTIFS=$5
OUTDIR=$6
THREADS=$7
REGIONS=${8:-} # can be not given
MODUS=$9 #Single or Differential
CONTROL=${10:-}
TREATMENT=${11:-}

# Create output directory for ATACorrect
mkdir -p "${OUTDIR}/ATACorrect_output"
mkdir -p "${OUTDIR}/Footprint_Scores"
mkdir -p "${OUTDIR}/BINDetect"

# Loop through BAM files
for bam_file in "${PATH_BAM_FILES}"/*.bam; do
    # Extract sample name (basename without extensions)
    filename=$(basename "${bam_file}")           # e.g. sample.mLb.clN.sorted.bam
    SAMPLE=${filename%%.*}                       # take everything before first '.'

    echo "Processing ${SAMPLE} with BAM: ${bam_file}"

    # --- ATACorrect from TOBIAS Pipeline ---
    TOBIAS ATACorrect \
        --bam "${bam_file}" \
        --genome "${GENOME}" \
        --peaks "${PEAKS}" \
        --blacklist "${BLACKLIST}" \
        --outdir "${OUTDIR}/ATACorrect_output" \
        --cores "${THREADS}"

    # --- Decide which region file to use to score bigwigs , can save computing time if less peaks ---
    if [[ -n "$REGIONS" ]]; then
        regions_file=$REGIONS
    else
        regions_file=$PEAKS
    fi

    # --- ScoreBigwig/ FootprintScores from TOBIAS pipeline ---
    TOBIAS FootprintScores \
    --signal "${OUTDIR}/ATACorrect_output/${SAMPLE}"*_corrected.bw \
    --regions "${regions_file}" \
    --output "${OUTDIR}/Footprint_Scores/${SAMPLE}_footprints_score.bw"  \
    --cores ${THREADS}

done


## Run Binddetect for single condition or differential binding

# Generate Peak Header file before running BINDetect
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' "${REGIONS}" > "${OUTDIR}/BINDetect/Header_ATAC_Peaks.txt"

# Check if user chose Single or Differential Analysis
if [[ "$MODUS" == "Single" ]]; then
    # Loop over all footprint bigWig files
    for bw_file in "${OUTDIR}/Footprint_Scores/"*_footprints_score.bw; do
        # Extract sample name from filename
        SAMPLE=$(basename "$bw_file" | sed 's/_footprints_score.bw//')

        echo "Running BINDetect for single sample: $SAMPLE"

        TOBIAS BINDetect \
            --motifs "${MOTIFS}" \
            --signals "$bw_file" \
            --genome "${GENOME}" \
            --peaks "${regions_file}" \
            --peak_header "${OUTDIR}/BINDetect/Header_ATAC_Peaks.txt" \
            --outdir "${OUTDIR}/BINDetect/${SAMPLE}" \
            --cond_names "${SAMPLE}"
    done

#if user chose Differential analysis
elif [[ "$MODUS" == "Differential" ]]; then
    # Make sure at least two samples exist
    if [[ -z "$SAMPLE_A" || -z "$SAMPLE_B" -lt 2 ]]; then
        echo "Error: For Differential mode, please provide SAMPLE_A and SAMPLE_B as arguments."
        exit 1
    fi

    echo "Running BINDetect for differential analysis: $SAMPLE_A vs $SAMPLE_B"

    TOBIAS BINDetect \
        --motifs "${MOTIFS}" \
        --signals "${OUTDIR}/Footprint_Scores/${SAMPLE_A}_footprints_score.bw" \
                  "${OUTDIR}/Footprint_Scores/${SAMPLE_B}_footprints_score.bw" \
        --genome "${GENOME}" \
        --peaks "${regions_file}" \
        --peak_header "${OUTDIR}/BINDetect/Header_ATAC_Peaks.txt" \
        --outdir "${OUTDIR}/BINDetect/${SAMPLE_A}_vs_${SAMPLE_B}" \
        --cond_names "${SAMPLE_A}" "${SAMPLE_B}"

else
    echo "Error: MODUS not given. Choose between 'Single' or 'Differential'."
    exit 1
fi


#Plottracks
    # ============================== CONFIG ==============================
REGIONS="/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/ATAC_peaks_overlapping_HiC_TSS_only_whole_promoter.bed"
GTF="gtf_no_chr.gtf"

# BigWig files
BIGWIGS=(
    "/media/rad/HDD1/ChipSeqThorstenWS6/TK_2022_Foxp1_pvalue005_consensus3/bwa/mergedLibrary/bigwig/PKF1OE_R4.bigWig"
    "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/bigwig/GFP0h.mRp.clN.bigWig"
    "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/GFP0h.mRp.clN.sorted_corrected.bw"
    "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/bigwig/KrFo24h.mRp.clN.bigWig"
    "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/KrFo24h.mRp.clN.sorted_corrected.bw"
)

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

# ============================== MAIN LOOP ==============================
Conditions=("GFP0h" "KrFo24h")

for Sample in "${Conditions[@]}"; do
    TFBS="merged_selected_TFs_list_${Sample}_bound_unique.bed"

    if [[ ! -f "$TFBS" ]]; then
        echo "?? Warning: TFBS file not found for $Sample: $TFBS"
        continue
    fi

    OUTDIR="TOBIAS_latest_samplespecific_final/${Sample}"
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

### GTF Filtering

#!/bin/bash

# ========= CONFIGURE THIS =========
INPUT_GTF="gencode.vM10.annotation.gtf.gz"                     # Path to compressed GTF
GENE_LIST=("Foxp1" "Frmd4b" "Lmod3" )                                    # Genes to keep
OUTPUT_GTF="filtered_TF_genes.gtf"                            # Output path
# ==================================

# Create a regex pattern like: (IRF1|IRF2|...)
PATTERN=$(IFS=\| ; echo "${GENE_LIST[*]}")

# Decompress and filter
zcat "$INPUT_GTF" | awk -v pat="gene_name \"(${PATTERN})\"" '
    $0 ~ pat { print }
' > "$OUTPUT_GTF"

echo "Filtered GTF written to: $OUTPUT_GTF"
