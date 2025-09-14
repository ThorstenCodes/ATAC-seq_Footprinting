#!/usr/bin/env bash
set -euo pipefail

# Arguments
PATH_BAM_FILES=$1   # Folder with BAM files
PEAKS=$2 # needs to be .bed format
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
# To generate Header automatically we first count columsn in region file
ncol=$(awk 'NR==1{print NF}' "${regions_file}")

# Generate the header for the number of columns, if no column name is given than just Col1, Col2, Col3 is written in the header.
awk -v n=$ncol 'BEGIN{
    OFS="\t";
    names[1]="Chr"; names[2]="Start"; names[3]="End"; names[4]="PeakID";
    names[5]="Strand"; names[6]="Score"; names[7]="Annotation";
    names[8]="GeneID"; names[9]="GeneName";
    for (i=1; i<=n; i++) {
        if (i>1) printf OFS;
        if (i in names) printf names[i]; else printf "Col" i;
    }
    print ""
}' > "${OUTDIR}/BINDetect/Header_ATAC_Peaks.txt"


# Check if user chose Single or Differential Analysis
if [[ "$MODUS" == "Single" ]]; then
    # Loop over all footprint bigWig files
    for bw_file in "${OUTDIR}/Footprint_Scores/"*_footprints_score.bw; do
        # Extract sample name from filename (strip path and suffix (replacement with nothing through sed command)).
        # This identifier will be used for folder structure and tracking.

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
    # Make sure CONTROL and TREATMENT Sample were provided
    if [[ -z "$CONTROL" || -z "$TREATMENT" ]]; then
        echo "Error: For Differential mode, please provide CONTROL and TREATMENT as arguments."
        exit 1
    fi

    echo "Running BINDetect for differential analysis: $CONTROL vs $TREATMENT"

    TOBIAS BINDetect \
        --motifs "${MOTIFS}" \
        --signals "${OUTDIR}/Footprint_Scores/${CONTROL}_footprints_score.bw" \
                  "${OUTDIR}/Footprint_Scores/${TREATMENT}_footprints_score.bw" \
        --genome "${GENOME}" \
        --peaks "${regions_file}" \
        --peak_header "${OUTDIR}/BINDetect/Header_ATAC_Peaks.txt" \
        --outdir "${OUTDIR}/BINDetect/${CONTROL}_vs_${TREATMENT}" \
        --cond_names "${CONTROL}" "${TREATMENT}"

else
    echo "Error: MODUS not given. Choose between 'Single' or 'Differential'."
    exit 1
fi


#### With Cmd + Shift + 7 one can toggle commands in Visual Studio Code!!!

#Testrun:
bash Pipeline.sh /media/rad/HDD1/TK_Test/ATAC_Pipeline \
      /media/rad/HDD1/TK_Test/ATAC_Pipeline/macs/narrowPeak/consensus/consensus_peaks.mRp.clN.bed \
      /media/rad/HDD1/2023_AnalysisTK/mmu_GRCm38_gencode.fa \
      /media/rad/HDD1/ChipSeqThorstenWS6/Blacklists/resources/mm10-blacklist.v2.bed \
      /media/rad/HDD1/motif_databases.12.23/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format_plus_Foxp1_GN.meme  \
      /media/rad/HDD1/TK_Test/ATAC_Pipeline \
      40 \
      /media/rad/HDD1/TK_Test/ATAC_Pipeline/ATAC_peaks_overlapping_HiC_TSS_including_promoter_region.bed \
      Differential \
      GFP0h \
      PKF1OE






# #Plottracks
#     # ============================== CONFIG ==============================
# REGIONS="/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/ATAC_peaks_overlapping_HiC_TSS_only_whole_promoter.bed"
# GTF="gtf_no_chr.gtf"

# # BigWig files
# BIGWIGS=(
#     "/media/rad/HDD1/ChipSeqThorstenWS6/TK_2022_Foxp1_pvalue005_consensus3/bwa/mergedLibrary/bigwig/PKF1OE_R4.bigWig"
#     "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/bigwig/GFP0h.mRp.clN.bigWig"
#     "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/GFP0h.mRp.clN.sorted_corrected.bw"
#     "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/bigwig/KrFo24h.mRp.clN.bigWig"
#     "/media/rad/HDD1/2023_AnalysisTK/2023_ACC_ATAC/Nextflow_Results_D_A_Comp/bwa/mergedReplicate/2025_Footprinting/KrFo24h.mRp.clN.sorted_corrected.bw"
# )

# LABELS=("Chip" "GFPATAC" "GFPFoot" "PKRFATAC" "PKRFFoot")

# COLORS=(
#     "green"
#     "orange"
#     "blue"
#     "orange"
#     "blue"
# )

# # Validate
# if [[ ${#BIGWIGS[@]} -ne ${#LABELS[@]} ]] || [[ ${#BIGWIGS[@]} -ne ${#COLORS[@]} ]]; then
#     echo "ERROR: Mismatch between number of bigwigs, labels, and colors"
#     exit 1
# fi

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
