#!/bin/bash

### This pipeline will follow the first Pipeline to detect Differential Footprints between conditions or just Footprint analysis.
### If you are interested in the certain TFs and where they bind based on the footprints detected in a certain non-coding region or promoter region
### you can give this region here and see on a Plottrack where they are. Furthermore you can analyse if there is different binding detected for this TFs in the given region!

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
