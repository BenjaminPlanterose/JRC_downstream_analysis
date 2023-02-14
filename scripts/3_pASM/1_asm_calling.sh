#!bin/bash

# WB
biscuit epiread -P -B pooled_blood_data_hg38_snps.bed /media/nw_disk2/bronte/cord_blood_data/2_samples/reference_genome/hg38.fa /media/nw_disk2/bronte/38_alignment/pooled_blood/output_38.bam | sort -k1,1 -k2,2n -k3,3n > pooled_blood_data_hg38.pairwise.epiread
biscuit asm pooled_blood_data_hg38.pairwise.epiread > pooled_blood_data_hg38.asm


# sperm
biscuit epiread -P -B pooled_sperm_data_hg38_snps.bed /media/nw_disk2/bronte/cord_blood_data/2_samples/reference_genome/hg38.fa /media/nw_disk2/bronte/38_alignment/pooled_sperm/output_38_sperm.bam
biscuit asm pooled_sperm_data_hg38.pairwise.epiread > pooled_sperm_data_hg38.asm

# MERGED
biscuit epiread -P -B merged_hg38_snps.bed /media/nw_disk2/bronte/cord_blood_data/2_samples/reference_genome/hg38.fa /media/nw_disk2/bronte/merged_hg38_analysis/merged_38.bam | sort -k1,1 -k2,2n -k3,3n > merged_hg38.pairwise.epiread
biscuit asm merged_hg38.pairwise.epiread > merged_hg38.asm
