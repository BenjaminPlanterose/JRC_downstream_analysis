#!bin/bash

# ERR2722068: WB_old
# ERR2722069: sperm_old
# ERR2722070: WB_young
# ERR2722071: sperm_young

#ERR2722068.bam
#ERR2722070.bam
#ERR2722069.bam
#ERR2722071.bam


#### BEFORE RUNNING:


genome='/media/ultron/2tb_disk1/ben/JRC_project/reference_genomes/hg38/hg38.fa'
cd /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/

bams=($( ls *.bam -1v))
l=${#bams[@]}
for ((i=0;i<=$l-1;i+=1))
do
 index="$(cut -d'.' -f1 <<< ${bams[$i]})" 
 echo $index
 biscuit pileup -@ 4 -o $index.vcf $genome ${bams[$i]}
 bgzip $index.vcf
 tabix -p vcf $index.vcf.gz
 biscuit vcf2bed -t cg $index.vcf.gz > $index.bed
done

mv *.bed /media/ultron/2tb_disk1/ben/JRC_project/bed_files/hg38/

cd /media/ultron/2tb_disk1/ben/JRC_project/bed_files/hg38/

beds=($( ls *.bed -1v))
l=${#beds[@]}
for ((i=0;i<=$l-1;i+=1))
do
 index="$(cut -d'.' -f1 <<< ${beds[$i]})" 
 echo $index
 biscuit mergecg $genome ${beds[$i]} > ${index}_mergecg.bed
done








