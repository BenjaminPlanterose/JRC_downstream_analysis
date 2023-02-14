#!bin/bash

# ERR2722068: WB_old
# ERR2722069: sperm_old
# ERR2722070: WB_young
# ERR2722071: sperm_young

# create test
#samtools view -H output_38.bam > test.sam
#samtools view output_38.bam | head -n 1000 >> test.sam
#samtools view -bS test.sam > test.bam

# test
#samtools view -H test.bam > ERR2722068.sam
#cp ERR2722068.sam ERR2722070.sam
#samtools view -@ 4 test.bam | grep ERR2722068 >> ERR2722068.sam
#samtools view -@ 4 test.bam | grep ERR2722070 >> ERR2722070.sam
#samtools view -@ 4 -bS ERR2722068.sam > ERR2722068.bam
#samtools view -@ 4 -bS ERR2722070.sam > ERR2722070.bam
# samtools view -h -@ 4 output_38.bam | grep ERR2722068 | samtools view -@ 4 -bS > ERR2722068.bam
# Looses header to grep




# WB
echo 'Whole blood'
cd /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/pooled_blood/
samtools view -H output_38.bam > ERR2722068.sam
cp ERR2722068.sam ERR2722070.sam

samtools view -@ 4 output_38.bam | grep ERR2722068 >> ERR2722068.sam
samtools view -@ 4 -bS ERR2722068.sam > ERR2722068.bam
mv ERR2722068.bam /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/
rm ERR2722068.sam

samtools view -@ 4 output_38.bam | grep ERR2722070 >> ERR2722070.sam
samtools view -@ 4 -bS ERR2722070.sam > ERR2722070.bam
mv ERR2722070.bam /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/
rm ERR2722070.sam

# Sperm
echo 'Sperm'
cd /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/pooled_sperm/
samtools view -H output_38_sperm.bam > ERR2722069.sam
cp ERR2722069.sam ERR2722071.sam

samtools view -@ 4 output_38_sperm.bam | grep ERR2722069 >> ERR2722069.sam
samtools view -@ 4 -bS ERR2722069.sam > ERR2722069.bam
mv ERR2722069.bam /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/
rm ERR2722069.sam

samtools view -@ 4 output_38_sperm.bam | grep ERR2722071 >> ERR2722071.sam
samtools view -@ 4 -bS ERR2722071.sam > ERR2722071.bam
mv ERR2722071.bam /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/
rm ERR2722071.sam

