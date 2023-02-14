#!bin/bash

# ERR2722068: WB_old
# ERR2722069: sperm_old
# ERR2722070: WB_young
# ERR2722071: sperm_young

#ERR2722068.bam
#ERR2722070.bam
#ERR2722069.bam
#ERR2722071.bam

cd /media/ultron/2tb_disk1/ben/JRC_project/bam_files/hg38/split/
bams=($( ls *.bam -1v))
l=${#bams[@]}
for ((i=0;i<=$l-1;i+=1))
do
 bami=${bams[$i]}
 echo $bami
 samtools index -@ 4 $bami
done





