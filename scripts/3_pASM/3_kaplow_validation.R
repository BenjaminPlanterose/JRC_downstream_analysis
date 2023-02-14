############################################################################
############################################################################
###########                                                      ###########
###########               Comparison with Kaplow et al           ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########           b.planterosejimenez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################
library(data.table)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/kaplow_2015/')
kaplow = fread('kaplow_mQTL.txt') # threshold = 0.001; precision = 1e-06; hg19
kaplow = kaplow[order(kaplow$`p-Value`),]
min(kaplow$`p-Value`[kaplow$`p-Value` != 0]) # 1e-06
dim(kaplow)

kaplow[which(kaplow$`p-Value` == 4e-05),] # rs6584097, rs6599131, rs10815832
kaplow[which(kaplow$`p-Value` == 5e-05),] # rs10777864, 


fisher.test(matrix(c(13,2,5,22), nrow = 2))



kaplow$`p-Value` = kaplow$`p-Value` + 1e-06
plot(kaplow$`Position of Variant`-kaplow$`Position of CpG`, -log10(kaplow$`p-Value`), pch = 19)

table(kaplow$`Chromosome of Variant`== kaplow$`Chromosome of CpG`) # No inter-chromosomal mQTL

head(kaplow)
name= 'kaplow_hg19'
bed_df <- data.frame(seqname = kaplow$`Chromosome of CpG`,
                     start = as.integer(kaplow$`Position of CpG`), end = as.integer(kaplow$`Position of CpG`+1),
                     name = paste(kaplow$`Chromosome of CpG`, paste(kaplow$`Position of CpG`, kaplow$`Position of Variant`, sep = '-'), sep = ':'),
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/pASM/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)

kaplow_hg38 = fread('tmp_kaplow_liftovertohg38.bed', header = F)
kaplow_hg38$chr = sapply(strsplit(kaplow_hg38$V1, split = ':'), function(x) x[1])
kaplow_hg38$partB = sapply(strsplit(kaplow_hg38$V1, split = ':'), function(x) x[2])
kaplow_hg38$start = sapply(strsplit(kaplow_hg38$partB, split = '-'), function(x) x[1])
kaplow_hg38$end = sapply(strsplit(kaplow_hg38$partB, split = '-'), function(x) x[2])


head(kaplow)
name= 'kaplow_hg38'
bed_df <- data.frame(seqname = kaplow_hg38$chr,
                     start = as.integer(kaplow_hg38$start)-1L, end = as.integer(kaplow_hg38$start),
                     name = paste(kaplow_hg38$chr, paste(as.integer(kaplow_hg38$start)-1L, as.integer(kaplow_hg38$start)-kaplow$`Position of CpG`+kaplow$`Position of Variant`-1L, sep = '-'), kaplow$`rsID of Variant`, sep = ':'),
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/pASM/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)

# (1) rs7366554 has extremely low MAF
# (2) There is something. Not sure if artifact
# (3) INDEL, let us begin with easy
# (4) rs60205880


bed_df[5,]
# Prep for IGV
c(181945159-500L, 181945159+500L)
# samtools view -b output_38.bam "chr1:181944659-181945659" > chr1_181944659_181945659.bam
# samtools index chr1_181944659_181945659.bam
# chr1:181945159-181945152:rs2262513

head(bed_df)

selected = c('rs10737680', 'rs9517668', 'rs12905925', 'rs617201', 'rs1113144', 
             'rs4809456', 'rs62223713', 'rs7705033')
bed_df[kaplow$`rsID of Variant` %in% selected,]
# seqname     start       end                                name score
# 199    chr15  64872654  64872655  chr15:64872654-64872646:rs12905925  1000
# 314    chr18  74107552  74107553   chr18:74107552-74107551:rs1113144  1000
# 481    chr17  32566127  32566128    chr17:32566127-32566117:rs617201  1000
# 524     chr1 196710327 196710328 chr1:196710327-196710325:rs10737680  1000
# 543     chr5 123439101 123439102  chr5:123439101-123439091:rs7705033  1000
# 840    chr20  63029518  63029519   chr20:63029518-63029497:rs4809456  1000
# 1033   chr21  28785724  28785725  chr21:28785724-28785719:rs62223713  1000
# 1100   chr15  64872623  64872624  chr15:64872623-64872646:rs12905925  1000
# 2055   chr13  99271535  99271536   chr13:99271535-99271586:rs9517668  1000

kaplow$`p-Value`[kaplow$`rsID of Variant` %in% selected]
# 0.000006 0.000015 0.000036 0.000044 0.000048 0.000126 0.000200 0.000224 0.000771
plot(density(kaplow$`p-Value`)); abline(v = kaplow$`p-Value`[kaplow$`rsID of Variant` %in% selected], lty = 2)

# 
# Next, we converted all HapMap Phase II (Frazer et al. 2007) and 1000 Genomes Phase I Integrated Version 
# 3 (The 1000 Genomes Project Consortium 2010) single nucleotide polymorphisms (SNPs) with MAF > 0.04 in 
# human genome version hg19 (The International Human Genome Sequencing Consortium 2001) into N’s in order 
# to eliminate sources of reference bias when mapping (Degner et al. 2009)
# The CpGs with called methylation statuses should not contain most SNPs because we masked SNPs with Ns, 
# so such CpGs would have become CNs.

