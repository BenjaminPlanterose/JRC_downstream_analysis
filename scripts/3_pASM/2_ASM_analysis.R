############################################################################
############################################################################
###########                                                      ###########
###########                   Pooled ASM analysis                ###########
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
library(hexbin)
library(GenomicRanges)
library(gplots)

# Chromosome name
# SNP location
# CpG location
# SNP_1/SNP_2 (SNPs for row_1/row_2 in contingency table)
# CpG_1/CpG_2 (CpGs for col_1/col_2 in contingency table)
# SNP_1-CpG_1 count (contingency table value in row_1-col_1)
# SNP_1-CpG_2 count (contingency table value in row_1-col_2)
# SNP_2-CpG_1 count (contingency table value in row_2-col_1)
# SNP_2-CpG_2 count (contingency table value in row_2-col_2)
# Two-tail Fisher’s exact p-value
# χ2 p-value

################################## Annotations ##################################

setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/blacklist_regions/')
blacklist = fread('hg38_blacklist_regions.bed')
blacklist = GRanges(seqnames = blacklist$V1, ranges = IRanges(start = blacklist$V2, end = blacklist$V3))

setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/mappability_files/')
map = fread('hg38_k100.bismap.bed', skip = 1) # skip descriptor
map = GRanges(seqnames = map$V1, ranges = IRanges(start = map$V2, end = map$V3))
map = reduce(map)

setwd('/media/ultron/2tb_disk1/ben/dbSNP/')
dbSNP = fread('dbSnp155Common.bed', nThread = 4)
dbSNP = dbSNP[,c(1:5, 7, 10, 14)]
dbSNP$V10 = as.numeric(sapply(strsplit(dbSNP$V10, split = ','), function(x) x[1]))
colnames(dbSNP) = c('chr', 'start', 'end', 'rs', 'ref', 'alt', 'maf', 'type')
#dbSNP = dbSNP[dbSNP$type == 'snv',]
dbSNP = dbSNP[dbSNP$maf >= 0.05, ]
dbSNP$SNP_ID = paste(dbSNP$chr, dbSNP$start+1L, sep = ':')
dup = duplicated(dbSNP[,1:3]); sum(dup) # 6251
dbSNP = dbSNP[!dup, ]
dbSNP = dbSNP[dbSNP$chr %in% paste('chr', c(1:22, 'X', 'Y'), sep = ''), ]
dim(dbSNP) # 8,163,050       9
sum(dbSNP$start == dbSNP$end) # 48536; reference = del
add = as.integer(dbSNP$start != dbSNP$end) # To avoid start > end
dbSNP_GR = GRanges(seqnames = dbSNP$chr, ranges = IRanges(start = dbSNP$start+add, dbSNP$end))
mean(width(dbSNP_GR) >= 1) # 1


library(BSgenome.Hsapiens.UCSC.hg38)  
chrs = names(Hsapiens)[1:24]
cgs = lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr_p1 = do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1))))
cpgr_p2 = cpgr_p1; start(cpgr_p2) = start(cpgr_p2) + 1L; end(cpgr_p2) = start(cpgr_p2)
length(cpgr_p1) # 29,401,360
length(cpgr_p2) # 29,401,360
cpginfo = paste(seqnames(cpgr_p1), start(cpgr_p1), sep = ':')

# Define CpG-SNP/indel set
CpG_SNP_annotA = subsetByOverlaps(cpgr_p1, dbSNP_GR); mean(width(CpG_SNP_annotA) == 1) # 1
cpgsnps.idCpG1 = paste(seqnames(CpG_SNP_annotA), start(CpG_SNP_annotA), sep = ':')
length(cpgsnps.idCpG1) # 612534
CpG_SNP_annotB = subsetByOverlaps(cpgr_p2, dbSNP_GR); mean(width(CpG_SNP_annotB) == 1) # 1
cpgsnps.idCpG2 = paste(seqnames(CpG_SNP_annotB), start(CpG_SNP_annotB)-1L, sep = ':')
length(cpgsnps.idCpG2) # 608643
cpgsnps.idCpG = unique(c(cpgsnps.idCpG1, cpgsnps.idCpG2)) 
length(cpgsnps.idCpG) # 1191339
612534+608643 # 1221177; not additive (CpGs with both C/D and G/H or indels)

################################## WB ##################################

# 1. Read all associations
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
WB = fread('pooled_blood_data_hg38.asm', nThread = 4)
colnames(WB) = c('chr', 'SNP_pos', 'CpG_pos', 'SNP_alleles', 'CpG_alleles', 'A1_E1', 'A1_E2', 'A2_E1', 'A2_E2', 'fisher_pvalue', 'chi2_pvalue')
WB$SNP_ID = paste(WB$chr, WB$SNP_pos, sep = ':')
WB$CpG_ID = paste(WB$chr, WB$CpG_pos, sep = ':')
dim(WB) # 53,709,032  13

# 2. Known SNP at SNP pos (removes CpG-SNPs type V) and CpG at known CpG pos
invented_variants = WB[!(WB$SNP_ID %in% dbSNP$SNP_ID),] # Or low freq
dim(invented_variants) # 49,056,058       13

quantile(-log10(invented_variants$fisher_pvalue)) # Most are non-significant
# 0%       25%      50%      75%          100% 
# 0.000000 0.000000 0.000000 0.252246      Inf 
sum(-log10(invented_variants$fisher_pvalue) >= 8) # 48900

WB = WB[WB$SNP_ID %in% dbSNP$SNP_ID,]
dim(WB) # 4,652,974      13
4652974/53709032 # 0.08663299
WB = WB[WB$CpG_ID %in% cpginfo,]
dim(WB) # 4,652,974      13

# 3. Filter CpG-SNPs/indels
WB = WB[!(WB$CpG_ID %in% cpgsnps.idCpG),]
dim(WB) # 3,943,855      13

# 4. Remove blacklist and low map
WB_GR = GRanges(seqnames = WB$chr, ranges = IRanges(start = WB$CpG_pos, end = WB$CpG_pos+1L))
length(WB_GR) # 3,943,855
WB_GR = subsetByOverlaps(WB_GR, map)
nonoverlaps = is.na(findOverlaps(WB_GR, blacklist, select="first"))
WB_GR = WB_GR[nonoverlaps]
length(WB_GR) # 3,903,362
ids = paste(as.character(seqnames(WB_GR)), start(WB_GR), sep = ':')
WB = WB[WB$CpG_ID %in% ids,]
dim(WB) # 3,903,362      13

# 5. Select significant associations
BT = 0.05/nrow(WB) # Bonferroni threshold: 1.280947e-08
WB_sig = WB[WB$fisher_pvalue < BT,]
dim(WB_sig) # 73 13

# 6. min diff of methylation
WB_sig$meth_diff = sapply(1:nrow(WB_sig), function(i) abs(WB_sig$A1_E1[i]/(WB_sig$A1_E1[i]+WB_sig$A1_E2[i]) - WB_sig$A2_E1[i]/(WB_sig$A2_E1[i]+WB_sig$A2_E2[i])))
WB_sig = WB_sig[WB_sig$meth_diff > 0.2, ]
WB_sig$odds_ratio = sapply(1:nrow(WB_sig), function(i) fisher.test(matrix(c(WB_sig$A1_E1[i], WB_sig$A1_E2[i], 
                                                                            WB_sig$A2_E1[i], WB_sig$A2_E2[i]), byrow = T, nrow = 2))$estimate)
WB_sig = WB_sig[order(WB_sig$meth_diff, decreasing = T),]
dim(WB_sig) # 68 15
head(WB_sig)

# Assign rs numbers
WB_sig$rs = dbSNP$rs[match(WB_sig$SNP_ID, dbSNP$SNP_ID)]
WB_sig$dbSNP_alleles = paste(dbSNP$ref[match(WB_sig$rs, dbSNP$rs)], dbSNP$alt[match(WB_sig$rs, dbSNP$rs)], sep = '/')
table(WB_sig$dbSNP_alleles, WB_sig$SNP_alleles) # Alleles change due to bisulfite

# 6. Export results
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
fwrite(WB_sig, 'WB_sig.txt', sep = '\t')

# 7. Export bed
df = data.frame(table(WB_sig$CpG_ID))
df$chr = sapply(strsplit(as.character(df$Var1), split = ':'), function(x) x[1])
df$pos = sapply(strsplit(as.character(df$Var1), split = ':'), function(x) x[2])
df$pos = as.integer(df$pos)
name = 'pasm_blood'
bed_df <- data.frame(seqname = df$chr,
                     start = df$pos, end = df$pos+1L,
                     name = paste(df$Var1, df$Freq, sep = '_'),
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/pASM/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''),
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4,
       sep = '\t', col.names = F, append = T)


############################### Verify by hand ###############################

library(data.table)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
WB_sig = fread('WB_sig.txt')
dim(WB_sig) # 68 17

head(WB_sig[,c(1:4, 10, 16)], 10)
# chr   SNP_pos   CpG_pos SNP_alleles fisher_pvalue         rs
# 1: chr10  12825148  12825134         T/C  1.181024e-08  rs3122548 # DpG/CpG
# 2: chr10  44358000  44357980         A/G  1.603614e-11  rs1107008 # CpH/CpG
# 3: chr13  25096713  25096675         A/G  2.854341e-09 rs75475407 # CpG-SNP rs3002212 (T=0.0431)
# 4: chr20  53864295  53864279         T/C  1.275853e-11 rs78093244 # CpG-SNP CpG/TpG
# 5: chr20  61920071  61920080         A/G  8.350044e-09 rs79594189 # CpG/CpH
# 6: chr21  10381114  10381078         A/C  5.678465e-13 rs79874033 # Probably OK
# 7: chr21  10381114  10381126         A/C  1.024330e-13 rs79874033 # Probably OK
# 8:  chr3  75694316  75694362         T/G  1.041052e-12 rs74281629 # Looks good!
# 9:  chr4 172546807 172546819         T/C  2.121180e-09 rs60860516 # CpG/TpG
# 10: chr5  39088559  39088543         A/G  3.078872e-10  rs1119493 # CpH/CpG



################################## Sperm ##################################

# 1. Read all associations
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
sperm = fread('pooled_sperm_data_hg38.asm', nThread = 4)
colnames(sperm) = c('chr', 'SNP_pos', 'CpG_pos', 'SNP_alleles', 'CpG_alleles', 'A1_E1', 'A1_E2', 'A2_E1', 'A2_E2', 'fisher_pvalue', 'chi2_pvalue')
sperm$SNP_ID = paste(sperm$chr, sperm$SNP_pos, sep = ':')
sperm$CpG_ID = paste(sperm$chr, sperm$CpG_pos, sep = ':')
dim(sperm) # 38,299,379  13

# 2. Known SNP at SNP pos
invented_variants = sperm[!(sperm$SNP_ID %in% dbSNP$SNP_ID),] # Or low freq
dim(invented_variants) # 34,661,536       13

quantile(-log10(invented_variants$fisher_pvalue)) # Most are non-significant
# 0%  25%  50%  75% 100% 
# 0    0    0    0  Inf 
sum(-log10(invented_variants$fisher_pvalue) >= 8) # 23890

sperm = sperm[sperm$SNP_ID %in% dbSNP$SNP_ID,]
dim(sperm) # 3,637,843      13
3637843/38299379 # 0.09498439
sperm = sperm[sperm$CpG_ID %in% cpginfo,]
dim(sperm) # 3,637,843      13

# 3. Filter CpG-SNPs/indels
sperm = sperm[!(sperm$CpG_ID %in% cpgsnps.idCpG),]
dim(sperm) # 2,986,335      13

# 4. Remove blacklist and low map
sperm_GR = GRanges(seqnames = sperm$chr, ranges = IRanges(start = sperm$CpG_pos, end = sperm$CpG_pos+1L))
length(sperm_GR) # 2986335
sperm_GR = subsetByOverlaps(sperm_GR, map)
nonoverlaps = is.na(findOverlaps(sperm_GR, blacklist, select="first"))
sperm_GR = sperm_GR[nonoverlaps]
length(sperm_GR) # 2949944
ids = paste(as.character(seqnames(sperm_GR)), start(sperm_GR), sep = ':')
sperm = sperm[sperm$CpG_ID %in% ids,]
dim(sperm) # 2,949,944      13

# 5. Select significant associations
BT = 0.05/nrow(sperm) # Bonferroni threshold: 1.694947e-08
sperm_sig = sperm[sperm$fisher_pvalue < BT,]
dim(sperm_sig) # 58 13

# 6. min diff of methylation
sperm_sig$meth_diff = sapply(1:nrow(sperm_sig), function(i) abs(sperm_sig$A1_E1[i]/(sperm_sig$A1_E1[i]+sperm_sig$A1_E2[i]) - sperm_sig$A2_E1[i]/(sperm_sig$A2_E1[i]+sperm_sig$A2_E2[i])))
sperm_sig = sperm_sig[sperm_sig$meth_diff > 0.2, ]
sperm_sig$odds_ratio = sapply(1:nrow(sperm_sig), function(i) fisher.test(matrix(c(sperm_sig$A1_E1[i], sperm_sig$A1_E2[i], 
                                                                                  sperm_sig$A2_E1[i], sperm_sig$A2_E2[i]), byrow = T, nrow = 2))$estimate)
sperm_sig = sperm_sig[order(sperm_sig$meth_diff, decreasing = T),]
dim(sperm_sig) # 56 15
# Assign rs numbers
sperm_sig$rs = dbSNP$rs[match(sperm_sig$SNP_ID, dbSNP$SNP_ID)]
sperm_sig$dbSNP_alleles = paste(dbSNP$ref[match(sperm_sig$rs, dbSNP$rs)], dbSNP$alt[match(sperm_sig$rs, dbSNP$rs)], sep = '/')

table(sperm_sig$dbSNP_alleles, sperm_sig$SNP_alleles) # Alleles change due to bisulfite

# 7. Export results
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
fwrite(sperm_sig, 'sperm_sig.txt', sep = '\t')

# 8. Export bed
df = data.frame(table(sperm_sig$CpG_ID))
df$chr = sapply(strsplit(as.character(df$Var1), split = ':'), function(x) x[1])
df$pos = sapply(strsplit(as.character(df$Var1), split = ':'), function(x) x[2])
df$pos = as.integer(df$pos)
name = 'pasm_sperm'
bed_df <- data.frame(seqname = df$chr,
                     start = df$pos, end = df$pos+1L,
                     name = paste(df$Var1, df$Freq, sep = '_'),
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/pASM/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''),
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4,
       sep = '\t', col.names = F, append = T)


############################### Verify by hand ###############################

library(data.table)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
sperm_sig = fread('sperm_sig.txt')
dim(sperm_sig) # 56 17

head(sperm_sig[,c(1:4, 10, 16)], 10)
# chr   SNP_pos   CpG_pos SNP_alleles fisher_pvalue         rs
# 1: chr10  12825148  12825134         T/C  3.495309e-09  rs3122548 # TpG/CpG, probably ok
# 2: chr12     60040     60057         A/T  7.086141e-09 rs57449118 # Looks good!
# 3: chr12     60041     60057         T/C  3.693478e-10 rs60417285 # Looks good!
# 4: chr15  97594529  97594521         A/G  7.750467e-09 rs36107004 # Looks good!
# 5: chr16  33623550  33623585         T/C  7.086141e-09  rs2625425 # CpG/TpG
# 6:  chr3  75694084  75694094         A/T  2.217403e-09 rs80109768 # Looks good!
# 7:  chr3  75694316  75694362         T/G  8.310326e-11 rs74281629 # Looks good!
# 8:  chr3  99501623  99501610         T/C  1.590040e-08  rs1513299 # Suspicious, probably ok
# 9:  chr4 189973491 189973474         T/C  1.111548e-10 rs10024866 # CpG/TpG
# 10: chr5   3016563   3016594         A/C  1.931223e-10  rs2961779 # SNP allele should be A/T instead of A/C due to BC





############################### Verify by hand ###############################

dim(WB_sig) # 68 17
dim(sperm_sig) # 56 17
length(unique(WB_sig$CpG_ID)) # 62
length(unique(sperm_sig$CpG_ID)) # 54
length(unique(WB_sig$SNP_ID)) # 61
length(unique(sperm_sig$SNP_ID)) # 53

# 1. Common between WB_sig and sperm_sig
venn = venn(list(SP = sperm_sig$SNP_ID, WB = WB_sig$SNP_ID))
attr(venn, 'intersection')$`SP:WB` # 21

############################### Follow-up ###############################

## Genetic artifacts

# (Type-IA: not possible with Biscuit)
# type-IB

jkl = WB[WB$CpG_ID %in% cpgsnps.idCpG,]
jkl = jkl[jkl$fisher_pvalue < 1e-08,]
jkl = jkl[jkl$SNP_pos != jkl$CpG_pos +1,]
head(jkl)

# chr   SNP_pos   CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue   chi2_pvalue          SNP_ID
# rs4660216 5: chr1 42126046 42126024         A/G         T/C     0    21    13     0  1.077605e-09
jkl[800:810,]


WB = WB[order(WB$fisher_pvalue), ]
head(WB, 50)

(WB[WB$SNP_pos != WB$CpG_pos+1,])[51:100,]

WB[110:130,]

CT = WB[WB$SNP_alleles == 'T/C',]
CT[341:360,]

# type-IIA
# chr   SNP_pos   CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue   chi2_pvalue          SNP_ID
# chr4  49092187  49092186         A/G         T/C   773     0     0    59  5.882620e-92 2.155239e-181   chr4:49092187

# type-IIB
# chr   SNP_pos   CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue   chi2_pvalue          SNP_ID
# chr1  73113554  73113571         A/G         T/C    21     0     0    32  3.144788e-15


jkl = WB[WB$fisher_pvalue < 1e-08,]
jkl = jkl[jkl$SNP_pos != jkl$CpG_pos +1,]
jkl = jkl[order(jkl$fisher_pvalue),]

jkl[120:130,]


# P1 (WB)
#       chr  SNP_pos  CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue  chi2_pvalue         SNP_ID
# 1: chr10 12825148 12825134         T/C         T/C    11     0     0    20  1.181024e-08 1.855391e-07 chr10:12825148
#             CpG_ID meth_diff odds_ratio        rs dbSNP_alleles
# 1: chr10:12825134         1        Inf rs3122548        T/C,G,
# (P2A)
# P2B (WB)

WB_sig[4, ] # INDEL: rs3046367
# chr SNP_pos CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue  chi2_pvalue        SNP_ID
# 1: chr12 9243167 9243178         T/C         T/C     0    13    36     0   3.80812e-12 2.289735e-11 chr12:9243167
# CpG_ID meth_diff odds_ratio       rs dbSNP_alleles
# 1: chr12:9243178         1          0 rs717181          T/C,



# samtools view -b output_38.bam "chr5:27403645-27404645" > typeIB_WB.bam
# samtools index typeIB_WB.bam
# samtools view -b output_38.bam "chr4:49091687-49092687" > typeIIA_WB.bam
# samtools index typeIIA_WB.bam
# samtools view -b output_38.bam "chr3:70468044-70469044" > typeIIB_WB.bam
# samtools index typeIIB_WB.bam
# samtools view -b output_38.bam "chr10:12824648-12825648" > P1_WB.bam
# samtools index P1_WB.bam
# samtools view -b output_38.bam "chr12:9242667-9243667" > P2B_WB.bam
# samtools index P2B_WB.bam


# TP
# chr   SNP_pos   CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue  chi2_pvalue         SNP_ID
# 1: chr5 170977913 170977961         A/T         T/C    19     0     0    12  7.086141e-09 1.855391e-07 chr5:170977913
# CpG_ID meth_diff odds_ratio         rs dbSNP_alleles
# 1: chr5:170977961         1        Inf rs10042357        T/A,G,

# chr5:170977961

# samtools view -b output_38.bam "chr5:170977413-170978413" > TP_WB.bam
# samtools index TP_WB.bam


# SP
# chr  SNP_pos  CpG_pos SNP_alleles CpG_alleles A1_E1 A1_E2 A2_E1 A2_E2 fisher_pvalue  chi2_pvalue        SNP_ID
# 1: chr3 75694084 75694094         A/T         T/C     0    41     8     0  2.217403e-09 2.289735e-11 chr3:75694084
# 2: chr3 75694316 75694362         T/G         T/C     0    13    27     0  8.310326e-11 2.061154e-09 chr3:75694316
# CpG_ID meth_diff odds_ratio         rs dbSNP_alleles
# 1: chr3:75694094         1          0 rs80109768          A/T,
# 2: chr3:75694362         1          0 rs74281629        C/G,T,

# samtools view -b output_38_sperm.bam "chr3:75693584-75694584" > TP_SP.bam
# samtools index TP_SP.bam


#############################################################################

library(data.table)
library(GenomicRanges)

datatable2grange = function(mat)
{
  GRanges(mat$chr, IRanges(start = mat$start, end = mat$end))
}

coordinates2datatable = function(coordinates)
{
  chr = sapply(strsplit(coordinates, split = ':'), function(x) x[1])
  pos = sapply(strsplit(coordinates, split = ':'), function(x) x[2])
  start = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[1]))
  end = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[2]))
  data.table(chr = chr, start = start, end = end, range = end - start)
}

setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/')
WB = fread('p_values.txt')
WB = WB[!is.na(WB$V2),]; dim(WB) # 630019      2
BT = 0.05/nrow(WB)
sig_WB = WB$V1[WB$V2 < BT]; length(sig_WB) # 4215
A = coordinates2datatable(sig_WB); dim(A) # 4215    4
A = datatable2grange(A)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binokulars_output/sperm_results/')
SP = fread('p_values.txt')
SP = SP[!is.na(SP$V2),]; dim(SP) # 510092      2
BT = 0.05/nrow(SP)
sig_SP = SP$V1[SP$V2 < BT]; length(sig_SP) # 7386
B = coordinates2datatable(sig_SP); dim(B) # 7386    4
B = datatable2grange(B)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
ASM_WB = fread('WB_sig.txt')
dim(ASM_WB) # 68  17
C = coordinates2datatable(unique(ASM_WB$CpG_ID))
C$end = C$start + 1
C = datatable2grange(C)
length(C) # 62

setwd('/media/ultron/2tb_disk1/ben/JRC_project/pASM_output/')
ASM_SP = fread('sperm_sig.txt')
dim(ASM_SP) # 56 17
D = coordinates2datatable(unique(ASM_SP$CpG_ID))
D$end = D$start + 1
D = datatable2grange(D)
length(D) # 54

length(intersect(A, C)) # 7
length(intersect(B, D)) # 12


length(intersect(A, C))/length(C) # 0.1129032
length(intersect(B, D))/length(D) # 0.2222222
