library(data.table)
library(GEOquery)
library(minfi)
library(UMtools)
library(gplots)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ENmix)

annotation = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Export bed
head(annotation)
name = 'EPIC_hg19'
bed_df <- data.frame(seqname = annotation$chr,
                     start = as.integer(annotation$pos), end = as.integer(annotation$pos)+1L,
                     name = rownames(annotation),
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/JRCs/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)



# JRC-SP
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binokulars_output/sperm_results/')
SP = fread('p_values.txt')
dim(SP) # 536691      2
SP = SP[order(SP$V2),]

#0.05/nrow(SP)

# chr8:143921000-143926800
SP[91:100,]

# chr21:10415086-10423761, BAGE
# chr1:4302200-4303003, STR?
# chrX:97843800-97848800
# chr8:57214200-57219400, LINC01606; could be false positive
# chr5:118971800-118974264, DTWD2; could be false positive
# chr15:65075820-65077400, KBTBD13; could be false positive
# chr2:158310288-158312400, CCDC148; could be false positive

## REAL
# chr22:19724400-19728600, GP1BB (age)
# chr18:8783000-8787400, MTCL1
# chr19:23415800-23418952, LINC01224
# chr17:1468000-1474000, MYO1C
# chr17:76007200-76010400, EVPL
# chr7:780400-785000, DNAAF5
# chr2:241550400-241554600, BOK-AS1
# chr13:52425800-52427416, VPS36
# chr19:2715600-2720000, DIRAS1
# chr1:16630166-16632400, CROCCP2
# chr7:157464400-157468000, NR_110157
# chr11:134138223-134166129, SNORD153
# chr5:661400-670000, TPPP
# chr1:1337200-1342200, TAS1R3
# chr21:43361200-43363600, LINC01679
# chr6:167961400-167967200, AFDN
# chr1:53091200-53103600, SLC1A7
# chr19:2673200-2680200, GNG7
# chr21:41834000-41837200, PRDM15
# chr7:157801000-157803200, PTPRN2
# chr2:2187600-2191000, MYT1L
# chr10:132963800-132966800, LINC01166, LINC01168
# chr4:55600-69400, ZNF595
# chr16:632200-634200, WFIKKN1
# chr10:15211200-15215000, FAM171A1
# chr4:1260600-1263400, CTBP1-AS2








# setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/sperm_EPIC/GSE185445_RAW/')
# rgSet = read.metharray.exp(getwd(), extended = F)
# beta = preprocessQuantile(rgSet)
# beta = getBeta(beta)
# raw = preprocessRaw(rgSet)
# M = getMeth(raw)
# U = getUnmeth(raw)
# 
# setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/sperm_EPIC/PROCESSED/')
# UMtools::export_bigmat(M, 'M.txt', 4)
# UMtools::export_bigmat(U, 'U.txt', 4)
# UMtools::export_bigmat(round(beta, 5), 'beta.txt', 4)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/sperm_EPIC/PROCESSED/')
M = UMtools::import_bigmat('2022-10-12_M.txt')
U = UMtools::import_bigmat('2022-10-12_U.txt')
beta = UMtools::import_bigmat('2022-10-12_beta.txt')

######## GEO
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/sperm_EPIC/')
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
phenotype <- getGEO('GSE185445', destdir=".")
pheno <- Reduce(rbind, lapply(phenotype, function(x) pData(phenoData(x))))
ID = sapply(strsplit(colnames(M), split = '_'), function(x) x[1])
pheno = pheno[ID,]
dim(pheno) # 379  65

####### Col
breaks=seq(-0.05, 1.05, 0.05)
my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)




#
'cg07172676' %in% rownames(beta)

UM_plot(M, U, 'cg07172676')
res = Visualize_cometh(beta_mat = beta, annotation = annotation, CpG = 'cg07172676', distance = 1000, L_bound = 5, R_bound = 4)

mdsPlot(beta[res,-c(1, length(res))])


heatmap.2(beta[res,], trace = 'n', Rowv = 'n', dendrogram = 'col')

UM_plot(M, U, 'cg04247530')
# CpG_rs   CpG_maf 
# rs62641755  0.014153  
annotation['cg04247530',]$SourceSeq
dim(beta)

UM_plot(M, U, 'cg19405177')



UM_plot(M, U, 'cg12044749')

UM_plot(M, U, 'cg01949498')
UM_plot(M, U, 'cg09143319')

UM_plot(M, U, 'cg26120132')
UM_plot(M, U, 'cg03354047')
UM_plot(M, U, 'cg00017203')


plec = c('cg01949498', 'cg09143319', 'cg26120132', 'cg03354047', 'cg00017203')


heatmap.2(beta[plec,], trace = 'n', Rowv = 'n', dendrogram = 'col', na.color = 'lightcoral', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          cexCol = 0.8)


plot(ecdf(-beta['cg04247530', ]), xlim = c(-1,0))

plot(beta[plec[1], ], beta['cg04247530', ], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2, col = 'red2')
plot(beta[plec[1], ], beta[plec[2], ], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2, col = 'red2')
plot(beta[plec[1], ], beta[plec[3], ], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2, col = 'red2')
plot(beta[plec[1], ], beta[plec[4], ], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2, col = 'red2')
plot(beta[plec[1], ], beta[plec[5], ], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2, col = 'red2')

which.min(beta[plec[1],])
pheno['GSM5615114',]


# cg04247530
col = UMtools::bGMM(M, U, 'cg04247530', 2, transform = T)
plot(beta[plec[1], ], beta['cg04247530', ], col  = col)
plot(density(beta['cg04247530',]))
plot(density(beta[plec[1],]))
plot(density(beta[plec[2],]))
plot(density(beta[plec[3],]))
plot(density(beta[plec[4],]))
plot(density(beta[plec[5],]))




table(col)
# blue  red 
# 106  273 
(106+273/2)/379 # 0.6398417

col = UMtools::bGMM(M, U, plec[1], 2, transform = F)
plot(beta[plec[1], ], beta['cg04247530', ], col  = col)
table(col)


# Good ageing CpGs: cg21064553, cg08663040, cg25093603, cg21843517
plot(beta['cg21064553',], beta['cg08663040',], xlim = c(0,1), ylim = c(0,1))
plot(beta['cg25093603',], beta['cg21843517',], xlim = c(0,1), ylim = c(0,1))
age = colMeans(beta[c('cg21064553', 'cg08663040',  'cg25093603', 'cg21843517'),])

plot(age, beta['cg04247530', ])

plot(pheno$`age_m:ch1`, beta['cg21064553',], ylim = c(0,1))
plot(pheno$`age_m:ch1`, beta['cg08663040',], ylim = c(0,1))
plot(pheno$`age_m:ch1`, beta['cg25093603',], ylim = c(0,1))
plot(pheno$`age_m:ch1`, beta['cg21843517',], ylim = c(0,1))

#  Smoking
res = Visualize_cometh(beta_mat = beta, annotation = annotation, CpG = 'cg05575921', distance = 1000, L_bound = 5, R_bound = 4)

plot(as.numeric(pheno$`age_m:ch1`), as.numeric(pheno$`cotinine_m:ch1`), ylim = c(0,1))

plot(pheno$`cotinine_m:ch1`, beta['cg05575921',], ylim = c(0,1))
plot(pheno$`cotinine_m:ch1`, beta['cg23576855',], ylim = c(0,1))


plot(pheno$`cotinine_m:ch1`, beta['cg04247530',])


head(pheno)

table(pheno$source_name_ch1) # sperm: 379
table(pheno$characteristics_ch1.1) # normal: 379


# age_f:ch1 age_m:ch1 bmi_f:ch1   bmi_m:ch1 cotinine_f:ch1 cotinine_m:ch1 sea_cpg:ch1 sea_dmr:ch1 status:ch1 ttp:ch1


plot(beta['cg04247530', ], pheno$`age_f:ch1`)
plot(beta['cg04247530', ], pheno$`age_m:ch1`)
plot(beta['cg04247530', ], pheno$`cotinine_f:ch1`)
plot(beta['cg04247530', ], pheno$`cotinine_m:ch1`)
plot(beta['cg04247530', ], pheno$`status:ch1`)


# characteristics_ch1.5
# characteristics_ch1.6
# characteristics_ch1.7
# characteristics_ch1.8

# - Epithelial contamination? There should be more CpGs correlated
# not somatic contamination; if would be obvious at imprinted genes (Laurentino)
# To verify our methylation data were not influenced by somatic contamination, we analyzed methylation at a maternally imprinted gene, 
# DLK1, previously shown (Jenkins et al., 2018) to be differentially methylated between sperm and somatic tissues as well as the 
# paternally imprinted locus, H19 (Supplementary Fig. S1). At the DLK1 and H19 loci, average methylation for all participants was 
# below 5%, and above 90%, respectively, suggesting negligible somatic cell contamination.

# - predict age; see if correlated. Maybe it is a missed age-DMR
# Not in  All chronological age-associated bonferroni signficant CpGs

# sperm cellular composition or internal parameter?
# 


plot(beta['cg04247530', ], beta['cg17412258',], xlim = c(0,1), ylim = c(0,1))



cor_vec = sapply(1:nrow(beta), function(x) cor(beta['cg04247530',], beta[x,]))
names(cor_vec) = rownames(beta)
cor_vec = cor_vec[order(cor_vec^2, decreasing = T)]

head(cor_vec, 20)



####################################################################################


## REAL
# chr22:19724400-19728600, GP1BB (age)
# chr18:8783000-8787400, MTCL1
# chr19:23415800-23418952, LINC01224
# chr17:1468000-1474000, MYO1C
# chr17:76007200-76010400, EVPL
# chr7:780400-785000, DNAAF5
UM_plot(M, U, 'cg22213837')
# chr2:241550400-241554600, BOK-AS1
UM_plot(M, U, 'cg23940744')
# chr13:52425800-52427416, VPS36
# chr19:2715600-2720000, DIRAS1
UM_plot(M, U, 'cg15694117')
# chr1:16630166-16632400, CROCCP2
UM_plot(M, U, 'cg01926717')
# chr7:157464400-157468000, NR_110157
UM_plot(M, U, 'cg16329658')
# chr11:134138223-134166129, SNORD153
UM_plot(M, U, 'cg03640071')
# chr5:661400-670000, TPPP
UM_plot(M, U, 'cg06219726')
# chr1:1337200-1342200, TAS1R3
UM_plot(M, U, 'cg13494355')
# chr21:43361200-43363600, LINC01679
UM_plot(M, U, 'cg10251070')
# chr6:167961400-167967200, AFDN; mQTL?
# chr1:53091200-53103600, SLC1A7
UM_plot(M, U, 'cg15526153')
# chr19:2673200-2680200, GNG7
# chr21:41834000-41837200, PRDM15
UM_plot(M, U, 'cg21063361')
# chr7:157801000-157803200, PTPRN2
UM_plot(M, U, 'cg20194947')
# chr2:2187600-2191000, MYT1L; mQTL
UM_plot(M, U, 'cg21069434')
# chr10:132963800-132966800, LINC01166, LINC01168
UM_plot(M, U, 'cg02276826')
# chr4:55600-69400, ZNF595
# chr16:632200-634200, WFIKKN1
UM_plot(M, U, 'cg06520088')
# chr10:15211200-15215000, FAM171A1
UM_plot(M, U, 'cg20534130')
# chr4:1260600-1263400, CTBP1-AS2
UM_plot(M, U, 'cg06466757')
UM_plot(M, U, 'cg06232075')


# PLEC
UM_plot(M, U, 'cg04247530')
plot(density(beta['cg04247530', ]))
plot(pheno$`age_m:ch1`, beta['cg04247530',], ylim = c(0,1))

# chr17:1468000-1474000, MYO1C
UM_plot(M, U, 'cg27546431')
plot(pheno$`age_m:ch1`, beta['cg27546431',], ylim = c(0,1))
plot(pheno$`age_m:ch1`, beta['cg26183265',], ylim = c(0,1))


# chr17:76007200-76010400, EVPL
UM_plot(M, U, 'cg20143393')
plot(pheno$`age_m:ch1`, beta['cg00714713',], ylim = c(0,1))
plot(pheno$`age_m:ch1`, beta['cg20143393',], ylim = c(0,1))




UM_plot(M, U, 'cg25649515')


maf = annotation[rownames(beta), c('CpG_maf', 'SBE_maf')]
pos = which(is.na(maf$CpG_maf) & is.na(maf$SBE_maf))
var = rowVars(beta[pos,])
names(var) = rownames(beta[pos,])
var = sort(var, decreasing = T)
i = sample(names(var[1:10000]), 1); UM_plot(M, U, i)
annotation[i,]
#



### umplot of the year cg19563671
### cg12292084, CR power; L2

#####################################################################################################

# Strange bimodal

# PLEC
UM_plot(M, U, 'cg04247530')
y = beta['cg04247530',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg04247530', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871

# DNAH1
UM_plot(M, U, 'cg26060971')
UM_plot(M, U, 'cg06068388')
y = beta['cg06068388',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg06068388', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871

# DNAAF5, chr7:780400-785000
UM_plot(M, U, 'cg18344811')
UM_plot(M, U, 'cg03275910')
UM_plot(M, U, 'cg16782493')

y = beta['cg18344811',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg18344811', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871
y = beta['cg03275910',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg03275910', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871
y = beta['cg16782493',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg16782493', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871



#####################################################################################################

###### mQTL in sperm
breaks=seq(-0.05, 1.05, 0.05)
my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)


# (1) chr1:12744911-12746800, C1orf158, mQTL; 1.152161e-09
i = 'chr1:12744911-12746800'
UM_plot(M, U, 'cg17878878')
res = Visualize_cometh(annotation = annotation, CpG = 'cg17878878', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 2)
markers = c("cg00766918", "cg00214548", "cg17878878")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(3,2,1)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 30       168       181 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))


# (2) chr15:101052200-101054800, LRRK1, mQTL; 1.96182e-18
i = 'chr15:101052200-101054800'
UM_plot(M, U, 'cg02262866')
res = Visualize_cometh(annotation = annotation, CpG = 'cg02262866', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 2)

markers = c("cg19786805", "cg19405229", "cg02262866", "cg14070203")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(2,1,3)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 132       201       46 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
 

# (3) chr16:49767200-49769200, ZNF423, mQTL; 4.702179e-10
i = 'chr16:49767200-49769200'
UM_plot(M, U, 'cg26900686')
res = Visualize_cometh(annotation = annotation, CpG = 'cg26900686', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 2)

markers = c("cg26900686", "cg27525020", "cg16383503")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(2,3,1)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 135       180       64 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))


# (4) chr5:10289600-10291555, CMBL, mQTL; 1.570009e-14
i = 'chr5:10289600-10291555'
UM_plot(M, U, 'cg18328428')
res = Visualize_cometh(annotation = annotation, CpG = 'cg18328428', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 2)

markers = c("cg18328428", "cg00123333", "cg04990989")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(2,1,3)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 250       118       11 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))




# (5) chr4:1206200-1210000, SPON2, mQTL; 5.198004e-35
i = 'chr4:1206200-1210000'
UM_plot(M, U, 'cg26130533')
UM_plot(M, U, 'cg26910825')
res = Visualize_cometh(annotation = annotation, CpG = 'cg26130533', distance = 1000, beta_mat = beta, 
                       L_bound = 3, R_bound = 2)

markers = c("cg27660248", "cg20756245", "cg23045716", "cg19056444", "cg26910825", "cg09865015", "cg05436939", "cg14527262", "cg04228083",
            "cg17227257", "cg16721321", "cg26130533", "cg14505741", "cg11104416", "cg18085660", "cg11888738")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(3,2,1)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 116       149       114 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))




# (6) chr4:7508800-7512400, SORCS2, mQTL; 4.079964e-09
i = 'chr4:7508800-7512400'
UM_plot(M, U, 'cg10703338')
res = Visualize_cometh(annotation = annotation, CpG = 'cg10703338', distance = 1000, beta_mat = beta, 
                       L_bound = 3, R_bound = 2)
markers = c("cg08554860", "cg03166265", "cg21013745", "cg10703338", "cg16133398")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(3,2,1)]
assign = as.character(assign)
table(assign)
# blue     purple    blue
# 156       141       82 
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))




# chr19:55531000-55534200, SBK2, mQTL; 6.985929e-08
# UM_plot(M, U, 'cg23226510')
# lab = bGMM(M = M, U = U, CpG = 'cg23226510', K = 3, transform = T)
# freq = table(lab); (freq['red'] + freq['blue']/2)/sum(freq) # 0.8139842
# 
# i = 'chr19:55531000-55534200'
# UM_plot(M, U, 'cg23226510')
# UM_plot(M, U, 'cg25161252')
# res = Visualize_cometh(annotation = annotation, CpG = 'cg23226510', distance = 1000, beta_mat = beta, 
#                        L_bound = 3, R_bound = 2)
# markers = c("cg23226510", "cg25161252", "cg00844791", "cg18391611")
# assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
# levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(3,2,1)]
# assign = as.character(assign)
# table(assign)
# # blue     purple    blue
# # 156       141       82 
# heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
#           na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
#           offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
#           RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
# 




################## DNAH1, chr3:52372400-52374400 (cannot validate age-dependence)

i = 'chr3:52372400-52374400'
UM_plot(M, U, 'cg06068388')
res = Visualize_cometh(annotation = annotation, CpG = 'cg06068388', distance = 1000, beta_mat = beta, 
                       L_bound = 4, R_bound = 4)
markers = c("cg06068388", "cg06637517", "cg26060971")
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5))

plot(as.numeric(pheno$`age_m:ch1`), colMeans(beta[markers,]))


# Positive control (DNAH1)
plot(pheno$`age_m:ch1`, beta['cg06068388',], ylim = c(0,1))
cor.test(as.numeric(pheno$`age_m:ch1`), beta['cg06068388',]) # cor = 0.1572734; p-value = 0.002135


################## Age-associated

# LMNB2
i = 'chr19:2429200-2431800'
# samtools view -b output_38_sperm.bam "chr19:2429200-2431800" > LMNB2_SP.bam
# samtools index LMNB2_SP.bam

UM_plot(M, U, 'cg11964578')
plot(pheno$`age_m:ch1`, beta['cg11964578',], ylim = c(0,1))
cor.test(as.numeric(pheno$`age_m:ch1`), beta['cg11964578',]) # cor = -0.2216836; p-value = 1.326e-05

res = Visualize_cometh(annotation = annotation, CpG = 'cg11964578', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 2)
markers = c("cg17151731", "cg11964578", "cg18610079", "cg11823328", "cg06131173",
            "cg04989970", "cg00775818")
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5))

lapply(1:length(markers), function(i) cor.test(as.numeric(pheno$`age_m:ch1`), beta[markers[i],]))
# [[1]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -4.875, df = 377, p-value = 1.606e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3360108 -0.1463749
# sample estimates:
#        cor 
# -0.2435188 
# [[2]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -4.4141, df = 377, p-value = 1.326e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3153755 -0.1237116
# sample estimates:
#   cor 
# -0.2216836 
# [[3]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -3.3775, df = 377, p-value = 0.0008076
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.26749543 -0.07188522
# sample estimates:
#        cor 
# -0.1713788 
#
# [[4]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -4.2202, df = 377, p-value = 3.061e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3065682 -0.1140993
# sample estimates:
#   cor 
# -0.2123927 
# 
# [[5]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -3.6582, df = 377, p-value = 0.0002902
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.28064811 -0.08601751
# sample estimates:
#        cor 
# -0.1851478 
# 
# [[6]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -3.036, df = 377, p-value = 0.002564
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2513066 -0.0545978
# sample estimates:
#   cor 
# -0.1544828 
# 
# [[7]] Pearson's product-moment correlation
# data:  as.numeric(pheno$`age_m:ch1`) and beta[markers[i], ]
# t = -3.308, df = 377, p-value = 0.00103
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2642166 -0.0683744
# sample estimates:
#        cor 
# -0.1679522 

y = beta['cg17151731',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg17151731', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); abline(lm(y ~ x), col = 'red2', lty = 2)


################## PLEC

# PLEC, chr8:143921000-143926800
i = 'chr8:143921000-143926800'
UM_plot(M, U, 'cg04247530')
res = Visualize_cometh(annotation = annotation, CpG = 'cg04247530', distance = 1000, beta_mat = beta, 
                       L_bound = 4, R_bound = 4)
markers = c('cg07597069', "cg10620577", "cg12861418", "cg12044749", "cg01949498", "cg09143319", "cg04247530", "cg26120132", "cg07172676", "cg03354047",
            "cg00017203", "cg11207081", 'cg22095445')
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5))

# Ageing?
plot(pheno$`age_m:ch1`, beta['cg04247530',], ylim = c(0,1))
cor.test(as.numeric(pheno$`age_m:ch1`), beta['cg04247530',]) # cor = -0.1062376; p-value = 0.03871

y = beta['cg04247530',]; x = as.numeric(pheno$`age_m:ch1`)
plot(x, y, ylim = c(0,1), xlim = c(18, 71), main = 'cg04247530', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x,y) # p-value = 0.03871

y = beta['cg04247530',]; x = as.numeric(pheno$`cotinine_m:ch1`)
plot(x, y, ylim = c(0,1), main = 'cg04247530', xlab = 'cotinine levels (ng/mL)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); cor.test(x, y) # p-value = 0.1241


################## WDR27 (the M4 from blood was also significant JRC_SP)

UM_plot(M, U, 'cg19089141')
UM_plot(M, U, 'cg11938672')
UM_plot(M, U, 'cg18322025')
col = bGMM(M = M, U = U, CpG = 'cg11938672', K = 3, transform = T)
res = Visualize_cometh(annotation = annotation, CpG = 'cg11938672', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 1, max_y = 5)
#markers = c("cg19089141", "cg11938672", "cg18322025")
markers = c("cg11938672", "cg18322025")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 3), levels = 1:3)
levels(assign) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)[c(1,2,3)]
assign = as.character(assign)
table(assign)
# blue purple    red
# 295    82       2 

(p1 = (2+82/2)/(295 + 82 + 2)) # 0.1134565
(n1 = (295 + 82 + 2))

p = (13+83/2)/(283+83+13)
freq_exp = c(`0` = (1-p)^2, `0.5` = 2*p*(1-p), `1` = p^2)
print(chisq.test(x = c(283, 83, 13), p = freq_exp, correct = T))
# Chi-squared test for given probabilities
# data:  c(283, 83, 13)
# X-squared = 4.6397, df = 2, p-value = 0.09829


i = 'chr6:169653600-169656200'
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
178/(178+ 549) # 0.2448418

#########
#########

p1 = (2+82/2)/(295 + 82 + 2)
n1 = 2*(295 + 82 + 2)
p2 = 11/(11+119) # 0.1134565
n2 = (11+119)

MAT = matrix(c(p1*n1, (1-p1)*n1,
               p2*n2, (1-p2)*n2), nrow = 2, byrow = T)
colnames(MAT) = c('M', 'U')
rownames(MAT) = c('WB', 'SP')
#     M   U
# WB 86 672
# SP 11 119
86/(86+672) # 0.1134565
11/(11+119) # 0.08461538
fisher.test(MAT)
# Fisher's Exact Test for Count Data
# data:  MAT
# p-value = 0.4458
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7086085 2.9632909
# sample estimates:
# odds ratio 
#   1.383994 



###### VTRNA2.1, chr5:136067400-136083800

# cg18678645, cg06536614
UM_plot(M, U, 'cg18678645')
res = Visualize_cometh(annotation = annotation, CpG = 'cg18678645', distance = 1000, beta_mat = beta, 
                       L_bound = 5, R_bound = 2, max_y = 12)
markers = c("cg11608150", "cg06478886", "cg04481923", "cg18678645", "cg06536614", "cg26328633", "cg25340688", "cg26896946", "cg00124993",
            "cg08745965", "cg16615357")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 2), levels = 1:2)
levels(assign) = scales::alpha(c('darkmagenta', 'dodgerblue3'), 0.5)[c(1, 2)]
assign = as.character(assign)
table(assign)
# purple    blue
# 537       190 
i = 'chr5:136067400-136083800'
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
190/(190+ 537) # 0.261348

