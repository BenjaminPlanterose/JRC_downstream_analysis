library(data.table)
library(GEOquery)
library(minfi)
library(UMtools)
library(gplots)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
breaks=seq(-0.05, 1.05, 0.05)
my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)


setwd('/media/ultron/2tb_disk2/Twins_project/RAW_DATA/2018/twins_project/population_study/GSE87571_RAW/')
rgSet = read.metharray.exp(getwd(), extended = F)
beta = preprocessQuantile(rgSet)
beta = getBeta(beta)
raw = preprocessRaw(rgSet)
M = getMeth(raw)
U = getUnmeth(raw)
raw = getBeta(raw) # 416

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/blood_450K/PROCESSED/')
UMtools::export_bigmat(M, 'M.txt', 4)
UMtools::export_bigmat(U, 'U.txt', 4)
UMtools::export_bigmat(round(beta, 5), 'beta.txt', 4)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/blood_450K/PROCESSED/')
M = import_bigmat('2022-10-18_M.txt')
U = import_bigmat('2022-10-18_U.txt')
beta = import_bigmat('2022-10-18_beta.txt')


setwd('/media/ultron/2tb_disk2/Twins_project/RAW_DATA/2018/twins_project/population_study/GSE87571_RAW/')
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
phenotype <- getGEO('GSE87571', destdir=".")
pheno <- Reduce(rbind, lapply(phenotype, function(x) pData(phenoData(x))))
ID = sapply(strsplit(colnames(M), split = '_'), function(x) x[1])
pheno = pheno[ID,]
dim(pheno) # 727  40
age = as.numeric(pheno$`age:ch1`)
sex = factor(pheno$`gender:ch1`, levels = c('Male', 'Female'))
levels(sex) = c('deepskyblue2',  'hotpink')
sex = as.character(sex)


#### EWAS age
library(CpGassoc)
model1 <- cpg.assoc(beta.val = beta,
                    indep = age, fdr.method = 'bonferroni')
res <- model1$results
pvals <- res$P.value
names(pvals) <- res$CPG.Labels
sig_age <- names(which(p.adjust(pvals, 'bonferroni') < 0.05))
length(sig_age) # 144695

res = res[order(abs(res$T.statistic), decreasing = T),]
res[111:120,]

# cg25389087; chr18:77108400-77118200 - celltype_pval = 8.552719e-30, Lymph/Gran
# cg03431918; chr17:79734200-79743200 - celltype_pval = 0.008754707; huge effects!
# cg09401099; chr3:156815400-156816800 - celltype_pval = 0.2116346; nice!
# cg07544187; chr19:19536600-19540800 - celltype_pval = 0.00148556; 
# cg10804656; chr10:22332200-22335800 - celltype_pval = 4.9094e-07; beautiful
# cg09099868, cg13327545, 
# cg08160331; chr11:75426200-75435800 - celltype_pval = 0.01831256; beautiful
# cg05093315; chr11:18104600-18109200 - celltype_pval = 0.6939493; 
# cg23744638; chr11:10301200-10303200 - celltype_pval = 1.65098e-16; 
# cg03399905; chr15:79283400-79285000 - celltype_pval = 0.009780998
# cg01557798; chr1:228210400-228215000 - celltype_pval = 0.477629

FlowSorted.Blood.450k.compTable['cg01557798',]


###### FUT4, chr11:94544400-94556400
# cg18023065, cg04757806, cg07235053, cg10283505
i = 'chr11:94544400-94556400'
UM_plot(M, U, 'cg07235053')
res = Visualize_cometh(annotation = annotation, CpG = 'cg07235053', distance = 1000, beta_mat = beta, 
                       L_bound = 0, R_bound = 3)


markers = c('cg03543495', "cg06526620", "cg05229803", "cg11412713", "cg08863777", "cg20533957", "cg07235053", "cg04757806", "cg18023065", "cg10283505", "cg13300301",
            "cg02856190", "cg24286220")
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5))

y = beta['cg11412713',]; x = age
plot(x, y, ylim = c(0,1), main = 'cg11412713', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); abline(lm(y ~ x), col = 'red2', lty = 2); cor.test(x,y)
lapply(1:length(markers), function(i) cor.test(age, beta[markers[i],])) # Not convinced!


###### HOXB3, chr11:94544400-94556400
i = 'chr17:48574800-48577600'

# cg19710451, cg27496615, cg20693334

UM_plot(M, U, 'cg27496615')
res = Visualize_cometh(annotation = annotation, CpG = 'cg27496615', distance = 1000, beta_mat = beta, 
                       L_bound = 0, R_bound = 3)


markers = c("cg19710451", "cg27496615", "cg20693334")
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5))

y = beta['cg27496615',]; x = age
plot(x, y, ylim = c(0,1), main = 'cg11412713', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); abline(lm(y ~ x), col = 'red2', lty = 2); cor.test(x,y)
lapply(1:length(markers), function(i) cor.test(age, beta[markers[i],])) # Not convinced!

###### HOXA4, chr7:27129600-27135000
i = 'chr7:27129600-27135000'

# cg11410718
UM_plot(M, U, 'cg11410718')
res = Visualize_cometh(annotation = annotation, CpG = 'cg11410718', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = -1, max_y = 12)

markers = res[c(-1, -length(res))]
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5)) #, RowSideColors = sex)

lapply(1:length(markers), function(i) cor.test(age, beta[markers[i],])) #
# [[1]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 1.829, df = 725, p-value = 0.06781
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.00496695  0.13979451
# sample estimates:
#   cor 
# 0.06777046 
# 
# 
# [[2]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 4.1639, df = 725, p-value = 3.505e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.08101482 0.22306162
# sample estimates:
#       cor 
# 0.1528275 
# 
# 
# [[3]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 5.9629, df = 725, p-value = 3.87e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1457973 0.2844588
# sample estimates:
#   cor 
# 0.2162181 
# 
# 
# [[4]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 6.3919, df = 725, p-value = 2.931e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1609622 0.2986687
# sample estimates:
#       cor 
# 0.2309718 
# 
# 
# [[5]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 6.9063, df = 725, p-value = 1.089e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1789707 0.3154643
# sample estimates:
#   cor 
# 0.2484504 
# 
# 
# [[6]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 7.336, df = 725, p-value = 5.907e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1938634 0.3292896
# sample estimates:
#       cor 
# 0.2628708 
# 
# 
# [[7]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 8.8898, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2464196 0.3776195
# sample estimates:
#   cor 
# 0.313515 
# 
# 
# [[8]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 6.2577, df = 725, p-value = 6.681e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1562307 0.2942417
# sample estimates:
#      cor 
# 0.226372 
# 
# 
# [[9]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 7.3307, df = 725, p-value = 6.127e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1936814 0.3291210
# sample estimates:
#   cor 
# 0.2626947 
# 
# 
# [[10]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 5.7474, df = 725, p-value = 1.335e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1381354 0.2772561
# sample estimates:
#       cor 
# 0.2087516 
# 
# 
# [[11]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 9.1065, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2535749 0.3841445
# sample estimates:
#   cor 
# 0.3203806 
# 
# 
# [[12]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 8.9014, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2468024 0.3779689
# sample estimates:
#       cor 
# 0.3138824 
# 
# 
# [[13]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 9.3169, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2604777 0.3904270
# sample estimates:
#   cor 
# 0.3269972 
# 
# 
# [[14]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 8.8234, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2442155 0.3756069
# sample estimates:
#       cor 
# 0.3113987 
# 
# 
# [[15]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 9.0191, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2506932 0.3815183
# sample estimates:
#   cor 
# 0.3176164 
# 
# 
# [[16]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 6.3966, df = 725, p-value = 2.848e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1611266 0.2988225
# sample estimates:
#       cor 
# 0.2311316 
# 
# 
# [[17]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 10.663, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.3036025 0.4294043
# sample estimates:
#   cor 
# 0.3681874 
# 
# 
# [[18]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 10.997, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3140247 0.4387544
# sample estimates:
#       cor 
# 0.3781042 
# 
# 
# [[19]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 6.7181, df = 725, p-value = 3.724e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1724043 0.3093500
# sample estimates:
#   cor 
# 0.2420825 
# 
# 
# [[20]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 7.0861, df = 725, p-value = 3.275e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1852190 0.3212718
# sample estimates:
#       cor 
# 0.2545043 
# 
# 
# [[21]]
# 
# 	Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 7.5051, df = 725, p-value = 1.803e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1996841 0.3346774
# sample estimates:
#   cor 
# 0.2684985 
# 
# 
# [[22]]
# 
# Pearson's product-moment correlation
# 
# data:  age and beta[markers[i], ]
# t = 15.022, df = 725, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4297247 0.5407691
# sample estimates:
#       cor 
# 0.4872139 


which.min(sapply(1:length(markers), function(i) cor.test(age, beta[markers[i],])$p.value)) # 22
UM_plot(M = M, U = U, CpG = 'cg19142026')
y = beta['cg19142026',]; x = age
plot(x, y, ylim = c(0,1), main = 'cg19142026', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); abline(lm(y ~ x), col = 'red2', lty = 2); cor.test(x,y)

# RAW
y = raw['cg19142026',]; x = age
plot(x, y, ylim = c(0,1), main = 'cg19142026', xlab = 'age (y)', ylab = 'methylation',
     col = scales::alpha('black', 0.3), pch = 19); abline(lm(y ~ x), col = 'red2', lty = 2); cor.test(x,y)



# samtools view -b output_38.bam "chr7:27129600-27135000" > HOXA4_WB.bam
# samtools index HOXA4_WB.bam
library(FlowSorted.Blood.450k)
data(FlowSorted.Blood.450k)

data("FlowSorted.Blood.450k")
beta_sorted = getBeta(preprocessQuantile(FlowSorted.Blood.450k))

df = data.frame(beta = beta_sorted['cg19142026',]) # 
df$type = sapply(strsplit(colnames(beta_sorted), split = '_'), function(x) x[1])
df$col = df$type
df$col[df$col %in% c('Eos', 'Gran', 'Neu', 'CD14+')] = 'myeloid'
df$col[df$col %in% c('CD8+', 'CD4+', 'CD19+', 'CD56+')] = 'lymphoid'
df$col[df$col %in% c('WB', 'PBMC')] = 'mixture'
df$type = factor(df$type, levels = c('Eos', 'Gran', 'Neu', 'CD14+', 'CD8+', 'CD4+', 'CD19+', 'CD56+', 'WB', 'PBMC'))
library(ggplot2)
p = ggplot(data = df, mapping = aes(x = type, y = beta, col = col, fill = col)) + geom_boxplot() + ylim(0,1)

data("FlowSorted.Blood.450k.compTable")
FlowSorted.Blood.450k.compTable['cg19142026',]
#              Fstat  p.value      CD8T      CD4T        NK   Bcell     Mono      Gran     low      high     range
# cg19142026 1.123532 0.369222 0.5545967 0.4975732 0.5245881 0.48203 0.478173 0.4545578 0.33479 0.6560412 0.3212513



# FOXK1
df = data.frame(beta = beta_sorted['cg25344401',]) # 
df$type = sapply(strsplit(colnames(beta_sorted), split = '_'), function(x) x[1])
df$col = df$type
df$col[df$col %in% c('Eos', 'Gran', 'Neu', 'CD14+')] = 'myeloid'
df$col[df$col %in% c('CD8+', 'CD4+', 'CD19+', 'CD56+')] = 'lymphoid'
df$col[df$col %in% c('WB', 'PBMC')] = 'mixture'
df$type = factor(df$type, levels = c('Eos', 'Gran', 'Neu', 'CD14+', 'CD8+', 'CD4+', 'CD19+', 'CD56+', 'WB', 'PBMC'))
library(ggplot2)
p = ggplot(data = df, mapping = aes(x = type, y = beta, col = col, fill = col)) + geom_boxplot() + ylim(0,1)


setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/')
ggsave(plot = p, filename = 'cell_type_comp_FOXK1.tiff', device = 'tiff', units = 'in', 
       width = 6.6, height = 1.1, dpi = 300)
ggsave(plot = p, filename = 'cell_type_comp_FOXK1_2.tiff', device = 'tiff', units = 'in', 
       width = 6.6, height = 2, dpi = 300)




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

# Trigeminal
setwd('/media/ultron/2tb_disk2/0_startallover/followup_meQTLs/diandric/GSE74738_RAW/triploid/')
rgSet = minfi::read.metharray.exp(getwd())
raw = preprocessRaw(rgSet)
beta_tri = getBeta(raw)
colnames(beta_tri) = sapply(strsplit(colnames(beta_tri), split = '_'), function(x) x[1])
dim(beta_tri) # 485512     10
setwd('/media/ultron/2tb_disk2/0_startallover/followup_meQTLs/diandric/')
pheno = fread('info', header = F)
type = sapply(strsplit(pheno$V2, split = ' '), function(x) x[5])
names(type) = pheno$V1
type = type[colnames(beta_tri)]
col = factor(type); levels(col) = c('blue', 'green'); col = as.character(col)

i = 'chr5:136067400-136083800'
heatmap.2(t(beta_tri[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5), RowSideColors = col)
# diandric: blue; digynic: green

# UUU: {d_U, d_U, v_U}, {d_U, v_U, v_U}, 
# UUM: {d_U, d_U, v_M}
# UMM: {d_U, v_M, v_M}



###### WDR27, chr6:169653600-169656200

# cg19089141, cg11938672, cg18322025
UM_plot(M, U, 'cg11938672')
res = Visualize_cometh(annotation = annotation, CpG = 'cg11938672', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 1, max_y = 5)
markers = c("cg11938672", "cg18322025")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 2), levels = 1:2)
levels(assign) = scales::alpha(c('darkmagenta', 'brown3'), 0.5)[c(1, 2)]
assign = as.character(assign)
table(assign)
# purple    red
# 615       112 

i = 'chr6:169653600-169656200'
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
112/(112+ 615) # 0.1540578

i = 'chr6:169653600-169656200'
heatmap.2(t(beta_tri[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5), RowSideColors = col)
# diandric: blue; digynic: green
# ??






# DUSP22
i = 'chr6:284200-300400'
UM_plot(M, U, 'cg07332563')
res = Visualize_cometh(annotation = annotation, CpG = 'cg07332563', distance = 1000, beta_mat = beta, 
                       L_bound = 0, R_bound = 4)
markers = c("cg07332563", "cg21548813", "cg03395511", "cg15383120", "cg18110333", "cg05064044", "cg11235426", "cg01516881", "cg26668828",
            "cg01171360")
heatmap.2(t(beta_tri[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5), RowSideColors = col)
# diandric: blue; digynic: green

# CYP2E1
i = 'chr10:133523800-133534200'
UM_plot(M, U, 'cg00321709')
res = Visualize_cometh(annotation = annotation, CpG = 'cg18984983', distance = 1000, beta_mat = beta, 
                       L_bound = 5, R_bound = 1)
markers = c("cg13315147", "cg19469447", "cg00321709", "cg10862468", "cg25330361", "cg23400446", "cg24530264", "cg18984983", "cg05194426",
            "cg11445109", "cg27214960")
heatmap.2(t(beta_tri[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5), RowSideColors = col)
# diandric: blue; digynic: green


# OR2L13
i = 'chr1:247936400-247939000'
UM_plot(M, U, 'cg20433858')
res = Visualize_cometh(annotation = annotation, CpG = 'cg20434529', distance = 1000, beta_mat = beta, 
                       L_bound = -2, R_bound = 2)
markers = c("cg20434529", "cg08260406", "cg04028570", "cg03748376", "cg00785941", "cg20507276", "cg08944170", "cg20433858")
heatmap.2(t(beta_tri[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          main = i, cexCol = 0.8, margins = c(6,5), RowSideColors = col)
# diandric: blue; digynic: green



#################################################################################################

# AFR

setwd("/media/ultron/2tb_disk2/Twins_project/RAW_DATA/2018/twins_project/gambia/GSE99863_RAW/")
length(list.files(pattern = '*.idat'))/2 # 240
rgSet = read.metharray.exp(getwd(), extended = F)
beta = preprocessQuantile(rgSet)
beta = getBeta(beta)
raw = preprocessRaw(rgSet)
M = getMeth(raw)
U = getUnmeth(raw)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/gambia/')
UMtools::export_bigmat(M, 'M.txt', 4)
UMtools::export_bigmat(U, 'U.txt', 4)
UMtools::export_bigmat(round(beta, 5), 'beta.txt', 4)

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/gambia/')
M = import_bigmat('2022-10-19_M.txt')
U = import_bigmat('2022-10-19_U.txt')
beta = import_bigmat('2022-10-19_beta.txt')


UM_plot(M, U, 'cg11938672')
res = Visualize_cometh(annotation = annotation, CpG = 'cg11938672', distance = 1000, beta_mat = beta, 
                       L_bound = 1, R_bound = 1, max_y = 5)
markers = c("cg11938672", "cg18322025")
assign = factor(cutree(hclust(d = dist(t(beta[markers,]))), 2), levels = 1:2)
levels(assign) = scales::alpha(c('darkmagenta', 'brown3'), 0.5)[c(1, 2)]
assign = as.character(assign)
table(assign)
# purple    red
# 215       25 

plot(beta['cg11938672',], beta['cg18322025',], xlim = c(0,1), ylim = c(0,1)); abline(0,1, lty = 2)


i = 'chr6:169653600-169656200'
heatmap.2(t(beta[markers,]), trace = 'none', Colv = 'none', breaks = breaks,
          na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
          RowSideColors = assign, main = i, cexCol = 0.8, margins = c(6,5))
25/(25+ 215) # 0.1041667


#####################################################

# chr6:169653600-169656200
# samtools view -b output_38.bam "chr6:169653600-169656200" > WDR27_WB.bam
# samtools index WDR27_WB.bam

# chr6:169653000-169656200
# samtools view -b output_38_sperm.bam "chr6:169653000-169656200" > WDR27_SP.bam
# samtools index WDR27_SP.bam

f = c(1, 111, 615)





#### VARIANTS

# AMERICAN
86/(86+672) # 0.1134565
11/(11+119) # 0.08461538

# EUROPEAN
112/(112+ 615) # 0.1540578

# AFRICAN
25/(25+ 215) # 0.1041667

p1 = 112/(112+ 615); n1 = (112+ 615)
p2 = 25/(25+ 215); n2 = (25+ 215)
MAT = matrix(c(p1*n1, (1-p1)*n1,
               p2*n2, (1-p2)*n2), nrow = 2, byrow = T)
colnames(MAT) = c('M', 'U')
rownames(MAT) = c('EUR', 'AFR')
#      M   U
# EUR 112 615
# AFR  25 215
fisher.test(MAT)
# Fisher's Exact Test for Count Data
# data:  MAT
# p-value = 0.05528
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9762233 2.5924837
# sample estimates:
#   odds ratio 
# 1.565512 

(p1*n1+p2*n2)/(n1+n2) # 0.1416753



library(data.table)
library(matrixStats)
setwd('/media/ultron/2tb_disk1/ben/dbSNP/')
dbSNP = fread('dbSnp155Common.bed', nThread = 4)
dim(dbSNP) # 15162080       17
# i = 'chr6:169653600-169656200'
dbSNP = dbSNP[(dbSNP$V1 == 'chr6') & (dbSNP$V2 > 169653600) & (dbSNP$V2 < 169656200), ]
dim(dbSNP) # 29       17

# SNP           MAF_AFR MAF_EUR  MAF_AMR
# rs115553893   A=0.0567  A=0.0050  A=0.019
# rs144460930 G=0.0567  G=0.0050  G=0.019
# rs75964125  A=0.1679   	A=0.0716  A=0.052
# rs4716337 C=0.3971  C=0.2097  C=0.532
# rs78828886  T=0.0552  T=0.0000  T=0.003
# rs73790083  G=0.1974  G=0.1252  G=0.115
# rs9383513 G=0.0098  G=0.0547  G=0.272
# rs9371121 C=0.0038  C=0.0547   	C=0.272
# rs75659875  C=0.1324  C=0.0716  C=0.050
# rs576807990 G=0.0212  G=0.0000  G=0.001
# rs541149154  	T=0.0212  T=0.0000  T=0.001
# rs529850726 C=0.0318   	C=0.0000  C=0.003
# rs192641460 T=0.0197  T=0.0000  T=0.000
# rs530751994 G=0.0204  G=0.0000  G=0.000
# rs568121327 A=0.1082  A=0.0726  A=0.058
# rs79371338  T=0.0953  T=0.0249  T=0.125
# rs73790085  T=0.1959  T=0.1252  T=0.114
# rs9383514 A=0.0098  A=0.0547   	A=0.272
# rs4716371 G=0.8752  G=0.3539  G=0.659
# rs9371159 A=0.0098  A=0.0547  A=0.272
# rs13362807  A=0.1944  A=0.1243  A=0.114
# rs9371160 C=0.0038  C=0.0547  C=0.272
# rs150200200 A=0.0250  A=0.0467  A=0.037
# rs61578090  T=0.1233  T=0.0249  T=0.125
# rs9371161 G=0.0098  G=0.0547  G=0.272
# rs6927257 A=0.1014   	A=0.0716  A=0.049
# rs7763094 A=0.1944   	A=0.1252  A=0.114
# rs34319658  C=0.0023  C=0.0626  C=0.027
# rs79824271  delTC=0.1014  delTC=0.0716  delTC=0.049


rs = matrix(c('rs115553893',   'A=0.0567',  'A=0.0050',  'A=0.019',
       'rs144460930', 'G=0.0567',  'G=0.0050',  'G=0.019',
       'rs75964125',  'A=0.1679',   	'A=0.0716',  'A=0.052',
       'rs4716337', 'C=0.3971',  'C=0.2097', 'C=0.532',
       'rs78828886',  'T=0.0552',  'T=0.0000',  'T=0.003',
       'rs73790083',  'G=0.1974',  'G=0.1252',  'G=0.115',
       'rs9383513', 'G=0.0098',  'G=0.0547',  'G=0.272',
       'rs9371121', 'C=0.0038',  'C=0.0547',   	'C=0.272',
       'rs75659875',  'C=0.1324',  'C=0.0716', 'C=0.050',
       'rs576807990', 'G=0.0212',  'G=0.0000',  'G=0.001',
       'rs541149154',  	'T=0.0212',  'T=0.0000',  'T=0.001',
       'rs529850726', 'C=0.0318',   	'C=0.0000',  'C=0.003',
       'rs192641460', 'T=0.0197', 'T=0.0000',  'T=0.000',
       'rs530751994', 'G=0.0204',  'G=0.0000',  'G=0.000',
       'rs568121327', 'A=0.1082', 'A=0.0726',  'A=0.058',
       'rs79371338',  'T=0.0953',  'T=0.0249',  'T=0.125',
       'rs73790085',  'T=0.1959',  'T=0.1252',  'T=0.114',
       'rs9383514', 'A=0.0098',  'A=0.0547',  	'A=0.272',
       'rs4716371', 'G=0.8752',  'G=0.3539', 'G=0.659',
       'rs9371159', 'A=0.0098',  'A=0.0547',  'A=0.272',
       'rs13362807',  'A=0.1944',  'A=0.1243',  'A=0.114',
       'rs9371160', 'C=0.0038',  'C=0.0547',  'C=0.272',
       'rs150200200', 'A=0.0250', 'A=0.0467',  'A=0.037',
       'rs61578090',  'T=0.1233',  'T=0.0249',  'T=0.125',
       'rs9371161', 'G=0.0098',  'G=0.0547',  'G=0.272',
       'rs6927257', 'A=0.1014',   	'A=0.0716', 'A=0.049',
       'rs7763094', 'A=0.1944',   	'A=0.1252',  'A=0.114',
       'rs34319658',  'C=0.0023',  'C=0.0626', 'C=0.027',
       'rs79824271',  'delTC=0.1014',  'delTC=0.0716',  'delTC=0.049'), ncol = 4, byrow = T)
dim(rs) # 29 4
rs = as.data.frame(rs)
colnames(rs) = c('rs', 'MAF_AFR', 'MAF_EUR', 'MAF_AMR')

dbSNP = cbind(dbSNP, rs)
mean(dbSNP$V4 == dbSNP$rs) # 1
dbSNP$ALT = sapply(strsplit(dbSNP$MAF_AFR, split = '='), function(x) x[1])
table(dbSNP$ALT, dbSNP$V7)
dbSNP$MAF_AFR = as.numeric(sapply(strsplit(dbSNP$MAF_AFR, split = '='), function(x) x[2]))
dbSNP$MAF_EUR = as.numeric(sapply(strsplit(dbSNP$MAF_EUR, split = '='), function(x) x[2]))
dbSNP$MAF_AMR = as.numeric(sapply(strsplit(dbSNP$MAF_AMR, split = '='), function(x) x[2]))

#dbSNP = dbSNP[rowMaxs(as.matrix(dbSNP[,c('MAF_AFR', 'MAF_EUR', 'MAF_AMR')])) > 0.05,]
dim(dbSNP) # 23 22

xlim = c(169653600, 169656200)
reg = c(169654291, 169655605)
cg_pos = c(169655060, 169655237)
cg_pos_out = c(169654635)




DO_PLOT <- function(reg, xlim, cg_pos, cg_pos_out, dbSNP)
{
  # Pre-computed
  EUR = 0.1540578; AFR = 0.1041667
  
  # Plot
  pos = dbSNP$V2
  plot(pos, rep(1, length(pos)), xlim = xlim, type = "n", ylab = "", 
       pch = 19, ylim = c(0, 1.05), main = "1000genomes") # i = 'chr6:169653600-169656200'
  polygon(x = c(reg[2], reg[2], reg[1], reg[1]), y = c(1, 1.08, 1.08, 1), border = F,
          col = "antiquewhite3")
  arrows(x0 = xlim[1], x1 = xlim[2], y0 = 1, y1 = 1, code = 3, angle = 90, col = "gray")
  points(pos, rep(1, times = length(pos)),  pch = 19)
  points(cg_pos, rep(1.05, times = length(cg_pos)), col = "red", pch = 19)
  points(cg_pos_out, rep(1.05, times = length(cg_pos_out)), col = "blue", pch = 19)
  abline(h = EUR, col = "antiquewhite3", lty = 2, lwd = 2)
  abline(h = AFR, col = "brown", lty = 2, lwd = 2)
  
  EUR_confint = prop.test(x=112, n=727, conf.level=.95, correct=T)$conf.int[1:2]
  AFR_confint = prop.test(x=25, n=240, conf.level=.95, correct=T)$conf.int[1:2]
  
  rect(xleft = min(xlim)-1000, xright = max(xlim)+1000, ybottom = EUR_confint[1], ytop = EUR_confint[2], col = scales::alpha('gray', 0.3), border = NA)
  rect(xleft = min(xlim)-1000, xright = max(xlim)+1000, ybottom = AFR_confint[1], ytop = AFR_confint[2], col = scales::alpha('brown', 0.3), border = NA)
  
  
  M = sapply(1:nrow(dbSNP), function(i) max(c(dbSNP$MAF_AFR[i], dbSNP$MAF_EUR[i])))
  m = sapply(1:nrow(dbSNP), function(i) min(c(dbSNP$MAF_AFR[i], dbSNP$MAF_EUR[i])))
  arrows(x0 = pos, x1 = pos, y0 = M, y1 = m, code = 0, col = "gray")
  points(pos, dbSNP$MAF_AFR, col = "brown", pch = 19, cex = 1.2)
  points(pos, dbSNP$MAF_EUR,  col = "antiquewhite3", pch = 19, cex = 1.2)
}

DO_PLOT(reg, xlim, cg_pos, cg_pos_out, dbSNP)

#score = (dbSNP$MAF_AFR-AFR)^2 + (dbSNP$MAF_EUR-EUR)^2 + (dbSNP$MAF_AMR-AMR)^2
score = ((dbSNP$MAF_AFR-AFR)^2 + (dbSNP$MAF_EUR-EUR)^2)/(AFR+EUR)
names(score)= dbSNP$V4


DO_PLOT(reg, xlim, cg_pos, cg_pos_out, dbSNP[order(score)[9],])

# rs73790083 	169,654,453 	G = 0.120008 			        AF_approved DOUBT_bypos
# rs75659875 	169,654,704 	C = 0.0880591 			EGR1
# rs79371338 	169,654,759 	T = 0.108626 			
# rs73790085 	169,654,875 	T = 0.117412 			        AF_approved
# rs9371159 	169,655,407 	A = 0.197085

abline(v = 169654453)
abline(v = 169654759)





abline(v = 169655521)


sort(score) # rs7763094  rs13362807  rs73790085 rs568121327  rs73790083
DO_PLOT(reg, xlim, cg_pos, cg_pos_out, dbSNP[order(score)[1:3],])


top10 = c('rs568121327', 'rs6927257', 'rs79824271', 'rs75659875', 'rs7763094', 'rs13362807', 'rs73790085',
          'rs73790083', 'rs75964125', 'rs79371338')
jkl = assoc[assoc$SNP %in% top10,]
jkl[order(jkl$`Effect Size`, decreasing = T),]

library(data.table)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/mQTL_WDR27/')
assoc = fread('cg11938672.csv')
assoc = assoc[assoc$Trans == 'N']
head(assoc)
table(assoc$Timepoint)

which.max(assoc$`Effect Size`)
assoc[25,]


assoc$SNP[order(assoc$`p-value`)] # "rs13362807"  "rs73790085"  "rs7763094"   "rs11463842"  "rs4440470"

pos = 169655060+(assoc$`SNP Pos` - assoc$`CpG Pos`)
plot(pos[order(pos)], -log10(assoc$`p-value`)[order(pos)], type = 'o', ylim = c(0,15), 
     xlim = c(169330111, 169848697+2000), xlab = 'coordinate (chr6)', ylab = '-log10(pval)')
polygon(x = c(reg[2], reg[2], reg[1], reg[1]), y = c(0, 15, 15, 0), border = F,
        col = scales::alpha("antiquewhite3", 0.5))


plot(pos[order(pos)], -log10(assoc$`p-value`)[order(pos)], type = 'o', ylim = c(0,15), 
     xlim = c(169620000, 169680000), xlab = 'coordinate (chr6)', ylab = '-log10(pval)')
polygon(x = c(reg[2], reg[2], reg[1], reg[1]), y = c(0, 15, 15, 0), border = F,
        col = scales::alpha("antiquewhite3", 0.5))

### Effect size
pos = 169655060+(assoc$`SNP Pos` - assoc$`CpG Pos`)
plot(pos[order(pos)], assoc$`Effect Size`[order(pos)], type = 'o', ylim = c(0,0.5), 
     xlim = c(169330111, 169848697+2000), xlab = 'coordinate (chr6)', ylab = '-log10(pval)')
polygon(x = c(reg[2], reg[2], reg[1], reg[1]), y = c(0, 15, 15, 0), border = F,
        col = scales::alpha("antiquewhite3", 0.5))

### AF
pos = 169655060+(assoc$`SNP Pos` - assoc$`CpG Pos`)
plot(pos[order(pos)], assoc$MAF[order(pos)], type = 'o', ylim = c(0,0.5), 
     xlim = c(169330111, 169848697+2000), xlab = 'coordinate (chr6)', ylab = '-log10(pval)')
polygon(x = c(reg[2], reg[2], reg[1], reg[1]), y = c(0, 15, 15, 0), border = F,
        col = scales::alpha("antiquewhite3", 0.5))



# chr6	170054800	A	C	intronic	WDR27	.	.	rs75659875	EGR1
# chr6	170056141	CT	-	intronic	WDR27	.	.	rs79824271	Meis1
# chr6	170055982	A	G	intronic	WDR27	.	.	rs9371161	Pdx1
# chr6	170055982	A	G	intronic	WDR27	.	.	rs9371161	Prrx2
# chr6	170055096	G	A	intronic	WDR27	.	.	rs9383514	REL
# chr6	170055618	G	A	intronic	WDR27	.	.	rs13362807	SOX10
# chr6	170055750	G	C	intronic	WDR27	.	.	rs9371160	Tcf3




abline(v = 169654348, col = 'red2')
abline(v = 169655497, col = 'red2')

'a'

head(dbSNP)
head(dbSNP)

dbSNP = dbSNP[,c(1:5, 7, 10, 14)]
dbSNP$V10 = as.numeric(sapply(strsplit(dbSNP$V10, split = ','), function(x) x[1]))
colnames(dbSNP) = c('chr', 'start', 'end', 'rs', 'ref', 'alt', 'maf', 'type')
dbSNP = dbSNP[dbSNP$type == 'snv',]
dbSNP = dbSNP[dbSNP$maf >= 0.05, ]






# chr6:169653600-169656200

SNPs = c('rs75964125', 'rs4716337', 'rs73790083', 'rs9383513', 'rs9371121', 'rs75659875', 'rs568121327', 
         'rs79371338', 'rs73790085', 'rs9383514', 'rs4716371', 'rs9371159', 'rs13362807', 'rs9371160', 
         'rs61578090', 'rs9371161', 'rs6927257', 'rs7763094', 'rs79824271')
pos = c(169654194, 169654271, 169654453, 169654634, 169654690, 169654704, 169654758,
        169654759, 169654875, 169655000, 169655220, 169655407, 169655522, 169655654,
        169655782, 169655886, 169655924, 169655958, 169656044)
plot(c(169653600, 169656200), c(1,1)) # 2.6 kb
abline(v = pos)

length(SNPs) # 19
c(0, 0, 0, 0, )

# rs75964125 	169,654,194 	single nucleotide variant 	WDR27 and 1 more 	intron variant 	Not-Provided 	A = 0.0978435 			
# rs4716337 	169,654,271 	single nucleotide variant 	WDR27 and 3 more 	intron variant 	Not-Provided 	C = 0.453474 			
# rs73790083 	169,654,453 	single nucleotide variant 	WDR27 and 1 more 	intron variant 	Not-Provided 	G = 0.120008 			
# rs9383513 	169,654,634 	single nucleotide variant 	WDR27 and 3 more 	3 prime UTR variant, intron variant 	Not-Provided 	G = 0.194888 			
# rs9371121 	169,654,690 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	C = 0.192093 			
# rs75659875 	169,654,704 	single nucleotide variant 	WDR27 and 1 more 	3 prime UTR variant, intron variant 	Not-Provided 	C = 0.0880591 			
# rs568121327 	169,654,758 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.0672923 			
# rs79371338 	169,654,759 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	T = 0.108626 			
# rs73790085 	169,654,875 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	T = 0.117412 			
# rs9383514 	169,655,000 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.195088 			
# rs4716371 	169,655,220 	single nucleotide variant 	WDR27 and 3 more 	3 prime UTR variant, intron variant 	Not-Provided 	C = 0.285543 			
# rs9371159 	169,655,407 	single nucleotide variant 	WDR27 and 1 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.197085 			
# rs13362807 	169,655,522 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.116813 			
# rs9371160 	169,655,654 	single nucleotide variant 	WDR27 and 3 more 	3 prime UTR variant, intron variant 	Not-Provided 	C = 0.19349 			
# rs61578090 	169,655,782 	single nucleotide variant 	WDR27 and 1 more 	3 prime UTR variant, intron variant 	Not-Provided 	T = 0.115415 			
# rs9371161 	169,655,886 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	G = 0.197085 			
# rs6927257 	169,655,924 	single nucleotide variant 	WDR27 and 1 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.0796725 			
# rs7763094 	169,655,958 	single nucleotide variant 	WDR27 and 2 more 	3 prime UTR variant, intron variant 	Not-Provided 	A = 0.117013 			
# rs79824271 	169,656,044 - 169,656,049 	indel 	WDR27 and 1 more 	3 prime UTR variant, intron variant 	Not-Provided 	CTC = 0.0796725


write.table(SNPs, sep = '\t', quote = F, row.names = F, col.names = F)
# chr6	170054800	A	C	intronic	WDR27	.	.	rs75659875	EGR1
# chr6	170056141	CT	-	intronic	WDR27	.	.	rs79824271	Meis1
# chr6	170055982	A	G	intronic	WDR27	.	.	rs9371161	Pdx1
# chr6	170055982	A	G	intronic	WDR27	.	.	rs9371161	Prrx2
# chr6	170055096	G	A	intronic	WDR27	.	.	rs9383514	REL
# chr6	170055618	G	A	intronic	WDR27	.	.	rs13362807	SOX10
# chr6	170055750	G	C	intronic	WDR27	.	.	rs9371160	Tcf3

#TF Name	#TF-SNP matches (genome-wide)	#TF-SNP hits (from query)	#Fraction of TF-SNP hits	#Enrichment	#P-value
# REL	6202	1	0.000161238310222509	92.2498656347415	0.0107912801867274
# Meis1	13726	1	7.28544368352033e-05	41.6824760794599	0.023752353968951
# Tcf3	23038	1	4.34065457070926e-05	24.8343461527332	0.0395972389305107
# EGR1	70916	1	1.41011901404479e-05	8.06776561941828	0.117722165711536
# Pdx1	102271	1	9.77794291636925e-06	5.5942903331997	0.165957166993889
# SOX10	139269	1	7.18034882134574e-06	4.10811929910222	0.220026965580507
# Prrx2	381402	1	2.6219054960383e-06	1.50008040510188	0.506709985738502


#################################################################################

# JRCs Vs EPIC/450K
library(data.table)
library(GenomicRanges)
library(ggplot2)

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


setwd('/media/ultron/2tb_disk1/ben/450K/')
bed_450K = fread('450K_hg38.bed', skip = 1) # liftover of IlluminaHumanMethylation450kanno.ilmn12.hg19
bed_450K = bed_450K[bed_450K$V1 != '*',]
bed_EPIC = fread('EPIC_hg38.bed', skip = 1) # liftover of IlluminaHumanMethylationEPICanno.ilm10b4.hg19 (more probes than IlluminaHumanMethylationEPICanno.ilm10b2.hg19?)
bed_EPIC = bed_EPIC[bed_EPIC$V1 != '*',]
dim(bed_450K) # 477344      5
dim(bed_EPIC) # 865762      5


coord_450K = paste(bed_450K$V1, paste(bed_450K$V2, bed_450K$V2+1L, sep = '-'), sep = ':')
coord_EPIC = paste(bed_EPIC$V1, paste(bed_EPIC$V2, bed_EPIC$V2+1L, sep = '-'), sep = ':')
set_450K = datatable2grange(coordinates2datatable(coord_450K))
set_EPIC = datatable2grange(coordinates2datatable(coord_EPIC))


WB_450K = countOverlaps(A, set_450K)
WB_EPIC = countOverlaps(A, set_EPIC)
SP_450K = countOverlaps(B, set_450K)
SP_EPIC = countOverlaps(B, set_EPIC)

A = coordinates2datatable(sig_WB); dim(A) # 4215    4
A$`450K` = WB_450K
A$EPIC = WB_EPIC

B = coordinates2datatable(sig_SP); dim(B) # 7386    4
B$`450K` = SP_450K
B$EPIC = SP_EPIC

df.a = melt(A, id.vars = c('chr', 'range'), measure.vars = c('450K', 'EPIC'))
df.b = melt(B, id.vars = c('chr', 'range'), measure.vars = c('450K', 'EPIC'))

ggplot(data = df.a, mapping = aes(x = log10(range), y = value, col = variable)) + geom_point(alpha = 0.1) + 
  geom_smooth(method = glm, method.args = list(family = "poisson")) + facet_wrap( ~ chr) + ylim(0,125) + xlim(2,5)
mod.a = glm(formula = value ~ log10(range), data = df.a, family = 'poisson')


summary(mod.a)
#   glm(formula = value ~ log10(range), family = "poisson", data = df.a)
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -7.5663  -1.9401  -0.4441   1.0184  18.8634  
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  -2.52499    0.04035  -62.57   <2e-16 ***
#   log10(range)  1.24709    0.01062  117.46   <2e-16 ***
# (Dispersion parameter for poisson family taken to be 1)
# Null deviance: 61569  on 8429  degrees of freedom
# Residual deviance: 46902  on 8428  degrees of freedom
# AIC: 74322
# Number of Fisher Scoring iterations: 5

ggplot(data = df.b, mapping = aes(x = log10(range), y = value, col = variable)) + geom_point(alpha = 0.1) + 
  geom_smooth(method = glm, method.args = list(family = "poisson")) + facet_wrap( ~ chr) + ylim(0,125) + xlim(2,5)
mod.b = glm(formula = value ~ log10(range), data = df.b, family = 'poisson')
summary(mod.b)
# Call:
# glm(formula = value ~ log10(range), family = "poisson", data = df.b)
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -6.3663  -1.5078  -1.0833   0.0155  16.7061  
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  -6.00520    0.08060  -74.50   <2e-16 ***
#   log10(range)  1.81450    0.02285   79.41   <2e-16 ***
# (Dispersion parameter for poisson family taken to be 1)
# Null deviance: 50118  on 14771  degrees of freedom
# Residual deviance: 43914  on 14770  degrees of freedom
# AIC: 58525
# Number of Fisher Scoring iterations: 6



df.a$tissue = 'WB'
df.b$tissue = 'SP'
df = rbind(df.a, df.b)
colnames(df) = c('chr', 'JRC_size', 'microarray', 'probe_count', 'tissue')
mod.b = glm(formula = probe_count ~ log10(JRC_size)*tissue + microarray, data = df, family = 'poisson')
summary(mod.b)
#   glm(formula = probe_count ~ log10(JRC_size) * tissue + microarray,
#       family = "poisson", data = df)
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -8.0927  -1.5625  -0.9512   0.4896  19.7539
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#   (Intercept)              -6.16067    0.08070  -76.34   <2e-16 ***
#   log10(JRC_size)           1.81450    0.02285   79.41   <2e-16 ***
#   tissueWB                  3.48021    0.09014   38.61   <2e-16 ***
#   microarrayEPIC            0.28998    0.00690   42.03   <2e-16 ***
#   log10(JRC_size):tissueWB -0.56741    0.02519  -22.52   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 177202  on 23201  degrees of freedom
# Residual deviance:  89032  on 23197  degrees of freedom
# AIC: 131064
# 
# Number of Fisher Scoring iterations: 6



