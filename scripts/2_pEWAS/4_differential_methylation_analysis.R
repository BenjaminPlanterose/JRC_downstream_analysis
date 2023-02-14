############################################################################
############################################################################
###########                                                      ###########
###########                      Pooled EWAS                     ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########           b.planterosejimenez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

library(DSS)
library(data.table)
library(GenomicRanges)

prep_format <- function(mat)
{
  data.frame(chr = mat$V1, pos = mat$V2, N = mat$V5, X = round(mat$V4*mat$V5))
}

showOneDML <- function (OneDMR, BSobj, ext = 500, ylim = c(0, 1)) 
{
  allchr = as.character(seqnames(BSobj))
  allpos = start(BSobj)
  X = getBSseq(BSobj, "M")
  N = getBSseq(BSobj, "Cov")
  chr = as.character(OneDMR$chr)
  ix.chr = which(allchr == chr)
  thispos = allpos[ix.chr]
  thisN = N[ix.chr, ]
  thisX = X[ix.chr, ]
  xlim = c(OneDMR$pos - ext, OneDMR$pos + ext)
  ix1 = which(thispos <= xlim[2] & thispos >= xlim[1])
  nSample = ncol(X)
  if (nSample > 2) {
    y.cex = 0.66
  }
  else y.cex = 1
  sNames = sampleNames(BSobj)
  par(mfrow = c(nSample, 1), mar = c(2.5, 2.5, 1.6, 2.5), mgp = c(1.5, 
                                                                  0.5, 0))
  thisP = thisX/thisN
  maxN = max(thisN[ix1, ])
  for (i in 1:ncol(X)) {
    plot(thispos[ix1], thisP[ix1, i], type = "h", col = "blue", 
         axes = F, lwd = 1.5, xlab = "", ylab = "", ylim = ylim, 
         xlim = xlim, main = sNames[i])
    abline(v = OneDMR$pos, lty = 2, col = 'yellow2')
    box(col = "black")
    axis(1, )
    axis(2, col = "blue", col.axis = "blue")
    mtext(chr, side = 1, line = 1.33, cex = y.cex)
    mtext("methyl%", side = 2, line = 1.33, col = "blue", 
          cex = y.cex)
    thisN.norm = thisN[ix1, i]/maxN * ylim[2]
    lines(thispos[ix1], thisN.norm, type = "l", col = "gray", 
          lwd = 1.5)
    if (maxN >= 4) {
      axis(side = 4, at = seq(0, ylim[2], length.out = 5), 
           labels = round(seq(0, maxN, length.out = 5)))
    }
    else {
      axis(side = 4, at = seq(0, ylim[2], length.out = maxN + 
                                1), labels = 0:maxN)
    }
    mtext("read depth", side = 4, line = 1.33, cex = y.cex)
    rect(OneDMR$pos, ylim[1], OneDMR$pos, ylim[2], col = "#FF00001A", 
         border = NA)
  }
}

count_cgs <- function(dmrs, files)
{
  counts = matrix(nrow = nrow(dmrs), ncol = length(files))
  rownames(counts) = paste(dmrs$chr, paste(dmrs$start, dmrs$end, sep = '-'), sep = ':')
  colnames(counts) = files
  for(j in 1:length(files))
  {
    print(paste(j, 'out of ', length(files)))
    for(i in 1:nrow(dmrs))
    {
      print(paste(i, 'out of ', nrow(dmrs)))
      CHR = dmrs$chr[i]
      START = dmrs$start[i]
      END = dmrs$end[i]
      command = paste("awk ", "'", "$1==", '"', CHR, '"', " && ", "$2>=", START, " && ", "$2<", END, "' ", file_i, ' | wc -l', sep = "")
      counts[i,j] = fread(cmd = noquote(command))$V1
    }
  }
  return(counts)
}

export_bed <- function(dmrs, name)
{
  bed_df <- data.frame(seqname = dmrs$chr,
                       start = as.integer(dmrs$start), end = as.integer(dmrs$end),
                       name = paste(dmrs$chr, paste(dmrs$start, dmrs$end, sep = '-'), sep = ':'),
                       score = 1000)
  write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
              file = paste(name, '.bed', sep = ''),
              quote = F, sep = '\t', row.names = F, col.names = F)
  fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
         sep = '\t', col.names = F, append = T)
}

# Sample Info
# ERR2722068: WB_old
# ERR2722069: sperm_old
# ERR2722070: WB_young
# ERR2722071: sperm_young

###################################### EWAS ######################################

# 1. Create BS object
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/hg38/')
WB_old = prep_format(fread('ERR2722068_mergecg.bed', nThread = 4))
sperm_old = prep_format(fread('ERR2722069_mergecg.bed', nThread = 4))
WB_young = prep_format(fread('ERR2722070_mergecg.bed', nThread = 4))
sperm_young = prep_format(fread('ERR2722071_mergecg.bed', nThread = 4))
BSobj = makeBSseqData( list(WB_old, sperm_old, WB_young, sperm_young),
                       c("WB_old", "sperm_old", "WB_young", "sperm_young"))
# 1: In min(dd) : no non-missing arguments to min; returning Inf (times 6)
BSobj # 26,112,563 methylation loci
rm(WB_old, sperm_old, WB_young, sperm_young); gc() # Free some RAM

# 2.1 Perform hypothesis testing - WB, old vs young
age_WB = DMLtest(BSobj, group1 = c("WB_old"), group2 = c("WB_young"), 
                     smoothing = TRUE, ncores = 4, equal.disp = FALSE, smoothing.span = 500)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(age_WB, 'age_WB.Rds')
rm(age_WB); gc()

# 2.2 Perform hypothesis testing - sperm, old vs young
age_sperm = DMLtest(BSobj, group1 = c("sperm_old"), group2 = c("sperm_young"), 
                     smoothing = TRUE, ncores = 4, equal.disp = FALSE, smoothing.span = 500)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(age_sperm, 'age_sperm.Rds')
rm(age_sperm); gc()

# 2.2 Perform hypothesis testing - sperm vs blood
tissue = DMLtest(BSobj, group1 = c("sperm_young", "sperm_old"), group2 = c("WB_old", "WB_young"), 
                 smoothing = TRUE, ncores = 4, equal.disp = FALSE, smoothing.span = 500)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(tissue, 'tissue.Rds')
#rm(tissue); gc()

###################################### Call DMRs ######################################

# read RDS
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
age_WB = readRDS('age_WB.Rds')
age_sperm = readRDS('age_sperm.Rds')
tissue = readRDS('tissue.Rds')

### Call sig DMPs (not so interesting to us)
#dmls.age_WB = callDML(age_WB, delta=0.1, p.threshold=1e-5)
#dmls.age_sperm = callDML(age_sperm, delta=0.1, p.threshold=1e-5)
#dmls.tissue = callDML(tissue, delta=0.1, p.threshold=1e-5)


######## WB - old Vs young ########

dmrs.age_WB = callDMR(age_WB, delta=0.1, p.threshold=1e-5, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(dmrs.age_WB) # 135      9
Bg.age_WB = callDMR(age_WB, delta=0, p.threshold=1, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(Bg.age_WB) # 2822    9
#DSS::showOneDMR(dmrs.age_WB[1,], BSobj[,c("WB_old", "WB_young")])
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(dmrs.age_WB, 'dmrs.age_WB.Rds')
saveRDS(Bg.age_WB, 'Bg.age_WB.Rds')


######## Sperm - old Vs young ########

dmrs.age_sperm = callDMR(age_sperm, delta=0.1, p.threshold=1e-5, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(dmrs.age_sperm) # 1826      9
Bg.age_sperm = callDMR(age_sperm, delta=0, p.threshold=1, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(Bg.age_sperm) # 2742      9
#DSS::showOneDMR(dmrs.age_sperm[1,], BSobj[,c("sperm_old", "sperm_young")])
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(dmrs.age_sperm, 'dmrs.age_sperm.Rds')
saveRDS(Bg.age_sperm, 'Bg.age_sperm.Rds')

######## WB Vs Sperm ########

dmrs.tissue = callDMR(tissue, delta=0.1, p.threshold=1e-5, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(dmrs.tissue) # 118864      9
Bg.tissue = callDMR(tissue, delta=0, p.threshold=1, minlen=200, minCG=5, dis.merge=200, pct.sig=0.5)
dim(Bg.tissue) # 2915      9

setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
saveRDS(dmrs.tissue, 'dmrs.tissue.Rds')
saveRDS(Bg.tissue, 'Bg.tissue.Rds')

###################################### EWAS visualizations ######################################

library(hexbin)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
dmrs.age_WB = readRDS('dmrs.age_WB.Rds')
dim(dmrs.age_WB) # 135   9
dmrs.age_sperm = readRDS('dmrs.age_sperm.Rds')
dim(dmrs.age_sperm) # 1826    9
dmrs.tissue = readRDS('dmrs.tissue.Rds')
dim(dmrs.tissue) # 118864      9


# WB_age
plot(hexbin(dmrs.age_WB$meanMethy1, dmrs.age_WB$meanMethy2, xbins = 30, xbnds = c(0,1), ybnds = c(0,1),
            xlab = 'Mean methylation - old', ylab = 'Mean methylation - young'), main = 'Age_WB')
plot(hexbin(log10(dmrs.age_WB$length), dmrs.age_WB$meanMethy1 - dmrs.age_WB$meanMethy2,
            xbins = 40, ybnds = c(-1,1),
            xlab = 'log10(size)', ylab = 'beta_old - beta_young'), main = 'Age_WB')

# WB_SP
plot(hexbin(dmrs.age_sperm$meanMethy1, dmrs.age_sperm$meanMethy2, xbins = 100, xbnds = c(0,1), ybnds = c(0,1),
            xlab = 'Mean methylation - old', ylab = 'Mean methylation - young'), main = 'Age_Sperm')
plot(hexbin(log10(dmrs.age_sperm$length), dmrs.age_sperm$meanMethy1 - dmrs.age_sperm$meanMethy2,
            xbins = 40, ybnds = c(-1,1),
            xlab = 'log10(size)', ylab = 'beta_old - beta_young'), main = 'Age_WB')

# Tissue

plot(hexbin(dmrs.tissue$meanMethy1, dmrs.tissue$meanMethy2, xbins = 100, xbnds = c(0,1), ybnds = c(0,1),
            xlab = 'Mean methylation in Sperm', ylab = 'Mean methylation in WB'), main = 'WB Vs Sperm')
plot(hexbin(log10(dmrs.tissue$length), dmrs.tissue$meanMethy1 - dmrs.tissue$meanMethy2,
            xbins = 100, ybnds = c(-1,1),
            xlab = 'log10(size)', ylab = 'beta_SP - beta_WB'), main = 'WB Vs Sperm')




###################################### Prep Genomic Ranges ######################################

dmrs2genomicranges <- function(mat)
{
  # remove non-standard assemblies
  LEV = paste('chr', c(1:22, 'X', 'Y'), sep = '')
  mat = mat[mat$chr %in% LEV,]
  # 
  GRanges(seqnames = mat$chr, ranges = IRanges(start = mat$start, end = mat$end))
}

#mat = dmrs.age_sperm
# filter_blacklist_target <- function(mat, blacklist, map, pct.unmap = 0.05, pct.black = 0.05)
# {
#   # Remove
#   GR = dmrs2genomicranges(mat)
#   w_vec = width(GR)
#   ratio_mapped = sapply(1:length(GR), function(x) sum(width(intersect(GR[x], map))/w_vec[x]))
#   GR = GR[ratio_mapped > 1-pct.unmap]
#   w_vec = width(GR)
#   ratio_notblack = sapply(1:length(GR), function(x) sum(width(setdiff(GR[x], blacklist))/w_vec[x]))
#   GR = GR[ratio_notblack > 1-pct.black]
#   
#   # Trim - select biggest chunk
#   GR_list = lapply(1:length(GR), function(x) intersect(GR[x], map))
#   GR_list = lapply(1:length(GR_list), function(x) setdiff(GR_list[[x]], blacklist))
#   indices = sapply(GR_list, function(x) which.max(width(x)))
#   GR_list = lapply(1:length(indices), function(x) GR_list[[x]][indices[x]])
#   return(Reduce(c, GR_list))
# }
# 
# filter_blacklist_bg <- function(mat, blacklist, map)
# {
#   GR = dmrs2genomicranges(mat)
#   return(setdiff(intersect(mat, map), blacklist))
# }


filter_blacklist <- function(mat, blacklist, map, l_filter = 200)
{
  GR = dmrs2genomicranges(mat)
  GR = setdiff(intersect(GR, map), blacklist)
  w_vec = width(GR)
  return(GR[w_vec > l_filter])
}


# Read blacklist
setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/blacklist_regions/')
blacklist = dmrs2genomicranges(fread('hg38_blacklist_regions.bed', 
                                     col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')))
# Read mappability
setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/mappability_files/')
map = dmrs2genomicranges(fread('hg38_k100.bismap.bed', skip = 1, # skip descriptor
            col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')))
map = reduce(map)

#
min(dmrs.age_WB$length) # 201
min(dmrs.age_sperm$length) # 201
min(dmrs.tissue$length) # 201

# Filter blacklist and unmappability.
GR_dmrs.age_WB = filter_blacklist(dmrs.age_WB, blacklist, map, l_filter = 200)
GR_dmrs.age_sperm = filter_blacklist(dmrs.age_sperm, blacklist, map, l_filter = 200)
GR_dmrs.tissue = filter_blacklist(dmrs.tissue, blacklist, map, l_filter = 200)
GR_Bg.age_WB = filter_blacklist(Bg.age_WB, blacklist, map, l_filter = 200)
GR_Bg.age_sperm = filter_blacklist(Bg.age_sperm, blacklist, map, l_filter = 200)
GR_Bg.tissue = filter_blacklist(Bg.tissue, blacklist, map, l_filter = 200)

# Export Granges for enrichment analysis
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/GenomicRanges/')
saveRDS(GR_dmrs.age_WB, 'GR_dmrs.age_WB.Rds')
saveRDS(GR_dmrs.age_sperm, 'GR_dmrs.age_sperm.Rds')
saveRDS(GR_dmrs.tissue, 'GR_dmrs.tissue.Rds')
saveRDS(GR_Bg.age_WB, 'GR_Bg.age_WB.Rds')
saveRDS(GR_Bg.age_sperm, 'GR_Bg.age_sperm.Rds')
saveRDS(GR_Bg.tissue, 'GR_Bg.tissue.Rds')

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/DMRS/')
export_bed(dmrs.age_WB, 'dmrs.age_WB')
export_bed(dmrs.age_sperm, 'dmrs.age_sperm')
export_bed(dmrs.tissue, 'dmrs.tissue')


###################################### Visualization and annotation ######################################


dmrs.tissue$ID = paste(dmrs.tissue$chr, paste(dmrs.tissue$start, dmrs.tissue$end, sep = '-'), sep = ':')
dmrs.age_WB$ID = paste(dmrs.age_WB$chr, paste(dmrs.age_WB$start, dmrs.age_WB$end, sep = '-'), sep = ':')
dmrs.age_sperm$ID = paste(dmrs.age_sperm$chr, paste(dmrs.age_sperm$start, dmrs.age_sperm$end, sep = '-'), sep = ':')

########### TISSUE
# DIP2C
dmrs.tissue[dmrs.tissue$ID == 'chr10:668932-673161',] # entry 148
# chr  start    end length nCG meanMethy1 meanMethy2 diff.Methy  areaStat                  ID
# 20788 chr10 668932 673161   4230 170  0.1987938  0.9177825 -0.7189887 -4146.781 chr10:668932-673161




dmrs.tissue[41:60,] # all others are also associated with age in WB or SP
# chr6:166898539-166902781 RPS6KA2-AS1
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr6:166898539-166902781',], BSobj, ylim = c(0,1))
# chr5:1461889-1468856 LPCAT1
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr5:1461889-1468856',], BSobj, ylim = c(0,1))
# chr4:150577542-150584585 MAB21L2
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr4:150577542-150584585',], BSobj, ylim = c(0,1))
# chr7:152126642-152132817 KMT2C
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr7:152126642-152132817',], BSobj, ylim = c(0,1))
# chr15:40485494-40491817 CHST13
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr15:40485494-40491817',], BSobj, ylim = c(0,1))
# chr4:56312841-56317557 KIAA1211
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr4:56312841-56317557',], BSobj, ylim = c(0,1))
# chr7:27141671-27147605 HOXA6
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr7:27141671-27147605',], BSobj, ylim = c(0,1))
# chr16:1735172-1740292 MAPK8IP3
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr16:1735172-1740292',], BSobj, ylim = c(0,1))
# chr11:64368079-64373448 RPS6KA4
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr11:64368079-64373448',], BSobj, ylim = c(0,1))
# chr12:7711053-7715744 DPPA3
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr12:7711053-7715744',], BSobj, ylim = c(0,1))
# chr19:12157125-12162595 ZNF136, unmap?
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr19:12157125-12162595',], BSobj, ylim = c(0,1))
# chr19:811039-814180 PLPPR3
DSS::showOneDMR(dmrs.tissue[dmrs.tissue$ID == 'chr19:811039-814180',], BSobj, ylim = c(0,1))




########### Age - WB
head(dmrs.age_WB, 20)
# chr7:27130540-27131135 HOXA4
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr7:27130540-27131135',], BSobj, ylim = c(0,1)) # hypermethylation
# chr11:94545379-94545729 FUT4
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr11:94545379-94545729',], BSobj, ylim = c(0,1)) # hypomethylation
# # chr18:26547287-26547647 KCDT1
# DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr18:26547287-26547647',], BSobj, ylim = c(0,1)) # 
# chr16:2798335-2798795 PRSS41
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr16:2798335-2798795',], BSobj, ylim = c(0,1)) # 
# chr1:19665766-19666491 HTR6
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr1:19665766-19666491',], BSobj, ylim = c(0,1)) # 
# chr8:144524645-144524944 LRRC14
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr8:144524645-144524944',], BSobj, ylim = c(0,1)) # 
# chr6:37648744-37649214 MDGA1
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr6:37648744-37649214',], BSobj, ylim = c(0,1)) # 
# chr17:48576620-48577121 HOXB3
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr17:48576620-48577121',], BSobj, ylim = c(0,1)) # 
# chr5:136080591-136080956 VTRNA2.1
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr5:136080591-136080956',], BSobj, ylim = c(0,1)) # 
# chr5:1594489-1594837 SDHAP3
#DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr5:1594489-1594837',], BSobj, ylim = c(0,1)) # 
# chr19:51098886-51099204 CTU1
DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr19:51098886-51099204',], BSobj, ylim = c(0,1)) # 
# chr2:24174821-24175182 FAM228A
#DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr2:24174821-24175182',], BSobj, ylim = c(0,1)) # 
# chr17:7436504-7436866 TMEM102
# DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr17:7436504-7436866',], BSobj, ylim = c(0,1)) # 
# chr22:50529598-50529827 TYMP
# DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr22:50529598-50529827',], BSobj, ylim = c(0,1)) # 
# chr8:29353062-29353332 DUSP4
# DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr8:29353062-29353332',], BSobj, ylim = c(0,1)) # 
# chr10:132964778-132964979 LINC01166 (it is everything, all JRC and all pEWAS)
# DSS::showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr10:132964778-132964979',], BSobj, ylim = c(0,1)) # 

########### Age - SP

head(dmrs.age_sperm, 20)

# chr22:19726824-19728111 GP1BB
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr22:19726824-19728111',], BSobj, ylim = c(0,1)) # 
# chr22:19725206-19725817 SEPT5-GP1BB
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr22:19725206-19725817',], BSobj, ylim = c(0,1)) # 
# chr19:2429943-2431172 LMNB2
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr19:2429943-2431172',], BSobj, ylim = c(0,1)) # 
# chr7:23415058-23415667 IGF2BP3
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr7:23415058-23415667',], BSobj, ylim = c(0,1)) # 
# chr7:159773-160459 IG
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr7:159773-160459',], BSobj, ylim = c(0,1)) # 
# chr19:3492666-3493596 DOHH
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr19:3492666-3493596',], BSobj, ylim = c(0,1)) # 
# chr13:19769175-19769866 PSPC1
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr13:19769175-19769866',], BSobj, ylim = c(0,1)) # 
# chr3:52373064-52373766 DNAH1
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr3:52373064-52373766',], BSobj, ylim = c(0,1)) # 
# chr12:132717766-132718290 PGAM5
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr12:132717766-132718290',], BSobj, ylim = c(0,1)) # 
# chr19:3710655-3711103 TJP3
#DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr19:3710655-3711103',], BSobj, ylim = c(0,1)) # 
# chr3:43363927-43365099 ANO10
# DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr3:43363927-43365099',], BSobj, ylim = c(0,1)) # 
# chr17:48059296-48059723 NFE2L1
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr17:48059296-48059723',], BSobj, ylim = c(0,1)) # 
# chr3:197728732-197729503 RUBCN
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr3:197728732-197729503',], BSobj, ylim = c(0,1)) # 
# chr7:157465930-157467042 NR_110157
DSS::showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr7:157465930-157467042',], BSobj, ylim = c(0,1)) # 
# 


showOneDMR_mod = function(OneDMR, BSobj, ext = 500, ylim = c(0, 1), ...)
{
  allchr = as.character(seqnames(BSobj))
  allpos = start(BSobj)
  X = getBSseq(BSobj, "M")
  N = getBSseq(BSobj, "Cov")
  chr = as.character(OneDMR$chr)
  ix.chr = which(allchr == chr)
  thispos = allpos[ix.chr]
  thisN = N[ix.chr, ]
  thisX = X[ix.chr, ]
  xlim = c(OneDMR$start - ext, OneDMR$end + ext)
  ix1 = which(thispos <= xlim[2] & thispos >= xlim[1])
  nSample = ncol(X)
  if (nSample > 2) {
    y.cex = 0.2
  }
  else y.cex = 1
  sNames = sampleNames(BSobj)
  par(mfrow = c(nSample, 1), mar=c(1, 1, 1, 1), mgp=c(0.5, 0.2, 0), las=0)
  thisP = thisX/thisN
  maxN = max(thisN[ix1, ])
  
  
  # plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'Non-JRC; p = 0.5', cex.main = 0.4, 
  #         cex.axis = 0.1, cex.point = 0.8, cex.lab = 0.4, xaxt = 'n', yaxt = 'n')
  # axis(side = 1, at = seq(0, 600, 100), tck=-0.01, cex.axis = 0.3, mgp = c(-0.7, -0.4, 0))
  # axis(side = 2, at = seq(0, 50, 10), tck=-0.01, cex.axis = 0.3, las = 2, mgp = c(0.5, 0.2, 0))
  
  for (i in 1:ncol(X)) {
    plot(thispos[ix1], thisP[ix1, i], type = "h", col = "blue", 
         axes = F, lwd = 1.5, xlab = '', ylab = "", ylim = ylim, 
         xlim = xlim, main = sNames[i], cex.main = 0.4, cex.axis = 0.1, cex.lab = 0.4, ...)
    box(col = "black")
    axis(1, tck=-0.01, cex.axis = 0.3, mgp = c(-0.7, -0.4, 0))
    axis(2, col = "blue", col.axis = "blue", tck=-0.01, cex.axis = 0.3, las = 2, mgp = c(0.5, 0.2, 0))
    #mtext(chr, side = 1, line = 1.33, cex = y.cex)
    #mtext("methyl%", side = 2, line = 1.33, col = "blue", cex = y.cex)
    thisN.norm = thisN[ix1, i]/maxN * ylim[2]
    lines(thispos[ix1], thisN.norm, type = "l", col = "gray", lwd = 1.5)
    if (maxN >= 4) 
    {
      axis(side = 4, at = seq(0, ylim[2], length.out = 5), labels = round(seq(0, maxN, length.out = 5)), 
           tck=-0.01, cex.axis = 0.3, las = 2, mgp = c(0.5, 0.2, 0))
    }
    else 
    {
      axis(side = 4, at = seq(0, ylim[2], length.out = maxN + 1), labels = 0:maxN)
    }
    mtext("read depth", side = 4, line = 1.33, cex = y.cex)
    rect(OneDMR$start, ylim[1], OneDMR$end, ylim[2], col = "#FF00001A", border = NA)
  }
}


### Export high res

# setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/')
# dmrs.age_WB = readRDS('dmrs.age_WB.Rds')
# dim(dmrs.age_WB) # 135   9
# dmrs.age_sperm = readRDS('dmrs.age_sperm.Rds')
# dim(dmrs.age_sperm) # 1826    9
# dmrs.tissue = readRDS('dmrs.tissue.Rds')
# dim(dmrs.tissue) # 118864      9
# dmrs.tissue$ID = paste(dmrs.tissue$chr, paste(dmrs.tissue$start, dmrs.tissue$end, sep = '-'), sep = ':')
# dmrs.age_WB$ID = paste(dmrs.age_WB$chr, paste(dmrs.age_WB$start, dmrs.age_WB$end, sep = '-'), sep = ':')
# dmrs.age_sperm$ID = paste(dmrs.age_sperm$chr, paste(dmrs.age_sperm$start, dmrs.age_sperm$end, sep = '-'), sep = ':')

## HOXA4
setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/')
tiff(filename = 'HOXA4_HD.tiff', width = 5.3, height = 4, units = 'in', res = 300)
showOneDMR(dmrs.age_WB[dmrs.age_WB$ID == 'chr7:27130540-27131135',], BSobj, ylim = c(0,1)) # hypomethylation
dev.off()

## LMNB2
setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/')
tiff(filename = 'LMNB2_HD.tiff', width = 5.3, height = 4, units = 'in', res = 300)
showOneDMR(dmrs.age_sperm[dmrs.age_sperm$ID == 'chr19:2429943-2431172',], BSobj, ylim = c(0,1)) # 
dev.off()



### Visualize BAM

# chr11:94545379-94545729 FUT4
# samtools view -b output_38.bam "chr11:94545379-94545729" > FUT4_WB.bam
# samtools index FUT4_WB.bam
# cg18023065


# chr3:52373064-52373766 DNAH1
# samtools view -b output_38_sperm.bam "chr3:52373064-52373766" > DNAH1_SP.bam
# samtools index DNAH1_SP.bam
# cg26060971





