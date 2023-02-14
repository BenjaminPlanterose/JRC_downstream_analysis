############################################################################
############################################################################
###########                                                      ###########
###########               Collocalization analysis               ###########
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
library(GenomicRanges)
library(parallel)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(regioneR)

coordinates2datatable = function(coordinates)
{
  chr = sapply(strsplit(coordinates, split = ':'), function(x) x[1])
  pos = sapply(strsplit(coordinates, split = ':'), function(x) x[2])
  start = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[1]))
  end = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[2]))
  data.table(chr = chr, start = start, end = end, range = end - start)
}

datatable2bed <- function(mat, name, where)
{
  mat = mat[,-4]
  mat$regions = paste(mat$chr, paste(mat$start, mat$end, sep = '-'), sep = ':')
  bed_df <- data.frame(seqname = mat$chr,
                       start = as.integer(mat$start), end = as.integer(mat$end),
                       name = mat$regions,
                       score = 1000)
  setwd(where)
  write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
              file = paste(name, '.bed', sep = ''),
              quote = F, sep = '\t', row.names = F, col.names = F)
  fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
         sep = '\t', col.names = F, append = T)
}

datatable2grange = function(mat)
{
  GRanges(mat$chr, IRanges(start = mat$start, end = mat$end))
}

grange2datatable = function(GRANGE)
{
  data.table(chr = as.character(seqnames(GRANGE)), 
             start = start(GRANGE), end = end(GRANGE), range = end(GRANGE) - start(GRANGE))
}

pairwise_distance = function(T_list)
{
  mat = matrix(nrow = length(T_list), ncol = length(T_list))
  for(i in 1:length(T_list))
  {
    print(i)
    for(j in 1:length(T_list))
    {
      mat[i,j] = d(T_list[[i]], T_list[[j]])
    }
  }
  return(mat)
}

d <- function(A, B)
{
  sum(width(intersect(A, B)))
}

##################################### Test #####################################

# JRC: WB
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/')
WB = fread('p_values.txt')
WB = WB[!is.na(WB$V2),]; dim(WB) # 630019      2
BT = 0.05/nrow(WB)
sig_WB = WB$V1[WB$V2 < BT]; length(sig_WB) # 4215
A = coordinates2datatable(sig_WB); dim(A) # 4215    4
U_A = coordinates2datatable(WB$V1); dim(U_A) # 630019      4
A = datatable2grange(A)
U_A = datatable2grange(U_A)

# JRC: SP
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binokulars_output/sperm_results/')
SP = fread('p_values.txt')
SP = SP[!is.na(SP$V2),]; dim(SP) # 510092      2
BT = 0.05/nrow(SP)
sig_SP = SP$V1[SP$V2 < BT]; length(sig_SP) # 7386
B = coordinates2datatable(sig_SP); dim(B) # 7386    4
U_B = coordinates2datatable(SP$V1); dim(U_B) # 510092      4
B = datatable2grange(B)
U_B = datatable2grange(U_B)


# JRC: SP_WB
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/binokulars_output/merged_results/')
SP_WB = fread('p_values.txt')
SP_WB = SP_WB[!is.na(SP_WB$V2),]; dim(SP_WB) # 646847      2
BT = 0.05/nrow(SP_WB)
sig_SP_WB = SP_WB$V1[SP_WB$V2 < BT]; length(sig_SP_WB) # 56896
C = coordinates2datatable(sig_SP_WB); dim(C) # 56896    4
U_C = coordinates2datatable(SP_WB$V1); dim(U_C) # 646847      4
C = datatable2grange(C)
U_C = datatable2grange(U_C)

# EWAS
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/DMRs/GenomicRanges/')
D = readRDS('GR_dmrs.age_WB.Rds')
E = readRDS('GR_dmrs.age_sperm.Rds')
G = readRDS('GR_dmrs.tissue.Rds')
U_D = readRDS('GR_Bg.age_WB.Rds')
U_E = readRDS('GR_Bg.age_sperm.Rds')
U_G = readRDS('GR_Bg.tissue.Rds')


################################## Visualization ##################################

# Background network
distance_matrix0 = pairwise_distance(list(U_A, U_B, U_C, U_D, U_E, U_G))
rownames(distance_matrix0) = c('JRC_WB', 'JRC_SP', 'JRC_WB+SP', 'EWAS_A_WB', 'EWAS_A_SP', 'EWAS_T')
colnames(distance_matrix0) = c('JRC_WB', 'JRC_SP', 'JRC_WB+SP', 'EWAS_A_WB', 'EWAS_A_SP', 'EWAS_T')

# Target network
distance_matrix = pairwise_distance(list(A, B, C, D, E, G))
rownames(distance_matrix) = c('JRC_WB', 'JRC_SP', 'JRC_WB+SP', 'EWAS_A_WB', 'EWAS_A_SP', 'EWAS_T')
colnames(distance_matrix) = c('JRC_WB', 'JRC_SP', 'JRC_WB+SP', 'EWAS_A_WB', 'EWAS_A_SP', 'EWAS_T')

# Visualize Bg matrix
M = distance_matrix0/(1*10^9) # Scale down 
net = network(M, ignore.eval = FALSE, names.eval = "weights")
net %v% "id" = rownames(M)
net %v% "size" = diag(M)
y = c("firebrick", "lightpink2", "coral2", "cyan2", "blue4", "darkslategray4")
names(y) = rownames(M)
set.seed(39)
p1 = ggnet2(net, size = 'size', color = 'id', max_size = 20, shape = 19, palette = y,
            edge.size = "weights") + 
  guides(id = T, size = F) +
  theme(panel.background = element_rect(fill = "antiquewhite2")) +
  ggtitle('Background') + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 10)) + 
  theme(legend.position = "top", panel.grid = element_blank(), axis.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5)))

setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/3_pEWAS/')
ggsave(filename = 'Bg_network.tiff', plot = p1, device = 'tiff', width = 2.9, 
       height = 3.3, dpi = 600, units = 'in', scale = 2.5)

# Visualize target matrix
M = distance_matrix/(5*10^6) # Scale down 
net = network(M, ignore.eval = FALSE, names.eval = "weights")
net %v% "id" = rownames(M)
net %v% "size" = diag(M)
y = c("firebrick", "lightpink2", "coral2", "cyan2", "blue4", "darkslategray4")
names(y) = rownames(M)
set.seed(39)
p2 = ggnet2(net, size = 'size', color = 'id', max_size = 20, shape = 19, palette = y,
       edge.size = "weights") + 
  guides(id = T, size = F) +
  theme(panel.background = element_rect(fill = "antiquewhite2")) +
  ggtitle('Target') + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 10)) + 
  theme(legend.position = "top", panel.grid = element_blank(), axis.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p2

setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/3_pEWAS/')
ggsave(filename = 'Tg_network.tiff', plot = p2, device = 'tiff', width = 2.9, 
       height = 3.3, dpi = 600, units = 'in', scale = 2.5)



# Assume: Male samples, thus no X-inactivation

# JRC_WB: cell-type, mQTL, imprinting, 
# JRC_WB+SP: imprinting + tissue-specific + JRC_WB + JRC_SP
# JRC_SP: mQTL, polymorphic imprinting, 

# EWAS SP_Vs_WB: imprinting + tissue-specific
# EWAS age-WB: CellComp-dependent changes with age or cellComp-independent
# EWAS age-SP: cellComp-independent


# library(NetworkDistance)
# library(BioNetStat)

################################## Testing for significance ##################################

# C Vs A+B+D+E+G
Joint = Reduce(union, list(A, B, D, E, G))
pt1 <- overlapPermTest(A=C, B=Joint, ntimes=10, alternative = 'greater', 
                       genome = 'BSgenome.Hsapiens.UCSC.hg38.masked', 
                       non.overlapping = F, verbose = T, parallel = 4); plot(pt1) # 39.282
lz1 <- localZScore(A=C, B=Joint, pt=pt1, window=100*10^3, 
                   step=1*10^3); plot(lz1)
z0 = pt1$numOverlaps$zscore # 697.6119


# Colocalization
# Statistic colocalized length
# Permutation test where windows are wiggled and check for colocalization statistic.
# https://academic.oup.com/bioinformatics/article/35/9/1615/5126923
# Binning restricts degrees of freedom!!!!!!

# 1. Can JRC_seeker substitute an EWAS for tissue?
# YES

# 2. Can JRC_seeker substitute an EWAS for age - WB?
# YES

# 3. Can JRC_seeker substitute an EWAS for age - Sperm? (does not include age-dependent changes in cell composition)
# YES

# 4. Can JRC_seeker substitute an pASM - WB? (blood mQTL)

# 5. Can JRC_seeker substitute an pASM - SP? (sperm mQTL)

#### OBS: a lot of sharing between pASM in blood and sperm. WEIRD

# 6. colocalization spermWB and WB





