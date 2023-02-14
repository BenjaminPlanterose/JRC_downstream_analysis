library(data.table)
library(scales)

######################################### Read data #########################################

# WB
setwd("/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/")
IM_WB = fread('im_regions.txt', header = F)$V1
WB = fread('p_values.txt')
dim(WB) # 662659      2

# SP
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binokulars_output/sperm_results/')
IM_SP = fread('im_regions.txt', header = F)$V1
SP = fread('p_values.txt')
dim(SP) # 536691      2


# WB+Sp
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/binokulars_output/merged_results/')
IM_WBSP = fread('im_regions.txt', header = F)$V1
WB_SP = fread('p_values.txt')
dim(WB_SP) # 667071      2

######################################### Total IM coverage #########################################

diff_WB = sapply(strsplit(sapply(strsplit(IM_WB, split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))
diff_SP = sapply(strsplit(sapply(strsplit(IM_SP, split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))
diff_WBSP = sapply(strsplit(sapply(strsplit(IM_WBSP, split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))

plot(density(log10(diff_WBSP)), col = alpha('black', 0.5), ylim = c(0,7), xlim = c(0,6), main = 'All IM-regions')
lines(density(log10(diff_WB)), col = alpha('red2', 0.5))
lines(density(log10(diff_SP)), col = alpha('blue2', 0.5))
legend('topright', bty = 'n', fill = c('black', 'red2', 'blue2'), legend = c('WBSP', 'WB', 'SP'), cex = 1.2)

######################################### Significant #########################################

# NA p-value
sum(is.na(WB$V2)) # 32640
sum(is.na(SP$V2)) # 26599
sum(is.na(WB_SP$V2)) # 20224

# Filter out NA p-value
WB = WB[!is.na(WB$V2),]
SP = SP[!is.na(SP$V2),]
WB_SP = WB_SP[!is.na(WB_SP$V2),]

BT_WB = 0.05/nrow(WB) # 
BT_SP = 0.05/nrow(SP) # 
BT_WBSP = 0.05/nrow(WB_SP) # 

# Counts
sum(WB$V2 < BT_WB) # 4215
sum(SP$V2 < BT_SP) # 7386
sum(WB_SP$V2 < BT_WBSP) # 56896

# Distribution of length
diff_WB_sig = sapply(strsplit(sapply(strsplit(WB$V1[WB$V2 < BT_WB], split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))
diff_SP_sig = sapply(strsplit(sapply(strsplit(SP$V1[SP$V2 < BT_SP], split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))
diff_WBSP_sig = sapply(strsplit(sapply(strsplit(WB_SP$V1[WB_SP$V2 < BT_WBSP], split = ':'), function(x) x[2]), split = '-'), function(y) diff(as.numeric(y)))

plot(density(log10(diff_WBSP_sig)), col = alpha('black', 0.5), ylim = c(0,7), xlim = c(0,6), main = 'JRCs')
lines(density(log10(diff_WB_sig)), col = alpha('red2', 0.5))
lines(density(log10(diff_SP_sig)), col = alpha('blue2', 0.5))
legend('topright', bty = 'n', fill = c('black', 'red2', 'blue2'), legend = c('WBSP', 'WB', 'SP'), cex = 1.2)

######################################### Create beds #########################################

# WB
name = 'JRC_blood'
df = data.frame(regions = WB$V1[WB$V2 < BT_WB])
df$chr = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[1])
df$pos = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[2])
df$start = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[1]))
df$end = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[2]))
bed_df <- data.frame(seqname = df$chr,
                     start = df$start, end = df$end,
                     name = df$regions,
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/JRCs/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)

# Sperm
name = 'JRC_sperm'
df = data.frame(regions = SP$V1[SP$V2 < BT_SP])
df$chr = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[1])
df$pos = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[2])
df$start = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[1]))
df$end = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[2]))
bed_df <- data.frame(seqname = df$chr,
                     start = df$start, end = df$end,
                     name = df$regions,
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/JRCs/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)

# WB+Sperm
name = 'JRC_sperm_and_WB'
df = data.frame(regions = WB_SP$V1[WB_SP$V2 < BT_WBSP])
df$chr = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[1])
df$pos = sapply(strsplit(as.character(df$regions), split = ':'), function(x) x[2])
df$start = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[1]))
df$end = as.integer(sapply(strsplit(as.character(df$pos), split = '-'), function(x) x[2]))
bed_df <- data.frame(seqname = df$chr,
                     start = df$start, end = df$end,
                     name = df$regions,
                     score = 1000)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/JRCs/')
write.table(x = paste('track name=', name, ' description=', name, ' useScore=1', sep = ''), 
            file = paste(name, '.bed', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = bed_df, file = paste(name, '.bed', sep = ''), nThread = 4, 
       sep = '\t', col.names = F, append = T)

######################################### Examples #########################################

head(WB)

WB[which.min(WB$V2),]
# samtools view -b output_38.bam "chr13:48317600-48324000" > RB1_WB.bam
# samtools index RB1_WB.bam

SP = SP[order(SP$V2),]
head(SP, 20)
# samtools view -b output_38_sperm.bam "chr16:34586868-34587884" > IG_SP.bam
# samtools index IG_SP.bam

# samtools view -b output_38_sperm.bam "chr8:143921000-143926800" > PLEC_SP.bam
# samtools index PLEC_SP.bam


WB_SP = WB_SP[order(WB_SP$V2),]
head(WB_SP, 20)
# samtools view -b output_38_sperm.bam "chr10:668200-675600" > DIP2C_SP.bam
# samtools index DIP2C_SP.bam
# samtools view -b output_38.bam "chr10:668200-675600" > DIP2C_WB.bam
# samtools index DIP2C_WB.bam



head(WB)
WB = WB[order(WB$V2),]
head(WB, 12)

# 2:    chr20:58831400-58865400  0.000000e+00
# 3:    chr20:58887600-58891000  0.000000e+00
# 4:    chr21:7915746-7959500  0.000000e+00
# 5:    chr21:8208939-8209533  0.000000e+00
# 6:    chr21:8213243-8213722  0.000000e+00
# 7:    chr21:8401619-8401913  0.000000e+00
# 8:    chr21:8443214-8443438  0.000000e+00
# 9:    chr16:34586868-34587884 9.524765e-305
# 10:   chr2:206244800-206268270 1.281595e-299
# 11:   chr16:46398376-46398588 7.515659e-274
# 12:   chr7:130482800-130495400 6.209928e-247









