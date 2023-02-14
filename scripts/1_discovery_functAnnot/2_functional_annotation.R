############################################################################
############################################################################
###########                                                      ###########
###########              JRC Genome Wide Discovery               ###########
###########             Authors: Bronte Kolar,                   ###########
###########                      Benjamin Planterose             ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########                                                      ###########
############################################################################
############################################################################

# Load libraries
library(data.table)
library(tidyr)
library(grid)
library(gridGraphics)
library(qqman)
library(ggplot2)
library(gridExtra)
library(gplots)
library(dplyr)
library(gprofiler2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(seqinr)

# Load functions
process.file <- function(filename)
{
  df = fread(filename)
  colnames(df) = c("region", "P")
  print(paste('number of NA p-values = ', sum(is.na(df$P))))
  df = na.omit(df)
  df$P = df$P + .Machine$double.xmin
  df$reg = df$region
  df = separate(data = df, col = region, into = c("chr", "range"), sep = ":")
  df = separate(data = df, col = range, into = c("region_start", "region_end"), sep = "\\-")
  df$region_start = as.integer(df$region_start)
  df$region_end = as.integer(df$region_end)
  return(df)
}

create_manhattan_plot = function(df, sample_name)
{
  df$SNP = df$reg
  p_val_cutoff = 0.05/nrow(df)
  print(paste('number of JRC = ', sum(df$P < p_val_cutoff)))
  df$BP = rowMeans(df[,c('region_start', 'region_end')], na.rm=TRUE)
  df$CHR = substring(df$chr, 4)
  df$CHR[df$CHR == "X"] = "23"
  df$CHR[df$CHR == "Y"] = "24"
  df$CHR = as.numeric(df$CHR)
  
  manhattan_data = df[, c("SNP", "CHR", "BP", "P")]
  
  chr_freq = as.data.frame(table(manhattan_data$CHR))
  
  tiff(filename = paste(sample_name, "manhattan_plot.tiff", sep = " "), width = 10, height = 10, units = 'in', res = 300)
  manhattan(manhattan_data, main = paste(sample_name, "Manhattan Plot - All Binpolish Regions", sep=" "), ylim=c(0,350), cex=0.6, cex.axis=0.7,
            col = c("cadetblue3", "coral1"), suggestiveline = -log10(p_val_cutoff), genomewideline = F, chrlabs = c(1:22, "X","Y"))
  dev.off()
}

create_qq_plot = function(pval, sample_name)
{
  tiff(filename = paste(sample_name, "qq_plot.tiff", sep = " "), width = 10, height = 10, units = 'in', res = 300)
  qq_plot = qq(pval, xlim = c(0, 10), ylim = c(0, 330), main = paste(sample_name, "Q-Q Plot - All Binpolish Regions", sep=" "))
  dev.off()
  
  tiff(filename = paste(sample_name, "pval_dist.tiff", sep = " "), width = 10, height = 10, units = 'in', res = 300)
  hist(x = pval, xlim = c(0,1), breaks = 100, main = paste(sample_name, "P-val distribution - All Binpolish Regions", sep=" "))
  abline(h = length(pval)/100, lty = 2)
  dev.off()
}

segmentation_qc = function(sample, state_labels, binpolish_segments, chromhmm_segments)
{
  states_df = fread(state_labels)
  chromhmm_df = fread(chromhmm_segments)
  binpolish_df = fread(binpolish_segments)
  
  chromhmm_df$V4 = factor(chromhmm_df$V4, levels = c('E1', 'E2', 'E3', 'E4'))
  levels(chromhmm_df$V4) = states_df$label[match(levels(chromhmm_df$V4), states_df$state)]
  chromhmm_df$V4 = as.character(chromhmm_df$V4)
  
  
  chromhmm_df$segment_length = chromhmm_df$V3 - chromhmm_df$V2
  binpolish_df$segment_length = binpolish_df$V3 - binpolish_df$V2
  binpolish_df = binpolish_df[binpolish_df$segment_length > 0,]
  
  states = unique(chromhmm_df$V4)
  print(states)
  max_log1 = max(log10(chromhmm_df$segment_length))
  max_log2 = max(log10(binpolish_df$segment_length))
  max_log = ceiling(max(max_log1, max_log2))
  
  for(state in states) 
  {
    #label = states_df$label[which(states_df$state == state)]
    c_state_df = chromhmm_df[which(chromhmm_df$V4 == state), ]
    b_state_df = binpolish_df[which(binpolish_df$V4 == state), ]
    p1 = ggplot(c_state_df, aes(x=log10(segment_length))) + xlab("log10(segment size)") + ggtitle(paste(sample, state, "Regions - ChromHMM - Segment Size Histogram", sep=" ")) + geom_histogram(bins = 50, color="darkblue", fill="lightblue") + scale_x_continuous(limits = c(0, max_log), oob = scales::squish)
    p2 = ggplot(b_state_df, aes(x=log10(segment_length))) + xlab("log10(segment size)") + ggtitle(paste(sample, state, "Regions - BinPolish - Segment Size Histogram", sep=" ")) + geom_histogram(bins = 50, color="darkgreen", fill="lightgreen") + scale_x_continuous(limits = c(0, max_log), oob = scales::squish)
    y_max = ceiling(max(max(ggplot_build(p1)$data[[1]]$count), max(ggplot_build(p2)$data[[1]]$count)))
    
    #ggsave(plot = p3, device = 'tiff', filename = paste(sample, state, "segmentation_plot.tiff", sep = " "), width = 10, height = 10, units = 'in', res = 300)
    tiff(filename = paste(sample, state, "segmentation_plot.tiff", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
    grid.arrange(p1, p2, nrow = 2)
    dev.off()
  }  
}

jrc_chr_dist = function(df, sample, colour1, colour2)
{
  df$CHR = substring(df$chr, 4)
  chr_freq = table(df$CHR)
  Names = names(chr_freq)
  Counts = unname(chr_freq)
  chr_freq = data.frame(CHR = factor(Names, levels = c(1:22, 'X', 'Y')), n = Counts)
  p1 = ggplot(chr_freq, aes(x=CHR, y=n.Freq)) + xlab("chromosome") + ylab("frequency") + geom_bar(stat="identity", color=colour1, fill=colour2) + ggtitle(paste(sample, "- JRC Chromosome Distribution", sep=" ")) 
  return(p1)
}

jrc_chr_dist_sig = function(df, sample, colour1, colour2)
{
  thrs = 0.05/nrow(df)
  df = df[df$P < thrs,]
  df$CHR = substring(df$chr, 4)
  chr_freq = table(df$CHR)
  Names = names(chr_freq)
  Counts = unname(chr_freq)
  chr_freq = data.frame(CHR = factor(Names, levels = c(1:22, 'X', 'Y')), n = Counts)
  p1 = ggplot(chr_freq, aes(x=CHR, y=n.Freq)) + xlab("chromosome") + ylab("frequency") + geom_bar(stat="identity", color=colour1, fill=colour2) + ggtitle(paste(sample, "- JRC Chromosome Distribution", sep=" ")) 
  return(p1)
}

width <- function(df)
{
  thrs = 0.05/nrow(df)
  width = df$region_end-df$region_start
  width_sig = width[df$P < thrs]
  return(list(width = width, width_sig = width_sig))
}

enrichment.per.category <- function(mat)
{
  pvals <- numeric(length = ncol(mat))
  odd_ratio <- numeric(length = ncol(mat))
  names(pvals) <- colnames(mat)
  names(odd_ratio) <- colnames(mat)
  
  for(i in 1:ncol(mat))
  {
    category <- colnames(mat)[i]
    data <- cbind(mat[,i], rowSums(mat[,-i]))
    rownames(data) <- c("unsig_jrc", "sig_jrc")
    colnames(data) <- c(category, paste('Not', category))
    message(category)
    print(data)
    res <- fisher.test(data)
    print(res)
    pvals[i] <- res$p.value
    odd_ratio[i] <-res$estimate
  }
  output <- list(pvals, odd_ratio)
  names(output) <- c('pvals', 'Odd_ratios')
  return(output)
}

get_types = function(x)
{
  if (is.na(x) | grepl('NA', x, fixed = TRUE))
  {
    "intergenic"
  } else {
    if (grepl('promoter', x, fixed = TRUE))
    {
      "promoter"
    } else if (grepl('exon', x, fixed = TRUE))
    {
      "exon"
    } else {
      "intron"
    }
  }
}

refseq_distribution = function(df, sample, refseq_annotation_file, seed = 1, B = 100000)
{
  p_val_cutoff = 0.05/nrow(df)
  
  #######################################################
  # Load annotation file
  #######################################################
  
  annotation = fread(refseq_annotation_file)
  colnames(annotation) = c("gene_name", "gene_name_2", "chrom", "strand", "start", "end", "type")
  annotation = annotation[with(annotation, order(chrom, start)),]
  
  #######################################################
  # Find intersection
  #######################################################
  
  setkey(annotation, chrom, start, end)
  print(sample)
  result_data = foverlaps(df, annotation, c("chr", "region_start", "region_end"), c("chrom", "start", "end"))
  result_indices = foverlaps(df, annotation, c("chr", "region_start", "region_end"), c("chrom", "start", "end"), which = TRUE, mult = "all")
  result_data = na.omit(result_data)
  
  #######################################################
  # Classify regions
  #######################################################
  
  # Hierarchy: promoter > exon > intron > intergenic
  
  result_indices$type = annotation[result_indices$yid, "type"]
  classify = aggregate(type ~ xid, data=result_indices, na.action=na.pass, FUN = function(x) paste0(x,collapse = ', '))
  types_list = lapply(1:(nrow(classify)), function(x) get_types(classify$type[x]))
  df$type = unlist(types_list)
  
  #######################################################
  # Enrichment Analysis
  #######################################################
  
  sig_jrc_regions = subset(df, P < p_val_cutoff)
  unsig_jrc_regions = subset(df, P >= p_val_cutoff)
  
  sig_region_dist = table(sig_jrc_regions$type)
  unsig_region_dist = table(unsig_jrc_regions$type)
  
  # Build contingency table
  merged_data = do.call("rbind", list(unsig_region_dist, sig_region_dist))
  rownames(merged_data) = c("unsig_jrc", "sig_jrc")
  merged_data = as.matrix(merged_data)
  
  # Calculate global enrichment
  message('global')
  set.seed(seed)
  print(fisher.test(merged_data, simulate.p.value = T, B = B))
  
  # Per category enrichment
  res1 <- enrichment.per.category(merged_data)
  
  #######################################################
  # Enrichment Analysis
  #######################################################
  
  plot_values = log10(1/res1$Odd_ratios)
  p_values = p.adjust(res1$pvals, method = 'bonferroni')
  
  colours = vector(mode="character", length=4)
  p_val_text = vector(mode="character", length=4)
  p_val_position = vector(mode = "numeric", length=4)
  epsilon = 0.02
  for (i in 1:length(plot_values))
  {
    if (plot_values[i] < 0) {
      colours[i] = 'brown3'
        p_val_position[i] = epsilon
    } else {
      colours[i] = 'chartreuse4'
        p_val_position[i] = -epsilon
    }
    
    if (p_values[i] < 0.001) {
      p_val_text = "***"
    } else if (p_values[i] < 0.01) {
      p_val_text = "**"
    } else if (p_values[i] < 0.05) {
      p_val_text = "*"
    } else {
      p_val_text = ""
    }
  }
  
  tiff(filename = paste(sample, "jrc_refseq_barplot.tiff", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
  #y = barplot(1/res1$Odd_ratios-1, horiz = T, xaxt = 'n', xlab = paste(sample, 'Target (JRC) to Background (IM region) Odds ratio', sep=" "), las = 1,
  #            col = colours, xlim = c(-1, 14))
  #axis(side = 1, labels = c(0.25, 0.5, 0.75, 1:8), at = c(0.25, 0.5, 0.75, 1:8)-1)
  par(mar=c(5.1, 6.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  y = barplot(log10(1/res1$Odd_ratios), horiz = T, xlab = paste(sample, 'Target (JRC) to Background (IM region) log10(OR)', sep=" "), las = 1,
             col = colours)
  abline(v = 0, lty = 1)
  text(x = p_val_position, y = y, p_val_text)
  dev.off()
  
  mat1 = do.call("rbind", list(unsig_region_dist, sig_region_dist))
  mat1[1,] <- mat1[1,]/sum(mat1[1,])
  mat1[2,] <- mat1[2,]/sum(mat1[2,])
  rownames(mat1) <- c('unsig_jrc', 'sig_jrc')
  mat1 <- mat1[, ncol(mat1):1]
  tiff(filename = paste(sample, "jrc_refseq_dist_balloon.tiff", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
  balloonplot(as.table(mat1), main = '', ylab = 'Region Status', xlab = 'Set')
  dev.off()
}

GO_enrichment <- function(df, sample, annotation_file, seed = 1)
{
  get_index = function(x)
  {
    if (is.element('promoter', x))
    {
      index = which(x %in% c("promoter"))
    } 
    else if (is.element('exon', x))
    {
      index = which(x %in% c("exon"))
    } 
    else 
    {
      index = which(x %in% c("intron"))
    }
    return(index)
  }
  
  get_types2 = function(x)
  {
    if(is.element('promoter', x))
    {
      type = "promoter"
    } 
    else if(is.element('exon', x))
    {
      type = "exon"
    } 
    else 
    {
      type = "intron"
    }
    return(type)
  }
  
  get_gene_names = function(x, y)
  {
    unique(x[[1]][y])
  }
  
  get_gene_2_names = function(x, y)
  {
    unique(y[[1]][x])
  }
  
  p_val_cutoff = 0.05/nrow(df)
  #get_types
  # Load annotation
  annotation = fread(annotation_file)
  colnames(annotation) = c("gene_name", "gene_name_2", "chrom", "strand", "start", "end", "type")
  annotation = annotation[order(annotation$chrom, annotation$start),]
  
  # Find intersection
  setkey(annotation, chrom, start, end)
  result_data = foverlaps(df, annotation, c("chr", "region_start", "region_end"), c("chrom", "start", "end"))
  result_indices = foverlaps(df, annotation, c("chr", "region_start", "region_end"), c("chrom", "start", "end"), which = TRUE, mult = "all")
  
  # Remove unmapped regions
  result_data = na.omit(result_data)
  result_indices = na.omit(result_indices)
  
  # Find annotation
  df$gene_name = NA
  df$gene_name_2 = NA
  non_intergenic = unique(result_indices$xid)
  
  # removed strand
  testing = result_data %>% dplyr::group_by(chr, region_start, region_end, P) %>% 
    dplyr::summarise(gene_name = list(gene_name), gene_name_2 = list(gene_name_2), type = list(type))
  
  index_list = lapply(1:(nrow(testing)), function(i) get_index(unlist(testing$type[i])))
  types_list = lapply(1:(nrow(testing)), function(i) get_types2(unlist(testing$type[i])))
  genes_list = lapply(1:(nrow(testing)), function(i) get_gene_names(testing$gene_name[i], unlist(testing$gene_name[i])))
  genes_2_list = lapply(1:(nrow(testing)), function(i) get_gene_2_names(unlist(index_list[[i]]),  testing$gene_name_2[i]))
  
  types_df = as.data.frame(do.call(rbind, types_list))
  colnames(types_df) = c("types")
  genes_df = data.frame(gene_name=sapply(genes_list, toString), stringsAsFactors=FALSE)
  genes_2_df = data.frame(gene_name_2=sapply(genes_2_list, toString), stringsAsFactors=FALSE)
  
  testing$type_fin = types_df$types
  testing$gene_name_fin = genes_df$gene_name
  testing$gene_name_2_fin = genes_2_df$gene_name_2
  final_merge = testing[, c("chr", "region_start", "region_end", "P", "type_fin", "gene_name_fin", "gene_name_2_fin")]
  result_dataframe = merge(x = df, y = final_merge, by = c("chr", "region_start", "region_end", "P"), all.x = TRUE)
  
  #######################################################
  # GO - gProfiler2
  #######################################################
  
  ## create input dataframe with candidate genes and background genes (all)
  
  sig_regions = subset(result_dataframe, P < p_val_cutoff)
  sig_genes = na.omit(sig_regions$gene_name_2_fin)
  test_bg_all = na.omit(result_dataframe$gene_name_2_fin)
  pos = match(sig_genes, test_bg_all)
  
  set.seed(seed) # If several genes associated to an IM region, select one randomly.
  test_bg_all = sapply(strsplit(test_bg_all, ', '), function(x) x[sample(1:length(x), 1)])
  sig_genes = test_bg_all[pos]
  
  res_gost = gost(query = sig_genes, correction_method = "gSCS", custom_bg = test_bg_all, 
                  organism = "hsapiens", domain_scope = 'custom_annotated') # Matters custom vs custom_annotated
  
  
  #######################################################
  # GO - clusterProfiler
  #######################################################
  
  # sig_regions = subset(result_dataframe, P < p_val_cutoff)
  # bg_regions = subset(result_dataframe, P >= p_val_cutoff)
  # sig_genes = na.omit(sig_regions$gene_name_2_fin)
  # sig_genes = unique(unlist(strsplit(sig_genes,", ")))
  # bg_genes = na.omit(bg_regions$gene_name_2_fin)
  # bg_genes = unique(unlist(strsplit(bg_genes,", ")))
  # test_bg_all = na.omit(result_dataframe$gene_name_2_fin)
  # test_bg_all = unique(unlist(strsplit(test_bg_all,", ")))
  # 
  # hs = org.Hs.eg.db
  # ego = enrichGO(gene          = sig_genes,
  #                universe      = test_bg_all,
  #                OrgDb         = org.Hs.eg.db,
  #                ont           = "all",
  #                pAdjustMethod = "BH",
  #                pvalueCutoff  = 0.01,
  #                qvalueCutoff  = 0.05,
  #                keyType = "SYMBOL")
  # 
  #tiff(filename = paste(sample, "clusterProfiler_barplot.tiff", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
  #barplot(wbresgos$ego, showCategory=20) # Can also do barplot: dotplot(ego, showCategory = 30)
  #dev.off()
  
  #return(list(result_dataframe = result_dataframe, res_gost = res_gost, ego = ego))
  return(list(result_dataframe = result_dataframe, res_gost = res_gost))
}

# 0. Binpolish QC
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/2_JRC_seeker/")
segmentation_qc('WB_hg38', 
                state_labels = '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binpolish/assets/state_labels.txt', 
                binpolish_segments = '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binpolish/polished_segmentation.bed', 
                chromhmm_segments = '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/chromhmm/output_files/pooled_WB_v2_4_segments.bed')

segmentation_qc('SP_hg38', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binpolish/assets/state_labels.txt', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binpolish/polished_segmentation.bed', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/chromhmm/output_files/pooled_sperm_v2_4_segments.bed')

segmentation_qc('WBSP_hg38', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/binpolish/assets/state_labels.txt', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/binpolish/polished_segmentation.bed', 
                '/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/chromhmm/output_files/merged_pool_v2_4_segments.bed')

# 1. Read p-values for all IM-regions
WB = process.file("/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/p_values.txt")
# number of NA p-values =  32640
dim(WB) # 630019      5

WB[WB$reg == 'chr12:629400-633600',]

SP = process.file("/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/sperm_pool/binokulars_output/sperm_results/p_values.txt")
# number of NA p-values =  26599
dim(SP) # 510092      5

SP[SP$reg == 'chr8:143921000-143926800',]


WB_SP = process.file("/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/merged_pool/binokulars_output/merged_results//p_values.txt")
# number of NA p-values =  20224
dim(WB_SP) # 646847      5

# 2. Create Manhattan plots
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/2_JRC_seeker/")
create_manhattan_plot(WB, "WB_hg38")
create_manhattan_plot(SP, "SP_hg38")
create_manhattan_plot(WB_SP, "WBSP_hg38")

# Tissue  n_JRC
# WB      4215
# SP      7386
# WB+SP   56896

# 3. QQ-plots+pvalDist
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/2_JRC_seeker/")
create_qq_plot(WB$P, 'WB_hg38')
create_qq_plot(SP$P, 'SP_hg38')
create_qq_plot(WB_SP$P, 'WBSP_hg38')

# 4. IM chromosome distribution
tiff(filename = paste('IM', "jrc_chr_distribution.png", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
p1 = jrc_chr_dist(WB, 'WB_hg38', "darkblue", "lightblue")
p2 = jrc_chr_dist(SP, 'SP_hg38', "darkgreen", "lightgreen")
p3 = jrc_chr_dist(WB_SP, 'WB+SP_hg38', "red4", "lightcoral")
grid.arrange(p1, p2, p3, nrow = 3)
dev.off()

# 5. JRC chromosome distribution
tiff(filename = paste('JRC', "jrc_chr_distribution.png", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
p1 = jrc_chr_dist_sig(WB, 'WB_hg38', "darkblue", "lightblue")
p2 = jrc_chr_dist_sig(SP, 'SP_hg38', "darkgreen", "lightgreen")
p3 = jrc_chr_dist_sig(WB_SP, 'WB+SP_hg38', "red4", "lightcoral")
grid.arrange(p1, p2, p3, nrow = 3)
dev.off()

# 6. Size distribution
width_list = lapply(X = list(WB, SP, WB_SP), FUN = width)
width_dfs = lapply(width_list, melt)
df = Reduce(rbind, width_dfs)
NROW = sapply(width_dfs, nrow); df$tissue = unlist(lapply(1:length(NROW), function(x) rep(c('WB', 'SP', 'WB_SP')[x], times = NROW[x])))
colnames(df) = c('Size', 'Set', 'Tissue')
df$Set = factor(df$Set); levels(df$Set) = c('IM', 'JRC')
tiff(filename = paste('size', "jrc_size_compare_dist.tiff", sep = "_"), width = 10, height = 10, units = 'in', res = 300)
ggplot(data = df, mapping = aes(x = Tissue, y = log10(Size), fill = Set)) + geom_violin(trim=FALSE) + ylim(0, 6) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
dev.off()

# 7. RefSeq distribution
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/2_JRC_seeker/")
refseq_distribution(WB, 'WB_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/5_functional_annotation/2_annotation/assets/hg38_refseq_annotation.tsv')
# global
# Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)
# data:  merged_data
# p-value = 1e-05
# alternative hypothesis: two.sided
# 
# exon
#            exon Not exon
# unsig_jrc 66379   559425
# sig_jrc     898     3317
#
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4068155 0.4725442
# sample estimates:
#   odds ratio 
# 0.4382926 
# 
# intergenic
# intergenic Not intergenic
# unsig_jrc     329399         296405
# sig_jrc         1092           3123
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.964940 3.409017
# sample estimates:
# odds ratio 
#   3.178229 
# 
# intron
#           intron Not intron
# unsig_jrc 203863     421941
# sig_jrc      717       3498
# 
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.174234 2.558309
# sample estimates:
#   odds ratio 
# 2.357147 
# 
# promoter
# promoter Not promoter
# unsig_jrc    26163       599641
# sig_jrc       1508         2707
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.07341021 0.08356482
# sample estimates:
# odds ratio 
#  0.0783253 

refseq_distribution(SP, 'SP_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/5_functional_annotation/2_annotation/assets/hg38_refseq_annotation.tsv')
# global
# Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)
# data:  merged_data
# p-value = 1e-05
# alternative hypothesis: two.sided
# 
# exon
#            exon Not exon
# unsig_jrc 43694   459012
# sig_jrc     854     6532
# 
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6773720 0.7833663
# sample estimates:
#   odds ratio 
# 0.7280956 
# 
# intergenic
# intergenic Not intergenic
# unsig_jrc     277388         225318
# sig_jrc         4085           3301
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value = 0.832
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9495938 1.0421565
# sample estimates:
# odds ratio 
#  0.9948216 
# 
# intron
#           intron Not intron
# unsig_jrc 160825     341881
# sig_jrc     2045       5341
# 
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value = 1.676e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.166899 1.293942
# sample estimates:
#   odds ratio 
# 1.228566 
# 
# promoter
# promoter Not promoter
# unsig_jrc    20799       481907
# sig_jrc        402         6984
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value = 8.678e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.6773286 0.8319710
# sample estimates:
# odds ratio 
#  0.7498221 

refseq_distribution(WB_SP, 'WBSP_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/5_functional_annotation/2_annotation/assets/hg38_refseq_annotation.tsv')
# global
# Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)
# data:  merged_data
# p-value = 1e-05
# alternative hypothesis: two.sided
# 
# exon
#            exon Not exon
# unsig_jrc 52373   537578
# sig_jrc    6877    50019
# 
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6898667 0.7279273
# sample estimates:
#   odds ratio 
# 0.7086035 
# 
# intergenic
# intergenic Not intergenic
# unsig_jrc     321676         268275
# sig_jrc        28177          28719
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.201223 1.243364
# sample estimates:
# odds ratio 
#   1.222112 
# 
# intron
#           intron Not intron
# unsig_jrc 197596     392355
# sig_jrc    12052      44844
# 
# 	Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.835244 1.913586
# sample estimates:
#   odds ratio 
# 1.873886 
# 
# promoter
# promoter Not promoter
# unsig_jrc    18306       571645
# sig_jrc       9790        47106
# 
# Fisher's Exact Test for Count Data
# data:  data
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.1500999 0.1581921
# sample estimates:
# odds ratio 
#  0.1540834 


# 7. db enrichment analysis (gprofiler2)

# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/")
wbresgos = GO_enrichment(WB, 'WB_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/1_discovery_functAnnot/0_assets/hg38_refseq_annotation.tsv')
100*mean(is.na(wbresgos$result_dataframe$gene_name_2_fin)) # 52.45731 % unmapped
gostplot(wbresgos$res_gost, capped = F, interactive = T) # SET INTERACTIVE TO TRUE IF YOU WANT TO 
dim(wbresgos$result_dataframe) # 630019     10
dim(wbresgos$res_gost$result) # 856  14
saveRDS(wbresgos, 'wbresgos.Rds')

spresgos = GO_enrichment(SP, 'SP_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/1_discovery_functAnnot/0_assets/hg38_refseq_annotation.tsv')
100*mean(is.na(spresgos$result_dataframe$gene_name_2_fin)) # 55.18083 % unmapped
gostplot(spresgos$res_gost, capped = F, interactive = T) # SET INTERACTIVE TO TRUE IF YOU WANT TO 
dim(spresgos$result_dataframe) # 510092     10
dim(spresgos$res_gost$result) # 94  14
saveRDS(spresgos, 'spresgos.Rds')

wbspresgos = GO_enrichment(WB_SP, 'WBSP_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/1_discovery_functAnnot/0_assets/hg38_refseq_annotation.tsv')
100*mean(is.na(wbspresgos$result_dataframe$gene_name_2_fin)) # 54.0859 % unmapped
gostplot(wbspresgos$res_gost, capped = FALSE, interactive = T) # SET INTERACTIVE TO TRUE IF YOU WANT TO 
dim(wbspresgos$result_dataframe) # 646847     10
dim(wbspresgos$res_gost$result) # 1010  14
saveRDS(wbspresgos, 'wbspresgos.Rds')

# Random
setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/2_JRC_seeker/")
random = WB; random$P = random$P[sample(1:nrow(WB), nrow(WB))]
randomgo = GO_enrichment(random, 'random_hg38', '/media/ultron/2tb_disk1/ben/JRC_project/scripts/1_discovery_functAnnot/0_assets/hg38_refseq_annotation.tsv')
gostplot(randomgo$res_gost, capped = F, interactive = T) # SET INTERACTIVE TO TRUE IF YOU WANT TO 
saveRDS(randomgo, 'randomgo.Rds')


# 8. Imprinting enrichment
library(data.table)
setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/imprinting/')
imprinted.genes = fread('imprinting.txt', header = F)
#imprinted.genes = imprinted.genes[imprinted.genes$V4 %in% c('Imprinted', 'Predicted'),]
imprinted.genes = imprinted.genes[imprinted.genes$V4 %in% c('Imprinted'),]
imprinted.genes = unique(c(imprinted.genes$V1, unlist(strsplit(imprinted.genes$V2, split = ', '))))

setwd("/media/ultron/2tb_disk1/ben/JRC_project/figures/tmp/2_JRC_seeker/GO/")
wbresgos = readRDS('wbresgos.Rds')$result_dataframe
dim(wbresgos) # 630019     10
spresgos = readRDS('spresgos.Rds')$result_dataframe
dim(spresgos) # 510092     10
wbspresgos = readRDS('wbspresgos.Rds')$result_dataframe
dim(wbspresgos) # 646847     10

enrichment_imprinting = function(annot, imprt.set)
{
  alpha_thr = 0.05/nrow(annot)
  sig = annot$P < alpha_thr
  split_genes = strsplit(annot$gene_name_2_fin, split =', ')
  impr = sapply(1:nrow(annot), function(x) length(intersect(split_genes[[x]], imprt.set)) >= 1)
  mat = table(sig, impr)[c('FALSE', 'TRUE'), c('FALSE', 'TRUE')]
  print(mat)
  print(fisher.test(mat))
  print(paste(annot$gene_name_2_fin[sig & impr], annot$reg[sig & impr]))
}



enrichment_imprinting(wbresgos, imprinted.genes)
#         impr
# sig      FALSE   TRUE
# FALSE 623559   2245
# TRUE    4159     56
# 
# Fisher's Exact Test for Count Data
# 
# data:  mat
# p-value = 7.089e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.811056 4.884266
# sample estimates:
# odds ratio 
#   3.739983 
# 
#  [1] "TP73 chr1:3655000-3670800"                          "TP73 chr1:3694200-3721400"                         
#  [3] "DIRAS3 chr1:68046200-68048600"                      "DIRAS3 chr1:68049200-68052600"                     
#  [5] "MIR675, H19 chr11:1997000-2005624"                  "MIR483, IGF2 chr11:2131000-2137600"                
#  [7] "IGF2, INS-IGF2, INS chr11:2145400-2163000"          "KCNQ1OT1 chr11:2696800-2706200"                    
#  [9] "KCNQ1 chr11:2777200-2783800"                        "KCNQ1 chr11:2790800-2793000"                       
# [11] "KCNQ1DN chr11:2867400-2874000"                      "OSBPL5 chr11:3147200-3164600"                      
# [13] "WT1 chr11:32425200-32427400"                        "WT1 chr11:32427600-32433800"                       
# [15] "ANO1 chr11:70140000-70151800"                       "NTM chr11:131684600-131700400"                     
# [17] "RB1 chr13:48317600-48324000"                        "DLK1 chr14:100722000-100731600"                    
# [19] "MEG3 chr14:100822800-100828600"                     "MIR4508, MKRN3 chr15:23551200-23568000"            
# [21] "MAGEL2 chr15:23646600-23655000"                     "NDN chr15:23683000-23697938"                       
# [23] "SNRPN, SNHG14 chr15:24837800-24849400"              "SNURF, SNRPN chr15:24953600-24959600"              
# [25] "ATP10A chr15:25858000-25862800"                     "ATP10A chr15:25863400-25864800"                    
# [27] "CHD2, RGMA chr15:93024600-93038800"                 "GNG13, PRR25 chr16:798400-809000"                  
# [29] "ZNF597 chr16:3430000-3433000"                       "PARD6G chr18:80157000-80162200"                    
# [31] "PEG3, MIMT1, ZIM2 chr19:56836400-56849400"          "NNAT, BLCAP chr20:37518400-37522800"               
# [33] "L3MBTL1 chr20:43512600-43516400"                    "GNAS, GNAS-AS1, GNAS, GNAS chr20:58831400-58865400"
# [35] "GNAS, LOC101927932 chr20:58887600-58891000"         "NAP1L5 chr4:88696400-88699083"                     
# [37] "VTRNA2-1 chr5:136067400-136083800"                  "FAM50B chr6:3846800-3852800"                       
# [39] "CRYBG1 chr6:106507800-106510400"                    "PLAGL1, HYMAI chr6:144006000-144009813"            
# [41] "HOXA4 chr7:27129600-27135000"                       "GRB10 chr7:50779200-50785800"                      
# [43] "PEG10, PEG10, PEG10, SGCE chr7:94655600-94658262"   "PEG10, PEG10, PEG10 chr7:94658297-94659400"        
# [45] "MEST, MESTIT1, MIR335 chr7:130482800-130495400"     "SVOPL chr7:138661400-138668600"                    
# [47] "DLGAP2 chr8:865800-871200"                          "DLGAP2 chr8:891200-895400"                         
# [49] "DLGAP2 chr8:943800-950800"                          "DLGAP2 chr8:1128200-1145200"                       
# [51] "DLGAP2 chr8:1157800-1167600"                        "DLGAP2 chr8:1371000-1375280"                       
# [53] "DLGAP2 chr8:1548200-1549600"                        "DLGAP2 chr8:1669600-1671600"                       
# [55] "DLGAP2 chr8:1701000-1702600"                        "PEG13 chr8:140094600-140102200"    



jkl = c("TP73", "TP73", "DIRAS3", "DIRAS3", "H19", "IGF2", "IGF2", "KCNQ1OT1",
  "KCNQ1", "KCNQ1", "KCNQ1DN", "OSBPL5", "WT1", "WT1", "ANO1", "NTM",                     
  "RB1", "DLK1", "MEG3", "MKRN3", "MAGEL2", "NDN", "SNRPN", "SNRPN", "ATP10A", "ATP10A",                    
  "CHD2", "GNG13, PRR25", "ZNF597", "PARD6G", "PEG3", "NNAT, BLCAP",
  "L3MBTL1", "GNAS", "GNAS", "NAP1L5", "VTRNA2-1", "FAM50B", "CRYBG1", "PLAGL1, HYMAI",            
  "HOXA4", "GRB10", "PEG10", "PEG10", "MEST", "SVOPL", "DLGAP2", "PEG13")    
length(unique(jkl)) # 39


# IGF2  chr11:2145400-2163000
# GNAS  chr20:58831400-58865400
# L3MBTL1 chr20:43512600-43516400
# NAP1L5  chr4:88696400-88699083
# RB1 chr13:48317600-48324000
# DIRAS3  chr1:68049200-68052600

# samtools view -b output_38.bam "chr11:2145400-2163000" > IGF2_WB.bam
# samtools index IGF2_WB.bam
# samtools view -b output_38.bam "chr20:58831400-58865400" > GNAS_WB.bam
# samtools index GNAS_WB.bam
# samtools view -b output_38.bam "chr20:43512600-43516400" > L3MBTL1_WB.bam
# samtools index L3MBTL1_WB.bam
# samtools view -b output_38.bam "chr4:88696400-88699083" > NAP1L5_WB.bam
# samtools index NAP1L5_WB.bam
# samtools view -b output_38.bam "chr1:68049200-68052600" > RB1_WB.bam
# samtools index RB1_WB.bam
# samtools view -b output_38.bam "chr13:48317600-48324000" > DIRAS3_WB.bam
# samtools index DIRAS3_WB.bam




enrichment_imprinting(spresgos, imprinted.genes)
# impr
# sig      FALSE   TRUE
# FALSE 500948   1758
# TRUE    7356     30
# 
# Fisher's Exact Test for Count Data
# 
# data:  mat
# p-value = 0.4261
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7811214 1.6665036
# sample estimates:
# odds ratio 
#   1.162124 
# 
#  [1] "KCNQ1 chr11:2468800-2470400"                 "NTM chr11:131579000-131586600"              
#  [3] "NTM chr11:131928000-131930600"               "ST8SIA1 chr12:22258400-22259634"            
#  [5] "ST8SIA1 chr12:22259672-22260400"             "ST8SIA1 chr12:22320971-22323000"            
#  [7] "ATP10A chr15:25712600-25715200"              "ATP10A chr15:25716800-25719400"             
#  [9] "ATP10A chr15:25725800-25728200"              "NLRP2 chr19:54962787-54964306"              
# [11] "GPR1, GPR1-AS chr2:206209784-206211000"      "GPR1 chr2:206212400-206218600"              
# [13] "DSCAM chr21:40574600-40582400"               "DSCAM chr21:40789765-40793600"              
# [15] "LOC102724770, DGCR6 chr22:18905438-18906000" "ERAP2 chr5:96915800-96918800"               
# [17] "ADTRP chr6:11769200-11772400"                "LIN28B chr6:104991686-104992600"            
# [19] "DDC chr7:50457800-50458996"                  "DDC chr7:50459189-50459491"                 
# [21] "DDC chr7:50459648-50460541"                  "MAGI2 chr7:78034200-78036600"               
# [23] "MAGI2 chr7:79390000-79393000"                "CPA4 chr7:130305400-130309400"              
# [25] "DLGAP2 chr8:747600-750600"                   "DLGAP2 chr8:973600-975000"                  
# [27] "DLGAP2 chr8:1071000-1087000"                 "DLGAP2 chr8:1358400-1362400"                
# [29] "DLGAP2 chr8:1668400-1671600"                 "DLGAP2 chr8:1673800-1677800"       

enrichment_imprinting(wbspresgos, imprinted.genes)
# impr
# sig      FALSE   TRUE
# FALSE 587787   2164
# TRUE   56640    256
# 
# Fisher's Exact Test for Count Data
# 
# data:  mat
# p-value = 0.002505
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.073955 1.398530
# sample estimates:
# odds ratio 
#   1.227662 
# [1] "TP73 chr1:3654800-3659200"                                  
# [2] "TP73 chr1:3676200-3680000"                                  
# [3] "TP73 chr1:3693400-3697400"                                  
# [4] "TP73 chr1:3699200-3705200"                                  
# [5] "TP73 chr1:3706800-3710800"                                  
# [6] "TP73 chr1:3715400-3721000"                                  
# [7] "RNU5F-1, RNU5D-1 chr1:44720400-44732800"                    
# [8] "DIRAS3 chr1:68041200-68048600"                              
# [9] "DIRAS3 chr1:68048800-68054200"                              
# [10] "MIR675, H19 chr11:1997000-2005624"                          
# [11] "IGF2 chr11:2144200-2153800"                                 
# [12] "KCNQ1 chr11:2445200-2448600"                                
# [13] "KCNQ1 chr11:2466600-2473800"                                
# [14] "KCNQ1 chr11:2479800-2483200"                                
# [15] "KCNQ1 chr11:2524400-2528200"                                
# [16] "KCNQ1 chr11:2530800-2535200"                                
# [17] "KCNQ1OT1 chr11:2697400-2702600"                             
# [18] "KCNQ1 chr11:2777400-2782600"                                
# [19] "KCNQ1 chr11:2789800-2795200"                                
# [20] "KCNQ1 chr11:2811200-2814400"                                
# [21] "KCNQ1 chr11:2822200-2829400"                                
# [22] "KCNQ1DN chr11:2867000-2868800"                              
# [23] "KCNQ1DN chr11:2869400-2872800"                              
# [24] "CDKN1C chr11:2881600-2884400"                               
# [25] "CDKN1C chr11:2885800-2889400"                               
# [26] "SLC22A18 chr11:2908000-2911200"                             
# [27] "SLC22A18 chr11:2913200-2915400"                             
# [28] "OSBPL5 chr11:3147600-3152400"                               
# [29] "OSBPL5 chr11:3152600-3162000"                               
# [30] "OSBPL5 chr11:3165400-3172200"                               
# [31] "WT1 chr11:32398800-32401600"                                
# [32] "WT1 chr11:32411000-32413800"                                
# [33] "WT1 chr11:32425200-32427200"                                
# [34] "WT1 chr11:32428600-32433600"                                
# [35] "WT1-AS chr11:32439400-32449200"                             
# [36] "ANO1 chr11:70110200-70112800"                               
# [37] "ANO1 chr11:70140000-70151000"                               
# [38] "ANO1 chr11:70164200-70168000"                               
# [39] "ANO1 chr11:70176000-70180800"                               
# [40] "ZC3H12C chr11:110144400-110149800"                          
# [41] "ZC3H12C chr11:110163400-110167000"                          
# [42] "NTM chr11:131369800-131377400"                              
# [43] "NTM chr11:131577400-131586600"                              
# [44] "NTM chr11:131625800-131629000"                              
# [45] "NTM chr11:131679000-131683800"                              
# [46] "NTM chr11:131685400-131697000"                              
# [47] "NTM chr11:131744600-131750800"                              
# [48] "NTM chr11:131759000-131766800"                              
# [49] "NTM chr11:131783200-131789498"                              
# [50] "NTM chr11:131789577-131791000"                              
# [51] "NTM chr11:131812600-131817600"                              
# [52] "NTM chr11:131836600-131839000"                              
# [53] "NTM chr11:131928000-131930600"                              
# [54] "NTM chr11:132122200-132129800"                              
# [55] "NTM chr11:132159600-132164000"                              
# [56] "NTM chr11:132235200-132239400"                              
# [57] "NTM chr11:132244600-132250200"                              
# [58] "ST8SIA1 chr12:22258400-22259634"                            
# [59] "ST8SIA1 chr12:22259672-22260400"                            
# [60] "ST8SIA1 chr12:22267600-22269400"                            
# [61] "HNF1A chr12:120977400-120984000"                            
# [62] "RB1 chr13:48314800-48324200"                                
# [63] "RB1 chr13:48456200-48460400"                                
# [64] "ESR2 chr14:64257000-64262200"                               
# [65] "ESR2 chr14:64308400-64311400"                               
# [66] "ESR2 chr14:64338800-64341600"                               
# [67] "SMOC1 chr14:69875600-69879600"                              
# [68] "DLK1 chr14:100724600-100731800"                             
# [69] "MEG3 chr14:100822800-100831000"                             
# [70] "RTL1 chr14:100903200-100907200"                             
# [71] "DIO3, DIO3OS, MIR1247 chr14:101561000-101563200"            
# [72] "MIR4508, MKRN3 chr15:23559400-23573200"                     
# [73] "MAGEL2 chr15:23644400-23655000"                             
# [74] "NPAP1 chr15:24662400-24679400"                              
# [75] "SNRPN chr15:24845000-24859600"                              
# [76] "SNRPN, SNHG14 chr15:24876000-24883200"                      
# [77] "SNHG14, SNRPN chr15:24915400-24923200"                      
# [78] "SNURF, SNRPN chr15:24953200-24959800"                       
# [79] "SNORD109B, SNORD109A chr15:25273000-25278373"               
# [80] "ATP10A chr15:25712600-25715200"                             
# [81] "ATP10A chr15:25716800-25719400"                             
# [82] "ATP10A chr15:25726000-25730000"                             
# [83] "ATP10A chr15:25778999-25779600"                             
# [84] "ATP10A chr15:25794600-25801400"                             
# [85] "ATP10A chr15:25809000-25816800"                             
# [86] "ATP10A chr15:25858400-25862800"                             
# [87] "ATP10A chr15:25863400-25865400"                             
# [88] "RASGRF1 chr15:79034600-79037800"                            
# [89] "CHD2 chr15:92910400-92924800"                               
# [90] "CHD2, RGMA chr15:93024600-93038400"                         
# [91] "GNG13, PRR25 chr16:798600-804800"                           
# [92] "PRR25 chr16:805600-809400"                                  
# [93] "ZNF597 chr16:3429600-3438600"                               
# [94] "NAA60, ZNF597 chr16:3441200-3445400"                        
# [95] "CMTM1 chr16:66564600-66567467"                              
# [96] "ZFP90 chr16:68537600-68539400"                              
# [97] "TP53 chr17:7679000-7683200"                                 
# [98] "ZNF396 chr18:35375000-35376800"                             
# [99] "ZNF396 chr18:35377000-35379400"                             
# [100] "PARD6G chr18:80156800-80162400"                             
# [101] "PARD6G chr18:80184400-80187200"                             
# [102] "PARD6G chr18:80247400-80250200"                             
# [103] "NLRP2 chr19:54962787-54964306"                              
# [104] "NLRP2 chr19:54964422-54964647"                              
# [105] "NLRP2 chr19:54965075-54965414"                              
# [106] "NLRP2 chr19:54965424-54965737"                              
# [107] "PEG3, MIMT1, ZIM2 chr19:56834600-56847400"                  
# [108] "GPR1, GPR1-AS chr2:206209784-206211000"                     
# [109] "GPR1 chr2:206212000-206221000"                              
# [110] "NNAT, BLCAP chr20:37517800-37523800"                        
# [111] "L3MBTL1 chr20:43508000-43510000"                            
# [112] "L3MBTL1 chr20:43512600-43517000"                            
# [113] "HNF4A chr20:44351200-44358600"                              
# [114] "HNF4A chr20:44368800-44371800"                              
# [115] "GNAS, GNAS-AS1, GNAS, GNAS chr20:58837000-58859200"         
# [116] "GNAS, LOC101927932 chr20:58886400-58890800"                 
# [117] "DSCAM chr21:40009200-40017200"                              
# [118] "DSCAM chr21:40030400-40033705"                              
# [119] "DSCAM chr21:40136600-40147600"                              
# [120] "DSCAM chr21:40166800-40169600"                              
# [121] "DSCAM chr21:40175200-40180200"                              
# [122] "DSCAM chr21:40335200-40340200"                              
# [123] "DSCAM chr21:40399000-40403200"                              
# [124] "DSCAM chr21:40485449-40486200"                              
# [125] "DSCAM chr21:40511000-40515200"                              
# [126] "DSCAM chr21:40553800-40560800"                              
# [127] "DSCAM chr21:40575000-40582400"                              
# [128] "DSCAM chr21:40703400-40707000"                              
# [129] "DSCAM chr21:40779000-40783400"                              
# [130] "DSCAM chr21:40789765-40799000"                              
# [131] "DSCAM chr21:40818034-40819800"                              
# [132] "DSCAM chr21:40833200-40843200"                              
# [133] "DSCAM chr21:40844600-40849200"                              
# [134] "NAP1L5 chr4:88695800-88699083"                              
# [135] "NAP1L5 chr4:88699138-88700600"                              
# [136] "NFKB1 chr4:102498400-102501200"                             
# [137] "ERAP2 chr5:96875200-96879000"                               
# [138] "ERAP2 chr5:96911000-96913400"                               
# [139] "ERAP2 chr5:96916000-96918600"                               
# [140] "VTRNA2-1 chr5:136078000-136083800"                          
# [141] "PXDC1 chr6:3731000-3733800"                                 
# [142] "PXDC1 chr6:3741200-3745000"                                 
# [143] "PXDC1 chr6:3745200-3751000"                                 
# [144] "FAM50B chr6:3847000-3852800"                                
# [145] "ADTRP chr6:11726200-11733200"                               
# [146] "ADTRP chr6:11769200-11772400"                               
# [147] "CRYBG1 chr6:106358200-106360800"                            
# [148] "CRYBG1 chr6:106361000-106363050"                            
# [149] "CRYBG1 chr6:106363173-106364000"                            
# [150] "CRYBG1 chr6:106435200-106443400"                            
# [151] "CRYBG1 chr6:106459800-106463600"                            
# [152] "CRYBG1 chr6:106481182-106483200"                            
# [153] "CRYBG1 chr6:106507800-106510200"                            
# [154] "PLAGL1, HYMAI chr6:144005200-144009813"                     
# [155] "RAC1 chr7:6372000-6374400"                                  
# [156] "RAC1 chr7:6375400-6379200"                                  
# [157] "RAC1 chr7:6395400-6397200"                                  
# [158] "HOXA4, HOXA-AS3, HOXA3, HOXA5, HOXA6 chr7:27128000-27148000"
# [159] "GLI3 chr7:41960800-41968000"                                
# [160] "HECW1 chr7:43193600-43200600"                               
# [161] "HECW1 chr7:43260400-43266000"                               
# [162] "HECW1 chr7:43272600-43275600"                               
# [163] "HECW1 chr7:43309800-43315800"                               
# [164] "HECW1 chr7:43443000-43447000"                               
# [165] "HECW1 chr7:43535600-43540400"                               
# [166] "DDC chr7:50459648-50460541"                                 
# [167] "DDC chr7:50467200-50469200"                                 
# [168] "DDC chr7:50490600-50498200"                                 
# [169] "DDC chr7:50526200-50530000"                                 
# [170] "GRB10 chr7:50659400-50668000"                               
# [171] "GRB10 chr7:50779200-50784400"                               
# [172] "GRB10 chr7:50787200-50792400"                               
# [173] "GRB10 chr7:50793400-50797000"                               
# [174] "MAGI2 chr7:78016200-78019600"                               
# [175] "MAGI2 chr7:78034200-78036800"                               
# [176] "MAGI2 chr7:78039000-78043200"                               
# [177] "MAGI2 chr7:78104360-78105600"                               
# [178] "MAGI2 chr7:78109191-78113000"                               
# [179] "MAGI2 chr7:78118000-78120213"                               
# [180] "MAGI2 chr7:78360200-78365400"                               
# [181] "MAGI2 chr7:78510000-78514000"                               
# [182] "MAGI2 chr7:78536321-78538000"                               
# [183] "MAGI2 chr7:78611400-78612848"                               
# [184] "MAGI2 chr7:78631400-78634600"                               
# [185] "MAGI2 chr7:78635800-78637800"                               
# [186] "MAGI2 chr7:78654200-78655786"                               
# [187] "MAGI2 chr7:78722400-78728655"                               
# [188] "MAGI2 chr7:78736000-78742000"                               
# [189] "MAGI2 chr7:78802570-78804630"                               
# [190] "MAGI2 chr7:78804673-78806800"                               
# [191] "MAGI2 chr7:78822400-78823390"                               
# [192] "MAGI2 chr7:78825200-78827600"                               
# [193] "MAGI2 chr7:78867400-78872400"                               
# [194] "MAGI2 chr7:78903406-78904600"                               
# [195] "MAGI2 chr7:78906800-78910800"                               
# [196] "MAGI2 chr7:78968200-78970600"                               
# [197] "MAGI2 chr7:79074400-79079000"                               
# [198] "MAGI2 chr7:79110600-79113000"                               
# [199] "MAGI2 chr7:79214200-79226200"                               
# [200] "MAGI2 chr7:79293000-79299800"                               
# [201] "MAGI2 chr7:79329200-79332800"                               
# [202] "MAGI2 chr7:79389800-79393000"                               
# [203] "MAGI2 chr7:79432000-79436600"                               
# [204] "TFPI2-DT, TFPI2 chr7:93890600-93894200"                     
# [205] "PEG10, PEG10, PEG10, SGCE chr7:94655800-94658262"           
# [206] "PEG10, PEG10, PEG10 chr7:94658297-94659519"                 
# [207] "PEG10, PEG10, PEG10 chr7:94662000-94669400"                 
# [208] "PPP1R9A chr7:94919600-94922800"                             
# [209] "PPP1R9A chr7:95044400-95048000"                             
# [210] "CPA4 chr7:130305400-130309200"                              
# [211] "MEST chr7:130482800-130486600"                              
# [212] "MEST, MESTIT1, MIR335 chr7:130486800-130495400"             
# [213] "KLF14 chr7:130730600-130734000"                             
# [214] "SVOPL chr7:138646200-138649800"                             
# [215] "SVOPL chr7:138661400-138668600"                             
# [216] "DLGAP2, LOC401442 chr8:735000-737200"                       
# [217] "DLGAP2 chr8:747600-750600"                                  
# [218] "DLGAP2 chr8:753000-755800"                                  
# [219] "DLGAP2 chr8:803000-807325"                                  
# [220] "DLGAP2 chr8:820200-823000"                                  
# [221] "DLGAP2 chr8:823600-826200"                                  
# [222] "DLGAP2 chr8:832200-840000"                                  
# [223] "DLGAP2 chr8:865800-871200"                                  
# [224] "DLGAP2 chr8:891800-894800"                                  
# [225] "DLGAP2 chr8:899800-904800"                                  
# [226] "DLGAP2 chr8:945200-952000"                                  
# [227] "DLGAP2 chr8:1053400-1055200"                                
# [228] "DLGAP2 chr8:1068200-1087600"                                
# [229] "DLGAP2 chr8:1104200-1109200"                                
# [230] "DLGAP2 chr8:1109934-1113800"                                
# [231] "DLGAP2 chr8:1141400-1145600"                                
# [232] "DLGAP2 chr8:1153200-1157800"                                
# [233] "DLGAP2 chr8:1159400-1161200"                                
# [234] "DLGAP2 chr8:1180988-1182600"                                
# [235] "DLGAP2 chr8:1198800-1202600"                                
# [236] "DLGAP2 chr8:1357200-1362200"                                
# [237] "DLGAP2 chr8:1371000-1375280"                                
# [238] "DLGAP2 chr8:1415400-1420600"                                
# [239] "DLGAP2 chr8:1435400-1440600"                                
# [240] "DLGAP2 chr8:1475600-1477400"                                
# [241] "DLGAP2 chr8:1481200-1484400"                                
# [242] "DLGAP2 chr8:1490400-1494800"                                
# [243] "DLGAP2 chr8:1503800-1510800"                                
# [244] "DLGAP2 chr8:1547200-1550800"                                
# [245] "DLGAP2 chr8:1668400-1671600"                                
# [246] "DLGAP2 chr8:1672000-1676200"                                
# [247] "DLGAP2 chr8:1687200-1691888"                                
# [248] "DLGAP2 chr8:1692378-1693400"                                
# [249] "DLGAP2 chr8:1699400-1704400"                                
# [250] "ZFAT chr8:134708600-134712800"                              
# [251] "KCNK9 chr8:139617000-139620800"                             
# [252] "KCNK9 chr8:139628800-139632600"                             
# [253] "KCNK9 chr8:139661200-139664600"                             
# [254] "PEG13 chr8:140095800-140102200"                             
# [255] "GLIS3 chr9:4115400-4121200"                                 
# [256] "GLIS3 chr9:4299600-4301400"        


## Polymorphic imprinting

enrichment_imprinting(wbresgos, c('DUSP22', 'VTRNA2-1', 'CYP2E1', 'PAX8-AS1', 'CFD'))
# [1] "CYP2E1 chr10:133523800-133534200"  "CFD chr19:857357-864400"           "PAX8-AS1 chr2:113235979-113236184"
# [4] "VTRNA2-1 chr5:136067400-136083800" "DUSP22 chr6:284200-300400"         "DUSP22 chr6:329600-342800"  

wbresgos[wbresgos$reg == 'chr6:284200-300400',]


# samtools view -b output_38.bam "chr10:133523800-133534200" > CYP2E1_WB.bam
# samtools index CYP2E1_WB.bam
# samtools view -b output_38.bam "chr19:857357-864400" > CFD_WB.bam
# samtools index CFD_WB.bam
# samtools view -b output_38.bam "chr2:113235979-113236184" > PAX8_WB.bam
# samtools index PAX8_WB.bam
# samtools view -b output_38.bam "chr5:136067400-136083800" > VTRNA2-1_WB.bam
# samtools index VTRNA2-1_WB.bam
# samtools view -b output_38.bam "chr6:284200-300400" > DUSP22_WB.bam
# samtools index DUSP22_WB.bam



# 9. Cell type composition
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
FlowSorted.Blood.450k.compTable$gene_name = annot[rownames(FlowSorted.Blood.450k.compTable),]$UCSC_RefGene_Name
FlowSorted.Blood.450k.compTable$chr = annot[rownames(FlowSorted.Blood.450k.compTable),]$chr
FlowSorted.Blood.450k.compTable = FlowSorted.Blood.450k.compTable[order(FlowSorted.Blood.450k.compTable$p.value),]

FlowSorted.Blood.450k.compTable[81:100,]


###
### Candidates
# cg04153882, WIPI2, chr7:5229800-5236200 (B-cell marker)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg04153882 2373.206 5.209077e-38 0.9499032 0.9563241 0.9048152 0.10240016 0.9184050 0.9507915 0.07339106 0.9711910

# cg25344401, FOXK1, chr7:4712000-4718200 (myeloid vs lymphoid)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg25344401 1494.073 5.255948e-35 0.9435842 0.9100685 0.9072202 0.92666774 0.11459038 0.16313355 0.08134981 0.9696060
data("FlowSorted.Blood.450k")
beta = getBeta(preprocessQuantile(FlowSorted.Blood.450k))

df = data.frame(beta = beta['cg25344401',]) # neighbours: cg18276112, cg15200418, cg04572930, cg26136772
df = data.frame(beta = beta['cg04572930',])
df$type = sapply(strsplit(colnames(beta), split = '_'), function(x) x[1])
df$col = df$type
df$col[df$col %in% c('Eos', 'Gran', 'Neu', 'CD14+')] = 'myeloid'
df$col[df$col %in% c('CD8+', 'CD4+', 'CD19+', 'CD56+')] = 'lymphoid'
df$col[df$col %in% c('WB', 'PBMC')] = 'mixture'
df$type = factor(df$type, levels = c('Eos', 'Gran', 'Neu', 'CD14+', 'CD8+', 'CD4+', 'CD19+', 'CD56+', 'WB', 'PBMC'))
library(ggplot2)
ggplot(data = df, mapping = aes(x = type, y = beta, col = col, fill = col)) + geom_boxplot() + ylim(0,1)


# cg17356733, IFNGR2, chr21:33399000-33403200 (myeloid vs lymphoid)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg17356733 1465.057 7.044402e-35 0.9360442 0.8860312 0.9009631 0.92735760 0.12817475 0.18032913 0.10036738 0.9475569

# cg03463948, UBASH3B, chr11:122740600-122747400 (myeloid vs lymphoid)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg03463948 1371.332 1.890326e-34 0.9139463 0.84209207 0.8794507 0.9028312 0.22603134 0.1974208 0.15508847 0.9253750

# cg15945333, SPI1, chr11:47371800-47378800 (myeloid vs lymphoid, B cell weird)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg15945333 1289.718 4.723794e-34 0.8375102 0.76445686 0.7717993 0.1417399 0.09675307 0.1498660 0.07962407 0.8632462

# cg13738327, LRP5, chr11:68369800-68375400 (B-cell marker)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg13738327 1214.981 1.151174e-33 0.9090035 0.9645673 0.8823188 0.09531174 0.97041521 0.9752966 0.07847517 0.9878111

# Targeted search for cell type specific markers
cell.type.specific <- function(FlowSorted.Blood.450k.compTable, cell_type, thr)
{
  types = c('CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran')
  meth = FlowSorted.Blood.450k.compTable[,types]
  meth$beta_a = rowMeans(meth[, !(types %in% cell_type)])
  
  if(length(cell_type) == 1)
  {
    meth$delta_square = (meth[,cell_type] - meth$beta_a)^2
  }
  else
  {
    meth$delta_square = (rowMeans(meth[,cell_type]) - meth$beta_a)^2
  }
  
  meth[meth$delta_square > thr,]
}

##### GRAN
Gran = cell.type.specific(FlowSorted.Blood.450k.compTable, 'Gran', 0.5)
dim(Gran)
head(Gran)
FlowSorted.Blood.450k.compTable['cg00439981',]

# cg04855678, SLC51A, chr3:196218200-196225400, 
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg04855678 913.2906 8.109486e-32 0.9324268 0.8932976 0.856004 0.8855972 0.4934109 0.07966842 0.05410017 0.9489897 0.8948896

# cg05398700, WDR20, chr14:102209600-102211800
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg05398700 896.6842 1.065796e-31 0.9796279 0.9798998 0.9813682 0.9651722 0.8785039 0.1080905 0.07900473 0.9886974 0.9096926

# cg25074794, MARCH8, chr10:45462830-45469673
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg25074794 785.3921 7.660586e-31 0.9407197 0.9273895 0.9182589 0.9133046 0.7667072 0.1757772 0.1270191 0.9535871 0.826568

# cg00660167, CSNK1D, chr17:82231600-82244000
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg00660167 763.2621 1.171894e-30 0.9733051 0.9680825 0.9664978 0.9746241 0.7256866 0.1086736 0.0774894 0.9935668 0.9160774

# cg16802439, GALNS, chr16:88836800-88842400
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg16802439 714.9393 3.098854e-30 0.9880443 0.9573347 0.964681 0.9741468 0.7092254 0.1546377 0.1127298 0.9948815 0.8821517

# cg08253808, WDR20, chr14:102209600-102211800
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg08253808 600.1794 4.165751e-29 0.9445466 0.9234356 0.9321285 0.9154405 0.7070534 0.1312166 0.08301834 0.9612943 0.8782759

# cg25605731, CALR, chr19:12939400-12948400
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg25605731 467.0178 1.713281e-27 0.953266 0.9439523 0.9381523 0.847237 0.4673408 0.07924584 0.03401039 0.9746165 0.9406061

# cg11070172, CEBPE, chr14:23108000-23122600
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg11070172 456.1831 2.424357e-27 0.9450691 0.9729939 0.9520934 0.9598095 0.8337401 0.1925944 0.1532681 0.982217 0.8289489

# cg13468144, MIR5096, chr17:4174600-4180000
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg13468144 417.9207 8.847334e-27 0.9344142 0.8722759 0.9227842 0.9313514 0.8865968 0.1863357 0.1293392 0.9645769 0.8352377

# cg09993145, RUNX3, chr1:24961600-24967800
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg09993145 409.4216 1.198304e-26 0.01529253 0.04134432 0.02914457 0.1073545 0.09327158 0.8073967 0.01124458 0.8522232

# cg00439981, FAM120B, chr6:170373400-170377065
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg00439981 404.7549 1.419257e-26 0.9613244 0.9256736 0.9234732 0.9548502 0.5778122 0.119916 0.05140111 0.9730139 0.9216128

##### Mono
Mono = cell.type.specific(FlowSorted.Blood.450k.compTable, 'Mono', 0.45)
FlowSorted.Blood.450k.compTable['cg25898577',]

# cg25898577, PPM1F, chr22:21935600-21954592
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg25898577 323.7457 3.809688e-25 0.8801805 0.8271134 0.8099808 0.8916123 0.1481245 0.7319802 0.08465875 0.9249095 0.8402507

##### Tcell
Tcell = cell.type.specific(FlowSorted.Blood.450k.compTable, c('CD4T', 'CD8T'), 0.45)
dim(Tcell)
FlowSorted.Blood.450k.compTable['cg07015803',]

# cg16452866, BCL11B, chr14:99185332-99246000
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg16452866 380.9034 3.477007e-26 0.1867661 0.1367108 0.7607488 0.9616786 0.9766029 0.9454502 0.09751938 0.9941536 0.8966342
# cg07015803, BCL11B, chr14:99185332-99246000
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg07015803 281.7582 2.927577e-24 0.1933961 0.145032 0.7125506 0.9481633 0.9589778 0.922007 0.1037235 0.9678475 0.864124
# cg02380585, SLFN13, chr17:35447000-35450800
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg02380585 336.6483 2.144235e-25 0.8216421 0.8154984 0.2860515 0.05555116 0.05643259 0.1033196 0.03548165 0.8851005


##### NK
NK = cell.type.specific(FlowSorted.Blood.450k.compTable, 'NK', 0.4)
dim(NK)
FlowSorted.Blood.450k.compTable['cg15956469',]

# cg15956469, KLRD1, chr12:10300600-10314800
#            Fstat      p.value      CD8T      CD4T       NK     Bcell      Mono       Gran        low      high     range
# cg15956469 102.8556 5.990472e-18 0.6829011 0.9646968 0.259083 0.9740762 0.9747581 0.9681489 0.1520337 0.9858352 0.8338015


##### Bcell
Bcell = cell.type.specific(FlowSorted.Blood.450k.compTable, 'Bcell', 0.6)
FlowSorted.Blood.450k.compTable['cg14482811',]
dim(Bcell)
# cg04153882, WIPI2, chr7:5229800-5236200 (B-cell marker)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg04153882 2373.206 5.209077e-38 0.9499032 0.9563241 0.9048152 0.10240016 0.9184050 0.9507915 0.07339106 0.9711910
# cg13738327, LRP5, chr11:68369800-68375400 (B-cell marker)
#             Fstat      p.value      CD8T      CD4T        NK      Bcell      Mono      Gran        low      high
# cg13738327 1214.981 1.151174e-33 0.9090035 0.9645673 0.8823188 0.09531174 0.97041521 0.9752966 0.07847517 0.9878111
# cg11699517, BAHCC1, chr17:81451000-81457800
# cg11699517 1008.147 1.860576e-32 0.9800161 0.9830573 0.9680402 0.1903934 0.954495 0.9427486 0.1313168 0.9878571 0.8565404
# cg14482811, LCN8, chr9:136757200-136765000
# cg14482811 422.1361 7.628536e-27 0.911173 0.9509423 0.7069898 0.1052258 0.9173946 0.929602 0.07817443 0.964382 0.8862075



# samtools view -b output_38.bam "chr7:5229800-5236200" > WIPI2_WB.bam
# samtools index WIPI2_WB.bam
# samtools view -b output_38.bam "chr7:4712000-4718200" > FOXK1_WB.bam
# samtools index FOXK1_WB.bam
# samtools view -b output_38.bam "chr11:68369800-68375400" > LRP5_WB.bam
# samtools index LRP5_WB.bam
# samtools view -b output_38.bam "chr3:196218200-196225400" > SLC51A_WB.bam
# samtools index SLC51A_WB.bam
# samtools view -b output_38.bam "chr22:21935600-21954592" > PPM1F_WB.bam
# samtools index PPM1F_WB.bam
# samtools view -b output_38.bam "chr14:99185332-99246000" > BCL11B_WB.bam
# samtools index BCL11B_WB.bam
# samtools view -b output_38.bam "chr12:10300600-10314800" > KLRD1_WB.bam
# samtools index KLRD1_WB.bam


# 10. mQTL



## DATABASE - confirmed ##
# cg05872129, cg24399712, cg11247378  chr22 39784982      SYNGR1  NOT_VERIFIED
# cg08704250, cg08109568, cg14829155  chr15 31115871      HERC2P10  NOT_VERIFIED
# cg21942184, cg17397159, cg10864200, cg00614959  chr4  720870     PCGF3  NOT_VERIFIED
# cg03295274, cg08779649, cg03651054  chr13 50194554      ARL11 NOT_VERIFIED
# cg07185983, cg06864789, cg25399239, cg18136963  chr6  139012860      FLJ46906 NOT_VERIFIED
# cg09408571, cg06223162, cg14855972, cg10298815 chr1  101003688    GPR88 NOT_VERIFIED
# cg16653991, cg08899523  chr10 123244591      FGFR2  NOT_VERIFIED
# cg03958363, cg18705808, cg04633141  chr7  1588319      TMEM184A NOT_VERIFIED
# cg27003165, cg01062020  chr7  156181990      SH2D1B NOT_VERIFIED
# cg08033130, cg19145607  chr3  45983597    FYCO1   NOT_VERIFIED
# cg08382737, cg10667338  chr19 49617042      LIN7B NOT_VERIFIED
# cg06791446, cg25052156, cg22633036, cg11430259, cg02210151, cg17681491, cg18566515  chr10 123356041      FGFR2  NOT_VERIFIED
# cg23558002, cg20101529, cg02645340, cg07911080, cg09781827, cg13156931, cg13847589, cg25591794, cg15332217, cg21950155  chr13 28555387      URAD  NOT_VERIFIED
# cg11440486, cg00901687, cg22571038, cg03611598, cg03168497  chr17 48585270      MYCBPAP NOT_VERIFIED
# cg22480875, cg07679219  chr12 77417738      E2F7  NOT_VERIFIED
# cg00071950, cg11266682  chr4  10021025      SLC2A9  NOT_VERIFIED
# cg08779649, cg03651054  chr13 50194643      ARL11 NOT_VERIFIED
# cg23649088, cg17644776  chr2  200775458      C2orf69  NOT_VERIFIED
# cg12454169, cg15652532, cg17749961  chr2  30669597      LCLAT1  NOT_VERIFIED
# cg09381666  chr16 88757816      SNAI3-AS1_RNF166  NOT_VERIFIED
# cg21138405, cg00255919  chr5  131827918      IRF1 NOT_VERIFIED
# cg03475776, cg25521439, cg27074582  chr2  240114406      HDAC4  NOT_VERIFIED
# cg16145915, cg05779406, cg18765753, cg12950816  chr7  1198977      ZFAND2A  NOT_VERIFIED
# cg07157030, cg18771300  chr14 63671737      RHOJ  NOT_VERIFIED
# cg20923605, cg19891452 chr13 112978683      LINC01044 NOT_VERIFIED
# cg11345693, cg14413466  chr17 79170920      CEP131  NOT_VERIFIED

# samtools view -b output_38.bam "chr22:39388416-39389373" > SYNGR1_WB.bam
# samtools index SYNGR1_WB.bam

# samtools view -b output_38.bam "chr13:27978527-27982836" > URAD_WB.bam
# samtools index URAD_WB.bam

# samtools view -b output_38.bam "chr10:121595171-121597691" > FGFR2_WB.bam
# samtools index FGFR2_WB.bam
# 


# FGFR2_WB - good_cov, good_rhocpg, good_IM - celltype_specific or environmental
# URAD_WB.bam - good_cov, good_rhocpg, good_IM - celltype_specific or environmental
# SYNGR1_WB.bam - low_cov




# chr20:62950800-62960600 cg12099423, cg17221813, cg19142181, cg00688402  chr20 61590751     SLC17A9  VERIFIED
# chr1:205849000-205851600  cg17178900, cg26354017, cg14159672, cg14893161, cg11965913, cg07167872, cg24503407, cg07157834  chr1  205819251      PM20D1  VERIFIED
# chr16:87645200-87655200 cg03834411, cg26764761  chr16 87682142      JPH3  VERIFIED
# chr2:208110200-208115600  cg24664470, cg10392614, cg20408181, cg04883911  chr2  208977395      NR_038437  VERIFIED
# chr19:13763400-13766400 cg16474696, cg25755428, cg14272688  chr19 13875014     CCDC130_MRI1 VERIFIED
# chr17:6654200-6655993 cg13207180, cg24686902, cg08103988, cg21358336, cg08750459  chr17 6558440      C17orf100_MIR4520-1  VERIFIED
# chr6:2950000-2953600  cg00101728, cg05082466  chr6  2953027      SERPINB6 VERIFIED
# chr10:132228600-132233800 cg14418756, cg27173819, cg08485649, cg01458105, cg19754622, cg16106427, cg25969878      STK32C  VERIFIED
# chr9:122225800-122228600  cg04282082, cg21213617, cg21469772, cg13571460, cg13862711, cg03363289, cg00142257  chr9  124989915      LHX6 VERIFIED
# chr12:629400-633600 cg14911689, ... chr12 630815  NINJ2-AS1 VERIFIED

# samtools view -b output_38.bam "chr20:62950800-62960600" > SLC17A9_WB.bam
# samtools index SLC17A9_WB.bam
# samtools view -b output_38.bam "chr1:205849000-205851600" > PM20D1_WB.bam
# samtools index PM20D1_WB.bam
# samtools view -b output_38.bam "chr16:87645200-87655200" > JPH3_WB.bam
# samtools index JPH3_WB.bam
# samtools view -b output_38.bam "chr2:208110200-208115600" > NR_038437_WB.bam
# samtools index NR_038437_WB.bam
# samtools view -b output_38.bam "chr19:13763400-13766400" > CCDC130_MRI1_WB.bam
# samtools index CCDC130_MRI1_WB.bam
# samtools view -b output_38.bam "chr17:6654200-6655993" > C17orf100_MIR4520-1_WB.bam
# samtools index C17orf100_MIR4520-1_WB.bam
# samtools view -b output_38.bam "chr6:2950000-2953600" > SERPINB6_WB.bam
# samtools index SERPINB6_WB.bam
# samtools view -b output_38.bam "chr10:132228600-132233800" > STK32C_WB.bam
# samtools index STK32C_WB.bam
# samtools view -b output_38.bam "chr9:122225800-122228600" > LHX6_WB.bam
# samtools index LHX6_WB.bam
# samtools view -b output_38.bam "chr12:629400-633600" > NINJ2-AS1_WB.bam
# samtools index NINJ2-AS1_WB.bam

# wc -l not_verified 
# 25 not_verified
# wc -l verified 
# 9 verified
9/(25+9) # 0.2647059


# # Trans-mQTL? NOPE
# # cg13592947  chr5  1111049 chr5:1106000-1111600  SLC12A7
# # cg10959984(27K?);cg21470387(450K) ALDH2     chr12:111763800-111768800 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4127343/table/T1/)
# wbresgos[wbresgos$reg == 'chr5:1106000-1111600',] # 1.767482e-57; SLC12A7
# wbresgos[wbresgos$reg == 'chr12:111763800-111768800',] # 2.222718e-08; ALDH2
# 
# # samtools view -b output_38.bam "chr5:1106000-1111600" > SLC12A7_WB.bam
# # samtools index SLC12A7_WB.bam
# # samtools view -b output_38.bam "chr12:111763800-111768800" > ALDH2_WB.bam
# # samtools index ALDH2_WB.bam



# 11. TFBS enrichment analysis

motif_enrichment = function(df, sample) 
{
  p_val_cutoff = 0.05/nrow(df)
  sig_regions = subset(df, P < p_val_cutoff)
  bg_regions = subset(df, P >= p_val_cutoff)
  sig_regions = sig_regions$reg
  bg_regions = bg_regions$reg
  fwrite(x = as.data.frame(sig_regions), file = paste("/home/bronte/opt/binokulars_analysis/2_annotation/sig_", sample, sep = ''), quote = F, row.names = F, col.names = F)
  fwrite(x = as.data.frame(bg_regions), file = paste("/home/bronte/opt/binokulars_analysis/2_annotation/bg_", sample, sep = ''), quote = F, row.names = F, col.names = F)
}

## To extract sequences for target/Bg we employed samtools:
# samtools faidx hg19.fasta.gz -r evCpGs500 > evCpGs500.fa
# samtools faidx hg19.fasta.gz -r bg500 > bg500.fa
#~/opt/samtools-1.14/samtools faidx ../../../disk3_data/ena_pooled_data/reference_genome/hg19.fa -r ~/opt/binokulars_analysis/2_annotation/sig_Pooled_Blood_hg19 > sig_pooled_blood_hg19.fa
## We then performed motif enrichment analysis with Homer by calling:
# findMotifs.pl evCpGs500.fa fasta . -fasta bg500.fa -p 4 -humanGO > log.txt

output_location = "/home/bronte/opt/binokulars_analysis/images/2_annotation"
setwd(output_location)

files = c("/media/nw_disk2/bronte/pooled_blood_hg19_analysis/binokulars_output/test_results/p_values.txt",
          "/media/nw_disk2/bronte/pooled_sperm_hg19_analysis/binokulars_output/test_results/p_values.txt",
          "/media/nw_disk2/bronte/pooled_blood_hg38_analysis/binokulars_output/test_results/p_values.txt",
          "/media/nw_disk2/bronte/pooled_sperm_hg38_analysis/binokulars_output/test_results/p_values.txt")

sample_names = c("Pooled_Blood_hg19",
                 "Pooled_Sperm_hg19",
                 "Pooled_Blood_hg38",
                 "Pooled_Sperm_hg38")

for(i in 1:length(files)) {
  motif_enrichment(files[i], sample_names[i])
}




get_version_info(organism = "hsapiens")
# $biomart
# [1] "Ensembl"
# $biomart_version
# [1] "106"
# $display_name
# [1] "Human"
# $genebuild
# [1] "GRCh38.p13"
# $gprofiler_version
# [1] "e106_eg53_p16_65fcd97"
# $organism
# [1] "hsapiens"
# $sources
# $sources$CORUM
# $sources$CORUM$name
# [1] "CORUM protein complexes"
# $sources$CORUM$version
# [1] "03.09.2018 Corum 3.0"
# $sources$`GO:BP`
# $sources$`GO:BP`$name
# [1] "biological process"
# $sources$`GO:BP`$version
# [1] "annotations: BioMart\nclasses: releases/2022-03-22"
# $sources$`GO:CC`
# $sources$`GO:CC`$name
# [1] "cellular component"
# $sources$`GO:CC`$version
# [1] "annotations: BioMart\nclasses: releases/2022-03-22"
# $sources$`GO:MF`
# $sources$`GO:MF`$name
# [1] "molecular function"
# $sources$`GO:MF`$version
# [1] "annotations: BioMart\nclasses: releases/2022-03-22"
# $sources$HP
# $sources$HP$name
# [1] "Human Phenotype Ontology"
# $sources$HP$version
# [1] "annotations: hpo.annotations #12\nclasses: None"
# $sources$HPA
# $sources$HPA$name
# [1] "Human Protein Atlas"
# $sources$HPA$version
# [1] "annotations: HPA website: 21-12-06\nclasses: script: 21-12-17"
# $sources$KEGG
# $sources$KEGG$name
# [1] "Kyoto Encyclopedia of Genes and Genomes"
# $sources$KEGG$version
# [1] "KEGG FTP Release 2022-05-16"
# $sources$MIRNA
# $sources$MIRNA$name
# [1] "miRTarBase"
# $sources$MIRNA$version
# [1] "Release 7.0"
# $sources$REAC
# $sources$REAC$name
# [1] "Reactome"
# $sources$REAC$version
# [1] "annotations: BioMart\nclasses: 2022-5-18"
# $sources$TF
# $sources$TF$name
# [1] "Transfac"
# $sources$TF$version
# [1] "annotations: TRANSFAC Release 2021.3\nclasses: v2"
# $sources$WP
# $sources$WP$name
# [1] "WikiPathways"
# $sources$WP$version
# [1] "20220510"
# $taxonomy_id
# [1] "9606"

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] seqinr_4.2-16         org.Hs.eg.db_3.15.0   AnnotationDbi_1.58.0  IRanges_2.30.0        S4Vectors_0.34.0     
# [6] Biobase_2.56.0        BiocGenerics_0.42.0   DOSE_3.22.1           clusterProfiler_4.4.4 gprofiler2_0.2.1     
# [11] dplyr_1.0.9           gplots_3.1.3          gridExtra_2.3         ggplot2_3.3.6         qqman_0.1.8          
# [16] gridGraphics_0.5-1    tidyr_1.2.0           data.table_1.14.2    
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.22.0           colorspace_2.0-3       ggtree_3.4.2           ellipsis_0.3.2         qvalue_2.28.0         
# [6] XVector_0.36.0         aplot_0.1.7            rstudioapi_0.13        farver_2.1.1           graphlayouts_0.8.1    
# [11] ggrepel_0.9.1          bit64_4.0.5            fansi_1.0.3            scatterpie_0.1.8       codetools_0.2-18      
# [16] splines_4.2.1          cachem_1.0.6           GOSemSim_2.22.0        polyclip_1.10-0        ade4_1.7-19           
# [21] jsonlite_1.8.0         GO.db_3.15.0           png_0.1-7              shiny_1.7.2            ggforce_0.3.4         
# [26] compiler_4.2.1         httr_1.4.3             assertthat_0.2.1       Matrix_1.4-1           fastmap_1.1.0         
# [31] lazyeval_0.2.2         cli_3.3.0              later_1.3.0            tweenr_2.0.2           htmltools_0.5.3       
# [36] tools_4.2.1            igraph_1.3.4           gtable_0.3.0           glue_1.6.2             GenomeInfoDbData_1.2.8
# [41] reshape2_1.4.4         DO.db_2.9              fastmatch_1.1-3        Rcpp_1.0.9             enrichplot_1.16.2     
# [46] vctrs_0.4.1            Biostrings_2.64.0      ape_5.6-2              nlme_3.1-159           crosstalk_1.2.0       
# [51] ggraph_2.0.6           stringr_1.4.0          mime_0.12              lifecycle_1.0.1        gtools_3.9.3          
# [56] zlibbioc_1.42.0        MASS_7.3-58.1          scales_1.2.0           tidygraph_1.2.2        promises_1.2.0.1      
# [61] parallel_4.2.1         RColorBrewer_1.1-3     yaml_2.3.5             memoise_2.0.1          downloader_0.4        
# [66] ggfun_0.0.7            yulab.utils_0.0.5      calibrate_1.7.7        stringi_1.7.8          RSQLite_2.2.15        
# [71] tidytree_0.4.0         caTools_1.18.2         BiocParallel_1.30.3    GenomeInfoDb_1.32.3    rlang_1.0.4           
# [76] pkgconfig_2.0.3        bitops_1.0-7           lattice_0.20-45        purrr_0.3.4            labeling_0.4.2        
# [81] treeio_1.20.2          patchwork_1.1.2        htmlwidgets_1.5.4      shadowtext_0.1.2       bit_4.0.4             
# [86] tidyselect_1.1.2       plyr_1.8.7             magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
# [91] DBI_1.1.3              pillar_1.8.0           withr_2.5.0            KEGGREST_1.36.3        RCurl_1.98-1.8        
# [96] tibble_3.1.8           crayon_1.5.1           KernSmooth_2.23-20     utf8_1.2.2             plotly_4.10.0         
# [101] viridis_0.6.2          blob_1.2.3             digest_0.6.29          xtable_1.8-4           httpuv_1.6.5          
# [106] munsell_0.5.0          viridisLite_0.4.0      ggplotify_0.1.0   
