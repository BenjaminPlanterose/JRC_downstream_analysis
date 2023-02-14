############################################################################
############################################################################
###########                                                      ###########
###########                 JRC Annotation                       ###########
###########              Author: Bronte Kolar                    ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########            Rotterdam, The Netherlands                ###########
###########                                                      ###########
###########                                                      ###########
############################################################################
############################################################################

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)

intron_function = function(x)
{
  list(x[[1]], x[[2]], x[[3]], x[[4]], as.integer(x[[6]] + 1), as.integer(x[[9]] - 1), "intron")
}

promoter_func = function(x)
{
  if(x[[2]] == "+")
  {
    if(x[[3]] != 0) {
      p_start = max(0, x[[3]] - 1500)
      p_end = x[[3]] - 1
      list(x[[8]], x[[9]], x[[1]], x[[2]], p_start, p_end, "promoter")
    }
  } else if (x[[2]] == "-") {
    p_start = x[[4]] + 1
    p_end = x[[4]] + 1500
    list(x[[8]], x[[9]], x[[1]], x[[2]], p_start, p_end, "promoter")
  } else {
    print("ERROR - strand notation incorrect.")
  }
}

refseq_annotation = function(genome_version, filename, where)
{
  setwd(where)
  refseq_genes = fread(filename)
  colnames(refseq_genes) = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  compr_refseq = refseq_genes[, c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "name2")]
  
  ##############################
  # Combine genes / long RNA at the same location
  compr_refseq = compr_refseq %>%
    group_by(chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds) %>%
    summarise(name = paste(name, collapse = ", "), name2 = paste(name2, collapse = ", "))
  est_exons = as.integer(sum(compr_refseq$exonCount))
  est_introns = est_exons
  est_promoters = as.integer(nrow(compr_refseq))
  est_size = est_introns + est_exons + est_promoters + 1000
  
  ##############################
  # Add EXONS to annotation
  
  message('exons')
  #compr_refseq = head(compr_refseq,100)
  compr_refseq_orig = copy(compr_refseq)
  compr_refseq = separate_rows(compr_refseq, exonStarts, exonEnds, convert = TRUE)
  compr_refseq = na.omit(compr_refseq)
  annotation = setNames(data.frame(matrix(ncol = 9, nrow = est_size)), c("chrom", "strand", "txStart", "txEnd", "exonCount", "exonStarts", "exonEnds", "name", "name2"))
  annotation[1:nrow(compr_refseq),] = compr_refseq[1:nrow(compr_refseq),]
  annotation = annotation[,c(8,9,1,2,6,7)]
  colnames(annotation) = c("gene_name", "gene_name_2", "chrom", "strand", "start", "end")
  annotation[1:nrow(compr_refseq), "type"] = "exon"
  
  ##############################
  # Add INTRONS to annotation
  # "gene_name", "gene_name_2", "chrom", "strand", "start", "end", "type"
  message('introns')
  annotation_copy = copy(annotation)
  annotation_copy = transform(annotation_copy, next_gene_name = c(gene_name[-1], gene_name[1]))
  annotation_copy = transform(annotation_copy, next_exon_start = c(start[-1], start[1]))
  
  # for hg38, regions were found where exon next exon start was 2 bp after previous exon end
  intron_index = which(annotation_copy$gene_name == annotation_copy$next_gene_name & annotation_copy$end + 2 < annotation_copy$next_exon_start)
  intron_df = lapply(intron_index, function(x) intron_function(annotation_copy[x,]))
  intron_df = as.data.frame(do.call(rbind, intron_df))
  colnames(intron_df) = c("gene_name", "gene_name_2", "chrom", "strand", "start", "end", "type")
  #intron_df = do.call(rbind, Map(data.frame, gene_name=intron_df$gene_name, gene_name_2=intron_df$gene_name_2, chrom=intron_df$chrom, strand=intron_df$strand, start=intron_df$start, end=intron_df$end, type=intron_df$type))
  for(i in 1:ncol(intron_df)) 
  {
    intron_df[[i]] = unlist(intron_df[[i]])
  }
  
  ##############################
  # Add PROMOTERS to annotation
  
  message('promoters')
  promoter_df = lapply(1:nrow(compr_refseq_orig), function(x) promoter_func(compr_refseq_orig[x,]))
  promoter_df = as.data.frame(do.call(rbind, promoter_df))
  colnames(promoter_df) = c("gene_name", "gene_name_2", "chrom", "strand", "start", "end", "type")
  for(i in 1:ncol(promoter_df)) 
  {
    promoter_df[[i]] = unlist(promoter_df[[i]])
  }
  annotation = na.omit(annotation)
  annotation_final = do.call("rbind", list(annotation, intron_df, promoter_df))
  write.table(annotation_final, file=paste(genome_version, '_refseq_annotation.tsv', sep=""), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)  
}

refseq_annotation('hg19', "hg19.refGene.txt", '/media/ultron/2tb_disk1/ben/JRC_project/scripts/5_functional_annotation/2_annotation/assets/')
refseq_annotation('hg38', "hg38.refGene.txt", '/media/ultron/2tb_disk1/ben/JRC_project/scripts/5_functional_annotation/2_annotation/assets/')

# hg19.refGene.txt and hg38.refGene.txt were downloaded from UCSC at:
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/
