############################################################################
############################################################################
###########                                                      ###########
###########                    Simulation study                  ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########           b.planterosejimenez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

# Load libraries
library(parallel)
library(ggplot2)

# Load functions
place_CpGs <- function(L, rho_cpg)
{
  n_cpgs = floor(L*rho_cpg)
  set = 1:(L-1)
  positions = integer(length = n_cpgs)
  for(i in 1:n_cpgs)
  {
    pos_i = sample(set, 1)
    positions[i] = pos_i
    set = set[!(set %in% c(pos_i-1, pos_i, pos_i+1))]
    if(length(set) %in% c(0,1))
    {
      stop("CpG density is too high")
    }
  }
  positions = sort(positions)
  return(positions)
}

process_read <- function(read_pos_i, JRC_indicator_i, p, cg_pos)
{
  n = sum(cg_pos > read_pos_i & cg_pos < read_pos_i+100)
  if(JRC_indicator_i)
  {
    p_hat = sample(c(0,1), 1, prob = c(1-p, p))
    M = rbinom(n = 1, size = n, prob = p_hat)
    U = n - M
  }
  else
  {
    M = rbinom(n = 1, size = n, prob = p)
    U = n - M
  }
  return(c(read_pos_i, U, M))
}

generate_UM_counts <- function(L, l, rho_cpg, phi, n_reads, p, cg_pos = NULL)
{
  if(is.null(cg_pos))
  {
    cg_pos = place_CpGs(L, rho_cpg)
  }
  read_pos = sample(1:(L-l), n_reads, replace = T)
  JRC_indicator = sample(c(0,1), n_reads, prob = c(1-phi, phi), replace = T)
  UM_counts = lapply(1:n_reads, function(i) process_read(read_pos[i], JRC_indicator[i], p, cg_pos))
  UM_counts = Reduce(f = rbind, UM_counts)
  colnames(UM_counts) = c("start", "U", "M")
  UM_counts = UM_counts[order(UM_counts[,1]),]
  return(list(cg_pos = cg_pos, UM = as.data.frame(UM_counts)))
}

plot_UM <- function(L, l, cg_pos, UM, cex.point = 1, ...)
{
  plot(0, 1, xlim = c(1,L), ylim = c(1, nrow(UM)), type = "n", xlab = 'Genomic position', ylab = 'Read',...)
  abline(v = cg_pos, col = "yellow2")
  rect(xleft = 0, xright = L, ybottom = 0-0.8, ytop = 0+0.05, border = NA, col = "cornflowerblue")
  abline(h = 1:nrow(UM), lty = 2, col = "gray")
  for(i in 1:nrow(UM))
  {
    arrows(x0 = UM$start[i], x1 = UM$start[i]+100, y0 = i, y1 = i, `length` = 0)
    if(UM[i,]$U + UM[i,]$M > 0)
    {
      where = cg_pos[cg_pos > UM$start[i] & cg_pos < UM$start[i] + 100]
      col = sample(c(rep("black", UM[i,]$M), rep("white", UM[i,]$U)), UM[i,]$M+UM[i,]$U, replace = F)
      points(x = where, y = rep(i, length(where)), col = col, pch = 19, cex = cex.point)
      points(x = where, y = rep(i, length(where)), col = "black", pch = 1, cex = cex.point)
    }
  }
}

permute <- function(beta_vec, n_vec)
{
  M = vapply(1:length(beta_vec), function(i) rbinom(n = 1, size = n_vec[i], prob = beta_vec[i]), as.integer(1L))
  -sum(stats::dbinom(x = M, size = n_vec, prob = beta_vec, log = T))
}

simulation <- function(L, l, b, rho_cpg, n_reads, p, phi, N_boot, N_iter, N_rep)
{
  pval_list = as.list(rep(c(NA),N_rep))
  for(i in 1:N_rep)
  {
    pval_vec = lapply(1:N_iter, function(i) tryCatch(iteration(L, l, b, rho_cpg, n_reads, p, phi, N_boot), error = function(x) NA))
    pval_list[[i]] = unlist(pval_vec)
  }
  return(pval_list)
}

summarize_simulation <- function(pval_list)
{
  sapply(pval_list, function(x) mean(x < 0.05, na.rm = T))
}

iteration <- function(L, l, b, rho_cpg, n_reads, p, phi, N_boot)
{
  # cg_pos is not specified so a different CpG distribution is generated in each iteration
  UM = generate_UM_counts(L = L, l = l, rho_cpg = rho_cpg, phi = phi, n_reads = n_reads, p = p)$UM
  
  # Aggregate counts to estimate methylation levels per bin
  breaks = seq(0, L+1, b)
  UM$bins = cut(UM$start, breaks = breaks, labels = F)
  M = aggregate(M ~ bins, data = UM, FUN = sum)$M
  n = aggregate(M+U ~ bins, data = UM, FUN = sum)$`M + U`
  beta = M/n
  l_bins = table(UM$bins)
  beta_vec = unlist(sapply(1:length(l_bins), function(i) rep(beta[i], l_bins[i])))
  n_vec = UM$M + UM$U
  M_vec = UM$M
  
  # Filter potential missing values in beta created by [0/0] and n_vec values of 0
  keep = (n_vec > 0) & !is.na(beta_vec)
  beta_vec = beta_vec[keep]
  n_vec = n_vec[keep]
  M_vec = M_vec[keep]
  
  # Compute logL under H0
  logL_H0 = vapply(1:N_boot, function(i) permute(beta_vec, n_vec), numeric(1))
  
  # Compute test statistic
  logL = -sum(stats::dbinom(x = M_vec, size = n_vec, prob = beta_vec, log = T))
  
  # Compute permutation pvalue
  pval = 1 - mean(logL > logL_H0)
  return(pval)
}


########################## Define Parameters ##########################

# N_boot: Number of bootstraps in Binokulars
# N_iter: Number of times each simulation conditions is tested
# L: IM region length
# l: read length
# b: bin size
# rho_cpg # CpG density
# n_reads # coverage
# p # Mhaplo_frequency
# phi # effect size (proportion of reads that display JRC behaviour)

########################## Visualize simulated data ##########################

# EWAS association - Non-JRC Vs JRC
par(mfrow = c(2,2))
set.seed(4)
cg_pos = place_CpGs(L = 600, rho_cpg = 0.02)
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.02, phi = 0, n_reads = 50, p = 0.4, cg_pos = cg_pos)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'non-JRC; Group 1; p = 0.4')
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.02, phi = 0, n_reads = 50, p = 0.8, cg_pos = cg_pos)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'non-JRC; Group 2; p = 0.8')
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.02, phi = 1, n_reads = 50, p = 0.4, cg_pos = cg_pos)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'JRC; Group 1; p = 0.4')
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.02, phi = 1, n_reads = 50, p = 0.8, cg_pos = cg_pos)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'JRC; Group 2; p = 0.8')


########################## Preliminary testing ##########################

# TN
set.seed(1)
iteration(L = 2000, l = 100, b = 200, rho_cpg = 0.035, n_reads = 100, 
          p = 0.5, phi = 0, N_boot = 500)

# TP
set.seed(1)
iteration(L = 2000, l = 100, b = 200, rho_cpg = 0.035, n_reads = 100, 
          p = 0.5, phi = 1, N_boot = 500)

# Reduced parameter space
L = 2000
l = 100
b = 200
rho_cpg = c(0.035, 0.055)
n_reads = c(100, 200)
p = c(0.4, 0.5)
phi = c(0, 1)
conditions = expand.grid(L, l, b, rho_cpg, n_reads, p, phi)
colnames(conditions) = c("L", "l", "b", "rho_cpg", "n_reads", "p", "phi")
N_boot = 1000
N_iter = 10
N_rep = 2

# Runtime expectations
start_time <- Sys.time()
jkl = iteration(L, l, b, rho_cpg[1], n_reads[1], p[1], phi[1], N_boot)
end_time <- Sys.time()
time_per_iter = difftime(end_time, start_time, units = "hours")
nCores = 4; nrow(conditions)*N_iter*N_rep*time_per_iter/nCores # 0.01740361 hours

# Observed runtime
start_time <- Sys.time()
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
pval_list = mclapply(1:nrow(conditions), function(i) simulation(L = conditions$L[i], l = conditions$l[i], b = conditions$b[i], 
                                                                rho_cpg = conditions$rho_cpg[i], phi = conditions$phi[i], 
                                                                n_reads = conditions$n_reads[i], p = conditions$p[i], N_boot = N_boot,
                                                                N_iter = N_iter, N_rep = N_rep),  mc.cores = 4)
end_time <- Sys.time()
difftime(end_time, start_time, units = "hours") # 0.02226207 hours

# Visualize results
results = conditions[rep(seq_len(nrow(conditions)), each = N_rep), ]
results$p_sig = unlist(lapply(1:length(pval_list), function(x) summarize_simulation(pval_list[[x]])))
results$rho_cpg = factor(results$rho_cpg)
results$p = factor(results$p)
results$phi = factor(results$phi)
ggplot(data = results, mapping = aes(x = n_reads, y = p_sig, col = phi)) +
  facet_wrap(p ~ rho_cpg, ncol = 2, labeller = labeller(p = label_both, rho_cpg = label_both)) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  geom_point(alpha = 0.8) + 
  geom_smooth(alpha = 0.5, method = "loess", level = 0.95) +
  theme(legend.text=element_text(size=8),
        strip.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste('L = ', L, '; l = ', l, '; b = ', b, '; N_iter = ', N_iter, '; N_boot = ', 
                N_boot, '; N_rep = ', N_rep, sep = '')) +
  ylab("1/N_iter * Sum_i I(pval_i < 0.05)") + ylim(c(0,1))




########################## Perform simulations ##########################

# Parameter space
L = 2000
l = 100
b = 200
#rho_cpg = seq(0.025, 0.065, 0.01)
rho_cpg = seq(0.002, 0.202, 0.05)
#n_reads = c(10, 25, 50, 100, 250, 500, 750, 1000)
n_reads = c(10, 25, 50, 100, 250, 500, 750, 1000, 1500)
p = c(0.1, 0.25, 0.5)
#phi = c(0, 0.25, 0.5, 0.75, 1)
phi = c(0, 0.1, 0.25, 0.5, 0.75, 1)
conditions = expand.grid(L, l, b, rho_cpg, n_reads, p, phi)
colnames(conditions) = c("L", "l", "b", "rho_cpg", "n_reads", "p", "phi")
head(conditions)
dim(conditions) # 600   7

# Simulation depth
N_boot = 1000
N_iter = 100
N_rep = 5
nCores = 4

# Runtime expectations (overestimate; low n_reads take very little time, expect half roughly)
start_time <- Sys.time()
jkl = iteration(L, l, b, rho_cpg[2], n_reads[6], p[3], phi[3], N_boot)
end_time <- Sys.time()
time_per_iter = difftime(end_time, start_time, units = "hours")
nCores = 4; nrow(conditions)*N_iter*N_rep*time_per_iter/nCores # 106.5571 hours

# Observed runtime
start_time <- Sys.time()
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
pval_list = mclapply(1:nrow(conditions), function(i) simulation(L = conditions$L[i], l = conditions$l[i], b = conditions$b[i], 
                                                                rho_cpg = conditions$rho_cpg[i], phi = conditions$phi[i], 
                                                                n_reads = conditions$n_reads[i], p = conditions$p[i], N_boot = N_boot,
                                                                N_iter = N_iter, N_rep = N_rep),  mc.cores = nCores)
end_time <- Sys.time()
difftime(end_time, start_time, units = "hours") # 51.17565 hours

# Save results
setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/simulations_binokulars/')
saveRDS(pval_list, paste(format(Sys.time(), '%Y%m%d'), 'simulation_results.Rds', sep = '_'))

# Check for missing values
sum(is.na(unlist(pval_list))) # 0

# Prep data.frame for visualization
results = conditions[rep(seq_len(nrow(conditions)), each = N_rep), ]
results$p_sig = unlist(lapply(1:length(pval_list), function(x) summarize_simulation(pval_list[[x]])))
results$rho_cpg = factor(results$rho_cpg)
results$p = factor(results$p)
results$phi = factor(results$phi)
results$coverage = results$n_reads*results$l/results$L # coverage = (read count * read length ) / total genome size.

# Visualize results
p1 = ggplot(data = results, mapping = aes(x = coverage, y = p_sig, col = phi)) +
  facet_wrap(p ~ rho_cpg, ncol = length(levels(results$rho_cpg)), labeller = labeller(p = label_both, rho_cpg = label_both)) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  geom_point(alpha = 0.6) + 
  geom_smooth(alpha = 0, method = "glm", method.args = list(family = 'quasibinomial', control = list(maxit = 500))) + # span = 5
  #geom_smooth(alpha = 0.5, method = "loess", level = 0.95, method.args = list(degree = 2, family = 'gaussian')) + # span = 5
  theme(legend.text=element_text(size=8),
        strip.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste('L = ', L, '; l = ', l, '; b = ', b, '; N_iter = ', N_iter, '; N_boot = ', 
                N_boot, '; N_rep = ', N_rep, sep = '')) +
  ylab("1/N_iter * Sum_i I(pval_i < 0.05)") + ylim(c(-0.05,1.05)) + xlim(c(0, max(results$coverage)))


########################## Improvements ##########################

### RUN at a different region size? why 2 kb?


########################## Limitations ##########################

# Limitations
# 1) Assumes equal size single-end reads
# 2) Assumes homogeneous coverage
# 3) Assumes single-sized region
# 4) Assumes infinite size of pool of equi-contributors - THE MOST LIMITING ASSUMPTION
# 5) Needs to be tested - therefore, identified by JRC_seeker

# Take 2n alleles in the pool from n individuals
# Probability p for a given allele. Probability of observing it is:
# P(X > 0 | p) = 1-P(X = 0 | p) = 1 - (1-p)^(2*n)
# 0.95 = 1 - (1-p)^(2*n)
n_individuals = 1:20
plot(n_individuals, 1-(0.05)^(1/(2*n_individuals)), ylim = c(0,1), type = 'o',
     ylab = 'Minimum p observed at 95% confidence', xlab = 'Number of individuals in the pool', cex = 0.5, cex.lab = 0.8)
abline(v = 12, lty = 2)
abline(h = 1-(0.05)^(1/(2*12)), lty = 2)
text(-0.4, 1-(0.05)^(1/(2*12)), round(1-(0.05)^(1/(2*12)), 3), xpd = NA, cex = 0.5)

# What's the coverage in the pool of blood?
# coverage = (read count * read length ) / total genome size.

(n_reads*100)/(3.3*10^9)


########################## Rho CpG density at hg38 ##########################

library(reshape2)
library(ggplot2)
binning_count <- function(cpgr, bin.length = 1000L)
{
  LEV = as.character(unique(seqnames(cpgr)))
  RES_list = rep(list(NA), length(LEV))
  for(i in 1:length(LEV))
  {
    print(LEV[i])
    pos = start(cpgr[seqnames(cpgr) == LEV[i]])
    m = as.integer(min(pos)-1L); M = as.integer(max(pos))
    n = as.integer(floor((M - m)/bin.length))
    RES_list[[i]] = table(cut(pos, breaks = seq(m, m + (n+1)*bin.length, bin.length)))
  }
  names(RES_list) = LEV
  return(RES_list)
}

library(BSgenome.Hsapiens.UCSC.hg38)  
chrs = names(Hsapiens)[1:24]
cgs = lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr = do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1))))
length(cpgr) # 29,401,360

bin.length = 1000L
count_vec = binning_count(cpgr = cpgr, bin.length = bin.length)
df = melt(count_vec)
colnames(df) = c('bin', 'count', 'chr')
df$rho = df$count/bin.length
  
p2 = ggplot(data = df, mapping = aes(x = rho, fill = chr)) + facet_wrap( ~ chr) + geom_density() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_vline(xintercept = c(0.002, 0.202), linetype = 'dashed')

# Min rho: to have a neighbour at 1 kb, density must be:
2/1001 # ~ 0.002

2/100 # read length limit

# Max rho: to make sure the CpG placing algorithm does not fail in a deterministic way:
floor(2000/3)/2000 ~ 0.333


### Export figures
setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/') # 12
set.seed(15)
cg_pos = place_CpGs(L = 600, rho_cpg = 0.02)
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.025, phi = 0, n_reads = 20, p = 0.5, cg_pos = cg_pos)
tiff(filename = 'nonJRC.tiff', width = 3.3, height = 2.4, units = 'in', res = 600)
par(mar=c(0.5, 0.5, 0.9, 0.5), mgp=c(0.5, 0.2, 0), las=0)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'Non-JRC; p = 0.5', cex.main = 0.4, 
        cex.axis = 0.1, cex.point = 0.8, cex.lab = 0.4, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = seq(0, 600, 100), tck=-0.01, cex.axis = 0.3, mgp = c(-0.7, -0.4, 0))
axis(side = 2, at = seq(0, 50, 10), tck=-0.01, cex.axis = 0.3, las = 2, mgp = c(0.5, 0.2, 0))
dev.off()
UM_list = generate_UM_counts(L = 600, l = 100, rho_cpg = 0.025, phi = 1, n_reads = 20, p = 0.5, cg_pos = cg_pos)
tiff(filename = 'JRC.tiff', width = 3.3, height = 2.4, units = 'in', res = 600)
par(mar=c(0.5, 0.5, 0.9, 0.5), mgp=c(0.5, 0.2, 0), las=0)
plot_UM(L = 600, l = 100, UM_list$cg_pos, UM_list$UM, main = 'JRC; p = 0.5', cex.main = 0.4, 
        cex.axis = 0.1, cex.point = 0.8, cex.lab = 0.4, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = seq(0, 600, 100), tck=-0.01, cex.axis = 0.3, mgp = c(-0.7, -0.4, 0))
axis(side = 2, at = seq(0, 50, 10), tck=-0.01, cex.axis = 0.3, las = 2, mgp = c(0.5, 0.2, 0))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
ggsave(filename = 'simulation_results.tiff', plot = p1, device = 'tiff', width = 6.7, height = 6.4, units = 'in', dpi = 600)
ggsave(filename = 'rho_cg_bin1kb.tiff', plot = p2, device = 'tiff', width = 6.7, height = 6.4, units = 'in', dpi = 300)
###

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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6
# 
# loaded via a namespace (and not attached):
# [1] fansi_1.0.3      withr_2.5.0      assertthat_0.2.1 dplyr_1.0.9      utf8_1.2.2       grid_4.2.1       R6_2.5.1        
# [8] DBI_1.1.3        lifecycle_1.0.1  gtable_0.3.0     magrittr_2.0.3   scales_1.2.0     pillar_1.8.0     rlang_1.0.4     
# [15] cli_3.3.0        rstudioapi_0.13  generics_0.1.3   vctrs_0.4.1      tools_4.2.1      glue_1.6.2       purrr_0.3.4     
# [22] munsell_0.5.0    compiler_4.2.1   pkgconfig_2.0.3  colorspace_2.0-3 tidyselect_1.1.2 tibble_3.1.8  


