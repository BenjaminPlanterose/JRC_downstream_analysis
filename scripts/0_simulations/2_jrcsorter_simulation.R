############################################################################
############################################################################
###########                                                      ###########
###########                  Simulation study 2                  ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########           b.planterosejimenez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

## Load libraries
library(gplots)
library(data.table)
library(parallel)
library(MASS)
library(MLmetrics)
library(ggplot2)

## Load functions
sample_mat = function(N, p)
{
  if(length(N) == 0)
  {
    return(NULL)
  }
  if(is.null(ncol(N)))
  {
    dim(N) = c(length(N), 1)
  }
  matrix(sapply(1:prod(dim(N)), function(x) rbinom(n = 1, size = c(N)[x], prob = p)), nrow = nrow(N))
}

simulate_UM <- function(bin_size, f, L, n, lambda, p, epsilon, model)
{
  n_bin = L/bin_size
  n_bin_flank = f/bin_size
  
  # Fill in the flanking bins
  N_mat = matrix(rpois(n = n_bin*n, lambda), nrow = n_bin)
  M_mat = matrix(NA, nrow = n_bin, ncol = n)
  M_mat[1:n_bin_flank, ] = sample_mat(N_mat[1:n_bin_flank, ], (1-epsilon))
  M_mat[(n_bin-n_bin_flank+1):n_bin, ] = sample_mat(N_mat[(n_bin-n_bin_flank+1):n_bin, ], (1-epsilon))

  if(model == 'M1')
  {
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), ] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), ], p)
  }
  if(model == 'M2')
  {
    classes = sample(c(0,1), n, prob = c(1-p, p), replace= T)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0], epsilon)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1], (1-epsilon))
  }
  if(model == 'M3')
  {
    p = p/2
    classes = sample(c(0,1), n, prob = c(1-2*p, 2*p), replace= T)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0], epsilon)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1], 0.5)
  }
  if(model == 'M4')
  {
    p = 1-p/2
    classes = sample(c(0,1), n, prob = c(2*(1-p), 1-2*(1-p)), replace= T)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0], 0.5)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1], (1-epsilon))
  }
  if(model == 'M5')
  {
    classes = sample(c(0,1,2), n, prob = c((1-p)^2, 2*p*(1-p), p^2), replace= T)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 0], epsilon)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 1], 0.5)
    M_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 2] = sample_mat(N_mat[(n_bin_flank+1):(n_bin-n_bin_flank), classes == 2], (1-epsilon))
  }
  U_mat = N_mat - M_mat
  return(list(M = M_mat, U = U_mat))
}

Visualize_JRC <- function(M, U)
{
  breaks=seq(-0.05, 1.05, 0.05)
  my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)
  BETA = M/(M+U)
  heatmap.2(t(BETA), trace = 'none', Colv = 'none', breaks = breaks,
            na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
            offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
            main = '', cexCol = 0.8) 
}

simulation = function(bin_size, f, L, n, lambda, p, epsilon, N_rep, bin_tol, mod_tol)
{
  sim.data = lapply(c('M1', 'M2', 'M3', 'M4', 'M5'), function(i) lapply(1:N_rep, function(j) simulate_UM(bin_size, f, L, n, lambda, p, epsilon, i)))
  likelihood_list = lapply(1:5, function(i) sapply(1:N_rep, function(j) test_region(sim.data[[i]][[j]], epsilon = epsilon, tol = bin_tol)))
  likelihood_list = lapply(likelihood_list, t)
  likelihoods = as.data.frame(Reduce(rbind, likelihood_list))
  likelihoods$real = rep(1:5, each = N_rep)
  likelihoods$predicted = sapply(1:nrow(likelihoods), function(i) classify(as.numeric(likelihoods[i, 1:5]), tol = mod_tol))
  likelihoods$predicted[!(likelihoods$predicted %in% c('L1', 'L2_L5', 'L3_L5', 'L4_L5', 'L5'))] = 'L1'
  likelihoods$predicted2 = factor(likelihoods$predicted, levels =c('L1', 'L2_L5', 'L3_L5', 'L4_L5', 'L5'))
  levels(likelihoods$predicted2) = 1:5
  #likelihoods = likelihoods[!is.na(likelihoods$predicted),]
  return(F1_Score_macro(likelihoods$predicted2, likelihoods$real))
}

rowMin <- function(mat)
{
  sapply(1:nrow(mat), function(x) which.min(mat[x,]))
}

JRC_sorter <- function(U, M, epsilon)
{
  total = U+M
  
  # Model 1
  p1 = sum(M)/sum(total)
  L1 = sum(dbinom(x = M, size = total, prob = p1, log = T))
  
  # Model 2
  D_mat = cbind((M-0)^2, (M-total)^2)
  cluster2 = rowMin(D_mat)-1
  indicator2 = cluster2 + (cluster2 == 0)*epsilon - (cluster2 == 1)*epsilon
  L2 = sum(dbinom(x = M, size = total, prob = indicator2, log = T))
  
  # Model 3
  D_mat = cbind((M-0.5*total)^2, (M-total)^2)
  cluster3 = (rowMin(D_mat))/2
  indicator3 = cluster3 + (cluster3 == 0)*epsilon - (cluster3 == 1)*epsilon
  L3 = sum(dbinom(x = M, size = total, prob = indicator3, log = T))
  
  # Model 4
  D_mat = cbind((M-0)^2, (M-0.5*total)^2)
  cluster4 = (rowMin(D_mat)-1)/2
  indicator4 = cluster4 + (cluster4 == 0)*epsilon - (cluster4 == 1)*epsilon
  L4 = sum(dbinom(x = M, size = total, prob = indicator4, log = T))
  
  # Model 5
  D_mat = Reduce(cbind, list((M-0)^2, (M-0.5*total)^2, (M-total)^2))
  cluster5 = (rowMin(D_mat)-1)/2
  #p5 = mean(cluster5)
  indicator5 = cluster5 + (cluster5 == 0)*epsilon - (cluster5 == 1)*epsilon
  L5 = sum(dbinom(x = M, size = total, prob = indicator5, log = T))
  
  L_vec = c(L1 = L1, L2 = L2, L3 = L3, L4 = L4, L5 = L5)
  L_vec
}

colMax_tol <- function(mat, tol)
{
  sapply(1:ncol(mat), function(i) paste(names(which(mat[,i] > max(mat[,i]) - tol*abs(max(mat[,i])))), collapse = '_'))
}

test_region <- function(RES_i, epsilon, tol)
{
  M = RES_i[[1]]; U = RES_i[[2]]
  # Likelihoods
  L_mat = sapply(1:nrow(M), function(i) JRC_sorter(M[i,], U[i,], epsilon = epsilon))
  colnames(L_mat) = 1:ncol(L_mat)
  state_bin = colMax_tol(L_mat, tol)
  state_bin[grepl('L1', state_bin)] = 'L1'
  
  s_vec = numeric(length = 5); names(s_vec) = paste('L', 1:5, sep='')
  s_vec['L1'] = sum(L_mat[1,])
  
  for(i in paste('L', 2:5, sep=''))
  {
    where = grepl(i, state_bin)
    if(sum(where) > 0)
    {
      values = L_mat[1,]
      values[where] = L_mat[i, where]
      s_vec[i] = sum(values)
    }
    else
    {
      s_vec[i] = -Inf
    }
  }
  #s_vec = rowSums(L_mat)
  return(likelihoods = s_vec)
}

classify <- function(s_vec, tol)
{
  names(s_vec) = paste('L', 1:length(s_vec), sep = '')
  M_L = max(s_vec)
  out_category = names(which(s_vec >= M_L - tol*abs(M_L)))
  out_label = paste(out_category, collapse = '_')
  return(out_label)
}

simulate_all_conditions <- function(bin_size, f, L, epsilon, conditions, N_rep, bin_tol, mod_tol, N_iter)
{
  for(j in 1:N_iter)
  {
    print(paste(j, 'of', N_iter))
    f1_score_j = sapply(1:nrow(conditions), function(i) simulation(bin_size, f, L, n = conditions$n[i], 
                                                                   lambda = conditions$lambda[i],
                                                                   p = conditions$p[i], epsilon, N_rep, bin_tol, mod_tol))
    # bin_size, f, L, n, lambda, p, epsilon, N_rep, bin_tol, mod_tol
    conditions = cbind(conditions, f1_score_j)
  }
  return(conditions)
}


####### Visualize simulator
sim.data = simulate_UM(bin_size = 200, f = 400, L = 3000, n = 130, lambda = 5, p = 0.5, epsilon = 0.05, model = 'M4')
rowMeans(sim.data$M/(sim.data$M + sim.data$U))
Visualize_JRC(M = sim.data$M, U = sim.data$U)

####### Define parameters
# Fix
bin_size = 200
f = 400
L = 3000
epsilon = 0.05
N_rep = 50
N_iter = 5
bin_tol = 0.15
mod_tol = 0.1
# Variable
n = c(50, 75, 100, 150, 250)
lambda = c(5, 7, 10, 50, 100)
p = seq(0.05, 0.45, 0.1)
conditions = expand.grid(n, lambda, p)
colnames(conditions) = c("n", "lambda", "p")
dim(conditions) # 125   3

####### Run Simulation
start_time <- Sys.time()
set.seed(123)
RES = simulate_all_conditions(bin_size, f, L, epsilon, conditions, N_rep, bin_tol, mod_tol, N_iter)
end_time <- Sys.time()
end_time - start_time #  Time difference of 2.488155 hours

####### Visualize results
colnames(RES) = c('n', 'lambda', 'p', paste('f1_score', 1:(ncol(RES)-3)))
RES.m = reshape2::melt(RES, id.vars = c('n', 'lambda', 'p'))
RES.m$p = factor(RES.m$p)
RES.m$lambda = factor(RES.m$lambda)
p1 = ggplot(data = RES.m, mapping = aes(x = n, y = value, col = lambda)) +facet_wrap( ~ p, ncol = length(unique(RES.m$p))) +
  geom_point(alpha = 0.6) + ylim(0,1) + xlim(0, max(RES.m$n)) + geom_smooth(alpha = 0, method = "glm", method.args = list(family = 'quasibinomial', control = list(maxit = 500))) +
  theme(legend.text=element_text(size=8),
        strip.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste('L = ', L, '; f = ', f, '; bin_size = ', bin_size, '; epsilon = ', epsilon, '; N_iter = ', N_iter, '; N_rep = ', 
                N_rep, '; bin_tol = ', bin_tol, '; mod_tol = ', mod_tol, sep = '')) + ylab("Macro F1-score")


setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/')
ggsave(filename = 'simulations_JRCsorter.tiff', plot = p1, device = 'tiff', dpi = 300, width = 6.7, height = 3, units = 'in')

setwd('/media/ultron/2tb_disk1/ben/JRC_project/results/simulation_JRC_sorter/')
saveRDS(RES.m, paste(format(Sys.time(), '%Y%m%d'), 'simulation_results.Rds', sep = '_'))


##### Limitations

# 1. Region (of constant length) is properly specified except for the first and last two bins. In reality, it is a bit more complicated.
#    Also, region of interest is quite large.
# 2. Epsilon set for generation is also epsilon set for likelihood. This is not necesarily the case.
# 3. Count distribution follows a Poisson distribution independent of bin and individual. This is unlikely the case.
# 4. Performance is measured on equal number of M1, M2, M3, M4, M5 cases. We don´t expect similar proportions in practice. M1 >> M5 > M3, M4 > M2

