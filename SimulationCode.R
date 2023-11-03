rm(list = ls())

####ICA Simulation
library(semTools)
library(ica)
library(pracma)
library(whitening)
library(progress)
library(lessR)
library(multimode)
library(parameters)
library(Rcpp)
library(RcppArmadillo)
library(rslurm)
library(fastICA)

rslurm_list_gen = function(mat_gen_conditions, condition_set, reps = 100){
  
  to_return = list()
  counter = 1
  for(i in 1:dim(mat_gen_conditions)){
    for(j in 1:reps){
      to_return[[counter]] = list(rep= j,
                                  seed = sample(10^6,1), 
                                  n_sources = mat_gen_conditions$n_sources[i],
                                  n_signals = mat_gen_conditions$n_signals[i],
                                  overlap = mat_gen_conditions$overlap[i],
                                  zero_out  = mat_gen_conditions$zero_out[i],
                                  model_iter = mat_gen_conditions$model_iter[i],
                                  start_iter = mat_gen_conditions$start_iter[i],
                                  aci_iter = mat_gen_conditions$aci_iter[i],
                                  condition_set = condition_set)
      counter = counter + 1
    }
    
  }
  
  return(to_return)
}  

run_condition_set_rslurm <- function(cond_list_elem){
  source("./helper_functions.R")
  sourceCpp("./l1_norm.cpp")
  condition_set = cond_list_elem[["condition_set"]]
  n_sources = cond_list_elem[["n_sources"]]
  n_signals = cond_list_elem[["n_signals"]]
  overlap = cond_list_elem[["overlap"]]
  zero_out = cond_list_elem[["zero_out"]]
  model_iter = cond_list_elem[["model_iter"]]
  start_iter = cond_list_elem[["start_iter"]]
  aci_iter = cond_list_elem[["aci_iter"]]
  
  res_list = list()
  condition_set$n_source = n_sources
  condition_set$n_signal = n_signals
  condition_set$overlap = overlap
  condition_set$zero_out = zero_out
  counter = 1
  for(i in 1:model_iter){
    # mixmat = gen_mixmat(n_sources, n_signals, zero_out)
    mixmat = gen_mixmat_norm(n_sources, n_signals, overlap)
    cat("Mixing Matrix ", i, " generated with ", n_sources, " sources, ", n_signals, " signals and ", zero_out, " knockouts\n", sep = "")
    for(j in 1:nrow(condition_set)){
      print(condition_set[j,])
      res_list[[counter]] = iterate_ica(condition_set[j,],
                                        mixmat = mixmat, 
                                        model_num = i, 
                                        start_iter = start_iter,
                                        aci_iter = aci_iter)
      counter = counter + 1
    }
  }
  return(do.call("rbind", res_list))
  
}


comm_kurt = expand.grid(c(0,.05,.1,.5), c(0,3))
design_mat = data.frame(n_sources = c(10, 50), # number of sources
                        n_signals = c(100, 144), # number of signals
                        overlap = c(0.3, 0.7), # overlap factor
                        zero_out = c(1,1), 
                        model_iter = 10, # number of runs of fastICA which generate a new mixing matrix per each run
                        start_iter = 100, # ICA returns a finite set of solutions which depend on where the starting value is. None of these solutions are correct, and the average of these solutions is also not correct.
                        aci_iter = 100) # adjusted concordance index, which evaluates how well our estimates of the mixing matrix match up with the original mixing matrix

# n_comm ?
cond_set1 = data.frame(n_comm = 1, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 100, zeroout = 1)
cond_set2 = data.frame(n_comm = 1, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 500, zeroout = 1)
cond_set3 = data.frame(n_comm = 3, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 500, zeroout = 1)
cond_set4 = data.frame(n_comm = 3, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 100, zeroout = 1)
cond_set5 = data.frame(n_comm = 1, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 1000, zeroout = 1)
cond_set6 = data.frame(n_comm = 3, comm_weight = comm_kurt[,1], kurt = comm_kurt[,2], n_timepoints = 1000, zeroout = 1)


cond_list1 = rslurm_list_gen(design_mat, cond_set1, 50)
cond_list2 = rslurm_list_gen(design_mat, cond_set2, 50)
cond_list3 = rslurm_list_gen(design_mat, cond_set3, 50)
cond_list4 = rslurm_list_gen(design_mat, cond_set4, 50)

cond_list5 = rslurm_list_gen(design_mat, cond_set5, 50)
cond_list6 = rslurm_list_gen(design_mat, cond_set6, 50)

rslurm_cond_list50n <- c(cond_list1[51:100], cond_list2[51:100], cond_list3[51:100], cond_list4[51:100])
rslurm_cond_list10n <- c(cond_list1[1:50], cond_list2[1:50], cond_list3[1:50], cond_list4[1:50])

rslurm_cond_list1000tp = c(cond_list5, cond_list6)


# syst[[1]] = system.time(test_10s100ts <- run_condition_set_rslurm(cond_list1[[1]]))
# syst[[2]] = system.time(test_50s100ts <- run_condition_set_rslurm(cond_list1[[200]]))
#  
# syst[[3]] = system.time(test_10s500ts <- run_condition_set_rslurm(cond_list2[[1]]))
# syst[[4]] = system.time(test_50s500ts <- run_condition_set_rslurm(cond_list2[[100]]))
# 
# syst[[5]] = system.time(test_10s500ts3comm <- run_condition_set_rslurm(cond_list3[[1]]))
# syst[[6]] = system.time(test_50s500ts3comm <- run_condition_set_rslurm(cond_list3[[100]]))


sjob1 <- slurm_map( rslurm_cond_list50n,run_condition_set_rslurm, jobname = "50s_run", nodes=length(rslurm_cond_list50n), cpus_per_node = 1, submit = TRUE,
                    slurm_options = c(account = "netlab", partition = "standard",time = "4:00:00"), preschedule_cores = F)

save(sjob1, file = "50s_run_sjob.Rdata")


sjob2 <- slurm_map( rslurm_cond_list10n,run_condition_set_rslurm, jobname = "10s_run", nodes=length(rslurm_cond_list10n),  cpus_per_node = 1, submit = TRUE,
                    slurm_options = c(account = "netlab", partition = "standard",time = "2:00:00"), preschedule_cores = F)

save(sjob2, file = "10s_run_sjob.Rdata")


sjob3 <- slurm_map( rslurm_cond_list1000tp, run_condition_set_rslurm, jobname = "1000tp_run", nodes=length(rslurm_cond_list10n),  cpus_per_node = 1, submit = TRUE,
                    slurm_options = c(account = "netlab", partition = "standard",time = "3:30:00"), preschedule_cores = F)

save(sjob3, file = "1000tp_run_sjob.Rdata")




results <- as.data.frame(apply(results, MARGIN = 2, FUN = unlist))
results$baselined_ACI = results$ACI

uniqs <- unique(results[,c("n_timepoints", "model_num")])

for(i in 1:nrow(uniqs)){
  
  targets = which(results$n_timepoints == uniqs$n_timepoints[i] & results$model_num == uniqs$model_num[i])
  mean_baseline_ACI = mean(results$ACI[which(results$n_timepoints == uniqs$n_timepoints[i] & results$model_num == uniqs$model_num[i] & results$comm_weight == 0.0)])
  print(mean_baseline_ACI)
  results[targets,"baselined_ACI"] = results[targets,"ACI"] - mean_baseline_ACI 
}


mm_test = aggregate(results$ACI, by = list(model_num = results$model_num, dependency= results$comm_weight, sample_size = results$n_timepoints), extract_multimodaltest)

library(ggplot2)

results_agg = aggregate(results$ACI, by = c(list(model_num = results$model_num, dependency= results$comm_weight, sample_size = results$n_timepoints)), median)
results_agg2 = aggregate(results$ACI, by = c(list(model_num = results$model_num, dependency= results$comm_weight, sample_size = results$n_timepoints)), quantile, probs =c(.25,.75))
results_agg$upper = results_agg2[,"x"][,2]
results_agg$lower = results_agg2[,"x"][,1]
results_agg$sample_size = factor(results_agg$sample_size, levels = c(100,500,1000), labels = c("n = 100", "n = 500", "n = 1000"))
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num))) + geom_line() + facet_grid( .~ sample_size) 
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num)))+ geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1) + geom_line(size = 1) + facet_grid( .~ sample_size) +
  theme_classic() +labs(y = "Median Adjusted Concordance Index", x = "Level of Dependency", color = "Mixing Matrix") + theme(legend.position = "none", axis.text =element_text(color = "black", size = 12), strip.text = element_text(size = 12),
                                                                                                                      axis.title = element_text(size = 14)) + ylim(0,1)
  
results_agg = aggregate(results2$ACI, by = c(list(model_num = results$model_num, dependency= results$comm_weight, sample_size = results$n_timepoints)), median)
results_agg2 = aggregate(results2$ACI, by = c(list(model_num = results$model_num, dependency= results$comm_weight, sample_size = results$n_timepoints)), quantile, probs =c(.25,.75))
results_agg$upper = results_agg2[,"x"][,2]
results_agg$lower = results_agg2[,"x"][,1]
results_agg$sample_size = factor(results_agg$sample_size, levels = c(100,500,1000), labels = c("n = 100", "n = 500", "n = 1000"))
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num))) + geom_line() + facet_grid( .~ sample_size) 
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num)))+ geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1) + geom_line(size = 1) + facet_grid( .~ sample_size) +
  theme_classic() +labs(y = "Median Adjusted Concordance Index", x = "Level of Dependency", color = "Mixing Matrix") + theme(legend.position = "none", axis.text =element_text(color = "black", size = 12), strip.text = element_text(size = 12),
                                                                                                                      axis.title = element_text(size = 14)) + ylim(0,1)
#Need to ensure nesting within mixing matrix generation. Rewrite the code to allow for that.




