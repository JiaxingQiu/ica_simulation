# continued from SimulationCode.R line 100
results <- run_condition_set_rslurm(rslurm_cond_list50n[[1]])
save(results, file="./tmp/results_rslurm_cond_list50n_1.R")

load("./tmp/results_rslurm_cond_list50n_1.R")
# look at ACI score from results
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
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num))) + geom_line() + facet_grid( .~ sample_size) +
  labs(y = "Median Adjusted Concordance Index", x = "Level of Dependency", color = "Mixing Matrix")
ggplot(results_agg, aes(x = dependency, y = x, color = as.factor(model_num)))+ geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1) + geom_line(size = 1) + facet_grid( .~ sample_size) +
  theme_classic() +labs(y = "Median Adjusted Concordance Index", x = "Level of Dependency", color = "Mixing Matrix") + theme(legend.position = "none", axis.text =element_text(color = "black", size = 12), strip.text = element_text(size = 12),
                                                                                                                             axis.title = element_text(size = 14)) #+ ylim(0,1)





library(dplyr)
source("./helper_functions.R")
over_df <- data.frame()
for(o in seq(0,1,0.1)){
  mix_mat <- gen_mixmat_norm(20,100,overlap=o)
  fl <- list()
  par(mfrow=c(6,6))
  for(i in 1:ncol(mix_mat)){
    fl[[i]] <- matrix(mix_mat[,i],ncol=sqrt(ncol(mix_mat)))
    library('plot.matrix')
    plot(fl[[i]], main=o)
  }
  .list <- lapply(fl,as.matrix)
  sum_matrix <- Reduce('+', .list)
  par(mfrow=c(1,1))
  plot(sum_matrix, main=o)
  
  library(datawizard)
  sum_vec <- c(sum_matrix)
  over_df <- bind_rows(over_df,
                       data.frame(o = o,
                                  m = mean(sum_vec,na.rm = TRUE), 
                                  sd = sd(sum_vec,na.rm = TRUE),
                                  skew = datawizard::skewness(sum_vec,na.rm = TRUE),
                                  kurt = datawizard::kurtosis(sum_vec,na.rm = TRUE)))
}
library(ggplot2)
p1<- ggplot(over_df, aes(x=o, y=m)) + 
  geom_line()
p2<- ggplot(over_df, aes(x=o, y=sd)) + 
  geom_line()
p3<- ggplot(over_df, aes(x=o, y=skew.Skewness)) + 
  geom_line()
p4 <- ggplot(over_df, aes(x=o, y=kurt.Kurtosis)) + 
  geom_line()
ggpubr::ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)


