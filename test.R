# continued from SimulationCode.R line 100
results <- run_condition_set_rslurm(rslurm_cond_list50n[[1]])

library(dplyr)
source("./helper_functions.R")
over_df <- data.frame()
for(o in seq(0,1,0.1)){
  mix_mat <- gen_mixmat_norm(100,36,overlap=o)
  fl <- list()
  par(mfrow=c(6,6))
  for(i in 1:nrow(mix_mat)){
    fl[[i]] <- matrix(mix_mat[i,],ncol=sqrt(ncol(mix_mat)))
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


