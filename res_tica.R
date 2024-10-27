# sim_id: indicator of each simulation condition
# iter_id: for each simulation condition (data generated 100 times), the indicator of each iteration
# ica_iter_id: for each iteration (data generation) of each simulation condition (ica ran repeated 100 times), ica_iter_id is the indicator of each run of ica 

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rm(list = ls())

library(dplyr)
library(rslurm)
library(ggplot2)
library(tidyr)
library(ggpubr)
source("./sim_conditions.R")
res_df <- readRDS("./res/output_tica_run.RDS")
agg_df <- res_df %>% 
  filter(n_source==5) %>%
  group_by(sim_id, iter_id) %>%
  summarise(aci = median(ACI, na.rm=T))
agg_df <- merge(agg_df, tica_simulation_conditions, all.x=T, by="sim_id")
ggplot(agg_df, aes(x = as.factor(w_corr), y = aci)) + 
  geom_boxplot() + 
  labs(x = "correlation", y = "ACI", title = "tICA", subtitle = "5 sources" )


res_df <- readRDS("./res/output_sica_run.RDS")
agg_df <- res_df %>% 
  filter(n_source==5, n_signal==100, zero_out==0) %>%
  group_by(sim_id, iter_id) %>%
  summarise(aci = median(ACI, na.rm=T))
agg_df <- merge(agg_df, tica_simulation_conditions, all.x=T, by="sim_id")
ggplot(agg_df, aes(x = as.factor(w_corr), y = aci)) + 
  geom_boxplot() + 
  labs(x = "degree of overlapping", y = "ACI", title = "sICA", subtitle = "5 sources" )



