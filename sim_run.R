# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

rm(list = ls())
library(semTools)
library(reshape2)
library(progress)
list.of.packages <- c("semTools",
                      "reshape2",
                      "progress")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib = "/sfs/qumulo/qhome/jq2uw/R/goolf/4.3")



path = paste0("./utils")
flst = list.files(path)
sapply(c(paste(path,flst,sep="/")), source, .GlobalEnv)
source("./sim_conditions.R")



# sim_condition <- sica_simulation_conditions[nrow(sica_simulation_conditions),]
run_wrapper_sica <- function(sim_condition) {
  results_list = list()
  for(i in 1:sim_condition$iter){
    
    tryCatch({
      
      # generate source
      smat <- gen_smat_sica(sim_condition$n_source,
                    sim_condition$dim_voxel,
                    sim_condition$overlap)
      # generate mixing matrix
      mixmat <- gen_mixmat(sim_condition$n_source,
                 sim_condition$n_signal,
                 sim_condition$zero_out)
      # generate signal
      xmat <- mixmat %*% smat
      # final generate data obj
      gen_data_obj <- list("smat" = smat,
                           "mixmat" = mixmat,
                           "xmat" = xmat)
      
      # iterate ica and return evalutation ACI
      results_list[[i]] <- iterate_ica(gen_data_obj, sim_condition)
      results_list[[i]]$iter_id <- i
      
    }, error = function(e){
      print(e)
      print(paste0("skip iteration ",i))
    })
  }
  results_list <- Filter(function(x) !is.null(x), results_list)
  toReturn = do.call("rbind", results_list)
  return(toReturn)
}







# sim_condition <- tica_simulation_conditions[nrow(tica_simulation_conditions),]
# sim_condition <- tica_simulation_conditions[1,]
run_wrapper_tica <- function(sim_condition) {
  results_list = list()
  for(i in 1:sim_condition$iter){
    
    tryCatch({
      
      # generate source
      smat <- gen_smat_tica(sim_condition$n_source,
                            sim_condition$n_time,
                            sim_condition$n_corr_block,
                            sim_condition$w_corr,
                            sim_condition$kurt)
      # generate mixing matrix
      mixmat <- gen_mixmat(sim_condition$n_source,
                           sim_condition$n_signal,
                           sim_condition$zero_out)
      # generate signal
      xmat <- mixmat %*% smat
      # final generate data obj
      gen_data_obj <- list("smat" = smat,
                           "mixmat" = mixmat,
                           "xmat" = xmat)
      
      # iterate ica and return evalutation ACI
      results_list[[i]] <- iterate_ica(gen_data_obj, sim_condition)
      results_list[[i]]$iter_id <- i
      
    }, error = function(e){
      print(e)
      print(paste0("skip iteration ",i))
    })
  }
  results_list <- Filter(function(x) !is.null(x), results_list)
  toReturn = do.call("rbind", results_list)
  return(toReturn)
}


