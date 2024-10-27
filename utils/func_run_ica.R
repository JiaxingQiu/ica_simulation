run_ica = function(gen_data_obj){
  
  ICA_res = fastICA::fastICA(t(gen_data_obj$xmat), # fastICA has this weird transposed setting
                             n.comp = dim(gen_data_obj$smat)[[1]])
  mixmat_hat = t(apply(abs(t(ICA_res$A)), MARGIN = 1, function(x){return(x/sum(x))}))
  
  return(list(smat_hat = t(ICA_res$S), 
              mixmat_hat = mixmat_hat, 
              ACI = ACI(gen_data_obj$mixmat, mixmat_hat), # adjusted concordance index (default 1000 iter), which evaluates how well our estimates of the mixing matrix match up with the original mixing matrix
              NDC = NDC(gen_data_obj$mixmat, mixmat_hat)))
  
}

iterate_ica <- function(gen_data_obj, sim_condition){
  ica_iter <- sim_condition$ica_iter
  
  # under the same simulation condition, run fastICA ica_iter times
  ACI_list <- list() 
  pb = progress_bar$new(total = ica_iter, format = "[:bar] :current/:total eta: :eta")
  for(j in 1:ica_iter ){
    ica_run <- run_ica(gen_data_obj)
    ACI_list[[j]] = cbind(data.frame(sim_condition), 
                          data.frame(ica_iter_id = j, 
                                     ACI = ica_run$ACI, NDC = ica_run$NDC))
    pb$tick()
  }
  to_return = as.data.frame(do.call("rbind", ACI_list))
  return(to_return)
}
