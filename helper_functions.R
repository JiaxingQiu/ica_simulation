
#Function to generate block diagonal correlation matrices

gen_mat = function(d, nb = 2, wb = 0){
  
  base = diag(d)
  
  rep_membs = sort(rep(1:nb, floor(d/nb)))
  
  b_membs = matrix(0, d, nb)
  
  for(i in 1:length(rep_membs)){
    b_membs[i,rep_membs[i]] = 1
  }
  
  corr_mat = pw_l1_norm(b_membs)
  corr_mat[which(corr_mat != 1, arr.ind = T)] = 0
  corr_mat = corr_mat*wb
  diag(corr_mat) = 1
  return(corr_mat)
  
}


#Function to generate mixing matrix

gen_mixmat = function(n_source, n_signal, zero_out = 1){
  
  mix_mat = matrix(runif(n_source*n_signal, 0, 1), n_signal, n_source)
  
  if(zero_out !=0){
    for(i in 1:n_signal){
      mix_mat[i,sample(1:n_source, zero_out, replace = F)] = 0
    }
  }
  mix_mat = t(apply(mix_mat, MARGIN = 1, function(x){return(x/sum(x))}))
  
  return(mix_mat)
}



#Function to generate mixing matrix using 2d Gaussian PDF

gen_mixmat_norm <- function(n_source, n_signal, overlap=0.3, tag01=F){
  #' Generate correlated loadings using guassian pdf within each 2D filter
  #' 
  #' @param n_signal A number of filters.
  #' @param n_source A number of filter, must be squared numbers
  #' @param overlap An overlapping factor, scalar between 0 - 0.5, higher value means higher level of overlapping
  #' @returns a list of a mixing matrix of stacked flattened filters and a list of 2d filters.
  
  n_source_sqrt = sqrt(n_source)
  stopifnot(n_source_sqrt == round(n_source_sqrt))
  print("using gen_mixmat_norm")
  library(SparseM)
  # ---- fix boundaries for the overlapping factor ----
  if(overlap<=0.05){overlap=0.05} 
  if(overlap>=0.95){overlap=0.95} 
  
  # ---- m ----
  # the coordinates of m will be sampled from a 2D uniform distribution
  m_w <- ifelse(round((1-overlap)*n_source_sqrt)>0, round((1-overlap)*n_source_sqrt), 1) # width of the uniform square
  m_w_h <- ceiling(m_w/2) # half of the square width
  if(m_w_h == n_source_sqrt-m_w_h){
    s <- rep(m_w_h,2)
  }else{
    s <- seq(m_w_h,n_source_sqrt-m_w_h,1)
  }
  stopifnot(length(s)>=2)
  #TRH: Why are we sampling the center of the uniform? Why not just have the center be at the very center of the grid?
  #i.e. m_m = rep(n_source_sqrt/2,2)
  #Joy: good point, I wanted to randomize the location of the uniform, which should be the same as setting the location at the very center. 
  # m_m <- sample(s, 2, replace = T) # center of the uniform square
  m_m <- rep(round(n_source_sqrt/2),2)
  m <- cbind( runif(n_signal, min = max(1,m_m[1]-(m_w_h-1)), max = min(n_source_sqrt,m_m[1]+m_w_h)),
              runif(n_signal, min = max(1,m_m[2]-(m_w_h-1)), max = min(n_source_sqrt,m_m[2]+m_w_h)) )
  
  # ---- varcov (sigma) ----
  vars <- cbind(runif(n_signal, min=0.01*overlap*n_source_sqrt, max=0.1*overlap*n_source_sqrt),#(exp(overlap)-0.9)
                  runif(n_signal, min=0.01*overlap*n_source_sqrt, max=0.1*overlap*n_source_sqrt))
  
  
  # ---- generate filters ----
  fl <- list() 
  for(i in 1:n_signal){
    sigma <- diag(2)
    sigma[1,1] <- vars[i,1]
    sigma[2,2] <- vars[i,2]
    sigma[1,2] <- sigma[2,1] <- sample(c(-1,1),1)*runif(1, 
                                                        min=0,
                                                        max=min(sigma[1,1], sigma[2,2]))
    get_p <- function(vec) {  
      p <- mvtnorm::pmvnorm(vec, vec+1, m[i,], sigma=sigma)
      stopifnot(all(!is.na(p)))
      return(p) 
    }
    #TRH: Let's make the simplifying assumption that we are evaluating the normal PDF at the center of each grid cell:
    # for example, for the cell: voxels[1,1], the center would be c(.5, .5)
    #This would change the generating statement to voxels <- as.data.frame(expand.grid(seq(.5,n_source_sqrt-.5,1), seq(.5,n_source_sqrt-.5,1)))
    #Joy: it's essentially the same, I am using the grid coordinates (the lower left point of a unit (1x1) square) directly, rather than the center of the unit square.
    voxels <- as.data.frame(expand.grid(seq(1,n_source_sqrt,1), seq(1,n_source_sqrt,1)))
    voxels$p <- apply(voxels,1,get_p) # it takes 3-4 sec to calculate 100*100 pdf 
    if(any(!is.finite(voxels$p))) stop("infinite or missing values in 'voxels$p'")
    
    colnames(voxels) <- c("x", "y", "p")
    #TRH: No need to use a sparse matrix here, in fact, the matrix is not going to be sparse at all.
    #Joy: thanks for pointing this out, I agree, we don't have to use sparse matrix here, 
    # however, there are 2 reasons why I consider using it:
    # 1) the normal PDF can produce many zeros (numerical) when it is far from the mean. 
    # sparse matrix can be used for memory efficiency for intermediate FL object, 
    # we eventually can generate the mixing matrix from fl easily and do surgery on mix_mat directly.
    # 2) in case of we want to use circle (set tag01 = T) instead of normal PDF.
    # another purpose of using sparse matrix is that I can set non-zero cells as 1, and keep zero cells as 0,
    # a circle tagged with 1s inside and 0s outside is returned for each filter.
    fl[[i]] <- as.matrix.coo(matrix(voxels$p, nrow = n_source_sqrt, byrow = TRUE)) # store in sparse matrix format to save storage
    if(tag01){fl[[i]]@ra=rep(1,length(fl[[i]]@ra))}
    #TRH: I wouldn't create the mixing matrix quite yet
    # Joy: Yes, I have moved this part outside the loop
  }
  #This is where you should normalizing the maps, as FL should be a list of square matrices
  #You don't need any of the adjustments to the mix_mat if you normalized the set of maps.
  #Joy: Please correct me if I understand you incorrectly!
  # I have 2 questions regarding normalizing -- let's say we have normal PDF maps all centered at the very center of the map.
  # 1) only cells in the middle are non-zero, fringing cells are numerically zero ( though in theroy pdf should be >0)
  # thus, when I normalize fringing cells, i have (0+0)/0 issue, which causes infinite and NA
  # 2) normalizing non-zero cells will remove the "normal PDF" spacial structure. 
  # i.e. here are some code to test see what normalization does to 2 overlapping maps.
  # ------------------------------------------------------------------
  # # step 1: set up parameters like this:
  # n_source = 196
  # n_signal = 2
  # overlap = 1
  # # step 2: run codes above step 1
  # # step 3: plot each step
  # par(mfrow=c(3,4))
  # library('plot.matrix')
  # plot(as.matrix(fl[[1]]), main="map 1")
  # hist(c(as.matrix(fl[[1]])), main="histogram map 1 loadings")
  # plot(as.matrix(fl[[2]]), main="map 2")
  # hist(c(as.matrix(fl[[2]])), main="histogram map 2 loadings")
  # .list <- lapply(fl,as.matrix)
  # mean_matrix <- Reduce('+', .list)/2
  # sd_matrix <- matrix(mapply(function(x,y) sd(c(x,y)),.list[[1]], .list[[2]]), ncol=ncol(.list[[2]]))
  # plot(mean_matrix, main="cell-wise mean")
  # hist(c(mean_matrix), main="histogram cell-wise mean")
  # plot(sd_matrix, main="cell-wise sd")
  # hist(c(sd_matrix), main="histogram cell-wise sd")
  # normalize <- function(map){(as.matrix(map) - mean_matrix)/sd_matrix }
  # norm_fl <- lapply(fl,normalize)
  # plot(as.matrix(norm_fl[[1]]), main="normalized map 1")
  # hist(c(as.matrix(norm_fl[[1]])), main="histogram normalized map 1 loadings")
  # plot(as.matrix(norm_fl[[2]]), main="normalized map 2")
  # hist(c(as.matrix(norm_fl[[2]])), main="histogram normalized map 2 loadings")
  # ------------------------------------------------------------------
  
  
  
  # construct mixing matrix
  mix_mat <- matrix(unlist(lapply(fl,as.matrix)), ncol = n_source, byrow = TRUE)
  
  mix_mat = mix_mat*10 # scale up first
  mix_mat_noise = matrix(runif(n_source*n_signal, 0, 0.1), n_signal, n_source)
  mix_mat = mix_mat + mix_mat_noise # cannot have too many zeros in mixing matrix
  mix_mat[mix_mat>2] <- 2 # set boundaries
  mix_mat[mix_mat<0] <- 0
  if(any(!is.finite(mix_mat))) stop("infinite or missing values in 'mix_mat'")
  
  return(mix_mat) # return( list(mix_mat = mix_mat, fl = fl) ) 
}


#Function for generating the sources and the signals

gen_data = function(source_corr, mixmat, tp, kurt){
  
  sources <- mvrnonnorm(tp, mu = rep(0,dim(source_corr)[[1]]), Sigma = source_corr, kurtosis = rep(kurt, dim(source_corr)[[1]]))
  if(all(source_corr == diag(dim(source_corr)[[1]]))){
    sources <- whiten(sources)
  }
  signals = t(mixmat %*% t(sources))
  
  return(list(sources = sources, signals = signals, source_corr = source_corr, mixmat = mixmat))
  
}




run_ica = function(gen_data_obj, iter = 25, start_vals = NA){
  
  if(is.na(start_vals)){
    ICA_res = fastICA::fastICA(gen_data_obj$signals, n.comp = dim(gen_data_obj$sources)[[2]])
  }else{
    
    ICA_res = fastICA::fastICA(gen_data_obj$signals, n.comp = dim(gen_data_obj$sources)[[2]])
    
  }
  norm_A = t(apply(abs(t(ICA_res$A)), MARGIN = 1, function(x){return(x/sum(x))}))
  return(list(org_data = gen_data_obj,
              mixmat = t(ICA_res$A), 
              est_sources = ICA_res$S, 
              norm_mixmat = norm_A, 
              ACI = ACI(gen_data_obj$mixmat, norm_A, iter = iter), 
              NDC = NDC(gen_data_obj$mixmat, norm_A)))
  
}



iterate_ica = function(condition_vector,mixmat, model_num, start_iter = 100, aci_iter = 100){
  
  ACI_list = list()
  pb = progress_bar$new(total = start_iter, format = "[:bar] :current/:total eta: :eta")
  counter = 1
  data = gen_data(gen_mat(as.numeric(condition_vector["n_source"]), as.numeric(condition_vector["n_comm"]), as.numeric(condition_vector["comm_weight"])),
                  mixmat, as.numeric(condition_vector["n_timepoints"]), kurt = as.numeric(condition_vector["kurt"]))
  
  for(j in 1:start_iter){
    ica_run = run_ica(data,iter = aci_iter)
    ACI_list[[counter]] = cbind(data.frame(condition_vector), data.frame(model_num = model_num , start_iter = j, ACI = ica_run$ACI, NDC = ica_run$NDC))
    counter = counter + 1
    pb$tick()
  }  
  
  
  to_return = as.data.frame(do.call("rbind", ACI_list))
  return(to_return)
}


extract_multimodaltest = function(x){
  
  return(modetest(x)$p.value)
  
}

run_condition_set <- function(n_sources, n_signals, zero_out = 1 ,condition_set, model_iter=10, start_iter = 1000, aci_iter = 1000){
  res_list = list()
  condition_set$n_source = n_sources
  condition_set$n_signal = n_signals
  condition_set$zero_out = zero_out
  counter = 1
  for(i in 1:model_iter){
    mixmat = run_condition_set(n_sources, n_signals, zero_out)
    cat("Mixing Matrix ", i, " generated with ", n_sources, " sources, ", n_signals, " signals and ", zero_out, " knockouts\n", sep = "")
    for(j in 1:nrow(condition_set)){
      print(condition_set[j,])
      res_list[[counter]] = iterate_ica(condition_set[j,],mixmat = mixmat, model_num = i, start_iter = start_iter, aci_iter = aci_iter)
      counter = counter + 1
    }
  }
  return(do.call("rbind", res_list))
  
}




                     
