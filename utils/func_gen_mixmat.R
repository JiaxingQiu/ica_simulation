# Function to generate mixing matrix
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
# gen_mixmat(5,50,0)


