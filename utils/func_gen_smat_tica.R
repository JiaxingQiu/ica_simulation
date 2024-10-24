# this script contains the functions to generate tica source matrix





#' Calculate Pairwise L1 (Manhattan) Distance Between Rows of a Matrix
#'
#' This function computes the pairwise L1 (Manhattan) distance between the rows of a matrix and returns a symmetric matrix of these distances.
#'
#' @param v1 A numeric matrix where each row represents a data point and each column a feature.
#'
#' @return A symmetric numeric matrix where each entry \[i, j\] represents the normalized pairwise L1 (Manhattan) distance between row `i` and row `j` of the input matrix `v1`.
#'
#' @details The function computes the L1 distance between every pair of rows in the input matrix. The L1 distance is normalized by the number of columns (features) in the matrix, and a symmetric matrix is returned.
#'
#' @examples
#' v1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#' pw_l1_norm(v1)
#'
#' @export
pw_l1_norm <- function(v1) {
  # Initialize a zero matrix for storing the L1 norms
  l1_norm <- matrix(0, nrow = nrow(v1), ncol = nrow(v1))
  
  # Loop through the rows of the matrix
  for(i in 1:nrow(v1)) {
    for(j in 1:i) {
      # Calculate the pairwise L1 (Manhattan) distance
      l1_norm[i,j] <- (1.0 / ncol(v1)) * sum(abs(v1[i,] - v1[j,]))
    }
  }
  
  # Symmetrize the matrix
  l1_norm <- l1_norm + t(l1_norm) - diag(diag(l1_norm))
  
  # Return the modified L1 distance matrix
  return(1.0 - l1_norm)
}







#' Generate block diagonal correlation matrix
#'
#' Function to generate block diagonal correlation matrices
#'
#' @param d total number of components/sources
#' @param nb number of blocks where components are correlated
#' @param wb weight /strength of the correlation
#'
#' @return corr_mat (a square matrix containing pairwise correlation of components).
#' 
#' @examples
#' # Example usage of the function
#' corr_mat <- gen_source_corr(10, 3, 0.5)
#'
#' @export
gen_source_corr <- function(d, nb = 2, wb = 0){
  
  # Step 1: Create a diagonal matrix of size d x d
  base <- diag(d)
  
  # Step 2: Generate membership for each element by splitting into blocks
  # 'rep_membs' assigns block numbers to each element, repeating block IDs evenly across rows
  rep_membs <- sort(rep(1:nb, floor(d/nb)))
  
  # Step 3: Initialize a membership matrix of size d x nb with all zeros
  b_membs <- matrix(0, d, nb)
  
  # Step 4: Assign 1 to the corresponding block membership for each row
  for(i in 1:length(rep_membs)){
    b_membs[i, rep_membs[i]] <- 1
  }
  
  # Step 5: Compute the pairwise L1 distance between rows of the membership matrix
  corr_mat <- pw_l1_norm(b_membs)
  
  # Step 6: Set all non-diagonal entries (where corr_mat != 1) to zero
  corr_mat[which(corr_mat != 1, arr.ind = TRUE)] <- 0
  
  # Step 7: Multiply the off-diagonal values by weight 'wb' to adjust correlations
  corr_mat <- corr_mat * wb
  
  # Step 8: Set diagonal elements to 1 (correlation of an element with itself is 1)
  diag(corr_mat) <- 1
  
  # Step 9: Return the final block diagonal correlation matrix
  return(corr_mat)
}









#' Generate independent/correlated tICA sources
#'
#' @param n_source number of components/sources
#' @param n_time number of time points for each source
#' @param n_corr_block number of blocks containing correlated components, lower number of blocks means more components are correlated
#' @param w_corr weight of correlation
#' @param kurt kurtosis of each source components for strength of non-gaussianity
#'
#' @return a source matrix (smat) is a n_source x n_time matrix
#' 
#' @examples
#' # Example usage of the function
#' smat <- gen_smat_tica(10, 100, 4, 0.5, 3)
#'
#' @export
gen_smat_tica <- function(n_source, n_time, n_corr_block, w_corr, kurt){
  
  corr_mat <- gen_source_corr(n_source, n_corr_block, w_corr)
  
  smat <- semTools::mvrnonnorm(n_time, 
                        mu = rep(0,dim(corr_mat)[[1]]), 
                        Sigma = corr_mat, 
                        kurtosis = rep(kurt, dim(corr_mat)[[1]]))
    
  smat <- t(smat)
  return(smat)
}
