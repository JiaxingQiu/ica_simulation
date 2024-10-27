# pw_l1_norm is defined in gen_tica

NDC <- function(v1, v2) {
  # The NDC function (Normalized Distance Consistency) computes a measure of similarity 
  # between two matrices based on pairwise L1 norms.
  
  # Compute the pairwise L1 norms for v1 and v2
  e1 <- pw_l1_norm(v1) 
  e2 <- pw_l1_norm(v2)
  
  # Extract the lower triangular part (excluding diagonal) for both matrices
  e1v <- e1[lower.tri(e1, diag = F)]
  e2v <- e2[lower.tri(e2, diag = F)]
  
  # Compute the normalized distance consistency
  ndc <- 1.0 - mean(abs(e1v - e2v))
  
  return(ndc)
}



ENDC <- function(v1, v2, iter = 1000) {
  # This function (Expected NDC) computes the expected NDC by shuffling the elements
  # of two pairwise L1 norms and taking the average over a specified number of iterations.
  
  # Compute the pairwise L1 norms for v1 and v2
  e1 <- pw_l1_norm(v1)
  e2 <- pw_l1_norm(v2)
  
  # Extract the lower triangular part (including diagonal) for both matrices
  e1v <- e1[lower.tri(e1, diag = F)]
  e2v <- e2[lower.tri(e2, diag = F)]
  
  # Initialize a vector to store results from each iteration
  temp <- numeric(iter)
  
  # Perform the iterations
  for (i in 1:iter) {
    # Shuffle the elements of the lower triangular vectors
    v1s <- sample(e1v)
    v2s <- sample(e2v)
    
    # Calculate the NDC for the shuffled vectors
    temp[i] <- 1.0 - mean(abs(v1s - v2s))
  }
  
  # Return the expected NDC (average over all iterations)
  return(mean(temp))
}


ACI <- function(v1, v2, iter = 1000) {
  # The ACI function calculates the adjusted concordance index based on NDC and ENDC.
  
  # Calculate NDC (Normalized Distance Consistency)
  ndc <- NDC(v1, v2)
  
  # Calculate ENDC (Expected NDC) with the specified number of iterations
  endc <- ENDC(v1, v2, iter = iter)
  
  # Calculate and return the adjusted concordance index (ACI)
  return((ndc - endc) / (1.0 - endc))
}

# Example usage:
# v1 and v2 should be your matrices
# result <- ACI(v1, v2, iter = 1000)

