# please refer to paper Human Brain Mapping 2001 Calhoun for definitions of sICA and tICA.
library(fastICA)
set.seed(123)



#  ---- sICA (more voxels than timepoints) ---- 
# S: Source matrix is a 2 by 100 matrix, 2 source components, each source is a flatten "image" with 100 voxels
S <- matrix(runif(100*2, 0, 1), 2, 100)  
# A: Mixing matrix is a 5 by 2 matrix, mapping 2 components to 5 observed images (100 voxels)
A <- matrix(c(1, 0.5, 
              1, 0.4, 
              1, 0.3, 
              1, 0.2, 
              1, 0.1), 5, 2, byrow=T)
A = A / rowSums(A) # the loadings for each signal image should sum up to 1
A
# X: Signal image is a 5 by 100 matrix, 5 observed signals, each signal is a flatten image with 100 voxels
X <- A %*% S 
# for fastICA, X need to be transposed:  bc its row is observations (voxels), column is linearly-combined variables (brain image across time)
ica_result <- fastICA(t(X), n.comp = 2) 
# estimated components
A_hat <- t(abs(ica_result$A))  # Estimated mixing matrix
A_hat <- A_hat / rowSums(A_hat) # the loadings for each signal image should sum up to 1
A_hat

# # if we do a transpose on X then apply ica it's wrong
# X <- t(A %*% S)
# # sICA
# # for fastICA, X need to be transposed: 
# # bc its row is observations (voxels), column is linearly-combined variables (brain image across time)
# ica_result <- fastICA(t(X), n.comp = 2) 
# # estimated components
# A_hat <- abs(ica_result$S)  # Estimated mixing matrix
# A_hat <- A_hat / rowSums(A_hat) # the loadings for each signal image should sum up to 1
# A_hat







#  ---- tICA (more timepoints than voxels) ---- 
# S: Source matrix is a 2 by 200 matrix, 2 source components, each source is a time series with 200 time points
S <- matrix(runif(200*2, 0, 1), 2, 200)  
# A: Mixing matrix is a 5 by 2 matrix, mapping 2 components to 5 observed time series per voxels
A <- matrix(c(1, 0.5, 
              1, 0.4, 
              1, 0.3, 
              1, 0.2, 
              1, 0.1), 5, 2, byrow=T)
A = A / rowSums(A) # the loadings for each signal image should sum up to 1
A
# X: Signal image is a 5 by 200 matrix, 5 observed signals, each signal is a time series with 200 time points
X <- A %*% S  
# for fastICA, X need to be transposed:  bc its row is observations (voxels), column is linearly-combined variables (brain image across time)
ica_result <- fastICA(t(X), n.comp = 2) 
# estimated components
A_hat <- t(abs(ica_result$A))  # Estimated mixing matrix
A_hat <- A_hat / rowSums(A_hat) # the loadings for each signal image should sum up to 1
A_hat


# # if time < voxel, tICA is wrong.
# # S: Source matrix is a 5 by 10 matrix, 5 source components, each source is a time series with 10 time points
# S <- matrix(runif(10*5, 0, 1), 5, 10)  
# # A: Mixing matrix is a 100 by 5 matrix, mapping 5 components to 100 observed time series per voxels
# A <- matrix(runif(100*5, 0, 1), 100, 5, byrow=T)
# A = A / rowSums(A) # the loadings for each signal image should sum up to 1
# A[1:10,]
# 
# # X: Signal image is a 100 by 10 matrix, 100 observed signals, each signal is a time series with 10 time points
# X <- A %*% S 
# # for fastICA, X need to be transposed:  bc its row is observations (voxels), column is linearly-combined variables (brain image across time)
# ica_result <- fastICA(t(X), n.comp = 5) 
# # estimated components
# A_hat <- t(abs(ica_result$A))  # Estimated mixing matrix
# A_hat <- A_hat / rowSums(A_hat) # the loadings for each signal image should sum up to 1
# A_hat[1:10,]






# To simulate a 3D Gaussian brain image, we can create a 3D grid of voxels and assign Gaussian-distributed intensities. In this simulation, we will model the brain as a combination of several Gaussian blobs with different parameters (e.g., centers and standard deviations), representing different anatomical regions or intensity variations.
# Load necessary libraries
library(reshape2)
library(plotly)

# Set dimensions for the 3D brain image
x_dim <- 64  # Voxel grid dimensions in x
y_dim <- 64  # Voxel grid dimensions in y
z_dim <- 64  # Voxel grid dimensions in z

# Initialize the voxel grid with zeros
voxel_grid <- array(0, dim = c(x_dim, y_dim, z_dim))

# Function to create a 3D Gaussian blob
create_gaussian_blob <- function(voxel_grid, x_center, y_center, z_center, sigma, amplitude) {
  for (x in 1:x_dim) {
    for (y in 1:y_dim) {
      for (z in 1:z_dim) {
        distance_squared <- (x - x_center)^2 + (y - y_center)^2 + (z - z_center)^2
        voxel_grid[x, y, z] <- voxel_grid[x, y, z] + amplitude * exp(-distance_squared / (2 * sigma^2))
      }
    }
  }
  return(voxel_grid)
}

# Set parameters for different brain regions modeled as Gaussian blobs
# Example: Two Gaussian blobs representing different brain regions
voxel_grid <- create_gaussian_blob(voxel_grid, x_center = 32, y_center = 32, z_center = 16, sigma = 10, amplitude = 100)  # Region 1
voxel_grid <- create_gaussian_blob(voxel_grid, x_center = 20, y_center = 45, z_center = 10, sigma = 8, amplitude = 80)    # Region 2
voxel_grid <- create_gaussian_blob(voxel_grid, x_center = 45, y_center = 20, z_center = 25, sigma = 12, amplitude = 90)   # Region 3

# Convert the 3D array into a data frame for 3D plotting
voxel_df <- melt(voxel_grid)
names(voxel_df) <- c("x", "y", "z", "intensity")
# Normalize the intensity values between 0 and 1 for use in opacity scaling
voxel_df$intensity <- (voxel_df$intensity - min(voxel_df$intensity)) /
  (max(voxel_df$intensity) - min(voxel_df$intensity))


voxel_df_plot <- voxel_df[voxel_df$x + voxel_df$y + voxel_df$z <= 64 ,]

# Create a color scale with transparency (using hex codes)
color_scale <- scales::col_numeric(
  palette = c("#0000FF00", "#FF0000FF"), # Transparent blue to opaque red
  domain = c(0, 1)  # Map intensity_scaled from 0 (transparent) to 1 (opaque)
)

# 3D scatter plot using plotly with transparency for lower intensities
fig <- plot_ly(
  voxel_df_plot, 
  x = ~x, 
  y = ~y, 
  z = ~z, 
  marker = list(size = 3, 
                color = ~intensity, 
                colorscale = color_scale(voxel_df$intensity)),
  type = "scatter3d",
  mode = "markers"
)
fig



