# this script contains the functions to generate sica source matrix
library(reshape2)
library(plotly)



# Function to generate x, y, z centers based on overlap
gen_xyz_centers <- function(n_source, overlap) {
  # function gen_xyz_centers that interprets the overlap parameter 
  # as defining the volume of a smaller cube, 
  # centered at (0.5,0.5,0.5) within the larger unit cube [0,1]×[0,1]×[0,1]. 
  # The centers of the 3D blobs will be sampled from this smaller cube.
  
  # n_source: Number of sources (or 3D blobs)
  # overlap: A scalar between 0 and 1, where a larger value indicates closer centers
  
  # The side length of the smaller cube, based on the overlap value
  # We interpret the 1-overlap as a volume fraction of the smaller cube inside the unit cube
  side_length <- (1-overlap)^(1/3)  # Side length of the small cube to achieve the overlap volume
  
  # Compute the bounds of the smaller cube, centered at (0.5, 0.5, 0.5)
  half_side <- side_length / 2
  lower_bound <- 0.5 - half_side
  upper_bound <- 0.5 + half_side
  
  # Initialize a matrix to store the x, y, z centers
  centers <- matrix(0, nrow = n_source, ncol = 3)
  
  # Generate random centers within the small cube
  centers <- matrix(runif(n_source * 3, min = lower_bound, max = upper_bound), ncol = 3)
  
  
  return(centers)
}
# # Example usage
# centers <- gen_xyz_centers(5, 0.9)
# print(centers)
# centers_df <- data.frame(x = centers[,1], y = centers[,2], z = centers[,3])
# tick_vals <- seq(0, 1, by = 0.2)
# fig <- plot_ly(centers_df, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
#                marker = list(size = 8, color = 'red')) %>%
#   layout(title = '3D Centers of Blobs', scene = list(
#     xaxis = list(title = 'X', range = c(0, 1), tickvals = tick_vals),
#     yaxis = list(title = 'Y', range = c(0, 1), tickvals = tick_vals),
#     zaxis = list(title = 'Z', range = c(0, 1), tickvals = tick_vals)
#   ))
# fig









# Function to create a 3D Gaussian blob
create_gaussian_blob <- function(dim_voxel, 
                                 x_center, y_center, z_center, 
                                 sigma=dim_voxel/15, 
                                 amplitude=100) {
  
  x_dim <- dim_voxel  # Voxel grid dimensions in x
  y_dim <- dim_voxel  # Voxel grid dimensions in y
  z_dim <- dim_voxel  # Voxel grid dimensions in z
  voxel_grid <- array(0, dim = c(x_dim, y_dim, z_dim)) # Initialize the voxel grid with zeros
  for (x in 1:x_dim) {
    for (y in 1:y_dim) {
      for (z in 1:z_dim) {
        distance_squared <- (x - x_center)^2 + (y - y_center)^2 + (z - z_center)^2
        voxel_grid[x, y, z] <- voxel_grid[x, y, z] + amplitude * exp(-distance_squared / (2 * sigma^2))
      }
    }
  }
  
  # Normalize the intensity values between 0 and 1 for use in opacity scaling
  voxel_grid <- (voxel_grid - min(voxel_grid)) / (max(voxel_grid) - min(voxel_grid))
  # Remove < 0.01 values
  voxel_grid[voxel_grid<0.01] <- 0
  return(voxel_grid)
}

# Function to plot a 3D Gaussian blob
plot_gaussian_blob <- function(voxel_grid, plotlower=T){
  # Convert the 3D array into a data frame for 3D plotting
  voxel_df <- melt(voxel_grid)
  names(voxel_df) <- c("x", "y", "z", "intensity")
  
  if(plotlower){
    voxel_df_plot <- voxel_df[voxel_df$x + voxel_df$y + voxel_df$z <= dim(voxel_grid)[[1]] ,]
  }else{
    voxel_df_plot <- voxel_df
  }
  
  # # Create a color scale with transparency (using hex codes)
  # color_scale <- scales::col_numeric(
  #   palette = c("#0000FF00", "#FF0000FF"), # Transparent white to opaque red
  #   domain = c(0, 1)  # Map intensity_scaled from 0 (transparent) to 1 (opaque)
  # )
  
  # 3D scatter plot using plotly with transparency for lower intensities
  fig <- plot_ly(
    voxel_df_plot, 
    x = ~x, 
    y = ~y, 
    z = ~z, 
    marker = list(size = 3),
    color = ~intensity, 
    # alpha = 0.5,
    # colorscale = color_scale(voxel_df$intensity),
    # colors = c('white', 'red'),
    type = "scatter3d",
    mode = "markers"
  )
  print(fig)
}
# # Set parameters for different brain regions modeled as Gaussian blobs
# # Example: Two Gaussian blobs representing different brain regions
# voxel_grid1 <- create_gaussian_blob(64, x_center = 32, y_center = 32, z_center = 16)  # Region 1
# plot_gaussian_blob(voxel_grid1)
# voxel_grid2 <- create_gaussian_blob(64, x_center = 20, y_center = 45, z_center = 10)    # Region 2
# plot_gaussian_blob(voxel_grid2)
# voxel_grid3 <- create_gaussian_blob(64, x_center = 45, y_center = 20, z_center = 12)   # Region 3
# plot_gaussian_blob(voxel_grid3)
# # combined sources
# voxel_grid <- voxel_grid1 + voxel_grid2 + voxel_grid3
# plot_gaussian_blob(voxel_grid, plotlower = T)
  








#' Generate independent/correlated sICA sources
#'
#'
#' To simulate a 3D Gaussian brain image, we can create a 3D grid of voxels and assign Gaussian-distributed intensities. 
#' In this simulation, we will model the brain as a combination of several Gaussian blobs with different parameters 
#' (e.g., centers and standard deviations), representing different anatomical regions or intensity variations.
#' 
#' @param n_source number of components/sources
#' @param dim_voxel the dim length (width/height) of the cube of 3d voxels
#' @param overlap the strength of overlapping (0-1)
#' @param w_corr weight of correlation
#' @param kurt kurtosis of each source components for strength of non-gaussianity
#'
#' @return a source matrix (smat) is a n_source x n_time matrix
#' 
#' @examples
#' # Example usage of the function
#' smat <- gen_smat_sica(5, 50, 0.9)
#'
#' @export
gen_smat_sica <- function(n_source, dim_voxel, overlap){
  # prepare the blob centers for the sources
  centers <- gen_xyz_centers(n_source, overlap)
  centers <- centers * dim_voxel
  # create the 3d blob signal
  voxel_grid_list <- lapply(1:nrow(centers), function(i) {
    create_gaussian_blob(dim_voxel,
                         x_center = centers[i, 1],
                         y_center = centers[i, 2],
                         z_center = centers[i, 3])
  })
  
  # Flatten each voxel grid into a vector
  voxel_vectors <- lapply(voxel_grid_list, as.vector)
  # Concatenate the vectors into an n_source by voxels matrix
  smat <- do.call(rbind, voxel_vectors)
  
  # voxel_grid <- Reduce(`+`, voxel_grid_list)
  # plot_gaussian_blob(voxel_grid, plotlower = T)
  
  return(smat)
}






jaccard_overlap <- function(smat, threshold = 0.01){
  binary_smat <- smat >= threshold
  voxel_contributions <- colSums(binary_smat)
  overlapping_voxel_count <- sum(voxel_contributions > 1)
  fired_voxel_count <- sum(voxel_contributions > 0)
  frac_overlap_voxels <- overlapping_voxel_count/fired_voxel_count
  return(frac_overlap_voxels)
}
# f <- c()
# for(i in 1:100){
#   smat <- gen_smat_sica(3, 50, 0.95)
#   f <- c(f, jaccard_overlap(smat))
# }
# hist(f)



