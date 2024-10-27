# refer to the sICA and tICA definitions in Human Brain Mapping 2001 Calhoun

# ----- sICA conditions ------
# Parameters
param_grid <- expand.grid(n_source = c(5),
                          dim_voxel = c(50), # number of voxels will be dim_voxel^3
                          n_signal = c(50),
                          overlap = c(0, 0.7, 0.9),
                          zero_out = c(0),
                          ica_iter = 10,
                          iter = 100 # simulation iteration
                          )
simulation_conditions <- as.data.frame(param_grid)
simulation_conditions$sim_id <- seq(1:nrow(simulation_conditions))
sica_simulation_conditions <- simulation_conditions

# ----- tICA conditions ------
param_grid <- expand.grid(n_source = c(5),
                          n_time = c(10000), # time point need to be more than signals
                          n_signal = c(50), 
                          n_corr_block = c(3),
                          w_corr = c(0, 0.5, 0.9),
                          kurt = c(3), # 0 is gaussian 3 is nongaussian
                          zero_out = c(0),
                          ica_iter = 10,
                          iter = 100 # simulation iteration
)
simulation_conditions <- as.data.frame(param_grid)
simulation_conditions$sim_id <- seq(1:nrow(simulation_conditions))
tica_simulation_conditions <- simulation_conditions


