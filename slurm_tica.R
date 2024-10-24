source("./sim_run.R")
source("./sim_conditions.R")


simulation_conditions <- tica_simulation_conditions
run_wrapper <- run_wrapper_tica

sjob_tica = slurm_map(
  split(simulation_conditions, simulation_conditions$sim_id),
  run_wrapper,
  nodes=nrow(simulation_conditions),
  cpus_per_node = 1,
  jobname = "tica_run",
  submit = TRUE,
  preschedule_cores = F,
  slurm_options = c(account = "netlab", partition = "standard", time = "5-00:00:00"), 
  global_objects = lsf.str()
)
save(sjob_tica, file = "tica_run.Rdata")
