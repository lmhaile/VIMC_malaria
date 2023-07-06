# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
ctx <- context::context_save(
  'pkgs',
  packages = packages,
  package_sources = src,
  sources = 'Q:/VIMC_malaria/run_malaria_model.R'
)

config <- didehpc::didehpc_config(
  use_rrq = FALSE,
  cores = 2,
  cluster = "wpia-hn" ,
  #"fi--dideclusthn", # , "fi--didemrchnb""fi--didemrchnb"
  template = "AllNodes",
  ## use for the wpia cluster
  parallel = FALSE
)

#load context into queue
obj <- didehpc::queue_didehpc(ctx, config)
didehpc::web_login()

# run jobs ---------------------------------------------
save_dir<- 'Q:/VIMC_files/central_estimates/'

bl <-
  list(
    'parameter_filepath' = paste0(save_dir, 'baseline/input_parameters_baseline_RFP.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/baseline/raw_model_output/'
  )

int <-
  list(
    'parameter_filepath' =  paste0(save_dir, 'intervention/input_parameters_vaccine_RFP.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/'
  )

filled<- list(bl, int)

test_case_jobs <- obj$lapply(filled, run_malaria_model_rfp)

