##################################################################
######
######    SCRIPT TO RUN ONE SCENARIO for HPC environment
######
######    NOTE:
######	  - 'scenario-builder.R' must be run _before_ this script
######    - 'multiscen-analyze.R' and 'multiscen-plot-only.R' can be run _after_ this script
######
######
##################################################################

# Arguments from the command line:
args <- commandArgs(trailingOnly = TRUE)

scen.id <- args[1]
print(paste('Running scenario #',scen.id,'...'))

master0 <- as.numeric(Sys.time())

# Run the simulations on the selected scenario:
source('multiscen-run-unit.R')
run.one.scenario(scenidx = scen.id)

master1 <- as.numeric(Sys.time()) 
msgt    <- paste("===> Time elapsed for scenario #",scen.id,' : ',round((master1-master0)/60,1),'minutes.')
print(msgt)
message(msgt)
