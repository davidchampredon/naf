##################################################################
######
######    SCRIPT TO RUN ONE SCENARIO for HPC environment
######
######    NOTE:
######	  - 'scenario-builder.R' must be run _before_ this script.
######    - 'multiscen-analyze.R' and 'multiscen-plot-only.R' 
######      can be run _after_ this script.
######
##################################################################

master0 <- as.numeric(Sys.time())

# Arguments from the command line:
args    <- commandArgs(trailingOnly = TRUE)
scenidx <- args[1]
print(paste('Running scenario #',scenidx,'...'))

# Run the simulations on the selected scenario:
source('utils-load-param.R')
source('utils-run.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]

# File defining parameters for various scenarios:
scen.list.file <- 'scenario-prm-list.csv'
# Overwrite parameter values
# associated with current scenario:
ov <- overwrite.selected.param(filename   = scen.list.file,
							   scen.id    = scenidx,
							   prm        = prm, 
							   interv.prm = interv.prm)
prm        <- ov[['prm']]
interv.prm <- ov[['interv.prm']]

# Run the simulation for that scenario:
run.simul(scen.id = scenidx, 
		  dir.save.rdata = dir.save.rdata,
		  force.light.output = TRUE)

master1 <- as.numeric(Sys.time()) 
msgt    <- paste("===> Time elapsed for scenario #",scenidx,' : ',round((master1-master0)/60,1),'minutes.')
print(msgt)
message(msgt)
