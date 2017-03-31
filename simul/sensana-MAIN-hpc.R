##################################################################
######
######    SCRIPT TO RUN ONE SCENARIO for HPC environment
######
######    NOTE:
######	  - 'scenario-builder-sensana.R' must be run _before_ this script.
######    - 'multiscen-analyze.R' and 'multiscen-plot-only.R' 
######      can be run _after_ this script.
######
##################################################################

master0 <- as.numeric(Sys.time())

# Arguments from the command line:
args    <- commandArgs(trailingOnly = TRUE)
scenidx <- as.numeric(args[1])
print(paste('Running sensitivity analysis scenario #',scenidx,'...'))

# Run the simulations on the selected scenario:
source('utils-load-param.R')
source('utils-run.R')
source('utils-misc.R')
source('utils-fit.R')

do.fit.R0 <- FALSE

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]

# File defining parameters for various scenarios:
scen.list.file <- 'scenario-prm-list-sensana.csv'
# Overwrite parameter values
# associated with current scenario:
ov <- overwrite.selected.param(filename   = scen.list.file,
							   scen.id    = scenidx,
							   prm        = prm, 
							   interv.prm = interv.prm)
prm        <- ov[['prm']]
interv.prm <- ov[['interv.prm']]

if(do.fit.R0){
# Re-fit mean contact rate to pre-specified R0:
    cr.mean.vec <- seq(1, 9, by = 0.5)
    fitted.MCR <- fit.cr.R0(R0.target = 1.4, 
                            cr.mean.vec = cr.mean.vec, 
                            scenidx = scenidx,
                            prm = prm, 
                            interv.prm.0 = interv.prm)
}
if(!do.fit.R0){
    # If no refit requested, just use the value from the file.:
    fitted.MCR <- prm[['contact_rate_mean']]
}

# overwrite the fitted value:
if(!is.na(fitted.MCR)){
    prm[['contact_rate_mean']] <- fitted.MCR
    # Run the simulation for that scenario:
    run.simul(scen.id = scenidx, 
              dir.save.rdata = dir.save.rdata,
              force.light.output = TRUE)
}
master1 <- as.numeric(Sys.time()) 
msgt    <- paste("===> Time elapsed for sensi. analysis scenario #",
                 scenidx,' : ',round((master1-master0)/60,1),'minutes.')
print(msgt)
message(msgt)
