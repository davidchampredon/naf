##################################################################
######
######    MASTER SCRIPT TO RUN MULTIPLE SCENARIOS
######
######
##################################################################

master0 <- as.numeric(Sys.time())

# Generate list of parameters for all scenarios:
source('scenario-builder.R')
# Run the simulations on all scenarios:
source('multiscen-run.R')
# Analyze and compare results from these simulations:
source('multiscen-analyze.R')

master1 <- as.numeric(Sys.time()) 
msgt    <- paste("===> Time elapsed for multi-scenarios:",round((master1-master0)/60,1),'minutes.')
print(msgt)
message(msgt)
