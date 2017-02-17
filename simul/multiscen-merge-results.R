###
###   ANALYZE & COMPARE MULTIPLE SCENARIOS PREVIOUSLY RUN
###

library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)

source('utils-analyze.R')
source('utils-compare.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]


# Define 'scen.id' by retrieving 
# all successful simulations:
fname <- 'result-scen-'
tmp     <- system(command = paste0('ls ',dir.save.rdata,fname,'*.RData'),
				  intern = TRUE)
tmp     <- gsub(x = tmp, pattern = paste0(dir.save.rdata,fname),replacement = '')
scen.id.vec <- gsub(x = tmp, pattern = '.RData',   replacement = '')

message(paste('Merging ',length(scen.id.vec),'RData files ...'))
# Merge results from all scenarios
# in file 'result-scen-all.RData':
res.all <- merge.result.scen(scen.id.vec, dir.save.rdata)

msg.ac <- 'Analysis and/or comparison of multi-scenarios done.'
print(msg.ac)
message(msg.ac)
