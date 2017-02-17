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

args <- commandArgs(trailingOnly = TRUE)
# Compare intervention of this scenario
# with common baseline:

compare.simul.scen(scen.id = args[1],
                   dir.save.rdata = dir.save.rdata,
                   dir.results = dir.results)


