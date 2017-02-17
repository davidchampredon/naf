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
wall.time <- args[1]

# Define 'scen.id' by retrieving 
# all successful simulations:
tmp     <- system(command = paste0('ls ',dir.save.rdata,'mc-simul-*.RData'),
				  intern = TRUE)
tmp     <- gsub(x = tmp, pattern = paste0(dir.save.rdata,'mc-simul-'),replacement = '')
scen.id.vec <- gsub(x = tmp, pattern = '.RData',   replacement = '')
ns <- length(scen.id.vec)

message(paste('Launching',ns,'jobs...'))

debug.local <- FALSE

# Generate simulation launch scripts:
script.name <- 'naf-analysis-'
for(i in 1:ns){
    scriptname <- paste0(script.name,i,'.sh')
    
    x <- character()
    x[1] <- '#!/bin/bash'
    x[2] <- paste0('#PBS -l walltime=',wall.time)
    x[3] <- '#PBS -l nodes=1:ppn=1
    #PBS -r n
    #PBS -m a
    #PBS -M david.champredon@gmail.com
    
    module load intel64/15.3.187 xz/5.2.2 gcc/5.2.0
    module load gcc bioinformatics/R/3.2.5 bioinformatics/Bioconductor/3.2
    
    cd /home/champrd/github/naf/simul
    '
    x[4] <- paste0('Rscript multiscen-analyze-unit.R ',i)
    
    if(debug.local) {
        x[2] <- ''
        x[3] <- ''
    }
    
    write(x,file = scriptname)
    system(paste0('chmod +x ',script.name,i,'.sh'))
}

# Send the jobs to the queuing system:
print(paste('Launching',ns,'job with wall.time.max =',wall.time,'...'))
for (i in 1:ns){
    cmd <- paste0('qsub ',script.name,i,'.sh')
    if(debug.local) cmd <- paste0('./',script.name,i,'.sh')
    system(cmd,intern = FALSE)
}



