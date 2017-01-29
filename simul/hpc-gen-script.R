###
### GENERATE SCRIPT TO LAUNCH IN HPC
###

args <- commandArgs(trailingOnly = TRUE)
wall.time <- args[1]
# wall.time <- '00:10:00'

# Clean up previous scripts:
system('rm -rf naf-*.sh')

# Generate scenarios:
source('scenario-builder.R')
fs <- read.csv('scenario-prm-list.csv')
ns <- nrow(fs)

debug.local <- FALSE

# Generate simulation launch scripts:
for(i in 1:ns){
	scriptname <- paste0('naf-',i,'.sh')
	
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
	x[4] <- paste0('Rscript multiscen-MAIN-hpc.R ',i)
	
	if(debug.local) {
		x[2] <- ''
		x[3] <- ''
	}
	
	write(x,file = scriptname)
	system(paste0('chmod +x naf-',i,'.sh'))
}

# Send the jobs to the queuing system:
print(paste('Launching',ns,'job with wall.time.max =',wall.time,'...'))
for (i in 1:ns){
	cmd <- paste0('qsub naf-',i,'.sh')
	if(debug.local) cmd <- paste0('./naf-',i,'.sh')
	system(cmd,intern = FALSE)
}