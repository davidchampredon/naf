##################################################
###
###  GENERATE PARAMETERS LIST FOR ALL SCENARIOS
###
##################################################


M <- list()

M[['interv_start']]       <- c(-28,0,28)
M[['interv_cvg_rate']]    <- c(0.0010, 0.0100)  #c(0.0010, 0.0030, 0.0100)
M[['interv_target']]      <- c('susceptible','young_old')
M[['interv_cvg_max_prop']]<- c(0.3, 0.8)
M[['vax_imm_hum_incr']]   <- c(0.7)  #c(0.2, 0.7)
M[['contact_rate_mean']]  <- 1.8   #c(1.8, 3.6)
M[['imm_hum_baseline']]   <- 0.1   #c(1.8, 3.6)


x <- as.data.frame(expand.grid(M))

# Parameters that are linked to other 
# (and do not enter the combination):
x$vax_imm_cell_incr <- x$vax_imm_hum_incr

x$contact_rate_stddev <- 0.4
x$contact_rate_stddev[x$contact_rate_mean==3.6] <- 0.7

# Export to CSV file:

x$scenario_id <- c(1:nrow(x))

write.csv(x = x, 
		  file = 'scenario-prm-list.csv',
		  quote = F,
		  row.names = F)

print(paste('Scenarios builder -->',nrow(x),'scenarios generated.'))


