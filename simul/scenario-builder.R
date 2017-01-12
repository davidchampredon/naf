##################################################
###
###  GENERATE PARAMETERS LIST FOR ALL SCENARIOS
###
##################################################


# Specify the values of the parameters
# that will change across the scenarios:

M <- list()

M[['interv_start']]        <- c(-56,-28,0,28)
M[['interv_cvg_rate']]     <- c(0.0010, 0.0030, 0.0060,0.0120)  #c(0.0010, 0.0030, 0.0100)
M[['interv_target']]       <- c('susceptible','young_old')
M[['interv_cvg_max_prop']] <- 0.40 #c(0.3, 0.8)
M[['interv_efficacy']]     <- 0.80 #c(0.2, 0.8)
M[['contact_rate_mean']]   <- c(3.0)   # c(1.8, 3.6)
M[['imm_hum_baseline']]    <- 0.1   

x <- as.data.frame(expand.grid(M))

# Parameters that are linked to other 
# (but do not enter in the combination):
# << add more here if necessary... >>

# Export to CSV file:
x$scenario_id <- c(1:nrow(x))
write.csv(x = x, 
		  file = 'scenario-prm-list.csv',
		  quote = F,
		  row.names = F)

print(paste('Scenarios builder -->',nrow(x),'scenarios generated.'))


