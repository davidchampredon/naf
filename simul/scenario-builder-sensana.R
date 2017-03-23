###
###  GENERATE PARAMETERS LIST FOR ALL SCENARIOS 
###  FOR SENSITIVITY ANALYSIS
###

# Specify the values of the parameters
# that will change across the scenarios:

# list 'M': Combination of all parameters
M <- list()
M[['interv_start']]        <- c(-28,14)
M[['interv_cvg_rate']]     <- c(0.003)   #c(0.0030, 0.0120)  
# M[['interv_efficacy']]     <- c(0.4, 0.9)
x <- as.data.frame(expand.grid(M))

# list 'A': Change one value at the time (NO combination!) 
#
# *** W A R N I N G *** 
# baseline value must be the first one for each vector!
#
A <- list()
A[['imm_hum_baseline']]    <- c(0.001, 0.2)
A[['contactAssort_lambda']]<- c(0.1, 0.2) #c(0.1, 0.05, 0.2)
# A[['frailty_sd']]          <- c(0.1, 0.02, 0.5)
# A[['imm_cell_max']]        <- c(0.45, 0.2, 0.9)
# A[['contact_rate_CV']]     <- c(0.75, 0.3, 1.5)
# A[['proba_move']]          <- c(0.98, 0.8)
# A[['asymptom_infectiousness_ratio']]<- c(0.1, 0.05, 0.2)

# # pubT prop
# proba_change_sp_other 

# ---- Create one-time-change matrices ----

# retrieve all baseline values
# (they are the first element of every vector)
baseline <- unlist(lapply(A, `[[`, 1))

# Create the "root" row vector (will be duplicated) 
tmp <- list()
n.A <- length(A)
cnt <- 1
for(i in 1:n.A){
    n.Ai <- length(A[[i]])
    if(n.Ai>1){
        for(j in 2:n.Ai){
            baseline.copy <- baseline
            baseline.copy[i] <- A[[i]][j]
            tmp[[cnt]] <- baseline.copy
            cnt <- cnt + 1
        }
    }
}
w <- do.call('rbind', tmp)
w <- rbind(baseline,w)

# Duplicate the "root" row vector,
# and merge it with the combination matrix "x"
tmp <- list()
for(i in 1:nrow(w)){
    tmp[[i]] <- t (matrix(w[i,], nrow = ncol(w), ncol=nrow(x)))
    tmp[[i]] <- cbind(x, tmp[[i]])
}
z <- do.call('rbind', tmp)


# ---- Export to CSV file ----
z$scenario_id <- c(1:nrow(z))
write.csv(x = z, 
		  file = 'scenario-prm-list-sensana.csv',
		  quote = F,
		  row.names = F)

print(paste('Scenarios builder (sensi. analysis) -->',
            nrow(z),'scenarios generated.'))


