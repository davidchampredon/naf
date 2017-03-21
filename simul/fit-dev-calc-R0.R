### 
###   CALCULATE R0 OF SIMULATIONS ALREADY RUN
###


# List simulation results 
z <- system('ls ./rdata/mc-simul-*.RData', intern = TRUE)

for(i in seq_along(z)){
    load(z[i])
    source('analysis_tools.R')  
    
    # Merge all simulations:
    n.cpu <- parallel::detectCores() - 1
    pop.all.mc <- merge.pop.mc(res.list.0, n.cpu = 1, doparallel = FALSE)
    
    # This is the time window (in days)
    # when R0 is estimated:
    # Window = [from] epidemic start [to] epidemic start + window.fit
    window.fit <- 10
    
    pdf(paste0('R0-fit-',i,'.pdf'), width = 16, height = 10)
    
    # Calculate R0 implied from SIR:
    R0.SIR <- calc.R0.SIR(pop.all.mc = pop.all.mc, 
                          res.list.0 = res.list.0,
                          t.max.fit = window.fit, 
                          do.plot = F)
    
    msgt <- paste0(i,': R0_sir = ', round(R0.SIR,4))
    print(msgt)
    message(msgt)
    dev.off()
}

