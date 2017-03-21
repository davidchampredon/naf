### 
###   CALCULATE R0 OF SIMULATIONS ALREADY RUN
###

R0.target <- 1.8

# Retrieve contact rate values used:
cr <- as.numeric(unlist(read.csv(file = 'fit-cr-mean-vec.csv')))

# List simulation results 
z <- system('ls ./rdata/mc-simul-*.RData', intern = TRUE)

R0.SIR <- numeric(length(z))

for(i in seq_along(z)){
    load(z[i])
    source('analysis_tools.R')  
    
    # Merge all simulations:
    n.cpu <- parallel::detectCores() - 1
    pop.all.mc <- merge.pop.mc(res.list.0, n.cpu = 1, doparallel = FALSE)
    
    # This is the time window (in days)
    # when R0 is estimated:
    # Window = [from] epidemic start [to] epidemic start + window.fit
    window.fit <- 15
    
    pdf(paste0('R0-fit-',i,'.pdf'), width = 16, height = 10)
    
    # Calculate R0 implied from SIR:
    R0.SIR[i] <- calc.R0.SIR(pop.all.mc = pop.all.mc, 
                             res.list.0 = res.list.0,
                             t.max.fit = window.fit, 
                             do.plot = TRUE)
    
    msgt <- paste0(i,': R0_sir = ', round(R0.SIR[i], 4))
    print(msgt)
    message(msgt)
    dev.off()
}

idx.best <- which.min((R0.SIR-R0.target)^2)

plot(x=cr, y=R0.SIR, typ='o',cex=0.8,pch=16, lwd=2, las=1)
grid()
abline(h=R0.target, col='red')
abline(h=R0.SIR[idx.best], lty=2)
abline(v=cr[idx.best], lty=2)

prm.name <- '../param-model/prm-epi.csv'
prm <- read.csv(prm.name, header = F, strip.white = T, as.is = T)
prm[prm[,1]=='contact_rate_mean', 2] <- cr[idx.best]
    
    