### 
###   ANALYZE SIMULATIONS (FOR SINGLE SCENARIO RUN)
###

t00 <- as.numeric(Sys.time())

### ==== Load simulation results ====


SIM <- c(1,100)

for(i in SIM){
    
    print(paste('Loading simulation results',i))
    load(paste0('./rdata/mc-simul-',i,'.RData'))
    t01 <- as.numeric(Sys.time()) ; msgt <- round((t01-t00)/60,2)
    print(paste('... simulation results loaded in',msgt,'minutes.'))
    
    source('analysis_tools.R')  
    dir.results <- '../results/'
    
    n.cpu <- parallel::detectCores() - 1
    res.select <- res.list
    
    # Merge populations of all MC iterations:
    pop.all.mc   <- merge.pop.mc(res.select,
                                 n.cpu = n.cpu, 
                                 doparallel = FALSE)
    
    # Filter out the MC iterations that produced fizzles:
    pop.nofizz <- filter.out.fizzle(pop.all.mc)
    
    # Merging populations of all MC iterations
    # can be very slow, and they very similar
    # across MC, so just display the first non-fizzled:
    idx.mc.no.fizz <- unique(pop.nofizz$mc)
    pop <- subset(pop.nofizz, mc==idx.mc.no.fizz[1])
    
    
    # Population:
    pdf(paste0(dir.results,'plot-agevax-',SIM,'.pdf'))
    try( plot.vax.age(pop.nofizz), silent = F)
    dev.off()
    
}
t1 <- as.numeric(Sys.time())
print(paste('Analysis completed in',round((t1-t00)/60,1),'min'))

