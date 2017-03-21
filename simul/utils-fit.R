

explore.cr <- function(cr.mean.vec) {
    ### ==== Libraries and sources ====
    R.library.dir    <- '../Rlibrary/lib'
    param.model.dir  <- '../param-model/'
    data.dir         <- '../data/'
    
    library(ggplot2)
    library(plyr)
    library(naf,lib.loc = R.library.dir)
    library(snowfall)
    
    source('analysis_tools.R')
    source('build-world-det-sizegen.R')
    source('utils-run.R')
    source(paste0(param.model.dir,'read-prm.R'))
    source(paste0(param.model.dir,'read-world-prm.R'))
    source(paste0(param.model.dir,'read-interv.R'))
    source(paste0(param.model.dir,'read-sched.R'))
    
    ### ==== Parameter file names ====
    fname.prm.epi      <- 'prm-epi.csv'
    fname.prm.simul    <- 'prm-simul.csv'
    fname.prm.au       <- 'prm-au-ontario.csv'
    fname.prm.interv.0 <- 'prm-interv-0.csv'
    fname.prm.interv   <- 'prm-interv.csv'
    fname.schedules    <- 'schedules.csv'
    
    save.plot.to.file <- TRUE
    
    prm        <- list()
    simul.prm  <- list()
    interv.prm <- list()
    world.prm  <- list()
    sched.prm  <- list()
    
    
    ### ==== Model Parameters ====
    
    prm[['debug_mode']] <- FALSE
    
    prm       <- read.prm(filename = paste0(param.model.dir,fname.prm.epi), 
                          prm = prm)
    simul.prm <- read.prm(filename = paste0(param.model.dir,fname.prm.simul), 
                          prm = simul.prm)
    
    stoch_build_world <- simul.prm[['build_world_stoch']]
    
    ### ==== World parameters ====
    
    world.prm <- load.world.prm(filename = paste0(param.model.dir,fname.prm.au), 
                                path.prm = param.model.dir)
    
    world.prm[['id_region']]  <- 0
    world.prm[['regionName']] <- "Canada"
    
    sf           <- as.numeric(simul.prm[['scale_factor']])
    world.prm    <- scale.world(1/sf, world.prm)
    print(paste('World size reduced by',sf))
    
    # age distributions, conditional on 
    # household composition:
    base.fname  <- paste0(data.dir,'households/hh_size_ad')
    world.prm   <- c(world.prm,
                     load.age.hh(base.fname, 
                                 max.hh.size = max(world.prm[['max_hh_size']])) )
    
    world.prm <- gen_world_ontario(world.prm)
    
    # schedule time slices:
    
    sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
                                  1.0/24, 2.0/24, 12.0/24)
    
    sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
    sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)
    
    ### ==== Intervention parameters ====
    
    # Baseline, no vax intervention:
    interv.prm.0 <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))
    # Vax intervention:
    interv.prm   <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))
    
    
    ### === Snowfall wrapper ===
    
    run.snow.wrap <- function(seedMC,
                              prm, 
                              simul.prm, 
                              interv.prm, 
                              world.prm, 
                              sched.prm,
                              stoch_build_world){
        
        if(stoch_build_world){
            res <- naf_run(prm, 
                           simul.prm, 
                           interv.prm, 
                           world.prm, 
                           sched.prm,
                           seedMC)
        }
        else{
            res <- naf_run_det(prm, 
                               simul.prm, 
                               interv.prm, 
                               world.prm, 
                               sched.prm,
                               seedMC)
        }
        return(res)
    }
    
    
    
    ### ==== Define parameter space ====
    
    # cr.mean.vec <- seq(2.1,6.1,by=2)
    write.csv(x = cr.mean.vec, file = 'fit-cr-mean-vec.csv', quote = F, row.names = F)
    
    # Force intervention start date (not used) to be 0
    # in order to get the simulation started at time t=0:
    for(k in seq_along(interv.prm)) 
        interv.prm[[k]][['interv_start']] <- 0
    
    ### ==== Run Simulation ====
    
    for(i in seq_along(cr.mean.vec)){
        print(paste('Simulate with CRmean =', cr.mean.vec[i]))
        # overwrite value:
        prm[['contact_rate_mean']] <- cr.mean.vec[i]
        
        run.simul(scen.id = i, 
                  dir.save.rdata = './rdata/',
                  baseonly = TRUE) 
    }
    
}

fit.cr.R0 <- function(R0.target){
    # R0.target <- 1.8
    
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
    
    
}