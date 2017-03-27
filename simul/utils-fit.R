#' Index that is unique to a scenario and 
#' to an iteration of the parameter explored.
sensana.id <- function(scenidx, i){
    return(scenidx*100000 + i)
}


#' Explore a range of values for the mean contact rate. 
#' The results will be used to fit to R0
#' @param cr.mean.vec Vector of mean contact rate values to explore.
#' @param scenidx ID number of the current scenario.
#' @return No object output. Save RData files (to be used by the fit).
explore.cr <- function(cr.mean.vec, scenidx) {
    ### ==== Libraries and sources ====
    R.library.dir    <- '../Rlibrary/lib'
    param.model.dir  <- '../param-model/'
    data.dir         <- '../data/'
    dirdef <- read.csv('dir-def.csv',header = F,strip.white = T,as.is = T)
    rdata.dir <- dirdef[dirdef[,1]=='dir.rdata',2]
    
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
    
    # Force intervention start date (not used) to be 0
    # in order to get the simulation started at time t=0:
    for(k in seq_along(interv.prm)) 
        interv.prm[[k]][['interv_start']] <- 0
    
    ### ==== Run Simulation ====
    
    for(i in seq_along(cr.mean.vec)){
        print(paste('Simulate with CRmean =', cr.mean.vec[i]))
        
        # overwrite value:
        prm[['contact_rate_mean']] <- cr.mean.vec[i]
        
        run.simul.fit.R(scen.id = sensana.id(scenidx, i), 
                        dir.save.rdata = rdata.dir, 
                        prm, simul.prm, 
                        interv.prm.0, world.prm, sched.prm, stoch_build_world)
    }
}

#' Fit the mean contact rate to a target R0.
fit.cr.R0 <- function(R0.target, cr.mean.vec, scenidx){
    
    # Explore parameter space:
    explore.cr(cr.mean.vec, scenidx)
    # Retrieve contact rate values used:
    cr <- cr.mean.vec
    # List simulation results 
    # for the associated scenario:
    dirdef    <- read.csv('dir-def.csv',header = F,strip.white = T,as.is = T)
    rdata.dir <- dirdef[dirdef[,1]=='dir.rdata',2]
    presimulated.rdata <- paste0('ls ',rdata.dir,'mc-simul-',
                                 sensana.id(scenidx,0)/1000,'*.RData')
    z <- system(presimulated.rdata, 
                intern = TRUE)
    
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
        
        # Calculate R0 implied from SIR:
        do.plot <- FALSE
        if(do.plot) pdf(paste0('R0-fit-',sensana.id(scenidx,i),'.pdf'), width = 16, height = 10)
        R0.SIR[i] <- calc.R0.SIR(pop.all.mc = pop.all.mc, 
                                 res.list.0 = res.list.0,
                                 t.max.fit = window.fit, 
                                 do.plot = do.plot)
        if(do.plot) dev.off()
        msgt <- paste0(i,': R0_sir = ', round(R0.SIR[i], 4))
        print(msgt)
        message(msgt)
    }
    
    idx.best <- which.min((R0.SIR - R0.target)^2)
    
    try({
        pdf(paste0('R0-fit-result-',scenidx,'.pdf'), width = 8, height = 8)
        plot(x=cr, y=R0.SIR, typ='o', cex=0.8, pch=16, lwd=2, las=1)
        grid()
        abline(h=R0.target, col='red')
        abline(h=R0.SIR[idx.best], lty=2)
        abline(v=cr[idx.best], lty=2)
        dev.off()
    },silent = FALSE)
    
    # Clean-up temporary RData files:
    delete.rdata <- paste0('rm ',rdata.dir,'mc-simul-',
                                 sensana.id(scenidx,0)/1000,'*.RData')
    system(command = delete.rdata, intern = F)
    
    # Save the fitted contact rate value in log file:
    a <- c(scenidx = scenidx, 
           fitted.cr = cr[idx.best], 
           R0.target = R0.target, 
           time = as.character(Sys.time()))
    write.table(x = t(a), file = 'fitted-cr.csv', 
                append = TRUE, sep=',', 
                row.names = F, col.names = F)
    
    return(cr[idx.best])
}


