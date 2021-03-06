### 
###   CREATE FIGURES FOR SUPPLEMENTARY INFORMATIONS
###


library(plyr)
library(ggplot2)

source('analysis_tools.R')  
source('utils-compare.R')
source('utils-misc.R')

dir.SI.fig <- '../results/figures/SI/'
file.scen.prm.list <- 'scenario-prm-list.csv'

# ---- Load simulation results  -----

load.simul.results <- function(rdataFile) {

    load(rdataFile)
    
    # Extract baseline simulations:
    res.select <- res.list.0
    
    # Merge populations of all MC iterations:
    pop.all.mc   <- merge.pop.mc(res.select,
                                 n.cpu = 1, 
                                 doparallel = FALSE)
    
    # Filter out the MC iterations that produced fizzles:
    pop.nofizz <- filter.out.fizzle(pop.all.mc)
    return(list(res.list.0 = res.list.0,
                pop.nofizz = pop.nofizz, 
                world.prm = world.prm))
}

# Single scenario for output on population:
try(expr = {
    rdataFile <- 'mc-simul.RData'
    RES <- load.simul.results(rdataFile)
    POP        <- RES[['pop.nofizz']]
    world.prm  <- RES[['world.prm']]
    res.list.0 <- RES[['res.list.0']]}
)

# Multi scenario
dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]
load(paste0(dir.save.rdata,'result-scen-all.RData'))
res.all <- list(main = result.scen.all,
                main.ageGroup = result.scen.all.ageGroup)

# ---- Internal functions ----

plot.sp.size.distrib <- function(res.list.0, world.prm, POP){
    
    par(mfrow=c(2,2))
    
    idx.mc.no.fizz <- unique(POP$mc)
    # Read target data:
    df <- read_input_sp_size_distribution(world.prm)
    # Retrieve all SocialPlaces types:
    sptype.vec <- seq_along(res.list.0[[1]][['sp_size_distrib']])
    
    
    for(k in seq_along(sptype.vec)) {
        
        sptype <- sp_type_string(k-1)
        data.df <- subset(df, sp_type_string==sptype)
        
        if(sptype != 'Hospital' & sptype != 'Other'){
            
            # Retrieve social place sizes of the kth type, 
            # for the mth MC iteration:
            p <- list()
            cnt <- 1
            for( m in idx.mc.no.fizz){
                sp.sz.d <- res.list.0[[m]][['sp_size_distrib']]
                h <- hist(x = sp.sz.d[[k]], plot = F, breaks = c(0,data.df$size))
                p[[cnt]] <- h$counts / sum(h$counts)
                cnt <- cnt + 1
            }
            # Summary stats across all MC iterations:
            pm <- matrix(unlist(p), ncol=cnt-1)
            p.avg <- apply(X = pm,MARGIN = 1,FUN = mean)
            p.max <- apply(X = pm,MARGIN = 1,FUN = max)
            p.min <- apply(X = pm,MARGIN = 1,FUN = min)
            
            # Plots:
            title <- sptype
            ylim  <- range(data.df$freq, p.avg)
            if(sptype == 'PubTransp'){
                title <- 'Public Transportation'
                ylim <- c(0,1)
            }
            plot(x = data.df$size, y = data.df$freq, 
                 ylim = ylim,
                 col='red',lwd=1, cex=2, 
                 typ='p', lty=2, las=1, 
                 main = title, xlab = 'size', ylab='proportion')
            grid()
            lines(x = data.df$size, y = p.avg, typ='p', 
                  pch=3, cex=1.7)
            legend(x='topright', legend = c('Data','Model'), 
                   col=c('red','black'), pch=c(1,3), lty=c(2,1))
        }
    }
}

plot.hh.size.distrib <- function(res.list.0, world.prm, POP){
    
    idx.mc.no.fizz <- unique(POP$mc)
    # Read target data:
    df <- read_input_sp_size_distribution(world.prm)
    # Retrieve all SocialPlaces types:
    sptype.vec <- seq_along(res.list.0[[1]][['sp_size_distrib']])
    
    # for(k in seq_along(sptype.vec)) {
        k <- 1
        sptype <- sp_type_string(k-1)
        data.df <- subset(df, sp_type_string==sptype)
        
        if(sptype == 'Household'){
            
            # Retrieve social place sizes of the kth type, 
            # for the mth MC iteration:
            p <- list()
            cnt <- 1
            for( m in idx.mc.no.fizz){
                sp.sz.d <- res.list.0[[m]][['sp_size_distrib']]
                h <- hist(x = sp.sz.d[[k]], plot = F, breaks = c(0,data.df$size))
                p[[cnt]] <- h$counts / sum(h$counts)
                cnt <- cnt + 1
            }
            # Summary stats across all MC iterations:
            pm <- matrix(unlist(p), ncol=cnt-1)
            p.avg <- apply(X = pm,MARGIN = 1,FUN = mean)
            p.max <- apply(X = pm,MARGIN = 1,FUN = max)
            p.min <- apply(X = pm,MARGIN = 1,FUN = min)
            
            # Plots:
            ylim  <- range(data.df$freq, p.avg)
            plot(x = data.df$size, y = data.df$freq, 
                 ylim = ylim,
                 col='red',lwd=2, cex=1.5, 
                 typ='p', lty=2, las=1, 
                 main = 'Households size distribution in Ontario', 
                 xlab = 'Household size', ylab='proportion')
            grid()
            points(x = data.df$size, y = p.avg, 
                  pch=3, cex=1.7, lwd =2)
            legend(x='topright', legend = c('Data','Model'), 
                   col=c('red','black'), pch=c(1,3), pt.lwd = 2, pt.cex=1.5)
        # }
    }
}

plot.age.distribution <- function(POP){
    # Age distributionfrom data(Stats Canada)
    dat <- read.csv('../data/ages/size-distrib-ages.csv')
    
    # Age distribution from simulations
    POP$ageround <- round(POP$age,0)
    z <- ddply(POP,'ageround', summarize, n = length(id_indiv))
    z$prop <- z$n / sum(z$n)
    
    # Cosmetics:
    
    cex <- 1.5
    dat.pch <- 1
    dat.lwd <- 2
    dat.col <- 'red'
    
    mod.pch <- 3
    mod.lwd <- 2
    mod.col <- 'black'
    
    # Plot:
    plot(dat$age, dat$prop, 
         cex = cex, 
         lwd = dat.lwd,
         pch = dat.pch,
         col = dat.col,
         las = 1,
         main = 'Ontario Age Distribution Fit',
         xlab = 'Age', ylab = 'Proportion')
    
    points(z$ageround, z$prop, 
           pch = mod.pch, 
           lwd = mod.lwd, 
           col = mod.col, 
           cex = cex)
    
    legend(x='topright', 
           legend = c('Data', 'Model'),
           pch = c(dat.pch, mod.pch),
           col = c(dat.col, mod.col),
           pt.cex = cex,
           pt.lwd = c(dat.lwd, mod.lwd))
    
}



# ---- SI Figures -----
figure.S1.calibration <- function(res.list.0, world.prm, POP){
    w <- 1000
    png(paste0(dir.SI.fig,'fig-S1.png'), width = w, height = w/2)
    par(mfrow = c(1,2))
    plot.age.distribution(POP)
    plot.hh.size.distrib(res.list.0, world.prm, POP)   
    dev.off()
}




figure.S2 <- function(zlist.full) {
    
    DAT <- do.call('rbind.data.frame', zlist.full)
    DAT <- reformat(DAT)
    
    # delete when sure: DAT <- subset(DAT, interv_target=='priority_age5_frailty')
    
    # Select a unique vaccination strategy for this plot:
    uvs <- as.character(unique(DAT$interv_target))
    sel.strategy.name <- 'priority_age'
    sel.vs <- uvs[grepl(sel.strategy.name, uvs)]
    
    if(length(sel.vs) != 1){
        warning('None, or more than one, vax strategy matched the strategy name... ABORTING!')
        warning(sel.strategy.name)
        warning(sel.vs)
        stop()
    }
    DAT <- subset(DAT, interv_target == sel.vs)
    
    DAT <- floor.results.0(DAT)
    ar <- unique(DAT$interv_start)
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, 
                         color=factor(interv_cvg_rate))) +
        geom_line(alpha=0.5, size=2) +
        geom_point(alpha=0.8, size=3) +
        scale_x_continuous(breaks=ar) +
        facet_grid( outcome ~ R0 + VE ) +
        guides(color=guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)")) +
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette2) +
        coord_cartesian(ylim=c(0,1)) +
        theme(legend.position=c(0.94,0.91)) +
        xlab("Vaccination time lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0(dir.SI.fig,'fig-S2.pdf'),
        width = 16, height = 10)    
    plot(g)
    dev.off()
    
    # Save data for this figure to a csv file:
    dat.sav <- ddply(DAT,c('R0','VE','outcome','interv_start','interv_cvg_rate'),
                     summarize, percentMean = 100 * mn)
    write.csv(x = dat.sav, file=paste0(dir.SI.fig,'fig-S2.csv'))
}

figure.S3 <- function(zlist.full) {
    
    DAT <- do.call('rbind.data.frame', zlist.full)
    DAT <- reformat(DAT)
    DAT <- floor.results.0(DAT)
    
    ar <- unique(DAT$interv_start)
    
    # Explicit name for plot:
    DAT$it <- NA
    DAT$it[DAT$interv_target=='never_sympt'] <- 'RVS'
    DAT$it[DAT$interv_target=='priority_age_frailty'] <- 'Priority'
    DAT$it[DAT$interv_target=='priority_age5_frailty'] <- 'PVS1'
    DAT$it[DAT$interv_target=='priority_age5_10_frailty'] <- 'Priority5-10'
    DAT$it[DAT$interv_target=='priority_age19_frailty'] <- 'PVS2'
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, color=factor(interv_cvg_rate))) +
        geom_line(alpha=0.5, size=2, aes(linetype=factor(it))) +
        geom_point(alpha=0.8, size=3) +
        scale_x_continuous(breaks=ar) +
        facet_grid( outcome ~ R0 + VE ) +
        guides(color=guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)")) +
        guides(linetype=guide_legend(title="Vaccination Strategy")) +
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette2) +
        coord_cartesian(ylim=c(0,1)) +
        xlab("Vaccination time lag (days)") + ylab("Mean relative reduction") +
        theme(legend.position=c(0.94,0.85)) +
        theme(legend.key.width=unit(2,"cm"))
    
    pdf(file = paste0(dir.SI.fig,'fig-S3.pdf'), width = 16, height = 10)    
    plot(g)
    dev.off()
}


figure.S4 <- function(zlist.ag){
    
    DAT <- do.call('rbind.data.frame', zlist.ag)
    DAT <- reformat(DAT)
    DAT <- floor.results.0(DAT)
    
    DAT$AG <- DAT$ageGroup
    DAT$AG[DAT$AG == '0_5'] <- '0 - 4'
    DAT$AG[DAT$AG == '5_18'] <- '5 - 18'
    DAT$AG[DAT$AG == '18_65'] <- '19 - 64'
    DAT$AG[DAT$AG == '65_over'] <- '65+'
    
    DAT$AG <- factor(DAT$AG, levels = c("0 - 4", "5 - 18", "19 - 64", "65+"))
    DAT$speed <- paste0('Vacc. admin. rate: ', DAT$interv_cvg_rate)
    lvl <- paste0('Vacc. admin. rate: ', sort(unique(DAT$interv_cvg_rate)))
    DAT$speed <- factor(DAT$speed,levels = lvl)
    
    ut    <- as.character(unique(DAT$interv_target))
    DAT.1 <- subset(DAT, interv_target==ut[1])
    DAT.2 <- subset(DAT, interv_target==ut[2])
    
    
    internal.plot <- function(DAT, title, label) {
        g <- ggplot(DAT, aes(x=AG, y=mn, 
                             fill=factor(interv_start),
                             color=factor(interv_start))) +
            geom_line(aes(group=factor(interv_start))) +
            geom_point(alpha=0.75, size=3, shape=22) +
            facet_grid( VE + R0 ~  speed) +
            guides(fill=guide_legend(title="Vacc. time lag (days)"), color=FALSE) +
            scale_fill_manual(values=mypalette) + 
            scale_color_manual(values=mypalette) + 
            coord_cartesian(ylim=c(0,1)) +
            ggtitle(label = paste(label,'vaccination strategy',title)) + 
            xlab("Age group") + ylab("Mean relative reduction") +
            theme(legend.position=c(0.1,0.9)) 
        return(g)
    }
    
    internal.translate <- function(interv_target) {
        res <- NA
        if(interv_target=='never_sympt') res <- 'RVS'
        if(interv_target=='priority_age_frailty') res  <- 'Priority'
        if(interv_target=='priority_age5_frailty') res  <- 'PVS1'
        if(interv_target=='priority_age5_10_frailty') res  <- 'Priority5-10'
        if(interv_target=='priority_age19_frailty') res  <- 'PVS2'
        return(res)
    }
    W <- 15
    pdf(file =  paste0(dir.SI.fig,'fig-S4.pdf'), width = W,height = 1.2*W)
    p1 <- internal.plot(DAT.1, title = internal.translate(ut[1]), label='B - ')
    p2 <- internal.plot(DAT.2, title = internal.translate(ut[2]), label='A - ')
    grid.arrange(p2,p1, nrow=2)
    dev.off()
}

#' Clinical (symptomatic infections) final size by age group
figure.S00 <- function(df) {

    df$CFS <- df$tot.sympt.baseline/df$popsize
    a <- ddply(df, c('ageGroup'), summarize, 
               m=mean(CFS),
               hi = max(CFS),
               lo = min(CFS))
    a$ageGroup <- factor(a$ageGroup, levels = c("0_5", "5_18", "18_65","65_over"))
    g <- ggplot(a)+ 
        geom_bar(aes(x=ageGroup, y=m), stat='identity')+
        geom_errorbar(aes(x=ageGroup, ymin=lo, ymax=hi),width=0.3)+
        ggtitle('Clinical attack rate (mean and range)') + ylab('')
    
    pdf(file =  paste0(dir.SI.fig,'fig-S00.pdf'), width = 12,height = 8)
    plot(g)
    dev.off()
}



# ---- Save all SI figures ----

df  <- res.all[['main.ageGroup']]
X <- process.outputs(df, dir.results, file.scen.prm.list)

zlist      <- X[['zlist']]
zlist.ag   <- X[['zlist.ag']]
zlist.full <- X[['zlist.full']]

try(figure.S1.calibration(res.list.0, world.prm, POP))
try(figure.S2(zlist.full) )
try(figure.S3(zlist.full) )
try(figure.S4(zlist.ag))
#figure.S00(df)

