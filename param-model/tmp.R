library(plyr)

load('mc-simul.RData')

mc <- 3

calc.inf.reduc <- function(mc) {
    
    extract.pop <- function(x, mc) {
        pop <- as.data.frame(x[[mc]]$world_final)
        pop$ageGroup <- round(pop$age/2)*2
        return(pop)    
    }
    
    pop.0 <- extract.pop(res.list.0, mc)
    pop   <- extract.pop(res.list, mc)
    
    calc.prop <- function(pop){
        z <- ddply(pop,'ageGroup', summarise, 
                   n = length(id_indiv),
                   pvax = sum(is_vaccinated)/length(is_vaccinated),
                   pinf = (sum(is_recovered)+sum(1-is_alive))/(length(is_recovered)+length(is_alive)))  # was_symptomatic
        return(z)
    }
    
    z.0 <- calc.prop(pop.0)
    z   <- calc.prop(pop)
    
    names(z.0) <- paste0(names(z.0),'.0')
    
    a <- cbind(z.0, z)
    a$diff.inf <- a$pinf - a$pinf.0
    plot(a$ageGroup, a$diff.inf, typ='h')
    return(a)  
}

par(mfrow=c(2,2))

tmp <- list()
for(i in seq_along(res.list)) { 
    tmp[[i]] <- calc.inf.reduc(i)
    tmp[[i]]$mc <- i
}
df <- do.call('rbind.data.frame',tmp)

res <- ddply(df,'ageGroup',summarize, diff = mean(diff.inf))

plot(res$ageGroup, res$diff, typ='h',lwd=5)
