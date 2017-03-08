library(plyr)
library(dplyr)
library(ggplot2)
load('../simul/mc-simul.RData')


# Retrieve time of infection for all simulations:
y <- list()
age.bucket <- 5
for(mc in 1:length(res.list.0)){
    x0 <- as.data.frame(res.list.0[[mc]]$world_final)
    fs <- sum(x0$is_recovered)/nrow(x0)
    # Filter out fizzles:
    if(fs>0.01){
        x <- subset(x0, t_infected<999)
        x$ageround <- round(x$age/age.bucket,0) * age.bucket
        y[[mc]] <- select(x, ageround, t_infected)
        y[[mc]]$mc <- mc
    }
}
df <- do.call('rbind.data.frame',y)
df$day_inf <- round(df$t_infected)

# Mean of infection times:
z <- ddply(df,'ageround',summarize, 
           m = mean(t_infected),
           qlo = quantile(t_infected, probs = 0.025),
           qhi = quantile(t_infected, probs = 0.975))

mm <- mean(df$t_infected)


# Calculate the midpoint of cumulative infections 
# See: Schanzer D, Vachon J, Pelletier L. Age-specific Differences in Influenza A Epidemic Curves: Do Children Drive the Spread of Influenza Epidemics? American Journal of Epidemiology 2011; 174: 109–17.
zz <- ddply(df,c('ageround','day_inf'),summarise, a=table(day_inf))
zz$b <- as.numeric(zz$a)
zz$csum <- ave(zz$b, zz$ageround, FUN=cumsum)
zz$csummax <- ave(zz$csum, zz$ageround, FUN=max)
zz$prop.inf <- zz$csum / zz$csummax

# Identify the time when 50% of infections occured:
zz$is.mid.val <- FALSE
zz$is.mid.val[abs(zz$prop.inf-0.5) < 0.05] <- TRUE

midinf0 <- select(zz[zz$is.mid.val,], ageround, day_inf)
midinf <- ddply(midinf0,'ageround',summarize, day_inf_m = mean(day_inf))

m.midinf <- mean(midinf$day_inf)

schanzer <- read.csv('../data/contact-rates/Schanzer-2011-Fig1-Canada.csv')

# plots:
pdf('fit-infectTime-age.pdf')
plot(z$ageround, z$m - mm, t='o', 
     main = 'Mean infection times',
     # ylim=range(z$qlo-mm,z$qhi-mm),
     xlab='age', ylab='Difference mean infection time (centred)')
abline(h=0, col='lightgrey', lty=3)
points(x=schanzer$agecat, y=schanzer$value, pch=15,col='red',lwd=2)
# lines(x = z$ageround, y=z$qlo - mm, lty=2)
# lines(x = z$ageround, y=z$qhi - mm, lty=2)


plot(midinf$ageround, midinf$day_inf_m-m.midinf, 
     main = 'Midpoint cumulative incidence',
     xlab = 'age', ylab = 'Centred Midpoints',
     typ='o')
points(x=schanzer$agecat, y=schanzer$value, pch=15,col='red',lwd=2)
dev.off()