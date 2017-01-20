library(ggplot2)
library(plyr)

# Data from Morrison KT, Buckeridge DL, Xiao Y, Moghadas SM. Health & Place. Health & Place 2014; 26: 53–9.
dat <- read.csv('time-to-hosp-Morrison-2014.txt')

is.data.frame(dat)
ggplot(dat)+geom_bar(aes(x=time_to_hosp, y=prop, fill=zone),
                     stat = 'identity',
                     position = 'dodge')


df <- ddply(dat,c('zone'), summarize, 
            m = sum(time_to_hosp * prop),
            v = sum((time_to_hosp^2 * prop)) - (sum(time_to_hosp * prop))^2)

df
mm <- 4#mean(df$m)
vv <- 0.7 #mean(df$v)

mu = log(mm^2/sqrt(vv+mm^2));
sigma = log(1.0 + vv/mm^2);


xx <- seq(0,10,by=0.01)
yy <- dlnorm(xx,x = mu, sdlog = sigma)

plot(xx,yy,typ='l')
    
