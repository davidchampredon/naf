library(ggplot2)

df <- list()

df[[1]] <- read.csv('prm-size-hh-ontario.csv')
df[[1]]$type <- 'households'

df[[2]] <- read.csv('prm-size-other-ontario.csv')
df[[2]]$type <- 'other social places'

df[[3]] <- read.csv('prm-size-pubt-ontario.csv')
df[[3]]$type <- 'public transportation'

df[[4]] <- read.csv('prm-size-school-ontario.csv')
df[[4]]$type <- 'schools'

df[[5]] <- read.csv('prm-size-wrk-ontario.csv')
df[[5]]$type <- 'workplace'

dfall <- do.call('rbind', df)

png('plot-sp-distributions.png', width = 1200, height = 1000)
ggplot(dfall,aes(x=size,y=prop))+
    geom_point() +
    geom_line() + 
    facet_wrap(~type, scales='free')+
    xlab('size of the social place')+
    ylab('proportion')+
    theme(text = element_text(size=20))
dev.off()
