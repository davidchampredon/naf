library(ggplot2)
x <- read.csv('ontario-influenza-hosp-death.csv')


g <- ggplot(x) + geom_boxplot(aes(x=age_group, y=rate_100000)) + facet_wrap(~type, scales = 'free')
plot(g)
