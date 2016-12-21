###
###   * * * WARNING * * * THIS DOES NOT FIT! IT SIMPLY DISPLAYS DATA 
###   THE VALUES OF THE PARAMETERS 'proba_hosp' AND 'proba_death_xxx'
###   ARE SET MANUALLY, SUCH THAT THE RESULTING PROBA FROM THE SIMULATION
###   MATCH THE DATA.
###


library(ggplot2)
x <- read.csv('../data/frailty/ontario/ontario-influenza-hosp-death.csv')

# reorder propoerly:
unique(x$age_group)
as.character(x$age_group)
x$age_group2 <- factor(x$age_group, levels=c('<1 ','1_4 ','5 _ 14 ','15 _ 24 ','25 _ 44 ','45 _ 64 ','65+ '))

pdf('proba-hosp-death-ontario.pdf', width = 10)
g <- ggplot(x) + geom_boxplot(aes(x=age_group2, y=rate_100000)) + facet_wrap(~type, scales = 'free')
g <- g + xlab('Age Group') + ylab('per 100,000 individuals per year') 
g <- g + ggtitle("Death and hospitalization rate for Ontario 2013-2016\n (source: Public Health Ontario)")
plot(g)
dev.off()