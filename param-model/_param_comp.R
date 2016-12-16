library(ggplot2)
library(tidyr)
x <- read.csv('_all_param_value_source.csv')
n <- ncol(x)

dat <- gather(x,'param','v',2:n)
tmp <- gsub(x = dat$param, pattern = '_lo',replacement = '')
tmp <- gsub(x = tmp, pattern = '_hi',replacement = '')

dat$p <- tmp
dat$val <- as.numeric(dat$v)

dat.plot <- subset(dat, !is.na(val))
pdf('_param_comp.pdf', width = 15)
g <- ggplot(dat.plot) + geom_point(aes(x=paper,y=val, colour=paper)) 
g <- g + facet_wrap(~p, scales = 'free')
g <- g + xlab('') + ylab('')
g <- g + theme(axis.text.x = element_text(angle = 25, hjust = 1))
plot(g)


g <- ggplot(dat.plot) + geom_histogram(aes(x=val), bins=10)+ facet_wrap(~p,scales='free_x')
g <- g + xlab('parameter value') 
plot(g)

dev.off()

# dat <- gather(x,'param','v',4:n)
# dat$val.lo <- dat$v
# dat$val.hi <- dat$v
# dat$val.lo[grepl('_hi',dat$param)] <- NA
