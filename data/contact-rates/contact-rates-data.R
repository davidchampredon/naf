library(plyr)
library(rgexf)
library(igraph)


### Workplace ###

# http://www.sociopatterns.org/datasets/contacts-in-a-workplace/
x <- read.table('tij_InVS.dat')


x$time <- (x$V1-min(x$V1))/86400
x$hour <- round((x$V1-min(x$V1))/60/60)
x$day <- round((x$V1-min(x$V1))/60/60/24)

x$key <- paste(x$V2,x$hour,sep='_')
summary(x$time)


z <- subset(x,subset = V2==492)


y <- ddply(x,c('V2','day', 'key'),summarise, n=length(key))

cr <- ddply(y,c('V2','day'), summarize, cd = length(key))

hist(cr$cd, breaks=0:15,
	 col = 'lightgrey',
	 main = 'Distribution of contact rate',
	 xlab = 'contacts / person / day')

summary(cr$cd)
sd(cr$cd)



### Primary School ###
summary.contact <- function(file=fname, do.plot = TRUE) {
	
	a <- igraph::read_graph(file = fname,format = 'gml')
	dd.x <- igraph::degree(a)
	deg.max <- max(dd.x)
	deg.vec <- 0:deg.max
	dd.y <- igraph::degree_distribution(graph = a)
	if(do.plot) plot(x=deg.vec, y=dd.y, typ='h',lwd=6)
	mean.deg <- sum( dd.y * deg.vec)
	stddev.deg <- sqrt( sum( dd.y * deg.vec^2) - mean.deg^2 )
	return(list(mean.deg=mean.deg, stddev.deg=stddev.deg))
}

# fname <- 'journal.pone.0023176.s003.GML'
par(mfrow=c(1,1))
summary.contact(file= 'journal.pone.0023176.s003.GML')
# SAME???! summary.contact(file= 'journal.pone.0023176.s004.GML')


