i <-1

library(plyr)

pop.list <- list()
for(i in seq_along(res.list)){
	print(i)
	res <- res.list[[i]]
	world0 <- res[['world']]
	z      <- lapply(world0, as.data.frame)
	pop    <- do.call('rbind.data.frame',z)
	pop$mc <- i
	pop.list[[i]] <- pop
}

pop.all <- do.call('rbind.data.frame',pop.list)


ts.list <- list()

for(i in seq_along(res.list)){
	print(i)
	res <- res.list[[i]]
	ts     <- as.data.frame(res[['time_series']])
	ts$mc <- i
	ts.list[[i]] <- ts
}
ts.all <-  do.call('rbind.data.frame',ts.list)

ggplot(ts.all) + geom_line(aes(x=time, y=nIa, colour=factor(mc))) + scale_y_log10()

ts.all$timeround <- round(ts.all$time)

df <- ddply(ts.all, c('timeround'),summarise,
			md=median(nIa),
			q.lo=quantile(nIa,probs = 0.025),
			q.hi=quantile(nIa,probs = 0.975))
		
g <- ggplot(df,aes(x=timeround))+geom_line(aes(y=md),size=1)
g <- g + geom_line(aes(y=q.lo)) + geom_line(aes(y=q.hi))
g <- g + scale_y_log10()
plot(g)
