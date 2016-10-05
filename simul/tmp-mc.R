i <-1
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

ggplot(ts.all) + geom_line(aes(x=time, y=nE, colour=factor(mc))) + scale_y_log10()
