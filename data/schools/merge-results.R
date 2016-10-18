x <- system('ls df-schools*RData',intern = TRUE)


dflist <- list()

for(i in seq_along(x)){
	load(x[i])
	dflist[[i]] <- df.school
}

dfall <- do.call('rbind',dflist)
save(dfall,file = 'dfall.RData')
