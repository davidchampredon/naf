###
###  READ SCHEDULE FROM FILES 
###


read.schedules <- function(filename, return.names) {
	
	x <- read.csv(filename, header = TRUE, as.is = TRUE)
	y <- list()
	for(i in seq_along(x)){
		y[[i]] <- paste0('SP_',trimws(x[,i]))
	}
	names(y) <- names(x)
	if(return.names) res <- names(y) 
	else res <- y
	return(res)
}


