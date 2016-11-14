dir.def <- function(fname) {
	x <- read.csv(fname, header = FALSE, as.is = TRUE)
	x[,2] <- trimws(x[,2])
	return(list(results = x[x[,1]=='dir.results',2],
				rdata   = x[x[,1]=='dir.rdata',2]))
}


