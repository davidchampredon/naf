###
###   READ PARAMETERS FOR WORLD CONSTRUCTION
###

read.prm.au <- function(filename) {
	
	x <- read.csv(file = filename, header = TRUE)
	
	big.list <- list()
	prm      <- list()
	
	for(i in seq_along(x$name)){
		# if 'id' changes, then store previous parameter
		# list in the big list, reset, and carry on scanning:
		
		num.prm <- ( sum(sapply(letters, grepl, x$value[i]))==0 )
		if(num.prm)  y <- as.numeric(as.character(x$value[i]))
		if(!num.prm) y <- trimws(as.character(x$value[i]))
		names(y) <- trimws(as.character(x$name[i]))
		prm <- append(prm,y)
		
		if(i<length(x$name)){
			if(x$id[i] != x$id[i+1]) {
				big.list[[ x$id[i] ]] <- prm
				prm <- list()
			}
		}
		if(i==length(x$name)){
			big.list[[x$id[i]]] <- prm
		}
	}
	return(big.list)
}

prm.w <- read.prm.au(filename = 'prm-au.csv')
