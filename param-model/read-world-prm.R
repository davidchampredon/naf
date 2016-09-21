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


read.size.distrib <- function(filename){
	x <- read.csv(filename)
	stopifnot(nrow(x)>0 | 
			  	names(x)==c('size','prop') |
			  	abs(sum(x$prop) - 1.0)<1E-8	)
	return(x)
}


get.distrib.type <- function(filename){
	d.type <- NA
	if(grepl('hh',filename))          d.type <- 'hh_size'
	else if(grepl('wrk',filename))    d.type <- 'wrk_size'
	else if(grepl('school',filename)) d.type <- 'school_size'
	else if(grepl('pubt',filename))   d.type <- 'pubt_size'
	else if(grepl('hosp',filename))   d.type <- 'hosp_size'
	else if(grepl('other',filename))  d.type <- 'other_size'
	if(is.na(d.type)){
		print(paste('Cannot determine size distribution type from file name',filename))
		stop()
	}
	return(d.type)
}



load.all.size.distrib <- function(filenames, prm.world) {
	for(i in seq_along(filenames)){
		x <- read.size.distrib(filenames[i])
		d.type <- get.distrib.type(filenames[i])
		n <- length(prm.world)
		prm.world[[n+1]] <- x$size
		names(prm.world)[n+1] <- d.type
		prm.world[[n+2]] <- x$prop
		names(prm.world)[n+2] <- paste0(d.type,'_proba')
	}
	return(prm.world)
}





filename <- 'prm-au.csv'
load.all.world.prm <- function(filename) {
	x <- read.prm.au(filename)
	
}



if(TRUE){
	prm.w <- read.prm.au(filename = 'prm-au.csv')
	filename <- 'prm-size-hh-ontario.csv'
	aa <- read.size.distrib(filename)

	filenames <- c('prm-size-hh-ontario.csv', 'prm-size-school-ontario.csv')
	prm.world <- list()
	prm.world <- load.all.size.distrib(filenames,prm.world)
	prm.world
}
