###
###   READ PARAMETERS FOR WORLD CONSTRUCTION
###

read.prm.au <- function(filename) {
	# Read the high level parameters
	# for all Area Units.
	
	x <- read.csv(file = filename, header = TRUE)
	big.list <- list()
	prm      <- list()
	
	for(i in seq_along(x$name)){
		# W A R N I N G : 'id' should start at 0 (C++ constraint downstream)
		
		# if 'id' changes, then store previous parameter
		# list in the big list, reset, and carry on scanning:
		num.prm <- ( sum(sapply(letters, grepl, x$value[i]))==0 )
		if(num.prm)  y <- as.numeric(as.character(x$value[i]))
		if(!num.prm) y <- trimws(as.character(x$value[i]))
		names(y) <- trimws(as.character(x$name[i]))
		prm <- append(prm,y)
		
		if(i<length(x$name)){
			if(x$id[i] != x$id[i+1]) {
				big.list[[ x$id[i] +1 ]] <- prm
				prm <- list()
			}
		}
		if(i==length(x$name)){
			big.list[[ x$id[i]+1 ]] <- prm
		}
	}
	return(big.list)
}


read.size.distrib <- function(filename){
	# Read size distribution from a file
	# and return a data frame.
	
	x <- read.csv(filename)
	stopifnot(nrow(x)>0 | 
			  	names(x)==c('size','prop') |
			  	abs(sum(x$prop) - 1.0)<1E-8	)
	return(x)
}


get.distrib.type <- function(filename){
	# From thename of the file,
	# determine the type of size distribution
	
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



load.all.size.distrib <- function(filenames, prm.list, path.prm) {
	# For given file names, return a 
	# named list of size distributions.
	for(i in seq_along(filenames)){
		full.path.fname <- paste0(path.prm,filenames[i])
		x <- read.size.distrib(full.path.fname)
		d.type <- get.distrib.type(full.path.fname)
		n <- length(prm.list)
		prm.list[[n+1]] <- x$size
		names(prm.list)[n+1] <- d.type
		prm.list[[n+2]] <- x$prop
		names(prm.list)[n+2] <- paste0(d.type,'_proba')
	}
	return(prm.list)
}


get.size.distrib.file <- function(x){
	# return the file names where the
	# size distribution are saved.
	# 'x' is the high-level parameters list
	# for a given Area Unit.
	idx <- which(grepl('size_distrib_file',names(x)))
	return( unlist(x[idx]) )
}

merge.lists.into.one.list <- function(L){
	# Merge several lists into a unique list.
	# All lists must have the same structure and element names.
	# If there are k lists, the ith element of the 
	# returned unique list is itself a list made of 
	# the ith element of the k input lists.
	
	N <- length(L)
	n <- length(L[[1]])
	for(k in 1:N) stopifnot(n==length(L[[k]]))
	
	z <- list()
	for(i in 1:n){
		tmp <- list(L[[1]][[i]])
		for(j in 2:N) tmp <- c(tmp,list(L[[j]][[i]]))
		z[[i]] <- tmp
	}
	names(z) <- names(L[[1]])
	return(z)
}


merge.lists.into.one.list.vec <- function(L){
	# Merge several lists into a unique list.
	# All lists must have the same structure and element names.
	# If there are k lists, the ith element of the 
	# returned unique list is itself a list made of 
	# the ith element of the k input lists.
	N <- length(L)
	n <- length(L[[1]])
	for(k in 1:N) stopifnot(n==length(L[[k]]))
	
	z <- list()
	for(i in 1:n){
		tmp <- L[[1]][[i]]
		for(j in 2:N) tmp <- c(tmp,L[[j]][[i]])
		z[[i]] <- tmp
	}
	names(z) <- names(L[[1]])
	return(z)
}

load.world.prm <- function(filename, path.prm){
	# Read from files and create list
	# of parameter values to build 
	# simulation world.
	
	prm.au <- read.prm.au(filename)
	n.au <- length(prm.au)
	
	sz.dst.list <- list()
	
	# Get size distributions
	for(i in 1:n.au){
		filenames <- get.size.distrib.file(prm.au[[i]])
		tmp <- list()
		tmp <- load.all.size.distrib(filenames,tmp,path.prm)
		sz.dst.list[[i]] <- tmp
	}
	
	prm.sz.dst <- merge.lists.into.one.list(sz.dst.list)
	prm.au.2 <- merge.lists.into.one.list.vec(prm.au)
	world.prm <- c(prm.sz.dst, prm.au.2)
	return(world.prm)
}


scale.world <- function(scale.factor, world.prm) {
	# Scale (up or down) the size of the world.
	
	idx <- which(substr(names(world.prm),1,2)=='n_')
	for (i in idx) world.prm[[i]] <- ceiling(world.prm[[i]]*scale.factor)
	return(world.prm)
}
