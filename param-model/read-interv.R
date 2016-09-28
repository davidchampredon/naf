###
###    READ INTERVENTION PARAMETERS
###

fname <- 'prm-interv.csv' # <-- DELETE

load.interv.prm <- function(fname) {

	x <- read.csv(file = fname, header = TRUE)	
	
	nam <- unique(x$interv_name)
	n <- length(nam)
	
	res <- list()
	char.test <- c(letters,LETTERS)
	
	for(i in 1:n){
		y <- subset(x, interv_name==nam[i])
		tmp <- list()
		for(k in 1:nrow(y)) {
			num.prm <- ( sum(sapply(char.test, grepl, y$value[k]))==0 )
			if(num.prm)  tmp[[k]] <- as.numeric(as.character(y$value[k]))
			if(!num.prm) tmp[[k]] <- trimws(as.character(y$value[k]))
			names(tmp)[k] <- trimws(as.character(y$variable[k]))
		}
		tmp[[nrow(y)+1]]      <- as.character(y$interv_name[1])
		names(tmp)[nrow(y)+1] <- 'interv_name'
		res[[i]] <- tmp
	}
	return(res)
}
