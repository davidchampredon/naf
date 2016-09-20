####
####   READ MODEL PARAMETERS FROM FILES
####


read.prm <- function(filename, prm) {
	
	x <- read.csv(file = filename, header = FALSE)
	names(x)<- c('name','value')
	
	for(i in seq_along(x$name)){
		
		num.prm <- ( sum(sapply(letters, grepl, x$value[i]))==0 )
		if(num.prm)  y <- as.numeric(as.character(x$value[i]))
		if(!num.prm) y <- trimws(as.character(x$value[i]))
		names(y) <- trimws(as.character(x$name[i]))
		prm <- append(prm,y)
	}
	return(prm)
}


