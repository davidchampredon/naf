library(RCurl)
library(XML)
library(stringr)
library(profvis)
# library(snowfall)


t0 <- as.numeric(Sys.time())
# profvis({

a <- commandArgs(trailingOnly = TRUE)

s.start <- as.numeric(a[1])
s.end   <- as.numeric(a[2])


# sfInit(parallel = TRUE, cpu = 4)
# sfLibrary(stringr) 
# sfLibrary(RCurl) 



read.value.from.url <- function(url) {
	html <- getURL(url, followlocation = TRUE)
	
	found <- is.na(str_locate(pattern='No schools were found',html)[,1])
	
	if(!found) return(NA)
	
	pos.a <- str_locate_all(pattern ='Address', html)[[1]][,2]
	pos.e <- str_locate_all(pattern ='Enrolment', html)[[1]][,2]
	
	val.a <- substr(html,pos.a+11,  pos.a+100)
	val.e <- substr(html,pos.e+11, pos.e+15)
	
	val.a <- trimws(substr(val.a,1,str_locate(val.a,"<")-1))
	val.e <- as.numeric(substr(val.e,1,str_locate(val.e,"<")-1))
	return(list(val.a,val.e))
}

read.one.url <- function(i,base.url) {
	ni <- nchar(as.character(i))
	school.num <- paste0(paste(rep('0',6-ni),collapse = ""),as.character(i))
	if(i%%10==0)  print(paste(school.num,'/',s.end))
	url <- paste0(base.url,school.num)
	return(read.value.from.url(url))
}


base.url <- 'https://www.app.edu.gov.on.ca/eng/sift/schoolProfile.asp?SCH_NUMBER='

# s.start <- 779
# s.end   <- 879

ns <- s.end-s.start+1
print(paste('Parsing',ns,'URLs ...', '(',s.start,s.end,')'))
navec <- rep(NA,ns)
df <- data.frame(id=navec, address=navec, enrolment=navec)


# idx.apply <- c(s.end:s.start)

### Parallel execution:
# sfExportAll()
# res<- sfSapply(idx.apply, read.one.url, base.url=base.url, simplify = FALSE)
# sfStop()

# print(res)


if(TRUE){
	cnt <- 1
	for (i in s.start:s.end){
		# ni <- nchar(as.character(i))
		# school.num <- paste0(paste(rep('0',6-ni),collapse = ""),as.character(i))
		# if(i%%10==0)  print(paste(school.num,'/',s.end))
		# url <- paste0(base.url,school.num)
		# x <- read.value.from.url(url)
		x <- read.one.url(i,base.url)
		if(!is.na(x[[1]])){
			df$id[cnt]        <- i
			df$address[cnt]   <- x[[1]][1]
			df$enrolment[cnt] <- x[[2]][1]
			cnt <- cnt + 1
		}
	}
	df.school <- na.omit(df)
	save(df.school,file=paste0('df-schools-ontario-',s.start,'.RData'))
}

# })

t1 <- as.numeric(Sys.time())
print(paste('Time elapsed:', round((t1-t0)/60,1),'minutes'))

