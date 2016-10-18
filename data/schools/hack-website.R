library(RCurl)
library(XML)
library(stringr)
library(profvis)


t0 <- as.numeric(Sys.time())
# profvis({

a <- commandArgs(trailingOnly = TRUE)
s.start <- as.numeric(a[1])
s.end   <- as.numeric(a[2])

base.url <- 'https://www.app.edu.gov.on.ca/eng/sift/schoolProfile.asp?SCH_NUMBER='

# Read the list of all school numbers in Ontario:
fname <- 'school-number-ontario.txt'
sn.all <- read.table(fname,stringsAsFactors = FALSE)[,1]

read.value.from.url <- function(url) {
	# Retrieve Address and Enrolment for a given URL (i.e one school)
	html <- getURL(url, followlocation = TRUE)
	
	found <- is.na(str_locate(pattern='No schools were found',html)[,1])
	
	if(!found) {
		print(paste('school not found.'))
		return(NA)
	}
	
	pos.a <- str_locate_all(pattern ='Address', html)[[1]][,2]
	pos.e <- str_locate_all(pattern ='Enrolment', html)[[1]][,2]
	
	if(length(pos.a)==0 | length(pos.e)==0) {
		print('information not found')
		return(NA)
	}
	
	val.a <- substr(html,pos.a+11, pos.a+100)
	val.e <- substr(html,pos.e+11, pos.e+17)
	
	val.a <- trimws(substr(val.a,1,str_locate(val.a,"<")-1))
	val.e <- as.numeric(substr(val.e,1,str_locate(val.e,"<")-1))
	return(list(val.a,val.e))
}

read.one.url <- function(i,base.url,sn.all) {
	ni <- nchar(as.character(sn.all[i]))
	school.num <- paste0(paste(rep('0',6-ni),collapse = ""),as.character(sn.all[i]))
	if(i%%5==0 | TRUE)  print(paste(s.start,'/',i,'/',s.end,'-->',school.num))
	url <- paste0(base.url,school.num)
	return(read.value.from.url(url))
}

# for tests:
# s.start <- 4700
# s.end   <- 4710

ns <- s.end-s.start+1
print(paste('Parsing',ns,'URLs ...', '(',s.start,s.end,')'))
navec <- rep(NA,ns)
df <- data.frame(id=navec, address=navec, enrolment=navec)

# Loop through the requested school numbers:
cnt <- 1
for (i in s.start:s.end){
	x <- read.one.url(i,base.url,sn.all)
	if(!is.na(x[[1]])){
		df$id[cnt]        <- i
		df$address[cnt]   <- x[[1]][1]
		df$enrolment[cnt] <- x[[2]][1]
		cnt <- cnt + 1
	}
}
df.school <- na.omit(df)
save(df.school,file=paste0('df-schools-ontario-',s.start,'-',s.end,'.RData'))

# }) # ((profvis))

t1 <- as.numeric(Sys.time())
print(paste('Time elapsed:', round((t1-t0)/60,1),'minutes'))

