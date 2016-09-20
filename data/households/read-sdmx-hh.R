library(rsdmx)
library(plyr)
library(dplyr)

# Source: Statistics Canada, 2011 Census of Population, Statistics Canada Catalogue no. 98-313-XCB2011022
# http://www12.statcan.gc.ca/census-recensement/2011/dp-pd/tbt-tt/Rp-eng.cfm?LANG=E&APATH=3&DETAIL=0&DIM=0&FL=A&FREE=0&GC=0&GID=0&GK=0&GRP=1&PID=102236&PRID=0&PTYPE=101955&S=0&SHOWALL=0&SUB=0&Temporal=2011&THEME=91&VID=0&VNAMEE=&VNAMEF=

sdmxobj.g <- readSDMX(file = 'test.xml', isURL = FALSE)
sdmxobj.s <- readSDMX(file = 'Structure_98-313-XCB2011022.xml', isURL = FALSE)

# Retrieve the codes from the RSDMX object:
a <- methods::slot(sdmxobj.s,'codelists')

data.header <- a@header@Name
print(data.header)

# Unpack all structural information:
b   <- a@codelists

geo.sdmx    <- b[[1]] # GEO
dtype.sdmx  <- b[[2]] # DTYPE
hhtype.sdmx <- b[[3]] # HHTYPE
nunits.sdmx <- b[[4]] # NUNITS
obs.sdmx    <- b[[5]] # OBS


extract.struct.label <- function(struct.sdmx, geo=FALSE) {
	# Extract the label of codes in the SDM structure:
	z <- struct.sdmx@Code
	n <- length(z)
	res <- data.frame(id = rep(NA,n), label=rep(NA,n))
	for(i in 1:n){
		res$id[i] <- z[[i]]@id
		if(geo) res$label[i] <- z[[i]]@label$default
		else res$label[i] <- z[[i]]@label$en
	}	
	return(res)
}

geo    <- extract.struct.label(geo.sdmx, geo=TRUE)
dtype  <- extract.struct.label(dtype.sdmx)
hhtype <- extract.struct.label(hhtype.sdmx)
nunits <- extract.struct.label(nunits.sdmx)


# WARNING: SLOW!
t0   <- as.numeric(Sys.time())
df.g <- as.data.frame(sdmxobj.g) ; save(df.g, file = 'dfg.RData')
t1   <- as.numeric(Sys.time())

print(paste('readSDMX:',format(object.size(sdmxobj.g), units='KB')))
print(paste('data.frame:',format(object.size(df.g), units='MB')))
print(paste("Data frame conversion took",round((t1-t0)/3600,2),"hours"))


# Ignore the type of dwelling, just take the total number:
df <- subset(df.g, DType=="1")
df <- select(df, -DType)

# Insert geographical information:
names(geo) <- c('GEO', 'GEO.label')
df        <- join(df, geo, by='GEO')

# Insert houshold type information:
hhtype.df <- data.frame(A11_HHtype_v1=hhtype$id, hhtype=hhtype$label)
df        <- join(df, hhtype.df, by='A11_HHtype_v1')
names(df)[names(df)=='A11_HHtype_v1'] <- 'hhtype.id'

# Insert Household size information:
hhsz <- data.frame(NUnits=nunits$id, size.label=nunits$label)
df   <- join(df,hhsz,by='NUnits')
df$todelete <- FALSE
df$todelete[grepl('Total',df$size.label)]   <- TRUE
df$todelete[grepl('Number',df$size.label)]  <- TRUE
df$todelete[grepl('Average',df$size.label)] <- TRUE
df <- subset(df, todelete==FALSE)
df <- select(df, -todelete)
df$size.label <- as.character(df$size.label)
df$hhsize <- as.numeric(substr(df$size.label,start = 3, stop = 3))

obj.to.save <- list(data.header, df)

save(obj.to.save, file='statCanada-census-households.RData')

print("Read SDMX completed.")