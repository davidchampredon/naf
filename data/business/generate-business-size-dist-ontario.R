library(tidyr)
library(stringr)

filename <- 'canada-business-size-prov.csv'
b.prov <- read.csv(filename) 

b.prov %>% gather(key=size, value=value,-province,-year) -> b.prov
b.prov %>% subset(size!='Total' ) -> b.prov


b.ont <- subset(b.prov, province=='Ontario' & year=='2012')
b.ont$prop <- b.ont$value / sum(b.ont$value)


filename2 <- 'canada-business-size-detail.csv'
b.canada <- read.csv(filename2)  
# b.canada$size <- factor(b.canada$size, levels = c('1_4','5_9','10_19','20_49','50_99','100_199','200_499','500+' ))
b.canada$size <- as.character(b.canada$size)
b.canada <- subset(b.canada, year=='2012')


b.ont
b.ont$size <- gsub(pattern = 'size_','',b.ont$size)
pos <- str_locate(b.ont$size, '_')
b.ont$uppersize <- as.numeric(substr(b.ont$size,pos[,1]+1, nchar(b.ont$size)))
b.ont$cum.prop <- cumsum(b.ont$prop)

b.canada
pos <- str_locate(b.canada$size, '_')
b.canada$uppersize <- as.numeric(substr(b.canada$size,pos[,1]+1, nchar(b.canada$size)))


# Assume distribution for Ontario is similar as for whole Canada

b.save <- data.frame(size = b.canada$uppersize, prop=b.canada$Total/sum(b.canada$Total))
b.save$size[nrow(b.save)] <- 1000
# business size saturate when considering human contacts:
saturate.size <- 199
b.save <- subset(b.save, size <=saturate.size)
n <- nrow(b.save)
b.save$prop[n] <- 1-sum(b.save$prop[1:(n-1)])

write.csv(b.save, file='size-distrib-workplace.csv', row.names = FALSE)

