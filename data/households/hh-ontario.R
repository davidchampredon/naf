###
###  Read and save household data for ONTARIO only
###

library(plyr)

# ==== Data manipulation ==== 

# Load data:
load('statCanada-census-households.RData')

if(exists('list.to.save')) obj.to.save <- list.to.save
df <- obj.to.save[[2]]

# Avoid too long names:
df$GEO.label <- substr(df$GEO.label,1,15)

# Select only total of hh types:
df <- subset(df, hhtype.id=='1')

# Province information:
df$prov.id <- substr(df$GEO,1,2)
prov <- subset(df, nchar(GEO)==2 & hhsize==1)
prov <- prov[,c('prov.id','GEO.label')]
names(prov)[2] <- 'province'
df <- join(df,prov,by='prov.id')

# Select only one province
geo.select <- '35'
df.select  <- subset(df, GEO==geo.select)

# Format before save:
n.hh <- sum(df.select$obsValue)
df.select$prop <- df.select$obsValue / n.hh
df.save <- df.select[,c('GEO.label','hhsize','prop')]

# Save to file:
write.csv(x = df.save, file = 'size-distrib-hh-ontario.csv',row.names = F)
write.table(n.hh,file = 'num-hh-ontario.csv', row.names = F,col.names = F)
