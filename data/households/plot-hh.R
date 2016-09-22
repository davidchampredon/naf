library(ggplot2)
library(plyr)

# ==== Data manipulation ==== 

# Load data:
load('statCanada-census-households.RData')

if(exists('list.to.save')) obj.to.save <- list.to.save
df <- obj.to.save[[2]]

# Avoid too long names:
df$GEO.label <- substr(df$GEO.label,1,15)

# Data at the provincial level only:
df.prov.only <- subset(df, nchar(GEO)==2)

# Select only total of hh types:
df <- subset(df, hhtype.id=='1')

# Province information:
df$prov.id <- substr(df$GEO,1,2)
prov <- subset(df, nchar(GEO)==2 & hhsize==1)
prov <- prov[,c('prov.id','GEO.label')]
names(prov)[2] <- 'province'
df <- join(df,prov,by='prov.id')

# ==== Distributio Ratios ====
df <- df[order(df$GEO,df$hhsize),]

df$ratio.to.2 <- NA
for(i in 1:nrow(df)){
	if(df$hhsize[i]==2){
		for( j in 0:5){
			df$ratio.to.2[i-1+j] <- df$obsValue[i-1+j]/df$obsValue[i]
		}
	}
}

ratio.x <- 1
ratio.y <- 3
df.rx <- subset(df, hhsize == ratio.x)
df.ry <- subset(df, hhsize == ratio.y)

# ==== PLOTS ==== 

pdf('plot_hh.pdf', height = 12, width = 20)

df.plot <- data.frame(rx = df.rx$ratio.to.2, 
					  ry = df.ry$ratio.to.2, 
					  GEO = df.rx$GEO,
					  GEO.label= df.rx$GEO.label,
					  province = df.rx$province,
					  prov.id = df.rx$prov.id)
df.plot$GEO <- as.character(df.plot$GEO)

tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")


g <- ggplot(df.plot) + geom_label(aes(x=rx, y=ry, label= GEO.label,
									  fill=province),
								  size = 3, alpha = 0.4)
g <- g + xlab(paste('ratio size',ratio.x,'/','size 2'))+ ylab(paste('ratio size',ratio.y,'/','size 2'))
g <- g + scale_fill_manual(values=tol14rainbow)
g <- g + ggtitle('Distribution size ratios - CANADA')
plot(g)

df.plot.prov <- subset(df.plot, nchar(GEO)==2)
g <- ggplot(df.plot.prov) + geom_label(aes(x=rx, y=ry, label= GEO.label,
									  fill=province),
								  size = 3, alpha = 0.8)
g <- g + xlab(paste('ratio size',ratio.x,'/','size 2'))+ ylab(paste('ratio size',ratio.y,'/','size 2'))
g <- g + scale_fill_manual(values=tol14rainbow)
g <- g + ggtitle('Distribution size ratios - CANADA')
plot(g)


g.prov <- ggplot(df.prov.only) + geom_bar(aes(x=hhsize,y=obsValue,fill=GEO.label), 
										  position='dodge',stat = 'identity')
g.prov <- g.prov + facet_wrap(~GEO.label, scale='free_y') + ggtitle('Households size distribution')
plot(g.prov)


df.ontario <- subset(df.plot, province=='Ontario')
g <- ggplot(df.ontario) + geom_label(aes(x=rx, y=ry, label= GEO.label),
								  size = 3, alpha = 0.4)
g <- g + xlab(paste('ratio size',ratio.x,'/','size 2'))+ ylab(paste('ratio size',ratio.y,'/','size 2'))
g <- g + scale_fill_manual(values=tol14rainbow)
g <- g + ggtitle('Distribution size ratios - Ontario')
plot(g)

dev.off()
