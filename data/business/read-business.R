library(tidyr)
library(ggplot2)

plot.biz.size.prov <- function(filename) {
	b.prov <- read.csv(filename) #'canada-business-size-prov.csv')
	
	b.prov %>% gather(key=size, value=value,-province,-year) -> b.prov
	b.prov %>% subset(size!='Total' ) -> b.prov
	
	g <- ggplot(b.prov) + geom_bar(aes(x=size, y=value, fill=factor(year)),stat='identity',position='dodge')
	g <- g +facet_wrap(~province, scales = 'free') + theme(axis.text.x = element_text(angle = 20, hjust = 1))
	g <- g + scale_y_log10()
	g <- g + ggtitle('Business sizes by Province')
	plot(g)
}


plot.biz.size.canada.finer <- function(filename) {
	b.det <- read.csv(filename)  #'canada-business-size-detail.csv')
	b.det$size <- factor(b.det$size, levels = c('1_4','5_9','10_19','20_49','50_99','100_199','200_499','500+' ))
	g <- ggplot(b.det) + geom_bar(aes(x=size, y=Cumulative.percentage, fill=factor(year)),stat='identity',position='dodge')
	g <- g + ggtitle('Business sizes for Canada') 
	g <- g + geom_text(aes(x=size,y=Cumulative.percentage+2,label=Cumulative.percentage,colour=factor(year),hjust=factor(year)))
	plot(g)
}


plot.biz.size.prov('canada-business-size-prov.csv')
plot.biz.size.canada.finer('canada-business-size-detail.csv')



