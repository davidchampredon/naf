load('dfall.RData')

print(paste(nrow(dfall),'schools retrieved.'))
print(summary(dfall))
write.csv(x = dfall, file = 'dfall.csv')
emax <- max(dfall$enrolment)+10

### ====  Histogram ====

size.bucket <- 100

h <- hist(dfall$enrolment, 
		  breaks = seq(0,emax, by=size.bucket),
		  las = 1,
		  freq = T, col='lightgrey')
m <- mean(dfall$enrolment)
CI <- 0.80
q <- quantile(dfall$enrolment, probs = c(0.5-CI/2, 0.5, 0.5+CI/2))

abline(v=q, lwd=4, col='orange')
abline(v=m, lwd=2, lty=3, col='orange')

df.save <- data.frame(school.size=h$mids, prop=h$counts/sum(h$counts))
write.table(df.save,file='size-distrib-schools-ontario.csv', 
			row.names = F, sep=',',quote = F)


h <- hist(dfall$enrolment, 
		  breaks = seq(0,emax, by=size.bucket),
		  las = 1,
		  freq = F, col='lightgrey')
lines(density(dfall$enrolment,adjust = 0.5))
xx <- seq(0,emax, by = 1)
yy <- dlnorm(xx,meanlog = mean(log(m)), sdlog = 0.6)
lines(xx,yy,col='red', lty=2)

