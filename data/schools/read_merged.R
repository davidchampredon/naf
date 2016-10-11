load('dfall.RData')
nrow(dfall)
h <- hist(dfall$enrolment, breaks=seq(0,1000,by=50),freq = T, col='grey')
m <- mean(dfall$enrolment)
CI <- 0.80
q <- quantile(dfall$enrolment, probs = c(0.5-CI/2, 0.5, 0.5+CI/2))

abline(v=q, lwd=4, col='orange')
abline(v=m, lwd=2, lty=3, col='orange')

df.save <- data.frame(school.size=h$mids, prop=h$counts/sum(h$counts))
write.table(df.save,file='size-distrib-schools-ontario.csv', 
			row.names = F, sep=',',quote = F)
