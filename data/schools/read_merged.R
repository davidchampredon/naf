load('./to-keep/dfall.RData')

h <- hist(dfall$enrolment, breaks=seq(0,1000,by=100),freq = T, col='grey')

df.save <- data.frame(school.size=h$mids, prop=h$counts/sum(h$counts))
write.table(df.save,file='size-distrib-schools-ontario.csv', row.names = F, sep=',',quote = F)
