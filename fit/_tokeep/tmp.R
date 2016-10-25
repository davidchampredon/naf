load('fit_hhsz_age_90K.RData')

print(xbest)
print(abmfit$fx.best)
par(mfrow=c(3,4))

for(i in 1:12){
	print(i)
	xx <- abmfit$x.all[,i]
	plot(xx, abmfit$fx.all, col=rgb(0,0,0,0.2), log='y', main=i)
	points(abmfit$x.best[i],abmfit$fx.best[1], pch=16,col='red',cex=2)
}
