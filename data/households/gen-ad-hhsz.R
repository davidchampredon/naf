###
###   GENERATE SYNTHETIC AGE DISTRIBUTIONS 
###   CONDITIONAL ON HOUSEHOLD SIZE
###

source('gen-ad-hhsz-fcts.R')

# Note: the values below for agemean.vec and a.vec
# were obtained by:
# 1/ running an ABC fit (90,000 iterations)
# 2/ manual adjusting the values to fine-tune fit
#
# The target data were ages distribution for Ontario only. 

agemean.vec <- c(49, 51, 23.2, 59.4, 52.1, 20.2)  
a.vec       <- c(1.57, 7.07, 4.54, 3.27, 4.37, 0.95)

pdf('plot-ad-hhsz.pdf', width=15, height = 10)
par(mfrow=c(6,6), mar=rep(0.9,4))
gen.all.ad.hhsz(agemean.vec, a.vec, doplot=TRUE)
dev.off()
