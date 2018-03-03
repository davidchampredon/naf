
# Retrieve Data 
dat <- read.csv('../../data/frailty/frailty-canada.csv')

# Segmented regression
library(segmented)

lin.mod <- lm(data = dat, formula = frailty ~age)
seglm <- segmented(lin.mod, seg.Z = ~age, psi=25)



#plot(seglm, add=T, lwd =8, col=rgb(0,0,0,0.3))

print(seglm)
cf <- seglm$coefficients
brk <- seglm$psi
print(paste('slope2=',cf[2]+cf[3]))


x1 <- seq(0,brk[2], by=0.1)
x2 <- seq(brk[2], 100, by=0.1)

y1 <- cf[1] + cf[2]*x1
y2 <- (cf[1] + cf[2]*brk[2]) + (cf[2]+cf[3])*(x2-brk[2])

col <- rgb(0,0,0,1)
lwd <- 8

# PLOT

pdf('../figures/frailty.pdf', width=10, height = 5)
par(cex.lab=1.5, cex.axis=1.5)
plot(dat$age, dat$frailty, 
     las=1,
     xlab = 'Age (years)',
     ylab = 'Frailty index',
     cex = 2, lwd = 6, col='red',
     xlim = range(c(0,dat$age)), ylim=c(0,1))
grid()
lines(x1,y1, lwd=lwd, col=col)
lines(x2,y2, lwd=lwd, col=col)
dev.off()

