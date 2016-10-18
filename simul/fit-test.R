agemean.vec <- c(65,59,59,35,35,12)    #
a.vec       <- c(1.5,2.5,2.5,1.8,1.8,2.1)
x0 <- c(agemean.vec, a.vec)

EF <- error.function(agemean.vec,a.vec)
plot(EF$t.ad,main='1')
lines(EF$target.reg)
points(EF$sim.age, EF$sim.age.prop, col='red',pch=2)
lines(EF$sim.reg, col = 'red')
print(EF$err)

### 2 

agemean.vec <- c(33, 43,31,25,22,9)   #c(65,59,59,35,35,12) 
a.vec       <- c(1.5,2.5,2.5,1.8,1.8,2.1)
x0 <- c(agemean.vec, a.vec)

EF <- error.function(agemean.vec,a.vec)
plot(EF$t.ad, main='2')
lines(EF$target.reg)
points(EF$sim.age, EF$sim.age.prop, col='red',pch=2)
lines(EF$sim.reg, col = 'red')
print(EF$err)