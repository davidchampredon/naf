### 
###   ANALYZE SIMULATIONS ALREADY RUN
###

library(ggplot2)
source('analysis_tools.R')

save.plot.to.file <- TRUE

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
print('... simulation results loaded.')


### ==== Merge all MC iterations ====

pop <- merge.pop.mc(res.list)
ts  <- merge.ts.mc(res.list)
tsc <- merge.ts.mc(res.list,is.contact = TRUE)
# acm <- average.age.contact(res.list)

### ==== Plots ====


# Population:
if (save.plot.to.file) pdf('plot_pop.pdf', width = 30, height = 18)

try( plot.population(pop, split.mc=T),  silent = T)
try( plot.n.contacts(tsc),  silent = T)
try( plot.age.contact.matrix.avg(res.list),  silent = T)
try( plot.sp.sz.distrib.new(pop,world.prm) , silent = T)
try( plot.share.same.hh(pop), silent = T)

if (save.plot.to.file) dev.off()
