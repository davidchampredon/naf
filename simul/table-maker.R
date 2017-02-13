library(plyr)
library(tidyr)

source('utils-compare.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]


# ==== Load saved data ====

print('Loading saved RData...')
load(paste0(dir.save.rdata,'result-scen-all.RData'))
res.all <- list(main = result.scen.all,
                main.ageGroup = result.scen.all.ageGroup)
print('Saved RData loaded.')
file.scen.prm.list <- 'scenario-prm-list.csv'
dir <- dir.results

df  <- res.all[['main.ageGroup']]
X <- process.outputs(df, dir, file.scen.prm.list)
DAT <- do.call('rbind.data.frame', X$zlist.full)

dh <- subset(DAT, outcome=='Hospitalized')

write.csv(x = dh, file = '../results/table_1.csv')
