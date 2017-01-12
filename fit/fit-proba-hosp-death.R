age <- read.csv('../data/ages/size-distrib-ages.csv')
frail <- read.csv('../data/frailty/frailty-canada.csv')

age.select <- subset(age, age%in%frail$age)

frail.wavg <- sum(frail$frailty * age.select$prop) / sum(age.select$prop)


TARGET.proba.hosp <- 75 / 1e5  # per 100,000

alpha.hosp.fit <- TARGET.proba.hosp / frail.wavg

TARGET.proba.death.cond.hosp <- 0.30
alpha.death.fit <- TARGET.proba.death.cond.hosp / frail.wavg


print('Suggested alpha values for a fit to:')
print(paste('Target proba hospitalization: ',TARGET.proba.hosp))
print(paste('Target proba death cond. on hosp: ',TARGET.proba.death.cond.hosp))
print('\n')
print(paste('--> alpha_hosp  =',alpha.hosp.fit))
print(paste('--> alpha_death =',alpha.death.fit))

