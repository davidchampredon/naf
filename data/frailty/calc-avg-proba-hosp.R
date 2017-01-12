age <- read.csv('../ages/size-distrib-ages.csv')
frail <- read.csv('frailty-canada.csv')

age.select <- subset(age, age%in%frail$age)

avg.proba.hosp.div.alpha <- sum(frail*age.select) / sum(age.select)
