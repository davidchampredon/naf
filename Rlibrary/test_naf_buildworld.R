library(ggplot2)
library(plyr)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

save.plot.to.file <- T

t0 <- as.numeric(Sys.time())

prm <- list()
simul.prm <- list()



prm <- list()

prm[["id_au"]] <- c(0,1,2,3)
prm[["name_au"]] <- c("AUzero","AUone","AUtwo","AUthree")

prm[["hh_size"]] <- c(1,2,3)
prm[["hh_size_proba"]] <- c(0.25, 0.50, 0.25)

prm[["age_hh_00"]] <- c(22,33,44)
prm[["age_hh_00_proba"]] <- c(0.2, 0.5, 0.3)
prm[["age_hh_10"]] <- c(22,33,44)
prm[["age_hh_10_proba"]] <- c(0.2, 0.5, 0.3)
prm[["age_hh_11"]] <- c(5,22,33,44)
prm[["age_hh_11_proba"]] <- c(0.2,0.1, 0.6, 0.1)
prm[["age_hh_20"]] <- c(22,33,44)
prm[["age_hh_20_proba"]] <- c(0.2, 0.5, 0.3)
prm[["age_hh_21"]] <- c(22,33,44)
prm[["age_hh_21_proba"]] <- c(0.2, 0.5, 0.3)
prm[["age_hh_22"]] <- c(5,22,33,44)
prm[["age_hh_22_proba"]] <- c(0.85, 0.05, 0.05, 0.05)

prm[["wrk_size"]] <- c(5,20,40,80)
prm[["wrk_size_proba"]] <- c(0.55, 0.3, 0.1, 0.05)
prm[["pubt_size"]] <- c(20,50,100)
prm[["pubt_size_proba"]] <- c(0.3,0.4,0.3)
prm[["school_size"]] <- c(30,60,90)
prm[["school_size_proba"]] <- c(0.7,0.2,0.1)

prm[["n_hh"]] <- c(5000,2500,1000,900)
prm[["n_wrk"]] <- c(10, 5, 3, 2)
prm[["n_pubt"]] <- c(7,2,1,0)
prm[["n_school"]] <- c(2,2,1,1)

x <- build_world(prm)
z <- lapply(x, as.data.frame)
world <- do.call('rbind',z)

ws <- ddply(world, c('id_au','sp_type','n_linked'), summarize, 
			n = length(id_indiv) ,
			ma = mean(age))
ws
ggplot(ws)+geom_bar(aes(x=n_linked, y=ma),stat = 'identity') + facet_wrap(~id_au)
