###
### Check various fits versus data
###
### WARNING: This will take the last run simulation results from the 'simul' folder.
###

# Folder paths:
dir.simul <- '../simul/'
dir.data  <- '../data/'

source(paste0(dir.simul,'analysis_tools.R'))


# ==== Load data ====

data.age        <- read.csv(paste0(dir.data,'ages/size-distrib-ages.csv'))
data.frailty    <- read.csv(paste0(dir.data,'frailty/frailty-canada.csv'))
data.hosp.death <- read.csv(paste0(dir.data,'frailty/ontario/ontario-influenza-hosp-death.csv'))

# clean up:
data.hosp.death$age_group <- as.character(data.hosp.death$age_group)


# ==== Load simulation ====

load(paste0(dir.simul,'mc-simul.RData'))
# Look at baseline only:
res.select <- res.list.0
# Merge populations of all MC iterations:
pop.all.mc   <- merge.pop.mc(res.select,
							 n.cpu = parallel::detectCores() - 1, 
							 doparallel = TRUE)
# Filter out the MC iterations that produced fizzles:
pop.nofizz <- filter.out.fizzle(pop.all.mc)
# For display, just the first non-fizzled:
idx.mc.no.fizz <- unique(pop.nofizz$mc)

pop <- subset(pop.nofizz, mc==idx.mc.no.fizz[1])

pop.nofizz$age.round <- round(pop.nofizz$age,0)


# ==== Age distribution ====

# Calculate proportion for each MC iteration:
d.num <- ddply(pop.nofizz,c('mc'), summarize, tot = length(id_indiv))
d.age <- ddply(pop.nofizz,c('age.round','mc'), summarize, n = length(id_indiv))
d.age <- join(d.age,d.num, by='mc')
d.age$p <- d.age$n/d.age$tot

g.age <- ggplot(data = d.age) + geom_boxplot(aes(x=age.round, y=p, group=age.round),fill='lightgrey', colour='grey')
g.age <- g.age + geom_point(data = data.age, aes(x=age,y=prop), colour='red')
g.age <- g.age + ggtitle('Age distribution') + xlab('age') + ylab('proportion')


# ==== Frailty/chronic disease ====

age.fra <- ddply(pop.nofizz,c('age.round', 'mc'),summarize, 
					 fra = mean(frailty))

g.age.fra <- ggplot(age.fra) + geom_boxplot(aes(x=age.round, y=fra, group=age.round),fill='lightgrey', colour='grey')
g.age.fra <- g.age.fra + coord_cartesian(ylim = c(0,1))
g.age.fra <- g.age.fra + geom_point(data=data.frailty, aes(x=age, y=frailty), colour='red', size=4)
g.age.fra <- g.age.fra + ggtitle('Frailty by age') + xlab('age') + ylab('frailty')


# ==== Hospitalization & death ====

# Summarize data:
pop.nofizz$is_dead <- 1-pop.nofizz$is_alive
pop.nofizz$age_group <- NA
pop.nofizz$age_group[pop.nofizz$age.round <1] <- '<1'
pop.nofizz$age_group[pop.nofizz$age.round >=1 & pop.nofizz$age.round <5] <- '1_4'
pop.nofizz$age_group[pop.nofizz$age.round >=5 & pop.nofizz$age.round <15] <- '5_14'
pop.nofizz$age_group[pop.nofizz$age.round >=15 & pop.nofizz$age.round <25] <- '15_24'
pop.nofizz$age_group[pop.nofizz$age.round >=25 & pop.nofizz$age.round <45] <- '25_44'
pop.nofizz$age_group[pop.nofizz$age.round >=45 & pop.nofizz$age.round <65] <- '45_64'
pop.nofizz$age_group[pop.nofizz$age.round >=65 ] <- '65+'

df.hd <- ddply(pop.nofizz, c('age_group','mc'), summarize,
			   nh = sum(was_hosp),
			   nd = sum(is_dead),
			   n  = length(id_indiv))

df.hd$prop.hosp <- df.hd$nh/df.hd$n * 100000
df.hd$prop.dead <- df.hd$nd/df.hd$n * 100000

data.hosp <- subset(data.hosp.death, type='hosp')
data.death <- subset(data.hosp.death, type='death')

g.age.hosp <- ggplot(df.hd) + geom_boxplot(aes(x=age_group, y=prop.hosp, group=age_group),fill='lightgrey', colour='black')
g.age.hosp <- g.age.hosp + ggtitle('Hospitalizations by age')
g.age.hosp <- g.age.hosp + xlab('age')+ylab('Hospitalization ratio per 100,000')
g.age.hosp <- g.age.hosp + geom_boxplot(data=data.hosp, aes(x=age_group, y=rate_100000), colour='red2',fill='red', alpha=0.15)
g.age.hosp <- g.age.hosp + scale_y_log10()

g.age.death <- ggplot(df.hd) + geom_boxplot(aes(x=age_group, y=prop.dead, group=age_group),fill='lightgrey', colour='black')
g.age.death <- g.age.death + ggtitle('Death by age')
g.age.death <- g.age.death + xlab('age')+ylab('Death ratio per 100,000')
g.age.death <- g.age.death + geom_boxplot(data=data.death, aes(x=age_group, y=rate_100000), colour='red2',fill='red', alpha=0.15)
g.age.death <- g.age.death + scale_y_log10()


# ==== Final plot ====

pdf('check-all-fits-data.pdf', width = 10, height = 10)

theme_set(theme_bw())
grid.arrange(g.age, 
			 g.age.fra,
			 g.age.hosp,
			 g.age.death)

# Social places size distribution 
plot.sp.sz.distrib(res.list.0, world.prm)

dev.off()