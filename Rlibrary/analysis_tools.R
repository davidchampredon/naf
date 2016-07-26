
library(plyr)
library(gridExtra)


plot.population <- function(pop) {
	
	g <- ggplot(pop)
	g.age <- g + geom_histogram(aes(x=age), binwidth=2, fill='darkgrey', colour = 'black')
	g.age <- g.age + ggtitle('Age distribution')
	
	linew = 1
	alpha = 0.8
	m_dol_drawn = mean(pop$dol_drawn)
	m_doi_drawn = mean(pop$doi_drawn)

	
	g.dol.drawn <- g + geom_histogram(aes(x=dol_drawn), fill='orange',colour='orange', size=linew, alpha=alpha, binwidth = 0.5)
	g.dol.drawn <- g.dol.drawn + geom_vline(xintercept = m_dol_drawn, colour='orange',linetype = 2)
	g.dol.drawn <- g.dol.drawn + ggtitle(paste0("DOL drawn distribution (mean = ", round(m_dol_drawn,3),")"))
	
	g.doi.drawn <- g + geom_histogram(aes(x=doi_drawn), fill='red',colour='red', size=linew, alpha=alpha, binwidth = 0.5)
	g.doi.drawn <- g.doi.drawn + geom_vline(xintercept = m_doi_drawn, colour='red',linetype = 2)
	g.doi.drawn <- g.doi.drawn + ggtitle(paste0("DOI drawn distribution (mean = ", round(m_doi_drawn,3),")"))
	
	pop.hosp <- subset(pop, (is_hosp>0 | is_discharged>0))
	m_dobh_drawn = mean(pop.hosp$dobh_drawn)
	m_doh_drawn  = mean(pop.hosp$doh_drawn)
	
	g.hosp <- ggplot(pop.hosp)
	g.dobh.drawn <- g.hosp + geom_histogram(aes(x=dobh_drawn), fill='purple',colour='purple', size=linew, alpha=alpha, binwidth = 0.5)
	g.dobh.drawn <- g.dobh.drawn + geom_vline(xintercept = m_dobh_drawn, colour='purple',linetype = 2)
	g.dobh.drawn <- g.dobh.drawn + ggtitle(paste0("DOBH drawn distribution, among hospitalized (mean = ", round(m_dobh_drawn,3),")"))
	
	g.doh.drawn <- g.hosp + geom_histogram(aes(x=doh_drawn), fill='green4',colour='green4', size=linew, alpha=alpha, binwidth = 0.5)
	g.doh.drawn <- g.doh.drawn + geom_vline(xintercept = m_doh_drawn, colour='green4',linetype = 2)
	g.doh.drawn <- g.doh.drawn + ggtitle(paste0("DOH drawn distribution, among hospitalized (mean = ", round(m_doh_drawn,3),")"))
	
	
	
	g.frail.hosp <- ggplot(pop)
	g.frail.hosp <- g.frail.hosp + geom_smooth(aes(x=frailty, y=is_discharged),method = "glm", 
						 method.args = list(family = "binomial"), 
						 se = FALSE, 
						 colour='red', size=2)
	g.frail.hosp <- g.frail.hosp + geom_point(aes(x=frailty, y=is_discharged), size=0.5)
	
	
	g.imm.hosp <- ggplot(pop)
	g.imm.hosp <- g.imm.hosp + geom_smooth(aes(x=immunity, y=is_discharged),method = "glm", 
											   method.args = list(family = "binomial"), 
											   se = FALSE, 
											   colour='red', size=2)
	g.imm.hosp <- g.imm.hosp + geom_point(aes(x=immunity, y=is_discharged), size=0.5)
	
	
	
	grid.arrange(g.age, 
				 g.dol.drawn, 
				 g.doi.drawn,
				 g.dobh.drawn,
				 g.doh.drawn,
				 g.imm.hosp,
				 g.frail.hosp)
	
}

plot.epi.timeseries <- function(ts){
	g.SR <- ggplot(ts, aes(x=time))
	g.SR <- g.SR + geom_step(aes(y=nS),colour='springgreen3') 
	g.SR <- g.SR + geom_step(aes(y=nR),colour='blue')
	g.SR <- g.SR + ggtitle("Susceptible and recovered") + ylab("")
	
	g.nE.nI <- ggplot(ts, aes(x=time))
	g.nE.nI <- g.nE.nI + geom_step(aes(y=nE),colour='orange') 
	g.nE.nI <- g.nE.nI + geom_step(aes(y=nIs),colour='red',size=2)
	g.nE.nI <- g.nE.nI + geom_step(aes(y=nIa),colour='red')
	g.nE.nI <- g.nE.nI + ggtitle("All latent (nE, orange) and Infectious (Ia & Is[bold])") + ylab("")
	
	g.prev <- ggplot(ts, aes(x=time))
	g.prev <- g.prev + geom_step(aes(y=prevalence))
	g.prev <- g.prev + ggtitle("Prevalence") + ylab("")
	
	grid.arrange(g.SR ,
				 g.nE.nI, 
				 g.prev)
}
