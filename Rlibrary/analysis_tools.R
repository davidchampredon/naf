
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
	
	g.dol.drawn <- g + geom_histogram(aes(x=dol_drawn), fill='orange',colour='orange', size=linew, alpha=alpha, bins = 50)
	g.dol.drawn <- g.dol.drawn + geom_vline(xintercept = m_dol_drawn, colour='orange',linetype = 2)
	g.dol.drawn <- g.dol.drawn + ggtitle(paste0("DOL drawn distribution (mean = ", round(m_dol_drawn,3),")"))
	
	g.doi.drawn <- g + geom_histogram(aes(x=doi_drawn), fill='red',colour='red', size=linew, alpha=alpha, bins = 50)
	g.doi.drawn <- g.doi.drawn + geom_vline(xintercept = m_doi_drawn, colour='red',linetype = 2)
	g.doi.drawn <- g.doi.drawn + ggtitle(paste0("DOI drawn distribution (mean = ", round(m_doi_drawn,3),")"))
	
	grid.arrange(g.age, 
				 g.dol.drawn, 
				 g.doi.drawn)
	
}

plot.epi.timeseries <- function(ts){
	g.SR <- ggplot(ts, aes(x=time))
	g.SR <- g.SR + geom_line(aes(y=nS),colour='springgreen3') 
	g.SR <- g.SR + geom_line(aes(y=nR),colour='blue')
	g.SR <- g.SR + ggtitle("Susceptible and recovered") + ylab("")
	
	g.nE.nI <- ggplot(ts, aes(x=time))
	g.nE.nI <- g.nE.nI + geom_line(aes(y=nE),colour='orange') 
	g.nE.nI <- g.nE.nI + geom_line(aes(y=nIs),colour='red')
	g.nE.nI <- g.nE.nI + ggtitle("All latent (nE) and Infectious (Ia+Is)") + ylab("")
	
	g.prev <- ggplot(ts, aes(x=time))
	g.prev <- g.prev + geom_line(aes(y=prevalence))
	g.prev <- g.prev + ggtitle("Prevalence") + ylab("")
	
	grid.arrange(g.SR ,
				 g.nE.nI, 
				 g.prev)
}
