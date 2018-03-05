### 
### merge results and calculate sensitivity for vaccine coverage
###


library(dplyr)

prm.fname <- '../param-model/prm-interv-sensi-vaxcvgmult-'
plist <- system(paste0('ls ',prm.fname,'*'), intern = TRUE)

fname <- '../results/hosp-death-ratio-sensi-vaxcvgmult-'
flist <- system(paste0('ls ',fname,'*'), intern = TRUE)

df <- list()

for(i in 1:length(flist)){
    
    tmp <- read.csv(plist[i], stringsAsFactors = FALSE, strip.white = TRUE)
    mult_val <- as.numeric(tmp$value[tmp$variable=='interv_cvg_age_mult'])
    
    df[[i]] <- read.csv(flist[i])
    df[[i]]$scen <- i
    df[[i]]$mult_val <- mult_val
}
df.all <- do.call('rbind', df)

z <- df.all %>%
    group_by(scen, mult_val) %>%
    summarize(hosp_100k_mean = mean(hosp.per.100000),
              death_100k_mean = mean(death.per.100000),
              n_mc = length(scen))

write.csv(x=z, file='out-sensi-vaxcvgmult.csv')

