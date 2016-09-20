load('size-distrib-ages-ontario.RData')

# Only data at the provincial level:
df <- subset(df.ontario, geo=='35' & sex==1)

# Format before saving:
df.sav <- data.frame(age=df$age, prop=df$val/sum(df$val))


write.csv(df.sav, file='size-distrib-ages.csv',row.names = F)
