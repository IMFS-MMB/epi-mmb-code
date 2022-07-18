currentdir = getwd()

tabval <- read.table(paste(currentdir , "\\initval.txt",sep = ""))

#iniI0 <- tabval[1,1]
dt_results = tabval

write.csv(dt_results,paste(currentdir , "\\dttest.csv",sep = ""), row.names = FALSE)