#read the file with climate data for all populations
clim<-read.table("timema_pops_pcvars.csv", sep=",", header=TRUE)

#subset climate data based on population for each species
bart<-c("BMCG3","BMP90","BMPCT","JL","PCTCR","PCT8000ft")
bart_clim<-clim[clim$Site %in% bart,]
cali<-c("LICK","LP","SM") #CORRECT
cali_clim<-clim[clim$Site %in% cali,]
chum<-c("BALD","BMT","BS","DZ243","GR104","GR603","GR806","HF4","HF6","HFRBP","HFRS","HFTP","PF243")
chum_clim<-clim[clim$Site %in% chum,]
cris<-c("BY","ECCAMP","OUT","R23","R9","VP") #CORRECT
cris_clim<-clim[clim$Site %in% cris,]
knul<-c("BCE","BCTUR","BCWP","H1M37","HB") #CORRECT
knul_clim<-clim[clim$Site %in% knul,]
land<-c("BCBOG","BCHC","BCOG","BCSUM") #correct
land_clim<-clim[clim$Site %in% land,]
podu<-c("BMCG3","BME","BMLC","BMOKC","BMPCT","BMT","BS","DZ243","PCT8000ft", "PF243","SRHWY")
podu_clim<-clim[clim$Site %in% podu,]
popp<-c("FROCK","LP","MM","SM","TBARN") #correct
popp_clim<-clim[clim$Site %in% popp,]

#transpose to get rows for writing out environment files
bart_t<-t(bart_clim[,-1])
cali_t<-t(cali_clim[,-1])
chum_t<-t(chum_clim[,-1])
cris_t<-t(cris_clim[,-1])
knul_t<-t(knul_clim[,-1])
land_t<-t(land_clim[,-1])
podu_t<-t(podu_clim[,-1])
popp_t<-t(popp_clim[,-1])

#for loops for writing out environment file for each climate layer
for (i in 1:nrow(bart_t)){
	write.table(t(bart_t[i,]), file = paste("bart",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(cali_t)){
        write.table(t(cali_t[i,]), file = paste("cali",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(cris_t)){
        write.table(t(cris_t[i,]), file = paste("cris",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(chum_t)){
        write.table(t(chum_t[i,]), file = paste("chum",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(knul_t)){
        write.table(t(knul_t[i,]), file = paste("knul",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(land_t)){
        write.table(t(land_t[i,]), file = paste("land",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(popp_t)){
        write.table(t(popp_t[i,]), file = paste("popp",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

for (i in 1:nrow(podu_t)){
        write.table(t(podu_t[i,]), file = paste("podu",as.character(i),".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
}

