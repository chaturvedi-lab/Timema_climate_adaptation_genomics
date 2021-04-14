###### climate layers ##########
pcvars<-list.files(path="../..", pattern = "bin_lgscaf_out.csv",full.names = TRUE)

##get top snps for each species and each PC
#function to do randomizations
funcRandomization<- function(species,out,pc){
                sdat<-cbind(datfile[,1:4], species)
                compdat<-sdat[complete.cases(sdat),]
                spstop<-which(compdat[,5] > quantile(compdat[,5], probs=0.9, na.rm=T))
                results<-compdat[spstop,]
                len<-dim(results)[1]
                spsc<-rep(paste(pc,out,sep = "_"), len)
                results<-cbind(spsc,results)
                colnames(results)<-c("comb","lg","scaffold","length","bin",out)
                results2<-cbind(results[,2:3], "TRUE")
                outname<-paste(pc,out,"topsnps0.9",sep = "_")
                write.table(as.data.frame(results),outname,col.names=T, row.names=F, sep=",",quote=F)
		outname2<-paste(pc,out,"lgscaf",sep = "_")
                write.table(as.data.frame(results2),outname2,col.names=T, row.names=F, sep=" ",quote=F)


}

datfile1<-read.csv(pcvars[1], header=T,sep=",")
datfile2<-read.csv(pcvars[2], header=T,sep=",")
datfile3<-read.csv(pcvars[3], header=T,sep=",")
#pc1
funcRandomization(datfile1$med_bart,"bart","pc1")
funcRandomization(datfile1$med_cali,"cali","pc1")
funcRandomization(datfile1$med_chum,"chum","pc1")
funcRandomization(datfile1$med_cris,"cris","pc1")
funcRandomization(datfile1$med_knul,"knul","pc1")
funcRandomization(datfile1$med_land,"land","pc1")
funcRandomization(datfile1$med_podu,"podu","pc1")
funcRandomization(datfile1$med_popp,"popp","pc1")
#pc2
funcRandomization(datfile2$med_bart,"bart","pc2")
funcRandomization(datfile2$med_cali,"cali","pc2")
funcRandomization(datfile2$med_chum,"chum","pc2")
funcRandomization(datfile2$med_cris,"cris","pc2")
funcRandomization(datfile2$med_knul,"knul","pc2")
funcRandomization(datfile2$med_land,"land","pc2")
funcRandomization(datfile2$med_podu,"podu","pc2")
funcRandomization(datfile2$med_popp,"popp","pc2")
#pc3
funcRandomization(datfile3$med_bart,"bart","pc3")
funcRandomization(datfile3$med_cali,"cali","pc3")
funcRandomization(datfile3$med_chum,"chum","pc3")
funcRandomization(datfile3$med_cris,"cris","pc3")
funcRandomization(datfile3$med_knul,"knul","pc3")
funcRandomization(datfile3$med_land,"land","pc3")
funcRandomization(datfile3$med_podu,"podu","pc3")
funcRandomization(datfile3$med_popp,"popp","pc3")



##get top shared SNPs
#function to do randomizations
funcRandomization<- function(species1,species2,out1,out2,pc){
		sdat<-cbind(datfile[,1:4], species1, species2)
		compdat<-sdat[complete.cases(sdat),]
		sps1top<-which(compdat[,5] > quantile(compdat[,5], probs=0.9, na.rm=T))
		sps2top<-which(compdat[,6] > quantile(compdat[,6], probs=0.9, na.rm=T))
		results<-compdat[which(sps1top %in% sps2top),]
		len<-dim(results)[1]
		spsc<-rep(paste(pc,out1,out2,sep = "_"), len)
		results<-cbind(spsc,results)
		colnames(results)<-c("comb","lg","scaffold","length","bin",out1,out2)
		#return(results)	
		outname<-paste(pc,out1,out2,"topsnps0.9",sep = "_")
		write.table(as.data.frame(results),outname,col.names=T, row.names=F, sep=",",quote=F)
}

datfile<-read.csv(pcvars[1], header=T,sep=",")
#bart
funcRandomization(datfile$med_bart,datfile$med_cris,"bart","cris","pc1")
funcRandomization(datfile$med_bart,datfile$med_podu,"bart","podu","pc1")
#cali
funcRandomization(datfile$med_cali,datfile$med_cris,"cali","cris","pc1")
funcRandomization(datfile$med_cali,datfile$med_land,"cali","land","pc1")
funcRandomization(datfile$med_cali,datfile$med_podu,"cali","podu","pc1")
#chum
funcRandomization(datfile$med_chum,datfile$med_cris,"chum","cris","pc1")
funcRandomization(datfile$med_chum,datfile$med_podu,"chum","podu","pc1")
#cris
funcRandomization(datfile$med_cris,datfile$med_knul,"cris","knul","pc1")
funcRandomization(datfile$med_cris,datfile$med_land,"cris","land","pc1")
funcRandomization(datfile$med_cris,datfile$med_podu,"cris","podu","pc1")
funcRandomization(datfile$med_cris,datfile$med_popp,"cris","popp","pc1")
#knul
funcRandomization(datfile$med_knul,datfile$med_land,"knul","land","pc1")
funcRandomization(datfile$med_knul,datfile$med_popp,"knul","popp","pc1")
#land
funcRandomization(datfile$med_land,datfile$med_podu,"land","podu","pc1")


datfile<-read.csv(pcvars[2], header=T,sep=",")
#bart
funcRandomization(datfile$med_bart,datfile$med_cris,"bart","cris","pc2")
funcRandomization(datfile$med_bart,datfile$med_podu,"bart","podu","pc2")
funcRandomization(datfile$med_bart,datfile$med_popp,"bart","popp","pc2")
#cali
funcRandomization(datfile$med_cali,datfile$med_popp,"cali","popp","pc2")
#knul
funcRandomization(datfile$med_knul,datfile$med_land,"knul","land","pc2")
funcRandomization(datfile$med_knul,datfile$med_podu,"knul","podu","pc2")
#land
funcRandomization(datfile$med_land,datfile$med_popp,"land","popp","pc2")


datfile<-read.csv(pcvars[3], header=T,sep=",")
#bart
funcRandomization(datfile$med_bart,datfile$med_cris,"bart","cris","pc3")
#cali
funcRandomization(datfile$med_cali,datfile$med_knul,"cali","knul","pc3")
funcRandomization(datfile$med_cali,datfile$med_land,"cali","land","pc3")
funcRandomization(datfile$med_cali,datfile$med_podu,"cali","podu","pc3")
#chum
funcRandomization(datfile$med_chum,datfile$med_cris,"chum","cris","pc3")
funcRandomization(datfile$med_chum,datfile$med_podu,"chum","podu","pc3")
#cris
funcRandomization(datfile$med_cris,datfile$med_knul,"cris","knul","pc3")
funcRandomization(datfile$med_cris,datfile$med_land,"cris","land","pc3")
funcRandomization(datfile$med_cris,datfile$med_popp,"cris","popp","pc3")
#knul
funcRandomization(datfile$med_knul,datfile$med_land,"knul","land","pc3")

##moved these files to the annotation folder
