###### climate layers ##########
pcvars<-list.files(path="../../../", pattern = "bin_lgscaf_out.csv",full.names = TRUE)


#function to do randomizations
funcRandomization<- function(species1,species2,quantile = NA){
		sdat<-cbind(datfile[,1:3], species1, species2)
		compdat<-sdat[complete.cases(sdat),]
		sps1top<-which(compdat[,4] > quantile(compdat[,4], probs=quantile, na.rm=T))
		sps2top<-which(compdat[,5] > quantile(compdat[,5], probs=quantile, na.rm=T))
		spslen<-dim(compdat)[1]

		null <- rep(NA,10000)
		for (j in 1:10000){
			x<-round((1-quantile)*spslen,0)
			y<-round((1-quantile)*spslen,0)
			bnull<-sample(1:spslen,x,replace=FALSE)
			cnull<-sample(1:spslen,y,replace=FALSE)
			null[j]<-sum(bnull %in% cnull)
		}

		obs<-sum(sps1top %in% sps2top)
		p<-mean(null >= obs)
		xfold<-obs / mean(null)
		results<-cbind(obs,xfold,p)
		#sharedsnps<-datfile[sps1top[which(sps1top %in% sps2top)],]
		return(results)	
}

quants<-c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)

for (var in pcvars[1:3]){
	#print(var)
	datfile<-read.csv(var, header=T,sep=",")

	quant_bart_cali<-vector("list",10)
	quant_bart_chum<-vector("list",10)
	quant_bart_cris<-vector("list",10)
	quant_bart_knul<-vector("list",10)
	quant_bart_land<-vector("list",10)
	quant_bart_podu<-vector("list",10)
	quant_bart_popp<-vector("list",10)
	#cali
	quant_cali_bart<-vector("list",10)
	quant_cali_chum<-vector("list",10)
	quant_cali_cris<-vector("list",10)
	quant_cali_knul<-vector("list",10)
	quant_cali_land<-vector("list",10)
	quant_cali_podu<-vector("list",10)
	quant_cali_popp<-vector("list",10)
	#chum
	quant_chum_bart<-vector("list",10)
	quant_chum_cali<-vector("list",10)
	quant_chum_cris<-vector("list",10)
	quant_chum_knul<-vector("list",10)
	quant_chum_land<-vector("list",10)
	quant_chum_podu<-vector("list",10)
	quant_chum_popp<-vector("list",10)
	#cris
	quant_cris_bart<-vector("list",10)
	quant_cris_cali<-vector("list",10)
	quant_cris_chum<-vector("list",10)
	quant_cris_knul<-vector("list",10)
	quant_cris_land<-vector("list",10)
	quant_cris_podu<-vector("list",10)
	quant_cris_popp<-vector("list",10)
	#knul
	quant_knul_bart<-vector("list",10)
	quant_knul_cali<-vector("list",10)
	quant_knul_chum<-vector("list",10)
	quant_knul_cris<-vector("list",10)
	quant_knul_land<-vector("list",10)
	quant_knul_podu<-vector("list",10)
	quant_knul_popp<-vector("list",10)
	#land
	quant_land_bart<-vector("list",10)
	quant_land_cali<-vector("list",10)
	quant_land_chum<-vector("list",10)
	quant_land_cris<-vector("list",10)
	quant_land_knul<-vector("list",10)
	quant_land_podu<-vector("list",10)
	quant_land_popp<-vector("list",10)
	#podu
	quant_podu_bart<-vector("list",10)
	quant_podu_cali<-vector("list",10)
	quant_podu_chum<-vector("list",10)
	quant_podu_cris<-vector("list",10)
	quant_podu_knul<-vector("list",10)
	quant_podu_land<-vector("list",10)
	quant_podu_popp<-vector("list",10)
	#popp
	quant_popp_bart<-vector("list",10)
	quant_popp_cali<-vector("list",10)
	quant_popp_chum<-vector("list",10)
	quant_popp_cris<-vector("list",10)
	quant_popp_knul<-vector("list",10)
	quant_popp_land<-vector("list",10)
	quant_popp_podu<-vector("list",10)

	for (i in 1:10){
		#bart
		quant_bart_cali[[i]]<-funcRandomization(datfile$med_bart,datfile$med_cali,quants[i])
		quant_bart_chum[[i]]<-funcRandomization(datfile$med_bart,datfile$med_chum,quants[i])
		quant_bart_cris[[i]]<-funcRandomization(datfile$med_bart,datfile$med_cris,quants[i])
		quant_bart_knul[[i]]<-funcRandomization(datfile$med_bart,datfile$med_knul,quants[i])
		quant_bart_land[[i]]<-funcRandomization(datfile$med_bart,datfile$med_land,quants[i])
		quant_bart_podu[[i]]<-funcRandomization(datfile$med_bart,datfile$med_podu,quants[i])
		quant_bart_popp[[i]]<-funcRandomization(datfile$med_bart,datfile$med_popp,quants[i])
		#cali
		quant_cali_bart[[i]]<-funcRandomization(datfile$med_cali,datfile$med_bart,quants[i])
		quant_cali_chum[[i]]<-funcRandomization(datfile$med_cali,datfile$med_chum,quants[i])
		quant_cali_cris[[i]]<-funcRandomization(datfile$med_cali,datfile$med_cris,quants[i])
		quant_cali_knul[[i]]<-funcRandomization(datfile$med_cali,datfile$med_knul,quants[i])
		quant_cali_land[[i]]<-funcRandomization(datfile$med_cali,datfile$med_land,quants[i])
		quant_cali_podu[[i]]<-funcRandomization(datfile$med_cali,datfile$med_podu,quants[i])
		quant_cali_popp[[i]]<-funcRandomization(datfile$med_cali,datfile$med_popp,quants[i])
		#chum
		quant_chum_bart[[i]]<-funcRandomization(datfile$med_chum,datfile$med_bart,quants[i])
		quant_chum_cali[[i]]<-funcRandomization(datfile$med_chum,datfile$med_cali,quants[i])
		quant_chum_cris[[i]]<-funcRandomization(datfile$med_chum,datfile$med_cris,quants[i])
		quant_chum_knul[[i]]<-funcRandomization(datfile$med_chum,datfile$med_knul,quants[i])
		quant_chum_land[[i]]<-funcRandomization(datfile$med_chum,datfile$med_land,quants[i])
		quant_chum_podu[[i]]<-funcRandomization(datfile$med_chum,datfile$med_podu,quants[i])
		quant_chum_popp[[i]]<-funcRandomization(datfile$med_chum,datfile$med_popp,quants[i])
		#cris
		quant_cris_bart[[i]]<-funcRandomization(datfile$med_cris,datfile$med_bart,quants[i])
		quant_cris_cali[[i]]<-funcRandomization(datfile$med_cris,datfile$med_cali,quants[i])
		quant_cris_chum[[i]]<-funcRandomization(datfile$med_cris,datfile$med_chum,quants[i])
		quant_cris_knul[[i]]<-funcRandomization(datfile$med_cris,datfile$med_knul,quants[i])
		quant_cris_land[[i]]<-funcRandomization(datfile$med_cris,datfile$med_land,quants[i])
		quant_cris_podu[[i]]<-funcRandomization(datfile$med_cris,datfile$med_podu,quants[i])
		quant_cris_popp[[i]]<-funcRandomization(datfile$med_cris,datfile$med_popp,quants[i])
		#knul
		quant_knul_bart[[i]]<-funcRandomization(datfile$med_knul,datfile$med_bart,quants[i])
		quant_knul_cali[[i]]<-funcRandomization(datfile$med_knul,datfile$med_cali,quants[i])
		quant_knul_chum[[i]]<-funcRandomization(datfile$med_knul,datfile$med_chum,quants[i])
		quant_knul_cris[[i]]<-funcRandomization(datfile$med_knul,datfile$med_cris,quants[i])
		quant_knul_land[[i]]<-funcRandomization(datfile$med_knul,datfile$med_land,quants[i])
		quant_knul_podu[[i]]<-funcRandomization(datfile$med_knul,datfile$med_podu,quants[i])
		quant_knul_popp[[i]]<-funcRandomization(datfile$med_knul,datfile$med_popp,quants[i])
		#land
		quant_land_bart[[i]]<-funcRandomization(datfile$med_land,datfile$med_bart,quants[i])
		quant_land_cali[[i]]<-funcRandomization(datfile$med_land,datfile$med_cali,quants[i])
		quant_land_chum[[i]]<-funcRandomization(datfile$med_land,datfile$med_chum,quants[i])
		quant_land_cris[[i]]<-funcRandomization(datfile$med_land,datfile$med_cris,quants[i])
		quant_land_knul[[i]]<-funcRandomization(datfile$med_land,datfile$med_knul,quants[i])
		quant_land_podu[[i]]<-funcRandomization(datfile$med_land,datfile$med_podu,quants[i])
		quant_land_popp[[i]]<-funcRandomization(datfile$med_land,datfile$med_popp,quants[i])
		#podu
		quant_podu_bart[[i]]<-funcRandomization(datfile$med_podu,datfile$med_bart,quants[i])
		quant_podu_cali[[i]]<-funcRandomization(datfile$med_podu,datfile$med_cali,quants[i])
		quant_podu_chum[[i]]<-funcRandomization(datfile$med_podu,datfile$med_chum,quants[i])
		quant_podu_cris[[i]]<-funcRandomization(datfile$med_podu,datfile$med_cris,quants[i])
		quant_podu_knul[[i]]<-funcRandomization(datfile$med_podu,datfile$med_knul,quants[i])
		quant_podu_land[[i]]<-funcRandomization(datfile$med_podu,datfile$med_land,quants[i])
		quant_podu_popp[[i]]<-funcRandomization(datfile$med_podu,datfile$med_popp,quants[i])
		#popp
		quant_popp_bart[[i]]<-funcRandomization(datfile$med_popp,datfile$med_bart,quants[i])
		quant_popp_cali[[i]]<-funcRandomization(datfile$med_popp,datfile$med_cali,quants[i])
		quant_popp_chum[[i]]<-funcRandomization(datfile$med_popp,datfile$med_chum,quants[i])
		quant_popp_cris[[i]]<-funcRandomization(datfile$med_popp,datfile$med_cris,quants[i])
		quant_popp_knul[[i]]<-funcRandomization(datfile$med_popp,datfile$med_knul,quants[i])
		quant_popp_land[[i]]<-funcRandomization(datfile$med_popp,datfile$med_land,quants[i])
		quant_popp_podu[[i]]<-funcRandomization(datfile$med_popp,datfile$med_podu,quants[i])

		}

	#write out to a file
	#bart
	bart_quants<- cbind(t(as.data.frame(quant_bart_cali)), t(as.data.frame(quant_bart_chum)), t(as.data.frame(quant_bart_cris)), t(as.data.frame(quant_bart_knul)), t(as.data.frame(quant_bart_land)), t(as.data.frame(quant_bart_podu)), t(as.data.frame(quant_bart_popp)))
	colnames(bart_quants)<-c("bart_cali","bart_chum","bart_cris","bart_knul","bart_land","bart_podu","bart_popp")
	outname<-paste(basename(var),"bart_quants",sep = "_")
	write.table(t(as.data.frame(bart_quants)),outname, col.names=T, row.names=T, quote=F, sep=" ")

	#cali
	cali_quants<- cbind(t(as.data.frame(quant_cali_bart)), t(as.data.frame(quant_cali_chum)), t(as.data.frame(quant_cali_cris)), t(as.data.frame(quant_cali_knul)), t(as.data.frame(quant_cali_land)), t(as.data.frame(quant_cali_podu)), t(as.data.frame(quant_cali_popp)))
	colnames(cali_quants)<-c("cali_bart","cali_chum","cali_cris","cali_knul","cali_land","cali_podu","cali_popp")
	outname<-paste(basename(var),"cali_quants",sep = "_")
	write.table(t(as.data.frame(cali_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#chum
	chum_quants<- cbind(t(as.data.frame(quant_chum_bart)), t(as.data.frame(quant_chum_cali)), t(as.data.frame(quant_chum_cris)), t(as.data.frame(quant_chum_knul)), t(as.data.frame(quant_chum_land)), t(as.data.frame(quant_chum_podu)), t(as.data.frame(quant_chum_popp)))
	colnames(chum_quants)<-c("chum_bart","chum_cali","chum_cris","chum_knul","chum_land","chum_podu","chum_popp")
	outname<-paste(basename(var),"chum_quants",sep = "_")
	write.table(t(as.data.frame(chum_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#cris
	cris_quants<- cbind(t(as.data.frame(quant_cris_bart)), t(as.data.frame(quant_cris_cali)), t(as.data.frame(quant_cris_chum)), t(as.data.frame(quant_cris_knul)), t(as.data.frame(quant_cris_land)), t(as.data.frame(quant_cris_podu)), t(as.data.frame(quant_cris_popp)))
	colnames(cris_quants)<-c("cris_bart","cris_cali","cris_chum","cris_knul","cris_land","cris_podu","cris_popp")
	outname<-paste(basename(var),"cris_quants",sep = "_")
	write.table(t(as.data.frame(cris_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#knul
	knul_quants<- cbind(t(as.data.frame(quant_knul_bart)), t(as.data.frame(quant_knul_cali)), t(as.data.frame(quant_knul_chum)), t(as.data.frame(quant_knul_cris)), t(as.data.frame(quant_knul_land)), t(as.data.frame(quant_knul_podu)), t(as.data.frame(quant_knul_popp)))
	colnames(knul_quants)<-c("knul_bart","knul_cali","knul_chum","knul_cris","knul_land","knul_podu","knul_popp")
	outname<-paste(basename(var),"knul_quants",sep = "_")
	write.table(t(as.data.frame(knul_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#land
	land_quants<- cbind(t(as.data.frame(quant_land_bart)), t(as.data.frame(quant_land_cali)), t(as.data.frame(quant_land_chum)), t(as.data.frame(quant_land_cris)), t(as.data.frame(quant_land_knul)), t(as.data.frame(quant_land_podu)), t(as.data.frame(quant_land_popp)))
	colnames(land_quants)<-c("land_bart","land_cali","land_chum","land_cris","land_knul","land_podu","land_popp")
	outname<-paste(basename(var),"land_quants",sep = "_")
	write.table(t(as.data.frame(land_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#podu
	podu_quants<- cbind(t(as.data.frame(quant_podu_bart)), t(as.data.frame(quant_podu_cali)), t(as.data.frame(quant_podu_chum)), t(as.data.frame(quant_podu_cris)), t(as.data.frame(quant_podu_knul)), t(as.data.frame(quant_podu_land)), t(as.data.frame(quant_podu_popp)))
	colnames(podu_quants)<-c("podu_bart","podu_cali","podu_chum","podu_cris","podu_knul","podu_land","podu_popp")
	outname<-paste(basename(var),"podu_quants",sep = "_")
	write.table(t(as.data.frame(podu_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

	#popp
	popp_quants<- cbind(t(as.data.frame(quant_popp_bart)), t(as.data.frame(quant_popp_cali)), t(as.data.frame(quant_popp_chum)), t(as.data.frame(quant_popp_cris)), t(as.data.frame(quant_popp_knul)), t(as.data.frame(quant_popp_land)), t(as.data.frame(quant_popp_podu)))
	colnames(popp_quants)<-c("popp_bart","popp_cali","popp_chum","popp_cris","popp_knul","popp_land","popp_podu")
	outname<-paste(basename(var),"popp_quants",sep = "_")
	write.table(t(as.data.frame(popp_quants)), outname, col.names=T, row.names=T, quote=F, sep=" ")

}
