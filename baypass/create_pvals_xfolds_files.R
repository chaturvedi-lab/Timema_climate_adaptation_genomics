library(DataCombine)
pcvars<-c("pc1","pc2","pc3")


for (var in pcvars){
#read in the files
        bart<-read.table(paste(var,"bin_lgscaf_out.csv","bart_quants",sep="_"), header=T)
        cali<-read.table(paste(var,"bin_lgscaf_out.csv","cali_quants",sep="_"), header=T)
        chum<-read.table(paste(var,"bin_lgscaf_out.csv","chum_quants",sep="_"), header=T)
        cris<-read.table(paste(var,"bin_lgscaf_out.csv","cris_quants",sep="_"), header=T)
        knul<-read.table(paste(var,"bin_lgscaf_out.csv","knul_quants",sep="_"), header=T)
        land<-read.table(paste(var,"bin_lgscaf_out.csv","land_quants",sep="_"), header=T)
        podu<-read.table(paste(var,"bin_lgscaf_out.csv","podu_quants",sep="_"), header=T)
        popp<-read.table(paste(var,"bin_lgscaf_out.csv","popp_quants",sep="_"), header=T)

        ### xfolds
        bartx<-InsertRow(data.frame(bart$xfold), NewRow=0, RowNum=1)
        calix<-InsertRow(data.frame(cali$xfold), NewRow=0, RowNum=2)
        chumx<-InsertRow(data.frame(chum$xfold), NewRow=0, RowNum=3)
        crisx<-InsertRow(data.frame(cris$xfold), NewRow=0, RowNum=4)
        knulx<-InsertRow(data.frame(knul$xfold), NewRow=0, RowNum=5)
        landx<-InsertRow(data.frame(land$xfold), NewRow=0, RowNum=6)
        podux<-InsertRow(data.frame(podu$xfold), NewRow=0, RowNum=7)
        poppx<-InsertRow(data.frame(popp$xfold), NewRow=0, RowNum=8)

        #xfold<-cbind(bartx,calix,chumx,crisx,knulx,landx,podux,poppx[-9,])
        #rownames(xfold)<-c("bart","cali","chum","cris","knul","land","podu","popp")
        #colnames(xfold)<-c("bart","cali","chum","cris","knul","land","podu","popp")
	#based on phylogeny
	xfold<-cbind(bartx,podux,chumx,crisx,knulx,poppx[-9,],landx,calix)
	colnames(xfold)<-c("bart","podu","chum","cris","knul","popp","land","cali")
	rownames(xfold)<-c("bart","cali","chum","cris","knul","land","podu","popp")
	xfoldr<-rbind(xfold[1,],xfold[7,],xfold[3,],xfold[4,],xfold[5,],xfold[8,],xfold[6,],xfold[2,])

        ##pvalues
        ###pvalues 
        bartp<-InsertRow(data.frame(bart$p), NewRow=1, RowNum=1)
        calip<-InsertRow(data.frame(cali$p), NewRow=1, RowNum=2)
        chump<-InsertRow(data.frame(chum$p), NewRow=1, RowNum=3)
        crisp<-InsertRow(data.frame(cris$p), NewRow=1, RowNum=4)
        knulp<-InsertRow(data.frame(knul$p), NewRow=1, RowNum=5)
        landp<-InsertRow(data.frame(land$p), NewRow=1, RowNum=6)
        podup<-InsertRow(data.frame(podu$p), NewRow=1, RowNum=7)
        poppp<-InsertRow(data.frame(popp$p), NewRow=1, RowNum=8)

        #pval<-cbind(bartp,calip,chump,crisp,knulp,landp,podup,poppp[-9,])
        #rownames(pval)<-c("bart","cali","chum","cris","knul","land","podu","popp")
        #colnames(pval)<-c("bart","cali","chum","cris","knul","land","podu","popp")
	#based on phylogeny
        pval<-cbind(bartp,podup,chump,crisp,knulp,poppp[-9,],landp,calip)
        colnames(pval)<-c("bart","podu","chum","cris","knul","popp","land","cali")
        rownames(pval)<-c("bart","cali","chum","cris","knul","land","podu","popp")
        pvalr<-rbind(pval[1,],pval[7,],pval[3,],pval[4,],pval[5,],pval[8,],pval[6,],pval[2,])
	#write out the files
        write.table(xfoldr,file=paste("pairwise_xfolds0.9",var,sep="_"), col.names=T, row.names=T, sep=" ", quote=F)
        write.table(pvalr,file=paste("pairwise_pvals0.9",var,sep="_"), col.names=T, row.names=T, sep=" ", quote=F)
}



