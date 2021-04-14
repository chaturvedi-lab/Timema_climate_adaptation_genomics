## updated on Jan 4th 2020
#plot bayesfactor scatterfplots for distribution of SNPs along the genome

#read in the matrix for pc1,2,3
#redoing this with new files with lg information
#folder: /uufs/chpc.utah.edu/common/home/gompert-group1/projects/timema_adaptation/analyses/baypass/pcvars_run/outfiles/bfmeans

#read the files
pc1<-read.table("pc1_bin_lgscaf_out.csv", header=T, sep=",")
pc2<-read.table("pc2_bin_lgscaf_out.csv", header=T, sep=",")
pc3<-read.table("pc3_bin_lgscaf_out.csv", header=T, sep=",")
species<-c("med_bart","med_podu","med_chum","med_cris","med_knul","med_popp","med_land","med_cali")

####### pc1 ##########
pdf("species_manhattanplot_pc1.pdf",width=25,height=20)
par(mfrow=c(4,2))
for (sps in species){
    dat<-cbind(pc1$lg,pc1$scaffold,pc1$bin,pc1[sps])
    datc<-dat[complete.cases(dat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc$lg == 0),]
    datco = data.frame(datc[order(datc[,1]),])
    lgo<-datco$lg
    bnds<-c(1,which(lgo[2:length(lgo)] != lgo[1:length(lgo)-1]),length(lgo))
    mids<-tapply(X=1:length(lgo),INDEX=lgo,mean)
    
    #plot
    par(mar=c(6,10,6,4))
    plot(datco$medbf, xaxt = 'n',yaxt = 'n', xlab = " ", ylab = "",cex.lab =3,cex=1.5, ylim=c(-20,20), type="n")
    mtext(text = expression("log"[10]("Bayes Factor")),side = 2,line = 5, cex =3)
    x = which(datco$medbf > 1)
    xtop = datco[x,]
    for(j in seq(1,13,2)){
        polygon(x=c(bnds[j],bnds[j],bnds[j+1],bnds[j+1]),y=c(-20,20,20,-20),col=alpha("gray",0.5),border=NA)
    }
    points(datco$medbf, col= "black",pch=16,cex=2)
    points(x=x, y = xtop$medbf, col="deeppink3", pch=16,cex=2)
    abline (h = 1, col = "black", lwd=4, lty = "dashed")
    axis(1,at=mids,c(1:13), cex.axis=3, las=3, font = 2)
    #axis(1,at=mids,c(1,"",3,"",5,"",7,"",9,"",11,"",13), cex.axis=4, las=2, font = 2)
    axis(2,at=c(-20,1,20), cex.axis=3)
    nx<-length(x)
    nxtop<-length(unique(xtop[,1]))
    title(main=paste(sps,nx,nxtop, sep=","),adj=0,cex.main=5)
    
}
dev.off()

#### pc2 #########
pdf("species_manhattanplot_pc2.pdf",width=25,height=20)
par(mfrow=c(4,2))
for (sps in species){
    dat<-cbind(pc2$lg,pc2$scaffold,pc2$bin,pc2[sps])
    datc<-dat[complete.cases(dat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc$lg == 0),]
    datco = data.frame(datc[order(datc[,1]),])
    lgo<-datco$lg
    bnds<-c(1,which(lgo[2:length(lgo)] != lgo[1:length(lgo)-1]),length(lgo))
    mids<-tapply(X=1:length(lgo),INDEX=lgo,mean)
    
    #plot
    par(mar=c(6,10,6,4))
    plot(datco$medbf, xaxt = 'n',yaxt = 'n', xlab = " ", ylab = "",cex.lab =3,cex=1.5, ylim=c(-20,20), type="n")
    mtext(text = expression("log"[10]("Bayes Factor")),side = 2,line = 5, cex =3)
    x = which(datco$medbf > 1)
    xtop = datco[x,]
    for(j in seq(1,13,2)){
        polygon(x=c(bnds[j],bnds[j],bnds[j+1],bnds[j+1]),y=c(-20,20,20,-20),col=alpha("gray",0.5),border=NA)
    }
    points(datco$medbf, col= "black",pch=16,cex=2)
    points(x=x, y = xtop$medbf, col="deeppink3", pch=16,cex=2)
    abline (h = 1, col = "black", lwd=4, lty = "dashed")
    axis(1,at=mids,c(1:13), cex.axis=3, las=3, font = 2)
    #axis(1,at=mids,c(1,"",3,"",5,"",7,"",9,"",11,"",13), cex.axis=4, las=2, font = 2)
    axis(2,at=c(-20,1,20), cex.axis=3)
    nx<-length(x)
    nxtop<-length(unique(xtop[,1]))
    title(main=paste(sps,nx,nxtop, sep=","),adj=0,cex.main=5)
    
}
dev.off()


### pc3 ###
####### pc1 ##########
pdf("species_manhattanplot_pc3.pdf",width=25,height=20)
par(mfrow=c(4,2))
for (sps in species){
    dat<-cbind(pc3$lg,pc3$scaffold,pc3$bin,pc3[sps])
    datc<-dat[complete.cases(dat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc$lg == 0),]
    datco = data.frame(datc[order(datc[,1]),])
    lgo<-datco$lg
    bnds<-c(1,which(lgo[2:length(lgo)] != lgo[1:length(lgo)-1]),length(lgo))
    mids<-tapply(X=1:length(lgo),INDEX=lgo,mean)
    
    #plot
    par(mar=c(6,10,6,4))
    plot(datco$medbf, xaxt = 'n',yaxt = 'n', xlab = " ", ylab = "",cex.lab =3,cex=1.5, ylim=c(-20,20), type="n")
    mtext(text = expression("log"[10]("Bayes Factor")),side = 2,line = 5, cex =3)
    x = which(datco$medbf > 1)
    xtop = datco[x,]
    for(j in seq(1,13,2)){
        polygon(x=c(bnds[j],bnds[j],bnds[j+1],bnds[j+1]),y=c(-20,20,20,-20),col=alpha("gray",0.5),border=NA)
    }
    points(datco$medbf, col= "black",pch=16,cex=2)
    points(x=x, y = xtop$medbf, col="deeppink3", pch=16,cex=2)
    abline (h = 1, col = "black", lwd=4, lty = "dashed")
    axis(1,at=mids,c(1:13), cex.axis=3, las=3, font = 2)
    #axis(1,at=mids,c(1,"",3,"",5,"",7,"",9,"",11,"",13), cex.axis=4, las=2, font = 2)
    axis(2,at=c(-20,1,20), cex.axis=3)
    nx<-length(x)
    nxtop<-length(unique(xtop[,1]))
    title(main=paste(sps,nx,nxtop, sep=","),adj=0,cex.main=5)
    
}
dev.off()

## get basic stats 
med_species<-c("med_bart","med_cali","med_chum","med_cris","med_knul","med_land","med_podu","med_popp")
getTop1<-function(species){
    datmat<-cbind(pc1[,1],pc1[,2],pc1[,3],pc1[[species]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    #bart1$lg[which(bart1$lg == 0)]<-14
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    #get number of SNPs in top quantile for each LG
    lg_nums<-data.frame(table(datco$lg))
    #get top number of SNPs
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    write.table(ntop_list, file=paste("pc1_topsnps",species,sep = "_"),sep=",",row.names=F, col.names=T, quote=F)
}

pc1_top<-c()
for (i in med_species){
    n<-getTop1(i)
    pc1_top[i]<-n
}
pc1_top

## SAM IS HERE
getToplg<-function(species){
    datmat<-cbind(pc1[,1],pc1[,2],pc1[,3],pc1[[species]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    #bart1$lg[which(bart1$lg == 0)]<-14
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    #get number of SNPs in top quantile for each LG
    lg_nums<-data.frame(table(datco$lg))
}

mat <- matrix(nrow = 8, ncol = 16)
for(column in 1:n.columns){
    mat[, column] <- vector
}

#med_bart med_cali med_chum med_cris med_knul med_land med_podu med_popp
#177      386      181      976      443      380      245      361

#pc2
getTop2<-function(species){
    datmat<-cbind(pc2[,1],pc2[,2],pc2[,3],pc2[[species]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    #bart1$lg[which(bart1$lg == 0)]<-14
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    write.table(ntop_list, file=paste("pc2_topsnps",species,sep = "_"),sep=",",row.names=F, col.names=T, quote=F)
}

for (i in med_species){
    n<-getTop2(i)
}

#pc3
getTop<-function(species){
    datmat<-cbind(pc3[,1],pc3[,2],pc3[,3],pc3[[species]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    #bart1$lg[which(bart1$lg == 0)]<-14
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    #write.table(ntop_list, file=paste("pc3_topsnps",species,sep = "_"),sep=",",row.names=F, col.names=T, quote=F)
}
pc3_top<-c()
for (i in med_species){
    n<-getTop(i)
    pc3_top[i]<-n
}

## LG size versus SNP associations
lg_scaf_len<-read.table("lg_scaf_length.txt", header = F)
#get length sum for each LG
lg_sumlen<-aggregate(lg_scaf_len[,3] ~ lg_scaf_len[,1], lg_scaf_len, sum)
lg_sumlen<-lg_sumlen[-1,]

#did pc1 separately for each species as podura had only 12 LGs
pdf("lg_topsnps_pc1.pdf", width=12, height = 20)
par(mfrow=c(4,2))
for (sps in species){
    datmat<-cbind(pc1[,1],pc1[,2],pc1[,3],pc1[[sps]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    ntop_snps<-as.data.frame(table(ntop_list[,1]))
    par(mar=c(5,10,4,4))
    plot(lg_sumlen[,2], ntop_snps[,2], xlab = " ", ylab = "No. of SNP bins", cex.lab=3, cex.axis=2, pch =16, cex=3, col=rainbow(13))
    mtext("Linkage group",side=1, line=3.5, cex=2)
    abline(lm(ntop_snps[,2] ~ lg_sumlen[,2]),col='black', lwd=2)
    legend("topleft",legend=c(1:13),fill=rainbow(13), cex = 2)
    title(sps,cex.main=4, adj=0)
}
dev.off()
#for podu
par(mar=c(5,10,4,4))
plot(lg_sumlen[-12,][,2], ntop_snps[,2], xlab = " ", ylab = "No. of SNP bins", cex.lab=3, cex.axis=2, pch =16, cex=3, col=rainbow(13))
mtext("Linkage group",side=1, line=3.5, cex=2)
abline(lm(ntop_snps[,2] ~ lg_sumlen[-12,][,2]),col='black', lwd=2)
legend("topleft",legend=c(1:11,13),fill=rainbow(13), cex = 2)
title(sps,cex.main=4, adj=0)

pdf("lg_topsnps_pc2.pdf", width=12, height = 20)
par(mfrow=c(4,2))
for (sps in species){
    datmat<-cbind(pc2[,1],pc2[,2],pc2[,3],pc2[[sps]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    ntop_snps<-as.data.frame(table(ntop_list[,1]))
    par(mar=c(5,10,4,4))
    plot(lg_sumlen[,2], ntop_snps[,2], xlab = " ", ylab = "No. of SNP bins", cex.lab=3, cex.axis=2, pch =16, cex=3, col=rainbow(13))
    mtext("Linkage group",side=1, line=3.5, cex=2)
    abline(lm(ntop_snps[,2] ~ lg_sumlen[,2]),col='black', lwd=2)
    legend("topleft",legend=c(1:13),fill=rainbow(13), cex = 2)
    title(sps,cex.main=4, adj=0)
}
dev.off()

pdf("lg_topsnps_pc3.pdf", width=12, height = 20)
par(mfrow=c(4,2))
for (sps in species){
    datmat<-cbind(pc3[,1],pc3[,2],pc3[,3],pc3[[sps]])
    datc<-datmat[complete.cases(datmat),]
    colnames(datc)<-c("lg","scaffold","bin","medbf")
    datc<-as.data.frame(datc)
    datc<-datc[-which(datc[,1] == 0),]
    datco = datc[order(datc[,1]),]
    datco = data.frame(datco)
    ntop<-length(which(datco$medbf > quantile(datco$medbf,probs=c(0.9))))
    ntop_list<-datco[which(datco$medbf > quantile(datco$medbf, probs=0.9, na.rm=T)),]
    ntop_snps<-as.data.frame(table(ntop_list[,1]))
    par(mar=c(5,10,4,4))
    plot(lg_sumlen[,2], ntop_snps[,2], xlab = " ", ylab = "No. of SNP bins", cex.lab=3, cex.axis=2, pch =16, cex=3, col=rainbow(13))
    mtext("Linkage group",side=1, line=3.5, cex=2)
    abline(lm(ntop_snps[,2] ~ lg_sumlen[,2]),col='black', lwd=2)
    legend("topleft",legend=c(1:13),fill=rainbow(13), cex = 2)
    title(sps,cex.main=4, adj=0)
}
dev.off()

######## EXTRA #############
########## this is all the plots for main text ###########
######### extra stuff here ################################

#getting simple summary stats
#get number of SNPs across scaffols
tail(data.frame(sort(table(bart$scaffold))),10)
tail(data.frame(sort(table(cali$scaffold))),10)
tail(data.frame(sort(table(chum$scaffold))),10)
tail(data.frame(sort(table(cris$scaffold))),10)
tail(data.frame(sort(table(knul$scaffold))),10)
tail(data.frame(sort(table(land$scaffold))),10)
tail(data.frame(sort(table(podu$scaffold))),10)
tail(data.frame(sort(table(popp$scaffold))),10)

#get the top snps
head(bart[order(bart$medbf, decreasing=T),],10)
head(cali[order(cali$medbf, decreasing=T),],10)
head(chum[order(chum$medbf, decreasing=T),],10)
head(cris[order(cris$medbf, decreasing=T),],10)
head(land[order(land$medbf, decreasing=T),],10)
head(knul[order(knul$medbf, decreasing=T),],10)
head(podu[order(podu$medbf, decreasing=T),],10)
head(popp[order(popp$medbf, decreasing=T),],10)

##top 0.9 SNPs
tail(data.frame(sort(table(bart_pc1_top$scaffold))),10)
tail(data.frame(sort(table(cali_pc1_top$scaffold))),10)
tail(data.frame(sort(table(chum_pc1_top$scaffold))),10)
tail(data.frame(sort(table(cris_pc1_top$scaffold))),10)
tail(data.frame(sort(table(land_pc1_top$scaffold))),10)
tail(data.frame(sort(table(knul_pc1_top$scaffold))),10)
tail(data.frame(sort(table(podu_pc1_top$scaffold))),10)
tail(data.frame(sort(table(popp_pc1_top$scaffold))),10)

#number of unique scaffolds
nrow(data.frame(sort(table(bart_pc1_top$scaffold))))
nrow(data.frame(sort(table(cali_pc1_top$scaffold))))
nrow(data.frame(sort(table(chum_pc1_top$scaffold))))
nrow(data.frame(sort(table(cris_pc1_top$scaffold))))
nrow(data.frame(sort(table(knul_pc1_top$scaffold))))
nrow(data.frame(sort(table(land_pc1_top$scaffold))))
nrow(data.frame(sort(table(podu_pc1_top$scaffold))))
nrow(data.frame(sort(table(popp_pc1_top$scaffold))))

####### pc2 ##########
#getting simple summary stats
#get number of SNPs across scaffols
tail(data.frame(sort(table(bart$scaffold))),10)
tail(data.frame(sort(table(cali$scaffold))),10)
tail(data.frame(sort(table(chum$scaffold))),10)
tail(data.frame(sort(table(cris$scaffold))),10)
tail(data.frame(sort(table(knul$scaffold))),10)
tail(data.frame(sort(table(land$scaffold))),10)
tail(data.frame(sort(table(podu$scaffold))),10)
tail(data.frame(sort(table(popp$scaffold))),10)

##top 0.9 SNPs
tail(data.frame(sort(table(bart_pc2_top$scaffold))),10)
tail(data.frame(sort(table(cali_pc2_top$scaffold))),10)
tail(data.frame(sort(table(chum_pc2_top$scaffold))),10)
tail(data.frame(sort(table(cris_pc2_top$scaffold))),10)
tail(data.frame(sort(table(land_pc2_top$scaffold))),10)
tail(data.frame(sort(table(knul_pc2_top$scaffold))),10)
tail(data.frame(sort(table(podu_pc2_top$scaffold))),10)
tail(data.frame(sort(table(popp_pc2_top$scaffold))),10)

#number of unique scaffolds
nrow(data.frame(sort(table(bart_pc2_top$scaffold))))
nrow(data.frame(sort(table(cali_pc2_top$scaffold))))
nrow(data.frame(sort(table(chum_pc2_top$scaffold))))
nrow(data.frame(sort(table(cris_pc2_top$scaffold))))
nrow(data.frame(sort(table(knul_pc2_top$scaffold))))
nrow(data.frame(sort(table(land_pc2_top$scaffold))))
nrow(data.frame(sort(table(podu_pc2_top$scaffold))))
nrow(data.frame(sort(table(popp_pc2_top$scaffold))))

##### pc3 ##########
#getting simple summary stats
#get number of SNPs across scaffols
tail(data.frame(sort(table(bart$scaffold))),10)
tail(data.frame(sort(table(cali$scaffold))),10)
tail(data.frame(sort(table(chum$scaffold))),10)
tail(data.frame(sort(table(cris$scaffold))),10)
tail(data.frame(sort(table(knul$scaffold))),10)
tail(data.frame(sort(table(land$scaffold))),10)
tail(data.frame(sort(table(podu$scaffold))),10)
tail(data.frame(sort(table(popp$scaffold))),10)

##top 0.9 SNPs
tail(data.frame(sort(table(bart_pc3_top$scaffold))),10)
tail(data.frame(sort(table(cali_pc3_top$scaffold))),10)
tail(data.frame(sort(table(chum_pc3_top$scaffold))),10)
tail(data.frame(sort(table(cris_pc3_top$scaffold))),10)
tail(data.frame(sort(table(land_pc3_top$scaffold))),10)
tail(data.frame(sort(table(knul_pc3_top$scaffold))),10)
tail(data.frame(sort(table(podu_pc3_top$scaffold))),10)
tail(data.frame(sort(table(popp_pc3_top$scaffold))),10)

#number of unique scaffolds
nrow(data.frame(sort(table(bart_pc3_top$scaffold))))
nrow(data.frame(sort(table(cali_pc3_top$scaffold))))
nrow(data.frame(sort(table(chum_pc3_top$scaffold))))
nrow(data.frame(sort(table(cris_pc3_top$scaffold))))
nrow(data.frame(sort(table(knul_pc3_top$scaffold))))
nrow(data.frame(sort(table(land_pc3_top$scaffold))))
nrow(data.frame(sort(table(podu_pc3_top$scaffold))))
nrow(data.frame(sort(table(popp_pc3_top$scaffold))))


#get top snps and write them out in annotation folder for annotation
##pc1 read the code from above for pc1 then subset the snps as follows and write them out
species<-c("bart","cali","chum","cris","knul","land","podu","popp")

pcvars<-c("pc1","pc2","pc3")
for (var in pcvars){
    pc<-read.table(paste(var,"bin_out.csv",sep="_"), header=T)
    dat<-cbind(pc$scaffold,pc$bin,pc$med_popp)
    popp<-ppdat[complete.cases(ppdat),]
    colnames(popp)<-c("scaffold","bin","medbf")
    popp<-as.data.frame(popp)
    pc1_top<-i[which(i$medbf > quantile(i$medbf, probs=0.9, na.rm=T)),]
    write.table(i, file=paste(i,"topsnps0.9", sep="_"),sep=",",row.names=F, col.names=T, quote=F)
}
