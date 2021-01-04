plot bayesfactor scatterfplots for distribution of SNPs along the genome

#read in the matrix for pc1,2,3

pc1<-read.table("pc1_bin_out.csv", header=T, sep=",")
pc2<-read.table("pc2_bin_out.csv", header=T, sep=",")
pc3<-read.table("pc3_bin_out.csv", header=T, sep=",")


####### pc1 ##########
#bart 1771
bdat<-cbind(pc1$scaffold,pc1$bin,pc1$med_bart)
bart<-bdat[complete.cases(bdat),]
colnames(bart)<-c("scaffold","bin","medbf")
bart<-as.data.frame(bart)
bart_pc1_top<-bart[which(bart$medbf > quantile(bart$medbf, probs=0.9, na.rm=T)),]
write.table(bart_pc1_top, file="../../../../annotation/pc1_bart_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cali 3852
cadat<-cbind(pc1$scaffold,pc1$bin,pc1$med_cali)
cali<-cadat[complete.cases(cadat),]
colnames(cali)<-c("scaffold","bin","medbf")
cali<-as.data.frame(cali)
cali_pc1_top<-cali[which(cali$medbf > quantile(cali$medbf, probs=0.9, na.rm=T)),]
write.table(cali_pc1_top, file="../../../../annotation/pc1_cali_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#chum 1806
chdat<-cbind(pc1$scaffold,pc1$bin,pc1$med_chum)
chum<-chdat[complete.cases(chdat),]
colnames(chum)<-c("scaffold","bin","medbf")
chum<-as.data.frame(chum)
chum_pc1_top<-chum[which(chum$medbf > quantile(chum$medbf, probs=0.9, na.rm=T)),]
write.table(chum_pc1_top, file="../../../../annotation/pc1_chum_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cris 9452
crdat<-cbind(pc1$scaffold,pc1$bin,pc1$med_cris)
cris<-crdat[complete.cases(crdat),]
colnames(cris)<-c("scaffold","bin","medbf")
cris<-as.data.frame(cris)
cris_pc1_top<-cris[which(cris$medbf > quantile(cris$medbf, probs=0.9, na.rm=T)),]
write.table(cris_pc1_top, file="../../../../annotation/pc1_cris_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#knul 4367
kdat<-cbind(pc1$scaffold,pc1$bin,pc1$med_knul)
knul<-kdat[complete.cases(kdat),]
colnames(knul)<-c("scaffold","bin","medbf")
knul<-as.data.frame(knul)
knul_pc1_top<-knul[which(knul$medbf > quantile(knul$medbf, probs=0.9, na.rm=T)),]
write.table(knul_pc1_top, file="../../../../annotation/pc1_knul_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#land 3799
ldat<-cbind(pc1$scaffold,pc1$bin,pc1$med_land)
land<-ldat[complete.cases(ldat),]
colnames(land)<-c("scaffold","bin","medbf")
land<-as.data.frame(land)
land_pc1_top<-land[which(land$medbf > quantile(land$medbf, probs=0.9, na.rm=T)),]
write.table(land_pc1_top, file="../../../../annotation/pc1_land_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#podu 2443
pddat<-cbind(pc1$scaffold,pc1$bin,pc1$med_podu)
podu<-pddat[complete.cases(pddat),]
colnames(podu)<-c("scaffold","bin","medbf")
podu<-as.data.frame(podu)
podu_pc1_top<-podu[which(podu$medbf > quantile(podu$medbf, probs=0.9, na.rm=T)),]
write.table(podu_pc1_top, file="../../../../annotation/pc1_podu_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#popp 3609
ppdat<-cbind(pc1$scaffold,pc1$bin,pc1$med_popp)
popp<-ppdat[complete.cases(ppdat),]
colnames(popp)<-c("scaffold","bin","medbf")
popp<-as.data.frame(popp)
popp_pc1_top<-popp[which(popp$medbf > quantile(popp$medbf, probs=0.9, na.rm=T)),]
write.table(popp_pc1_top, file="../../../../annotation/pc1_popp_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

library(RColorBrewer)
cols<-brewer.pal(8, name="Dark2")

pdf("species_bfmed_scatterplot_pc1.pdf", height=10, width=12)
par(mfrow=c(4,2))
#bart
par(mar=c(5,5,3,5))
plot(bart$scaffold,bart$medbf,col="black",bg=cols[4],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. bartmani")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cali
par(mar=c(5,5,3,5))
plot(cali$scaffold,cali$medbf,col="black",bg=cols[8],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. californicum")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#chum
par(mar=c(5,5,3,5))
plot(chum$scaffold,chum$medbf,col="black",bg=cols[1],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. chumash")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cris
par(mar=c(5,5,3,5))
plot(cris$scaffold,cris$medbf,col="black",bg=cols[6],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. cristinae")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#knul
par(mar=c(5,5,3,5))
plot(knul$scaffold,knul$medbf,col="black",bg=cols[3],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. knulli")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#land
par(mar=c(5,5,3,5))
plot(land$scaffold,land$medbf,col="black",bg=cols[2],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. landelsensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#podu
par(mar=c(5,5,3,5))
plot(podu$scaffold,podu$medbf,col="black",bg=cols[5],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1,main = expression(italic("T. podura")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#popp
par(mar=c(5,5,3,5))
plot(popp$scaffold,popp$medbf,col="black",bg=cols[7],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. poppensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)

dev.off()

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
#bart 1771
bdat<-cbind(pc2$scaffold,pc2$bin,pc2$med_bart)
bart<-bdat[complete.cases(bdat),]
colnames(bart)<-c("scaffold","bin","medbf")
bart<-as.data.frame(bart)
bart_pc2_top<-bart[which(bart$medbf > quantile(bart$medbf, probs=0.9, na.rm=T)),]
write.table(bart_pc2_top, file="../../../../annotation/pc2_bart_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cali 3852
cadat<-cbind(pc2$scaffold,pc2$bin,pc2$med_cali)
cali<-cadat[complete.cases(cadat),]
colnames(cali)<-c("scaffold","bin","medbf")
cali<-as.data.frame(cali)
cali_pc2_top<-cali[which(cali$medbf > quantile(cali$medbf, probs=0.9, na.rm=T)),]
write.table(cali_pc2_top, file="../../../../annotation/pc2_cali_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#chum 1806
chdat<-cbind(pc2$scaffold,pc2$bin,pc2$med_chum)
chum<-chdat[complete.cases(chdat),]
colnames(chum)<-c("scaffold","bin","medbf")
chum<-as.data.frame(chum)
chum_pc2_top<-chum[which(chum$medbf > quantile(chum$medbf, probs=0.9, na.rm=T)),]
write.table(chum_pc2_top, file="../../../../annotation/pc2_chum_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cris 9452
crdat<-cbind(pc2$scaffold,pc2$bin,pc2$med_cris)
cris<-crdat[complete.cases(crdat),]
colnames(cris)<-c("scaffold","bin","medbf")
cris<-as.data.frame(cris)
cris_pc2_top<-cris[which(cris$medbf > quantile(cris$medbf, probs=0.9, na.rm=T)),]
write.table(cris_pc2_top, file="../../../../annotation/pc2_cris_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#knul 4367
kdat<-cbind(pc2$scaffold,pc2$bin,pc2$med_knul)
knul<-kdat[complete.cases(kdat),]
colnames(knul)<-c("scaffold","bin","medbf")
knul<-as.data.frame(knul)
knul_pc2_top<-knul[which(knul$medbf > quantile(knul$medbf, probs=0.9, na.rm=T)),]
write.table(knul_pc2_top, file="../../../../annotation/pc2_knul_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#land 3799
ldat<-cbind(pc2$scaffold,pc2$bin,pc2$med_land)
land<-ldat[complete.cases(ldat),]
colnames(land)<-c("scaffold","bin","medbf")
land<-as.data.frame(land)
land_pc2_top<-land[which(land$medbf > quantile(land$medbf, probs=0.9, na.rm=T)),]
write.table(land_pc2_top, file="../../../../annotation/pc2_land_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#podu 2443
pddat<-cbind(pc2$scaffold,pc2$bin,pc2$med_podu)
podu<-pddat[complete.cases(pddat),]
colnames(podu)<-c("scaffold","bin","medbf")
podu<-as.data.frame(podu)
podu_pc2_top<-podu[which(podu$medbf > quantile(podu$medbf, probs=0.9, na.rm=T)),]
write.table(podu_pc2_top, file="../../../../annotation/pc2_podu_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#popp 3609
ppdat<-cbind(pc2$scaffold,pc2$bin,pc2$med_popp)
popp<-ppdat[complete.cases(ppdat),]
colnames(popp)<-c("scaffold","bin","medbf")
popp<-as.data.frame(popp)
popp_pc2_top<-popp[which(popp$medbf > quantile(popp$medbf, probs=0.9, na.rm=T)),]
write.table(popp_pc2_top, file="../../../../annotation/pc2_popp_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

pdf("species_bfmed_scatterplot_pc2.pdf", height=10, width=12)
par(mfrow=c(4,2))
#bart
par(mar=c(5,5,3,5))
plot(bart$scaffold,bart$medbf,col="black",bg=cols[4],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. bartmani")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cali
par(mar=c(5,5,3,5))
plot(cali$scaffold,cali$medbf,col="black",bg=cols[8],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. californicum")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#chum
par(mar=c(5,5,3,5))
plot(chum$scaffold,chum$medbf,col="black",bg=cols[1],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. chumash")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cris
par(mar=c(5,5,3,5))
plot(cris$scaffold,cris$medbf,col="black",bg=cols[6],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. cristinae")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#knul
par(mar=c(5,5,3,5))
plot(knul$scaffold,knul$medbf,col="black",bg=cols[3],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. knulli")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#land
par(mar=c(5,5,3,5))
plot(land$scaffold,land$medbf,col="black",bg=cols[2],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. landelsensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#podu
par(mar=c(5,5,3,5))
plot(podu$scaffold,podu$medbf,col="black",bg=cols[5],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1,main = expression(italic("T. podura")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#popp
par(mar=c(5,5,3,5))
plot(popp$scaffold,popp$medbf,col="black",bg=cols[7],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. poppensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)

dev.off()

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


####### pc3 ###########
#bart 1771
bdat<-cbind(pc3$scaffold,pc3$bin,pc3$med_bart)
bart<-bdat[complete.cases(bdat),]
colnames(bart)<-c("scaffold","bin","medbf")
bart<-as.data.frame(bart)
bart_pc3_top<-bart[which(bart$medbf > quantile(bart$medbf, probs=0.9, na.rm=T)),]
write.table(bart_pc3_top, file="../../../../annotation/pc3_bart_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cali 3852
cadat<-cbind(pc3$scaffold,pc3$bin,pc3$med_cali)
cali<-cadat[complete.cases(cadat),]
colnames(cali)<-c("scaffold","bin","medbf")
cali<-as.data.frame(cali)
cali_pc3_top<-cali[which(cali$medbf > quantile(cali$medbf, probs=0.9, na.rm=T)),]
write.table(cali_pc3_top, file="../../../../annotation/pc3_cali_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#chum 1806
chdat<-cbind(pc3$scaffold,pc3$bin,pc3$med_chum)
chum<-chdat[complete.cases(chdat),]
colnames(chum)<-c("scaffold","bin","medbf")
chum<-as.data.frame(chum)
chum_pc3_top<-chum[which(chum$medbf > quantile(chum$medbf, probs=0.9, na.rm=T)),]
write.table(chum_pc3_top, file="../../../../annotation/pc3_chum_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#cris 9452
crdat<-cbind(pc3$scaffold,pc3$bin,pc3$med_cris)
cris<-crdat[complete.cases(crdat),]
colnames(cris)<-c("scaffold","bin","medbf")
cris<-as.data.frame(cris)
cris_pc3_top<-cris[which(cris$medbf > quantile(cris$medbf, probs=0.9, na.rm=T)),]
write.table(cris_pc3_top, file="../../../../annotation/pc3_cris_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#knul 4367
kdat<-cbind(pc3$scaffold,pc3$bin,pc3$med_knul)
knul<-kdat[complete.cases(kdat),]
colnames(knul)<-c("scaffold","bin","medbf")
knul<-as.data.frame(knul)
knul_pc3_top<-knul[which(knul$medbf > quantile(knul$medbf, probs=0.9, na.rm=T)),]
write.table(knul_pc3_top, file="../../../../annotation/pc3_knul_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#land 3799
ldat<-cbind(pc3$scaffold,pc3$bin,pc3$med_land)
land<-ldat[complete.cases(ldat),]
colnames(land)<-c("scaffold","bin","medbf")
land<-as.data.frame(land)
land_pc3_top<-land[which(land$medbf > quantile(land$medbf, probs=0.9, na.rm=T)),]
write.table(land_pc3_top, file="../../../../annotation/pc3_land_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#podu 2443
pddat<-cbind(pc3$scaffold,pc3$bin,pc3$med_podu)
podu<-pddat[complete.cases(pddat),]
colnames(podu)<-c("scaffold","bin","medbf")
podu<-as.data.frame(podu)
podu_pc3_top<-podu[which(podu$medbf > quantile(podu$medbf, probs=0.9, na.rm=T)),]
write.table(podu_pc3_top, file="../../../../annotation/pc3_podu_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)

#popp 3609
ppdat<-cbind(pc3$scaffold,pc3$bin,pc3$med_popp)
popp<-ppdat[complete.cases(ppdat),]
colnames(popp)<-c("scaffold","bin","medbf")
popp<-as.data.frame(popp)
popp_pc3_top<-popp[which(popp$medbf > quantile(popp$medbf, probs=0.9, na.rm=T)),]
write.table(popp_pc3_top, file="../../../../annotation/pc3_popp_topsnps.txt",sep=",",row.names=F, col.names=T, quote=F)


pdf("species_bfmed_scatterplot_pc3.pdf", height=10, width=12)
par(mfrow=c(4,2))
#bart
par(mar=c(5,5,3,5))
plot(bart$scaffold,bart$medbf,col="black",bg=cols[4],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. bartmani")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cali
par(mar=c(5,5,3,5))
plot(cali$scaffold,cali$medbf,col="black",bg=cols[8],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. californicum")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#chum
par(mar=c(5,5,3,5))
plot(chum$scaffold,chum$medbf,col="black",bg=cols[1],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. chumash")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#cris
par(mar=c(5,5,3,5))
plot(cris$scaffold,cris$medbf,col="black",bg=cols[6],pch=21,xlab="",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. cristinae")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#knul
par(mar=c(5,5,3,5))
plot(knul$scaffold,knul$medbf,col="black",bg=cols[3],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. knulli")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#land
par(mar=c(5,5,3,5))
plot(land$scaffold,land$medbf,col="black",bg=cols[2],pch=21,xlab="",ylab="log10(Bayesfactors)",cex.lab=1.8,cex.axis=1, cex=1.5, main = expression(italic("T. landelsensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#podu
par(mar=c(5,5,3,5))
plot(podu$scaffold,podu$medbf,col="black",bg=cols[5],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1,main = expression(italic("T. podura")), cex.main = 2)
abline(h=2, col="black", lwd=2)

#popp
par(mar=c(5,5,3,5))
plot(popp$scaffold,popp$medbf,col="black",bg=cols[7],pch=21,xlab="Physical distance(bp)",ylab="log10(Bayesfactors)",cex=1.5,cex.lab=1.8,cex.axis=1, main = expression(italic("T. poppensis")), cex.main = 2)
abline(h=2, col="black", lwd=2)


dev.off()

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
