## read in the tree distances
trdist<-read.csv("species_pairwise_distances_edit.csv", header=T, sep=",")
pc1file<-read.table("pc1_obs_exp_pdist.txt",sep=",", header=T)
pc2file<-read.table("pc2_obs_exp_pdist.txt",sep=",", header=T)
pc3file<-read.table("pc3_obs_exp_pdist.txt",sep=",", header=T)

library("ggpubr")
par(mfrow=c(1,3))
#correlations
x1<-ggscatter(pc1file, x = "pdist", y = "xfold",title = "(A) PC1",
          size = 3,
          color = "#0072B2", 
          main.plot.size = 3,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "#0072B2", fill = "lightgray"),
          cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
          xlab = "Pairwise distances ", ylab = "X-fold", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text") +  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 25), axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) 

x2<-ggscatter(pc2file, x = "pdist", y = "xfold",title = "(B) PC2",
              size = 3,
              color = "#0072B2", 
              main.plot.size = 3,
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "#0072B2", fill = "lightgray"),
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
              xlab = "Pairwise distances ", ylab = "X-fold", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text") +  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 25), axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

x3<-ggscatter(pc3file, x = "pdist", y = "xfold",title = "(B) PC3",
              size = 3,
              color = "#0072B2", 
              main.plot.size = 3,
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "#0072B2", fill = "lightgray"),
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
              xlab = "Pairwise distances ", ylab = "X-fold", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text") +  theme(plot.title = element_text(color="black", size=40, face="bold"),axis.title = element_text(size = 30), axis.text = element_text(size = 30),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

#pvalues
p1<-ggscatter(pc1file, x = "pdist", y = "p",title = "(A) PC1",
              size = 3,
              color = "orangered", 
              main.plot.size = 3,
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "orangered", fill = "lightgray"),
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
              xlab = "Pairwise distances ", ylab = "P-value", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text", size=50) +  theme(plot.title = element_text(color="black", size=23, face="bold"),axis.title = element_text(size = 20), axis.text = element_text(size = 12),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

p2<-ggscatter(pc2file, x = "pdist", y = "p",title = "(B) PC2",
              size = 3,
              color = "orangered", 
              main.plot.size = 3,
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "orangered", fill = "lightgray"),
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
              xlab = "Pairwise distances ", ylab = "P-value", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text", size=50) +  theme(plot.title = element_text(color="black", size=23, face="bold"),axis.title = element_text(size = 20), axis.text = element_text(size = 12),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

p3<-ggscatter(pc1file, x = "pdist", y = "p",title = "(B) PC3",
              size = 3,
              color = "orangered", 
              main.plot.size = 3,
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "orangered", fill = "lightgray"),
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman"),
              xlab = "Pairwise distances ", ylab = "P-value", ggtheme = theme_bw()) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
  font("legend.text", size=50) +  theme(plot.title = element_text(color="black", size=23, face="bold"),axis.title = element_text(size = 20), axis.text = element_text(size = 12),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

pdf("pairdist_xfolds_spearmancor_pc3.pdf", width = 9, height = 9)
print(x3)
dev.off()

pdf("pairdist_xfolds_spearmancor_pc12.pdf", width = 6, height = 11)
ggarrange(x1,x2, ncol = 1)
dev.off()

pdf("pairdist_pvals_spearmancor_plot.pdf", width = 6, height = 15)
ggarrange(p1,p2,p3,ncol =1)
dev.off()

### Zach ecology vs genetics plots
library(plotrix)
genes_e<-c(-0.77,-0.17,-0.55)
genes_ul<-c(-0.28,0.22,-0.21)
genes_ll<-c(-1.20,-0.57,-0.87)

eco_e<-c(0.29,-0.10,-0.47)
eco_ul<-c(0.77,0.30,-0.14)
eco_ll<-c(-0.20,-0.50,-0.80)

pdf("eco_gen_barplots_pc3.pdf",width=10,height=8)
par(mar=c(5,8,6,8))
plotCI(x=c(0.4,0.6),y=c(genes_e[3],eco_e[3]),li=c(genes_ll[3],eco_ll[3]),ui=c(genes_ul[3],eco_ul[3]), col=c("#E7B800","#00AFBB"), xaxt="n",yaxt="n",xlab="", ylab="", cex.lab=3,ylim = c(-1,0),xlim=c(0.3,0.7), pch=20,lwd=4, cex=6, cex.axis=3)
mtext("Estimate",side=2, line=4, cex = 3)
axis(1,at=c(0.4,0.6),c("Genetics","Ecology"), cex.axis=3)
axis(2,at=c(-1,-0.5,0), cex.axis=3)
dev.off()

pdf("eco_gen_barplots_pc12.pdf",width=10,height=15)
par(mfrow=c(2,1))
par(mar=c(5,8,6,8))
plotCI(x=c(0.4,0.6),y=c(genes_e[1],eco_e[1]),li=c(genes_ll[1],eco_ll[1]),ui=c(genes_ul[1],eco_ul[1]), col=c("#E7B800","#00AFBB"), xaxt="n",yaxt="n",xlab="", ylab="", cex.lab=3,ylim = c(-1.5,1),xlim=c(0.3,0.7), pch=20,lwd=4, cex=6, cex.axis=3)
mtext("Estimate",side=2, line=4, cex = 3)
axis(1,at=c(0.4,0.6),c("Genetics","Ecology"), cex.axis=3)
axis(2,at=c(-1.5,0,1), cex.axis=3)
title("(A) PC1", adj=0, cex.main=5)

par(mar=c(5,8,6,8))
plotCI(x=c(0.4,0.6),y=c(genes_e[2],eco_e[2]),li=c(genes_ll[2],eco_ll[2]),ui=c(genes_ul[2],eco_ul[2]), col=c("#E7B800","#00AFBB"), xaxt="n",yaxt="n",xlab="", ylab="", cex.lab=3,ylim = c(-0.6,0.5),xlim=c(0.3,0.7), pch=20,lwd=4, cex=6, cex.axis=3)
mtext("Estimate",side=2, line=4, cex = 3)
axis(1,at=c(0.4,0.6),c("Genetics","Ecology"), cex.axis=3)
axis(2,at=c(-0.6,0,0.5), cex.axis=3)
title("(B) PC2", adj=0, cex.main=5)

dev.off()


##make it into a list
library(reshape2)

pcvar<-c("pc1","pc2","pc3")
quants<-c("0.90",0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99)
library(RColorBrewer)
cols<-brewer.pal("Set1",n=3)
#did this on the desktop

pdf("phylo_xfold_pval.pdf", width = 18, height = 11)
par(mfrow=c(2,3))
xfold<-read.table("../baypass/median/quantfiles/pairwise_xfolds0.9_pc1",sep="")
xfold_l<-melt(as.matrix(xfold))
#drop same species line
xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
xfold_log<-log(xfold_l[,3])
xfold_log[which(xfold_log == "-Inf")] <- -1000000
par(mar=c(6,5,5,2))
plot(pcfile$pdist,log(pcfile$xfold),pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
myline.fit <- lm(log(pcfile$xfold) ~ pcfile$pdist)
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit, lwd=2)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n',cex=2)
title(main="(A) PC1", cex.main=2, adj=0 )

#pc2
xfold<-read.table("../baypass/median/quantfiles/pairwise_xfolds0.9_pc2",sep="")
xfold_l<-melt(as.matrix(xfold))
#drop same species line
xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
xfold_log<-log(xfold_l[,3])
xfold_log[which(xfold_log == "-Inf")] <- -1000000
par(mar=c(6,5,5,2))
plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
myline.fit <- lm(xfold_log ~ trdist[,3])
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit, lwd=2)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n',cex=2)
title(main="(B) PC2", cex.main=2, adj=0 )

#pc3
xfold<-read.table("../baypass/median/quantfiles/pairwise_xfolds0.9_pc3",sep="")
xfold_l<-melt(as.matrix(xfold))
#drop same species line
xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
xfold_log<-log(xfold_l[,3])
xfold_log[which(xfold_log == "-Inf")] <- -1000000
par(mar=c(6,5,5,2))
plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
myline.fit <- lm(xfold_log ~ trdist[,3])
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit, lwd=2)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n',cex=2)
title(main="(C) PC3", cex.main=2, adj=0 )

############# pvals ##############################
#pc1
pval<-read.table("../baypass/median/quantfiles/pairwise_pvals0.9_pc1",sep="")
pval_l<-melt(as.matrix(pval))
pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
par(mar=c(6,5,5,2))
plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
myline.fit <- lm(pval_l[,3] ~ trdist[,3])
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n',cex=2)
title(main="(A) PC1", cex.main=2, adj=0 )

#pc2
pval<-read.table("../baypass/median/quantfiles/pairwise_pvals0.9_pc2",sep="")
pval_l<-melt(as.matrix(pval))
pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
par(mar=c(6,5,5,2))
plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
myline.fit <- lm(pval_l[,3] ~ trdist[,3])
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n',cex=2)
title(main="(B) PC2", cex.main=2, adj=0 )

#pc3
pval<-read.table("../baypass/median/quantfiles/pairwise_pvals0.9_pc3",sep="")
pval_l<-melt(as.matrix(pval))
pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
par(mar=c(6,5,5,2))
plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
myline.fit <- lm(pval_l[,3] ~ trdist[,3])
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot
abline(myline.fit)
r2 = summary(myline.fit)$r.squared
my.p = summary(myline.fit)$coefficients[,4][2]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n',cex=2)
title(main="(C) PC3", cex.main=2, adj=0 )
dev.off()

##### trying just simple linear model ########
pdf("pc1_pairdist_logxfolds_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="../baypass/median/quantfiles", pattern="pc1", all.files=T, full.names=T)
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc1",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  myline.fit <- lm(xfold_log ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit, lwd=2)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC1",q,sep=" - "), cex.main=2)
}
dev.off()

#pc2
pdf("pc2_pairdist_logxfolds_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc2", all.files=T, full.names=T)
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc2",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  myline.fit <- lm(xfold_log ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit, lwd=2)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC2",q,sep=" - "), cex.main=2)
}
dev.off()

#pc3
pdf("pc3_pairdist_logxfolds_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc3", all.files=T, full.names=T)
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc3",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  myline.fit <- lm(xfold_log ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit, lwd=2)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC3",q,sep=" - "), cex.main=2)
}
dev.off()

pdf("pc1_pairdist_pvals_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc1", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc1",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  myline.fit <- lm(pval_l[,3] ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC1",q,sep="-"), cex.main=2)
}
dev.off()

#pc2
pdf("pc2_pairdist_pvals_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc2", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc2",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  myline.fit <- lm(pval_l[,3] ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC2",q,sep="-"), cex.main=2)
}
dev.off()

#pc3
pdf("pc3_pairdist_pvals_median_lm.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc3", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc3",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  myline.fit <- lm(pval_l[,3] ~ trdist[,3])
  # get information about the fit
  summary(myline.fit)
  # draw the fit line on the plot
  abline(myline.fit)
  r2 = summary(myline.fit)$r.squared
  my.p = summary(myline.fit)$coefficients[,4][2]
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC3",q,sep="-"), cex.main=2)
}
dev.off()

## plot pairwise correlations ###
library("ggpubr")
library("ggplot2")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(reshape)
pcvar<-c("pc1","pc2","pc3")
quants<-c("0.9",0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99)

######## Spearmans correlation (non-parametric) ###########
pdf("pc1_pairdist_logxfolds_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/analyses/baypass/pcvars_run/outfiles/bfmeans/quants/median", pattern="pc1", all.files=T, full.names=T)
  xfold<-read.table(paste("/uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/analyses/baypass/pcvars_run/outfiles/bfmeans/quants/median/pairwise_xfolds",q,"_","pc1",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  cr<-cor.test(trdist[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC1",q,sep=" - "), cex.main=2)
}
dev.off()

#pc2
pdf("pc2_pairdist_logxfolds_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc2", all.files=T, full.names=T)
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc2",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  cr<-cor.test(trdist[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC2",q,sep=" - "), cex.main=2)
}
dev.off()

#pc3
pdf("pc3_pairdist_logxfolds_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc3", all.files=T, full.names=T)
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc3",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  par(mar=c(6,5,5,2))
  plot(trdist[,3],xfold_log,pch=21, ylab="X-fold", xlab="Pairwise distance", cex.lab=2, cex=2, bg="#31a354", col="black")
  cr<-cor.test(trdist[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC3",q,sep=" - "), cex.main=2)
}
dev.off()

###### pvalues ##############
pdf("pc1_pairdist_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc1", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc1",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  # get information about the fit
  cr<-cor.test(trdist[,3],pval_l[,3], method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC1",q,sep=" - "), cex.main=2)
}
dev.off()

pdf("pc2_pairdist_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc2", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc2",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  # get information about the fit
  cr<-cor.test(trdist[,3],pval_l[,3], method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC2",q,sep=" - "), cex.main=2)
}
dev.off()

pdf("pc3_pairdist_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc3", all.files=T, full.names=T)
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc3",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  par(mar=c(6,5,5,2))
  plot(trdist[,3],pval_l[,3],pch=21, xlab="Pairwise distance", ylab="P-value", cex.lab=2, cex=2, bg="#de2d26", col="black")
  # get information about the fit
  cr<-cor.test(trdist[,3],pval_l[,3], method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC3",q,sep=" - "), cex.main=2)
}
dev.off()

##### correlation of xfolds and pvalues #################
pdf("pc1_logxfolds_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc1", all.files=T, full.names=T)
  #xfold
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc1",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  #pval
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc1",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  #plot
  par(mar=c(6,5,5,2))
  plot(pval_l[,3],xfold_log,pch=21, xlab="X-fold", ylab="P-value", cex.lab=2, cex=2, bg="#2980B9", col="black")
  cr<-cor.test(pval_l[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC1",q,sep=" - "), cex.main=2)
}
dev.off()

#pc2
pdf("pc2_logxfolds_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc2", all.files=T, full.names=T)
  #xfold
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc2",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  #pval
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc2",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  #plot
  par(mar=c(6,5,5,2))
  plot(pval_l[,3],xfold_log,pch=21, xlab="X-fold", ylab="P-value", cex.lab=2, cex=2, bg="#2980B9", col="black")
  cr<-cor.test(pval_l[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC2",q,sep=" - "), cex.main=2)
}
dev.off()

#pc3
pdf("pc3_logxfolds_pvals_median_spearmanscor.pdf", width = 22, height=10)
#par(mfrow=c(5,2))
par(mfrow=c(2,5))
for (q in quants){
  files<-list.files(path="/Volumes/biksu/timema/baypass/median/quantfiles", pattern="pc3", all.files=T, full.names=T)
  #xfold
  xfold<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_xfolds",q,"_","pc3",sep=""))
  xfold_l<-melt(as.matrix(xfold))
  #drop same species line
  xfold_l<-xfold_l[-c(1,10,19,28,37,46,55,64),]
  xfold_log<-log(xfold_l[,3])
  xfold_log[which(xfold_log == "-Inf")] <- -1000000
  #pval
  pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals",q,"_","pc3",sep=""))
  pval_l<-melt(as.matrix(pval))
  pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
  #plot
  par(mar=c(6,5,5,2))
  plot(pval_l[,3],xfold_log,pch=21, xlab="X-fold", ylab="P-value", cex.lab=2, cex=2, bg="#2980B9", col="black")
  cr<-cor.test(pval_l[,3],xfold_log, method="spearman")
  pcor = cr$estimate
  my.p = cr$p.value
  mylabel = bquote(italic(R) == .(format(pcor, digits = 2)))
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R) == MYVALUE),
                     list(MYVALUE = format(pcor,dig=2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('bottomleft', legend = rp, bty = 'n',cex=2)
  title(main=paste("PC3",q,sep=" - "), cex.main=2)
}
dev.off()


## beta regression for pvalues
library (betareg)
pval<-read.table(paste("/Volumes/biksu/timema/baypass/median/quantfiles/pairwise_pvals0.9_pc1",sep=""))
pval_l<-melt(as.matrix(pval))
pval_l<-pval_l[-c(1,10,19,28,37,46,55,64),]
betaMod <- betareg(pval_l[,3] ~ trdist[,3]) # train model. Tune var names.
summary (betaMod) # model summary
predict (betaMod, testData)
