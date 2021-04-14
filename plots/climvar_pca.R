#map timema

library(ggplot2)
library(ggmap)
library(ggpubr)
library(maps)
library(mapdata)
library(factoextra)
library(viridis)
library(rcartocolor)
library(gridExtra)



#color scheme
cols<-carto_pal(8,"Bold")

#read in file
datfile<-read.csv("Timema_climatedata.csv", header=T,as.is = TRUE)
head(datfile)
datfile<-data.frame(datfile,stringsAsFactors = FALSE)
datfile<-datfile[-c(1,19,45),]

##pca
pca.clim<-prcomp(datfile[,-c(1,2)],center = T,scale=T)
#rename clim variables for plot
cv<-c("Tann_mean","Dmean_range", "Isothermality", "Tseason", "Tmax","Tmin","Tann_range","Twet","Tdry","Twarm", "Tcold","Pannual","Pwet","Pdry","Pseason","Pwet_quart","Pdry_quart", "Pwarm", "Pcold", "Altitude","Latitude", "Longitude")
cv<-c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19","Elev","Lat","Long")
rownames(pca.clim$rotation)<-cv

#plots
scree.pc12<-fviz_pca_biplot(pca.clim,axes = c(1,2),title = "(A) PC1 vs. PC2", 
                            # Individuals
                            geom = "point",repel = TRUE,mean.point=FALSE,
                            fill.ind = datfile$species, col.ind = "black",
                            pointshape = 21, pointsize = 3.5,labelsize=4,
                            palette = cols,
                            # Variables
                            col.var = "black",arrowsize=0.3,
                            legend.title = list(fill = "Species")) + xlim(-7.5,7.5) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  xlab("PC1 (51.7%)") + ylab("PC2 (24.4%)")


scree.pc23<-fviz_pca_biplot(pca.clim,axes = c(2,3),title = "(B) PC2 vs. PC3", 
                            # Individuals
                            geom = "point",repel = TRUE,mean.point=FALSE,
                            fill.ind = datfile$species, col.ind = "black",
                            pointshape = 21, pointsize = 3.5,labelsize=4,
                            palette = cols,
                            # Variables
                            col.var = "grey18",arrowsize = 0.3,
                            legend.title = list(fill = "Species")) + xlim(-7.5,7.5) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  xlab("PC2 (24.4%)") + ylab("PC3 (16.1%)")

pdf("fviz_pca_climvars_23.pdf", width = 8, height = 12, bg = "white")
grid.arrange(scree.pc12, scree.pc23, ncol = 1)
dev.off()


scree.pc13<-fviz_pca_biplot(pca.clim,axes = c(1,3),title = "(A) PC1 vs. PC3", 
                            # Individuals
                            geom = "point",repel = TRUE,mean.point=FALSE,
                            fill.ind = datfile$species, col.ind = "black",
                            pointshape = 21, pointsize = 3.5,labelsize=6,
                            palette = cols,
                            # Variables
                            col.var = "grey18",arrowsize = 0.3,
                            legend.title = list(fill = "Species")) + xlim(-7.5,7.5) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  xlab("PC1 (51.7%)") + ylab("PC3 (16.1%)")
#save plots
pdf("fviz_pca_climvars_13.pdf", width = 10, height = 8, bg = "white")
print(scree.pc13)
dev.off()





#analyses of variance of PC scores withing and among species
vardat<-cbind(datfile[,c(1:2)], pca.clim$x[,1], pca.clim$x[,2], pca.clim$x[,3])
vardat<-vardat[-c(1,19,45),]
colnames(vardat)[3]<-"pc1"
colnames(vardat)[4]<-"pc2"
colnames(vardat)[5]<-"pc3"

#boxplots
spnames<-c(expression(italic("T. bartmani")),expression(italic("T. podura")),expression(italic("T. chumash")),expression(italic("T. cristinae")),expression(italic("T. knulli")),expression(italic("T. poppensis")),expression(italic("T. landelsensis")),expression(italic("T. californicum")))

#reorder in phylogeny order
vardat$species <- factor(vardat$species, levels=c("bart","podu","chum","cris","knul","popp","land","cali"))
#rearrange colors
colsr<-cols[c(1,7,3,4,5,8,6,2)]
pdf("pcscores12_perspecies_boxplots.pdf", width=13, height=12)
par(mfrow=c(2,1))
par(mar=c(9.5,5,4,5))
boxplot(vardat$pc1 ~ vardat$species, col=colsr, names=spnames,las=2, xlab=" ", ylab="PC1 scores", main ="(C) PC1", cex.main=2, cex.lab=2, cex.axis=1.5)
par(mar=c(9.5,5,4,5))
boxplot(vardat$pc2 ~ vardat$species, col=colsr, names=spnames, las=2, xlab=" ", ylab="PC2 scores", main ="(D) PC2", cex.main=2, cex.lab=2, cex.axis=1.5)
dev.off()

pdf("pcscores3_perspecies_boxplots.pdf", width=15, height=8)
par(mar=c(11,5,4,5))
boxplot(vardat$pc3 ~ vardat$species, col=colsr, names=spnames, las=2, xlab=" ", ylab="PC3 scores", main ="PC3", cex.main=2, cex.lab=2, cex.axis=1.5)
dev.off()

pdf("pcscores_perspecies_boxplots.pdf", width=10, height=12, bg="white")
par(mfrow=c(3,1))
par(mar=c(1,10,6,4))
#par(mar=c(9.5,5,4,5))
boxplot(vardat$pc1 ~ vardat$species, col=colsr, las=2, xlab=" ", ylab="", main ="", cex.axis=3, xaxt='n')
title(main = "(C) PC1", cex.main=5, adj=0)
mtext(text = "PC1 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 3)

par(mar=c(1,10,6,4))
#par(mar=c(9.5,5,4,5))
boxplot(vardat$pc2 ~ vardat$species, col=colsr, las=2, xlab=" ", ylab="", main ="",cex.axis=3, xaxt='n')
title(main = "(D) PC2", cex.main=5, adj=0)
mtext(text = "PC2 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 3)

par(mar=c(1,10,6,4))
#par(mar=c(9.5,5,4,5))
boxplot(vardat$pc3 ~ vardat$species, col=colsr, names=spnames, las=2, ylab="", xlab="", main ="",cex.axis=3)
title(main = "(E) PC3", cex.main=5, adj=0)
mtext(text = "PC3 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 3)

dev.off()

##empty plot for multipanel figure

pdf("figure1_outline.pdf",width = 15, height = 10)
par(mfrow=c(3,2))
plot(0, type="n",main="")
title(main="(A) Map of species ranges", adj=0, cex.main=3)
plot(0, type="n")
title(main="(B) PCA biplot of climate variables",adj=0, cex.main=3)
plot(0, type="n")
title(main="(C) PC1", cex.main=3,adj=0)
plot.new()
plot(0, type="n")
title(main="(D) PC2", cex.main=3,adj=0)
plot.new()
dev.off()

pdf("figure3_outline.pdf",width = 22, height = 10)
par(mfrow=c(2,2))
plot(0, type="n",main="")
title(main="(A) PC1 - Shared SNP windows", adj=0, cex.main=3)
plot(0, type="n")
title(main="(B) PC1 - Pairwise comparison ",adj=0, cex.main=3)
plot(0, type="n")
title(main="(C) PC1 - Multi-species comparison", cex.main=3,adj=0)
plot(0, type="n")
title(main="(D) PC1 - Field experiment comparison", cex.main=3,adj=0)
dev.off()

#manova summaries
res.man <- manova(cbind(vardat$pc1, vardat$pc2,vardat$pc3) ~ vardat$species)
summary(res.man)

#save the pca variables to run baypass
###############################################
#create output file for Zach for Bayesian model
pcdat<-cbind(as.character(datfile[,1]),as.character(datfile[,2]), pca.clim$x[,1], pca.clim$x[,2], pca.clim$x[,3])
colnames(pcdat)<-c("site","species","PC1","PC2","PC3")
pcdat<-as.data.frame(pcdat)
write.table(pcdat,file="timema_pops_pcvars.csv", sep=",",quote = F, row.names = F, col.names = T)

#subset by species to get means
bart_pc<-pcdat[which(pcdat[,2] == "bart"),]
cali_pc<-pcdat[which(pcdat[,2] == "cali"),]
chum_pc<-pcdat[which(pcdat[,2] == "chum"),]
cris_pc<-pcdat[which(pcdat[,2] == "cris"),]
knul_pc<-pcdat[which(pcdat[,2] == "knul"),]
land_pc<-pcdat[which(pcdat[,2] == "land"),]
podu_pc<-pcdat[which(pcdat[,2] == "podu"),]
popp_pc<-pcdat[which(pcdat[,2] == "popp"),]

#get means
bart_m1<-mean(as.numeric(bart_pc$PC1))
bart_m2<-mean(as.numeric(bart_pc$PC2))
bart_m3<-mean(as.numeric(bart_pc$PC3))
cali_m1<-mean(as.numeric(cali_pc$PC1))
cali_m2<-mean(as.numeric(cali_pc$PC2))
cali_m3<-mean(as.numeric(cali_pc$PC3))
chum_m1<-mean(as.numeric(chum_pc$PC1))
chum_m2<-mean(as.numeric(chum_pc$PC2))
chum_m3<-mean(as.numeric(chum_pc$PC3))
cris_m1<-mean(as.numeric(cris_pc$PC1))
cris_m2<-mean(as.numeric(cris_pc$PC2))
cris_m3<-mean(as.numeric(cris_pc$PC3))
knul_m1<-mean(as.numeric(knul_pc$PC1))
knul_m2<-mean(as.numeric(knul_pc$PC2))
knul_m3<-mean(as.numeric(knul_pc$PC3))
land_m1<-mean(as.numeric(land_pc$PC1))
land_m2<-mean(as.numeric(land_pc$PC2))
land_m3<-mean(as.numeric(land_pc$PC3))
podu_m1<-mean(as.numeric(podu_pc$PC1))
podu_m2<-mean(as.numeric(podu_pc$PC2))
podu_m3<-mean(as.numeric(podu_pc$PC3))
popp_m1<-mean(as.numeric(popp_pc$PC1))
popp_m2<-mean(as.numeric(popp_pc$PC2))
popp_m3<-mean(as.numeric(popp_pc$PC3))

#get pairwise distance
library(spaa)

m1<-c(bart_m1,cali_m1,chum_m1,cris_m1,land_m1,knul_m1,podu_m1,popp_m1)
m1_dist<-dist(m1,method = "euclidean",diag = F,upper = F)
m1_dist_list<-dist2list(m1_dist)
#get absolute differences
m1_diff<-outer(m1,m1,"-")
m1_diff_list<-m1_diff[lower.tri((m1_diff))]

m2<-c(bart_m2,cali_m2,chum_m2,cris_m2,land_m2,knul_m2,podu_m2,popp_m2)
m2_dist<-dist(m2,method = "euclidean",diag = F,upper = F)
m2_dist_list<-dist2list(m2_dist)
#get absolute differences
m2_diff<-outer(m2,m2,"-")
m2_diff_list<-m2_diff[lower.tri((m2_diff))]

m3<-c(bart_m3,cali_m3,chum_m3,cris_m3,land_m3,knul_m3,podu_m3,popp_m3)
m3_dist<-dist(m3,method = "euclidean",diag = F,upper = F)
m3_dist_list<-dist2list(m3_dist)
#get absolute differences
m3_diff<-outer(m3,m3,"-")
m3_diff_list<-m3_diff[lower.tri((m3_diff))]

clim_dist<-cbind(m1_dist_list,m2_dist_list[,3],m3_dist_list[,3])
clim_dist<-clim_dist[-c(1,9,10,17,18,19,25,26,27,28,33:37,41:46,49:55,57:64),]

pc1file<-read.table("~/Desktop/pc1_obs_exp_pdist.txt",sep=",", header=T)
pc2file<-read.table("~/Desktop/pc2_obs_exp_pdist.txt",sep=",", header=T)
pc3file<-read.table("~/Desktop/pc3_obs_exp_pdist.txt",sep=",", header=T)

final_dat<-cbind(rownames(pc1file),pc1file$pdist,pc1file$xfold,pc2file$xfold,pc3file$xfold,clim_dist[,c(3:5)],abs(m1_diff_list),abs(m2_diff_list),abs(m3_diff_list))
colnames(final_dat)<-c("species_pair","pairwise_dist","pc1_xfold","pc2_xfold","pc3_xfold","pc1_climdist","pc2_climdist","pc3_climdist","pc1_climdiff","pc2_climdiff","pc3_climdiff")
colnames(final_dat)[7]<-"pc2_climdist"
colnames(final_dat)[8]<-"pc3_climdist"

write.table(final_dat, "Timema_pdist_xfolds_climdist.csv", col.names = T, row.names = F, quote = F, sep=",")

##################################################################

####### EXTRA #####################
#correlation matrix
datfile2<-read.csv("Timema_climatedata.csv", header=T)
head(datfile2)
res<-cor(datfile2[,-c(1,2)])
round(res, 2)
cv<-c("Ann. mean temperature", "Mean diurnal range"," Isothermality", " Temp. seasonality", " Max. temperature","Min. temperature","Temperature ann. range","Mean temperature wet","Mean temperature dry","Mean temperature warm", "Mean temperature cold","Annual precipitation","Precipitation wet","Precipitation dry","Precipitation seasonal","Precipitation wet quarter","Precipitation dry quarter", "Precipitation warm quarter", "Precipitation cold quarter", "Altitude","Latitude", "Longitude")
rownames(res)<-cv
colnames(res)<-cv

datfile2<-data.frame(datfile2,stringsAsFactors = FALSE)
##pca
pca.clim2<-prcomp(datfile2[,-c(1,2)],center = T,scale=T)
pc_drop2<-cbind(datfile2[,2],pca.clim2$x[,1], pca.clim2$x[,2],pca.clim2$x[,3])
head(pc_drop2)
pc_drop2<-pc_drop2[-c(1,19,45),]
head(pc_drop2)
pc_drop2<-data.frame(pc_drop2)
#get the maximum rotation value for each PC
#pc1 = altitude (0.28)
#pc2 = meantempwarm (0.338) BIO10
#pc3 = precipwetquart (0.359) BIO16

library(viridis)
# Build a Pannel of 100 colors with Rcolor Brewer
my_colors<-plasma(10)
library(corrplot)
pdf("climatevar_corrplot.pdf", width = 15, height=15)
par(mar=c(3,3,3,4))
corrplot(res, order = "hclust", tl.col = "black", tl.srt = 45,is.corr=FALSE,method="color", type="lower",tl.cex=1, title="",mar=c(1,1,1,1),cl.cex=1, cl.pos="b", diag=FALSE)
dev.off()

pdf("climatevar_heatmap.pdf", width = 15, height=15)
#par(mar=c(7,5,4,7))
#col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(res,margins= c(15,15),cexRow=2,cexCol=2,trace="none",scale="row", key=F)

dev.off()



                #Df Pillai approx F num Df den Df    Pr(>F)    
#vardat$species  7 2.3975   19.328     21    102 < 2.2e-16 ***
  #Residuals      34  

#one way anova for effect of species on pcs
pc1.aov<-aov(vardat$pc1 ~ vardat$species)
pc2.aov<-aov(vardat$pc2 ~ vardat$species)
pc3.aov<-aov(vardat$pc3 ~ vardat$species)

#anova summaries
summary(pc1.aov)
#Df Sum Sq Mean Sq F value Pr(>F)
#vardat$species  7  465.3   66.47   99.58 <2e-16 ***
#Residuals      36   24.0    0.67

summary(pc2.aov)
#Df Sum Sq Mean Sq F value   Pr(>F)
#vardat$species  7 129.28   18.47   7.996 7.62e-06 ***
#Residuals      36  83.15    2.31

summary(pc3.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)
#vardat$species  7 119.66   17.09   29.98 3.43e-13 ***
#Residuals      36  20.53    0.57

#As the ANOVA test is significant, we can compute Tukey HSD (Tukey Honest Significant Differences, R function: TukeyHSD()) for performing multiple pairwise-comparison between the means of groups. The function TukeyHD() takes the fitted ANOVA as an argument.

TukeyHSD(pc1.aov)
TukeyHSD(pc2.aov)
TukeyHSD(pc3.aov)

#results
#pc1
Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = vardat$pc1 ~ vardat$species)

$`vardat$species`
                 diff         lwr        upr     p adj
cali-bart -7.50847624  -9.3660671 -5.6508854 0.0000000
chum-bart -4.09096523  -5.3875310 -2.7943994 0.0000000
cris-bart -8.77431068 -10.2910273 -7.2575941 0.0000000
knul-bart -9.75299299 -11.3437388 -8.1622472 0.0000000
land-bart -8.61966023 -10.3154009 -6.9239195 0.0000000
podu-bart -2.29691080  -3.9926515 -0.6011701 0.0024234
popp-bart -7.42670631  -9.2842972 -5.5691154 0.0000000
chum-cali  3.41751101   1.7348650  5.1001570 0.0000036
cris-cali -1.26583444  -3.1234253  0.5917564 0.3805791
knul-cali -2.24451675  -4.1630284 -0.3260051 0.0125404
land-cali -1.11118399  -3.1176115  0.8952435 0.6363757
podu-cali  5.21156544   3.2051380  7.2179929 0.0000000
popp-cali  0.08176993  -2.0631913  2.2267311 1.0000000
cris-chum -4.68334544  -5.9799113 -3.3867796 0.0000000
knul-chum -5.66202776  -7.0444625 -4.2795930 0.0000000
land-chum -4.52869499  -6.0307570 -3.0266330 0.0000000
podu-chum  1.79405443   0.2919924  3.2961164 0.0101573
popp-chum -3.33574108  -5.0183871 -1.6530951 0.0000058
knul-cris -0.97868232  -2.5694281  0.6120635 0.5096434
land-cris  0.15465045  -1.5410903  1.8503912 0.9999888
podu-cris  6.47739987   4.7816592  8.1731406 0.0000000
popp-cris  1.34760437  -0.5099865  3.2051952 0.3046379
land-knul  1.13333277  -0.6289327  2.8955982 0.4536035
podu-knul  7.45608219   5.6938168  9.2183476 0.0000000
popp-knul  2.32628668   0.4077751  4.2447983 0.0086732
podu-land  6.32274942   4.4651585  8.1803403 0.0000000
popp-land  1.19295392  -0.8134735  3.1993814 0.5521495
popp-podu -5.12979551  -7.1362230 -3.1233680 0.0000000

#pc2
Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = vardat$pc2 ~ vardat$species)

$`vardat$species`
                diff        lwr        upr     p adj
cali-bart  1.6625597 -1.7928409  5.1179604 0.7770607
chum-bart  3.4128070  1.0009980  5.8246160 0.0013828
cris-bart  2.6267415 -0.1945813  5.4480643 0.0834330
knul-bart  0.5466551 -2.4123733  3.5056834 0.9987479
land-bart  0.9126895 -2.2416453  4.0670243 0.9808885
podu-bart  2.8659519 -0.2883829  6.0202867 0.0978427
popp-bart -2.5459401 -6.0013408  0.9094606 0.2866610
chum-cali  1.7502473 -1.3797294  4.8802239 0.6253349
cris-cali  0.9641817 -2.4912189  4.4195824 0.9844491
knul-cali -1.1159047 -4.6846271  2.4528178 0.9706226
land-cali -0.7498702 -4.4821295  2.9823891 0.9978644
podu-cali  1.2033922 -2.5288671  4.9356514 0.9653362
popp-cali -4.2084998 -8.1984528 -0.2185468 0.0325690
cris-chum -0.7860655 -3.1978745  1.6257435 0.9632846
knul-chum -2.8661519 -5.4376902 -0.2946137 0.0200056
land-chum -2.5001175 -5.2941805  0.2939456 0.1076823
podu-chum -0.5468551 -3.3409181  2.2472079 0.9981906
popp-chum -5.9587471 -9.0887237 -2.8287704 0.0000124
knul-cris -2.0800864 -5.0391148  0.8789419 0.3422668
land-cris -1.7140520 -4.8683868  1.4402829 0.6576216
podu-cris  0.2392104 -2.9151244  3.3935452 0.9999968
popp-cris -5.1726816 -8.6280822 -1.7172809 0.0006389
land-knul  0.3660345 -2.9120464  3.6441154 0.9999554
podu-knul  2.3192968 -0.9587841  5.5973777 0.3344432
popp-knul -3.0925952 -6.6613176  0.4761273 0.1303981
podu-land  1.9532624 -1.5021383  5.4086630 0.6127834
popp-land -3.4586296 -7.1908889  0.2736297 0.0860717
popp-podu -5.4118920 -9.1441513 -1.6796327 0.0009959

#pc3
Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = vardat$pc3 ~ vardat$species)

$`vardat$species`
                diff        lwr        upr     p adj
cali-bart  3.7475497  2.0307203  5.4643792 0.0000008
chum-bart  3.5431008  2.3447840  4.7414176 0.0000000
cris-bart  2.5768529  1.1750675  3.9786383 0.0000235
knul-bart  0.1305659 -1.3396391  1.6007708 0.9999906
land-bart  0.5261861 -1.0410576  2.0934298 0.9570632
podu-bart  1.5059033 -0.0613404  3.0731470 0.0671310
popp-bart  5.0628417  3.3460122  6.7796711 0.0000000
chum-cali -0.2044489 -1.7595902  1.3506923 0.9998660
cris-cali -1.1706969 -2.8875263  0.5461326 0.3797457
knul-cali -3.6169839 -5.3901177 -1.8438500 0.0000033
land-cali -3.2213637 -5.0757514 -1.3669759 0.0000630
podu-cali -2.2416464 -4.0960342 -0.3872586 0.0089590
popp-cali  1.3152919 -0.6671320  3.2977158 0.4141095
cris-chum -0.9662479 -2.1645647  0.2320689 0.1917621
knul-chum -3.4125349 -4.6902139 -2.1348560 0.0000000
land-chum -3.0169147 -4.4051560 -1.6286734 0.0000009
podu-chum -2.0371975 -3.4254387 -0.6489562 0.0008450
popp-chum  1.5197409 -0.0354004  3.0748821 0.0594199
knul-cris -2.4462870 -3.9164920 -0.9760821 0.0001285
land-cris -2.0506668 -3.6179105 -0.4834231 0.0036908
podu-cris -1.0709496 -2.6381933  0.4962942 0.3771221
popp-cris  2.4859888  0.7691593  4.2028182 0.0010151
land-knul  0.3956202 -1.2331072  2.0243477 0.9931148
podu-knul  1.3753375 -0.2533900  3.0040649 0.1508388
popp-knul  4.9322758  3.1591419  6.7054097 0.0000000
podu-land  0.9797172 -0.7371122  2.6965467 0.6016706
popp-land  4.5366556  2.6822678  6.3910434 0.0000001
popp-podu  3.5569383  1.7025506  5.4113261 0.0000108

#plot anova residuals
pdf("pc_anova_res_fitvals.pdf", width = 15, height = 4)
par(mfrow=c(1,3))
par(mar=c(5,5,4,4))
plot(pc1.aov, 1, main="PC1", cex.lab=1.5, cex.main=2)
par(mar=c(5,5,4,4))
plot(pc2.aov, 1, main="PC2", cex.lab=1.5, cex.main=2)
par(mar=c(5,5,4,4))
plot(pc3.aov, 1, main="PC3", cex.lab=1.5, cex.main=2)
dev.off()

