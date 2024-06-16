#PACKAGES----

library("dplyr")
library("ggpubr")
library("patchwork")
library("AICcmodavg")
library("MASS")
library("corrplot")
library('lmerTest')
library("scales")
library("lme4")
library("arm")
library('glmmADMB')
library("plotrix")
library(("R2admb"))


#DATA EXPLORATION----

##Sorting and filtering camera data----
rawdata<-read.csv("00. Raw Data/Raw.Image.Data.csv")

names(rawdata)[1]<-"Image_ID"

nofalse<-rawdata[rawdata$Pres_Abs==1,]

uniqueanimID<-nofalse %>% distinct(Animal_ID,.keep_all = TRUE)

onlyspecies<-filter(uniqueanimID,Taxinomic_Rank=="Species")

specieslist<-onlyspecies %>% distinct(Identification, .keep_all=TRUE) %>% 
  dplyr::select(Identification)

mammals<-filter(onlyspecies,Bird_Mammal=="Mammal")

#adding site data
rawsite<-read.csv("00. Raw Data/Site.Data.csv")

names(rawsite)[1]<-"Site"
names(rawsite)[12]<-"Understory"

rawsite$Understory <- factor(rawsite$Understory,levels = c("Sparse","Moderate","Dense"))
rawsite$Midstory <- factor(rawsite$Midstory,levels = c("Open","Sparse","Moderate","Dense"))
rawsite$Canopy_Cover <- factor(rawsite$Canopy_Cover,levels = c("Open","Sparse","Moderate","Dense"))
rawsite$Fire_habitat_catergory <- factor(rawsite$Fire_habitat_catergory,levels = c("Unburnt_Rainforest","Burnt_Rainforest","Burnt_Sclerophyll"))


image_site<-merge(rawsite,mammals,by="Site",all.x=T,all.y=F)

##Creating species richness data frame----

MAD<-data.frame(mammals)

us_ind<-unlist(gregexpr("_",MAD$Identification))

MAD$species<-NA
MAD$species[which(us_ind>0)]<-substr(MAD$Identification,us_ind+1,nchar(MAD$Identification))[which(us_ind>1)]

speciesdata<-aggregate(species~Site,data=MAD,FUN = function(x)length(unique(x)))

modelrich<-merge(rawsite,speciesdata,by="Site",all.x=T,all.y=F) %>% 
  arrange(desc(Fire_habitat_catergory))

##Summary stats ----

culled<-image_site %>% dplyr::select(-Latitude,-Longitude,-Vertical_Cam,-Horizontal_Cam,To_Slope_CM,-Height_CM,-Confidence,-Photo_Series,-Trigger_Group,-Trigger,-File,-Cam_ID, -Pres_Abs,-No_Animals,-Animal_Image,-Image_ID, -Cam_Type,-Old_ID,-Old_Rank,-Date,-Time,-Temp,-Taxinomic_Rank,-Bird_Mammal) %>% 
  na.omit()

##Species richness at each site and orientation----

richness<-data.frame(culled)

richness1<-aggregate(data=richness,
                     Identification~Site+Orientation,
                     function(x)length(unique(x)))

ggboxplot(richness1,x= "Orientation" ,y="Identification",ylab="Species Richness")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

summary(richness1$Identification)

richness2<-filter(richness1,Orientation=="Vertical")
summary(richness2$Identification)

##Comparing abundances between sites and orientation----

Abundances<-culled %>% group_by(Site,Orientation) %>% summarise(abundance=n())
Abundances1<-culled %>% group_by(Site,Orientation,Identification) %>% summarise(abun=n())


##Correlation----

modelrich$Understory <- as.numeric(modelrich$Understory)
modelrich$Midstory <- as.numeric(modelrich$Midstory)
modelrich$Canopy_Cover <- as.numeric(modelrich$Canopy_Cover)


cordata1<-modelrich %>% dplyr::select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Understory,Midstory,Canopy_Cover,Fire_habitat_catergory)
cordata1$Fire_habitat_catergory <- as.numeric(cordata1$Fire_habitat_catergory)


cr1<-cor(cordata1, method="pearson")

corrplot::corrplot(cr1,method="color",  
         type="upper", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE) 
#Fire habitat category correlated with canopy cover so this was excluded from modelling 


#PCA----

##Preparing data----

pcadata<-modelrich

speciescount<-MAD

count<-MAD %>% group_by(Site) %>% count(Identification)
count<-as.data.frame(count)

speciescount$antechinus<-ifelse(speciescount$Identification=="Antechinus_stuartii",1,0)
speciescount$melomys<-ifelse(speciescount$Identification=="Melomys_cervinipes",1,0)
speciescount$bandicoot<-ifelse(speciescount$Identification=="Perameles_nasuta",1,0)
speciescount$koala<-ifelse(speciescount$Identification=="Phascolarctos_cinereus",1,0)
speciescount$bushrat<-ifelse(speciescount$Identification=="Rattus_fuscipes",1,0)
speciescount$pig<-ifelse(speciescount$Identification=="Sus_scrofa",1,0)
speciescount$stigmatica<-ifelse(speciescount$Identification=="Thylogale_stigmatica",1,0)
speciescount$thetis<-ifelse(speciescount$Identification=="Thylogale_thetis",1,0)
speciescount$possum<-ifelse(speciescount$Identification=="Trichosurus_caninus",1,0)

antechinus<-aggregate(speciescount$antechinus,by=list(Site=speciescount$Site),FUN=sum)
melomys<-aggregate(speciescount$melomys,by=list(Site=speciescount$Site),FUN=sum)
bandicoot<-aggregate(speciescount$bandicoot,by=list(Site=speciescount$Site),FUN=sum)
koala<-aggregate(speciescount$koala,by=list(Site=speciescount$Site),FUN=sum)
bushrat<-aggregate(speciescount$bushrat,by=list(Site=speciescount$Site),FUN=sum)
pig<-aggregate(speciescount$pig,by=list(Site=speciescount$Site),FUN=sum)
stigmatica<-aggregate(speciescount$stigmatica,by=list(Site=speciescount$Site),FUN=sum)
thetis<-aggregate(speciescount$thetis,by=list(Site=speciescount$Site),FUN=sum)
possum<-aggregate(speciescount$possum,by=list(Site=speciescount$Site),FUN=sum)

pcadata<-modelrich

pcadata<-full_join(pcadata,antechinus)
colnames(pcadata)[24]<-c("antechinus")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,melomys)
colnames(pcadata)[25]<-c("melomys")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,bandicoot)
colnames(pcadata)[26]<-c("bandicoot")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,koala)
colnames(pcadata)[27]<-c("koala")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,bushrat)
colnames(pcadata)[28]<-c("bushrat")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,pig)
colnames(pcadata)[29]<-c("pig")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,stigmatica)
colnames(pcadata)[30]<-c("stigmatica")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,thetis)
colnames(pcadata)[31]<-c("thetis")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,possum)
colnames(pcadata)[32]<-c("possum")
pcadata[is.na(pcadata)] = 0


pcadata<- rename(pcadata,A.sti=antechinus,M.cer=melomys,P.nas=bandicoot,P.cin=koala,R.fus=bushrat,S.scr=pig,T.sti=stigmatica,T.the=thetis,T.can=possum)

pca2<-prcomp(~.,data=pcadata[,which(colnames(pcadata)=="A.sti"):which(colnames(pcadata)=="T.can"),])
 
PoV2<-summary(pca2)$importance[2,]

pcadata$pca.comp1<-pca2$x[,1]
pcadata$pca.comp2<-pca2$x[,2]


shapes<-c(15,17,16)
shapes<-shapes[as.factor(pcadata$Fire_habitat_catergory)]
col.1<-c("grey20","grey40","grey60")
col.1<-col.1[as.factor(pcadata$Fire_habitat_catergory)]


dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0))
plot(x=1:length(PoV2),y=PoV2,ylab="Propotion Varience Explained",xlab="Components",type="p")
lines(x=1:length(PoV2),y=PoV2)
mtext("(a)",3,0.7,F,0)
mtext("58.17",3,-1,F,2,cex=0.7)
mtext("37.60",3,-4.1,F,3,cex=0.7)

plot(pcadata$pca.comp1,pcadata$pca.comp2,pch=shapes, xlab="",ylab="",cex=2,col=alpha(col.1,1))
  legend("topleft",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),pch=c(15,17,16),pt.cex=2,col=c("grey20","grey40","grey60"))
  mtext("(b)",3,0.4,F,adj=0,at=-140)
  title(ylab="PC2",adj=0.1,cex=1.2,line=2)
  mtext(expression(italic("Rattus fuscipes")),side=2,adj=1,line=3,font=2,cex=0.7)
  mtext(expression(italic("Melomys crevinipes")),side=2,adj=0.4,line=2.5,font=2,cex=0.7)
  par(xpd=T)
  arrows(-150,0,-150,45,length=0.1)
  par(xpd=NA)
  title(xlab="PC1",adj=0.75,cex=1.2,line=1.9)
  mtext(expression(italic("Thylogale stigmatica")),side=1,adj=0.3,line=2.6,font=2,cex=0.70)
  par(xpd=T)
  arrows(-35,-32,-110,-32,length=0.1)
  par(xpd=NA)

  dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
  biplot(pca2,xlab="Component 1",ylab="Component 2",col=c("grey40","black"),var.axes=T,arrow.len=0.1,ylim=c(-0.3,0.6))

##comp 1 modeling----
  
  head(pcadata,5);dim(pcadata)
  
  summary(comp1null<-lm(pca.comp1~1,data=pcadata))
  
  summary(comp1a<-lm(pca.comp1~Fire_habitat_catergory+Elevation,data=pcadata))
  summary(comp1b<-lm(pca.comp1~Fire_habitat_catergory*Elevation,data=pcadata))
  AICc(comp1a);AICc(comp1b)
  
  summary(comp1c<-lm(pca.comp1~Fire_habitat_catergory+Dis_to_road,data=pcadata))
  summary(comp1d<-lm(pca.comp1~Fire_habitat_catergory*Dis_to_road,data=pcadata))
  AICc(comp1c);AICc(comp1d)
  
  summary(comp1e<-lm(pca.comp1~Fire_habitat_catergory+Dis_to_path,data=pcadata))
  summary(comp1f<-lm(pca.comp1~Fire_habitat_catergory*Dis_to_path,data=pcadata))
  AICc(comp1e);AICc(comp1f)
  
  summary(comp1g<-lm(pca.comp1~Fire_habitat_catergory+Dis_to_rainforest_boundry,data=pcadata))
  summary(comp1h<-lm(pca.comp1~Fire_habitat_catergory*Dis_to_rainforest_boundry,data=pcadata))
  AICc(comp1g);AICc(comp1h)
  
  summary(comp1i<-lm(pca.comp1~Fire_habitat_catergory+Understory,data=pcadata))
  summary(comp1j<-lm(pca.comp1~Fire_habitat_catergory*Understory,data=pcadata))
  AICc(comp1i);AICc(comp1j)
  
  summary(comp1k<-lm(pca.comp1~Fire_habitat_catergory+Midstory,data=pcadata))
  summary(comp1l<-lm(pca.comp1~Fire_habitat_catergory*Midstory,data=pcadata))
  AICc(comp1k);AICc(comp1l)
  
  summary(comp1m<-lm(pca.comp1~Fire_habitat_catergory,data=pcadata))
  
  comp1final<-list("Null Model"=comp1null,"FHC + Elevation"=comp1a,"FHC + Distance To Road"=comp1c,"FHC + Distance To Path"=comp1e,"FHC + Distance To Rainforest Boundary"=comp1g,"FHC + Understory"=comp1i,"FHC + Midstory"=comp1k,"FHC"=comp1m)
  
  aictab(comp1final)
  
##comp2 modelling----
  
  summary(comp2null<-lm(pca.comp2~1,data=pcadata))
  
  summary(comp2fhc<-lm(pca.comp2~Fire_habitat_catergory,data=pcadata))
  
  summary(comp2a<-lm(pca.comp2~Fire_habitat_catergory+Elevation,data=pcadata))
  summary(comp2b<-lm(pca.comp2~Fire_habitat_catergory*Elevation,data=pcadata))
  AICc(comp2a);AICc(comp2b)
  
  summary(comp2c<-lm(pca.comp2~Fire_habitat_catergory+Dis_to_road,data=pcadata))
  summary(comp2d<-lm(pca.comp2~Fire_habitat_catergory*Dis_to_road,data=pcadata))
  AICc(comp2c);AICc(comp2d)
  
  summary(comp2e<-lm(pca.comp2~Fire_habitat_catergory+Dis_to_path,data=pcadata))
  summary(comp2f<-lm(pca.comp2~Fire_habitat_catergory*Dis_to_path,data=pcadata))
  AICc(comp2e);AICc(comp2f)
  
  summary(comp2g<-lm(pca.comp2~Fire_habitat_catergory+Dis_to_rainforest_boundry,data=pcadata))
  summary(comp2h<-lm(pca.comp2~Fire_habitat_catergory*Dis_to_rainforest_boundry,data=pcadata))
  AICc(comp2g);AICc(comp2h)
  
  summary(comp2i<-lm(pca.comp2~Fire_habitat_catergory+Understory,data=pcadata))
  summary(comp2j<-lm(pca.comp2~Fire_habitat_catergory*Understory,data=pcadata))
  AICc(comp2i);AICc(comp2j)
  
  summary(comp2k<-lm(pca.comp2~Fire_habitat_catergory+Midstory,data=pcadata))
  summary(comp2l<-lm(pca.comp2~Fire_habitat_catergory*Midstory,data=pcadata))
  AICc(comp2k);AICc(comp2l)
  
  comp2final<-list("Null"=comp2null,"FHC + Elevation"=comp2a,"FHC + Distance To Road"=comp2c,"FHC + Distance To Path"=comp2e,"FHC + Distance To Rainforest Boundary"=comp2g,"FHC + Level Of Underbrush"=comp2i,"FHC X Level Of Midstory"=comp2l,"FHC"=comp2fhc)
  
  aictab(comp2final)
  
###Component 2 predictions----

  predictions1<-data.frame(Fire_habitat_catergory=factor(c(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),4)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Midstory=c(0,0,0,1,1,1,2,2,2,3,3,3))
  
  pr1<-predict(object=comp2l,newdata=predictions1,se.fit =T,type="response")
  
  midstory<-data.frame(predictions1,fit=pr1$fit,se=pr1$se.fit)
  midstory$lci<-midstory$fit-(1.96*midstory$se)
  midstory$uci<-midstory$fit+(1.96*midstory$se)
  
  ###component 2 graphing----
  
  dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
  par(mar=c(4,4,1,1))
  
  dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
  plot(1:4,midstory$fit[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(midstory$lci[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"]),max(midstory$uci[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"])),xlim=c(min(0.5),max(4.5)),xlab="Midstory Index",ylab="Principal Component 2",las=1,cex=3,pch=15,xaxt="n",col="grey20")
  axis(side=1,at=1:4,labels=0:3,xlab="midstory")
  lines(1:4,midstory$fit[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"],col="grey20")
  arrows(x0=1:4, y0=midstory$lci[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=midstory$uci[midstory$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="grey20", lwd=2)
  
  points((1:4)-0.2,midstory$fit[midstory$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=17, col="grey40")
  lines((1:4)-0.2,midstory$fit[midstory$Fire_habitat_catergory=="Burnt_Rainforest"],col="grey40")
  arrows(x0=(1:4)-0.2, y0=midstory$lci[midstory$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=midstory$uci[midstory$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="grey40", lwd=2)
  
  points((1:4)+0.2,midstory$fit[midstory$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=16,col="grey60")
  lines((1:4)+0.2,midstory$fit[midstory$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='grey60')
  arrows((1:4)+0.2, y0=midstory$lci[midstory$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=midstory$uci[midstory$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="grey60", lwd=2)
  
  legend("bottomright",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=c("grey20","grey40","grey60"),pch=c(15,17,16),pt.cex=1.5,title="Fire Habitat Catergory",cex=0.9)

#RAREFACTION AND EXTRAPOLATION----
  head(singledata)
  singledata[,c(1,18:26)]
  
  iN_data <- list(
    UBR = t(singledata[singledata$Fire_habitat_catergory == "Unburnt_Rainforest", c(18:26)]),
    BR = t(singledata[singledata$Fire_habitat_catergory == "Burnt_Rainforest", c(18:26)]),
    BS = t(singledata[singledata$Fire_habitat_catergory == "Burnt_Sclerophyll", c(18:26)]))
  
  iN_results <- iNEXT(x = iN_data, q = c(0,1,2), datatype = 'abundance')
  
  str(iN_results)
  
  iN_results$iNextEst$size_based$Assemblage <- factor(iN_results$iNextEst$size_based$Assemblage, levels = c("UBR","BR","BS"))
  
  iN_results$iNextEst$coverage_based$Assemblage <- factor(iN_results$iNextEst$coverage_based$Assemblage, levels = c("UBR","BR","BS"))
  
  levels(iN_results$AsyEst$Assemblage)
  
  ##Sample-Size-Based R/E Curve----
  
  dev.new(height=7,width=12,dpi=80,pointsize=14,noRStudioGD = T)
  p1 <- ggiNEXT(iN_results, type=1, se=TRUE, facet.var="Order.q", color.var="Assemblage", grey=FALSE) +
    theme(legend.position="right")
  
  #ggiNEXT(iN_results, type=1, se=TRUE, facet.var="Assemblage", grey=FALSE)  
  
  
  ##Sample Completeness Curve----
  
  dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
  ggiNEXT(iN_results, type=2, se=TRUE, facet.var="Order.q", color.var="Assemblage", grey=FALSE)  
  
  #Coverage-Based R/E Curve----
  
  dev.new(height=7,width=12,dpi=80,pointsize=14,noRStudioGD = T)
  p2 <- ggiNEXT(iN_results, type=3, se=TRUE, facet.var="Order.q", color.var="Assemblage", grey=FALSE) +
    theme(legend.position="none")
  
  #ggiNEXT(iN_results, type=3, se=TRUE, facet.var="Assemblage", grey=FALSE)
  
  
  g1 <- list (p1,p2)
  
  dev.new(height=14,width=22,dpi=80,pointsize=14,noRStudioGD = T)
  grid.arrange(
    grobs = g1,
    widths = c(2, 2, 1),
    layout_matrix = rbind(c(1, 1, 1),
                          c(2, 2, NA)))
  
  
  
#SPECIES PRESENCE/ABSENCE DATA----


propzero<-data.frame(species=names(apply(pcadata[,which(colnames(pcadata)=="A.sti"):which(colnames(pcadata)=="T.can")],2,FUN = function(x)table(x==0)[2]/sum(table(x==0)))),propzero=apply(pcadata[,which(colnames(pcadata)=="A.sti"):which(colnames(pcadata)=="T.can")],2,FUN = function(x)table(x==0)[2]/sum(table(x==0))))
#only use species with prop zero less then 80

propzero$abundance<-apply(pcadata[,which(colnames(pcadata)=="A.sti"):which(colnames(pcadata)=="T.can")],2,FUN = function(x)sum(x))
propzero$binomial<-ifelse(propzero$propzero<0.9,1,0)
binomsp<-propzero$species[propzero$binomial==1]

singledata<-pcadata
colnames(singledata)[24]<-c("antechinus")
colnames(singledata)[25]<-c("melomys")
colnames(singledata)[26]<-c("bandicoot")
colnames(singledata)[27]<-c("koala")
colnames(singledata)[28]<-c("bushrat")
colnames(singledata)[29]<-c("pig")
colnames(singledata)[30]<-c("stigmatica")
colnames(singledata)[31]<-c("thetis")
colnames(singledata)[32]<-c("possum")
singledata<-dplyr::select(singledata,-koala,-pig,-pca.comp1,-pca.comp2)

##Melomys----

mel<-glm(ifelse(singledata[,which(colnames(singledata)=='melomys')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

melnull<-glm(ifelse(singledata[,which(colnames(singledata)=='melomys')]==0,0,1)~1,data=singledata,family=binomial )

mellrt<-anova(melnull,mel,test="LRT")
melp<-round(mellrt$`Pr(>Chi)`[2],3)

meltab<-list("FHC"=mel,"Null"=melnull)
aictab(meltab)

###Predictions

mel.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

mel.4<-predict(object=mel,newdata=mel.1,se.fit =T,type="link")

mel.5<-data.frame(mel.1,fit.link=mel.4$fit,se.link=mel.4$se.fit)
mel.5$lci.link<-mel.5$fit.link-(1.96*mel.5$se.link)
mel.5$uci.link<-mel.5$fit.link+(1.96*mel.5$se.link)

mel.5$fit<-invlogit(mel.5$fit.link)
mel.5$se<-invlogit(mel.5$se.link)
mel.5$lci<-invlogit(mel.5$lci.link)
mel.5$uci<-invlogit(mel.5$uci.link)

##Bandicoot----

band<-glm(ifelse(singledata[,which(colnames(singledata)=='bandicoot')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

bandnull<-glm(ifelse(singledata[,which(colnames(singledata)=='bandicoot')]==0,0,1)~1,data=singledata,family=binomial )

bandlrt<-anova(bandnull,band,test="LRT")
bandp<-round(bandlrt$`Pr(>Chi)`[2],3)

bandtab<-list("FHC"=band,"Null"=bandnull)
aictab(bandtab)

###Predictions

band.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

band.4<-predict(object=band,newdata=band.1,se.fit =T,type="link")

band.5<-data.frame(band.1,fit.link=band.4$fit,se.link=band.4$se.fit)
band.5$lci.link<-band.5$fit.link-(1.96*band.5$se.link)
band.5$uci.link<-band.5$fit.link+(1.96*band.5$se.link)

band.5$fit<-invlogit(band.5$fit.link)
band.5$se<-invlogit(band.5$se.link)
band.5$lci<-invlogit(band.5$lci.link)
band.5$uci<-invlogit(band.5$uci.link)


###Melomys and Bandicoot Figure----

dev.new(height=4,width=8,dpi=70,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,3,3),mfrow=c(1,2))
plot(1:3,mel.5$fit,type="p", ylim=c(min(0,min(mel.5$lci)),max(1,max(mel.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Habitat Use",las=1,cex=2,pch=19,xaxt="n",main= expression(paste("(a) ", italic("Melomys cervinipes"), sep = " ")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=2:3, y0=mel.5$lci[2:3],x1=2:3, y1=mel.5$uci[2:3],angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",melp,sep=""))
#mtext("(a)",3,line = 0.8,cex=1.5,adj=0)
points(1.2,0.92,pch="*",cex=2)

plot(1:3,band.5$fit,type="p", ylim=c(min(0,min(band.5$lci)),max(1,max(band.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Habitat Use",las=1,cex=2,pch=19,xaxt="n",main=expression(paste("(b) ",italic("Perameles nasuta")), sep = " "),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1, y0=band.5$lci[1],x1=1, y1=band.5$uci[1],angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",bandp,sep=""))
#mtext("(b)",3,0.4,F,0)
points(2.2,0.1,pch="*",cex=2)
points(3.2,0.1,pch="*",cex=2)

##Antechinus----

ante<-glm(ifelse(singledata[,which(colnames(singledata)=='antechinus')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

antenull<-glm(ifelse(singledata[,which(colnames(singledata)=='antechinus')]==0,0,1)~1,data=singledata,family=binomial )

antelrt<-anova(antenull,ante,test="LRT")
antep<-round(antelrt$`Pr(>Chi)`[2],3)

antetab<-list("FHC"=ante,"Null"=antenull)
aictab(antetab)

###Predictions

ante.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

ante.4<-predict(object=ante,newdata=ante.1,se.fit =T,type="link")

ante.6<-data.frame(ante.1,fit.link=ante.4$fit,se.link=ante.4$se.fit)
ante.6$lci.link<-ante.6$fit.link-(1.96*ante.6$se.link)
ante.6$uci.link<-ante.6$fit.link+(1.96*ante.6$se.link)

ante.6$fit<-invlogit(ante.6$fit.link)
ante.6$se<-invlogit(ante.6$se.link)
ante.6$lci<-invlogit(ante.6$lci.link)
ante.6$uci<-invlogit(ante.6$uci.link)

##Bushrat----

bush<-glm(ifelse(singledata[,which(colnames(singledata)=='bushrat')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

bushnull<-glm(ifelse(singledata[,which(colnames(singledata)=='bushrat')]==0,0,1)~1,data=singledata,family=binomial )

bushlrt<-anova(bushnull,bush,test="LRT")
bushp<-round(bushlrt$`Pr(>Chi)`[2],3)

bushtab<-list("FHC"=bush,"Null"=bushnull)
aictab(bushtab)

###Predictions

bush.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

bush.4<-predict(object=bush,newdata=bush.1,se.fit =T,type="link")

bush.5<-data.frame(bush.1,fit.link=bush.4$fit,se.link=bush.4$se.fit)
bush.5$lci.link<-bush.5$fit.link-(1.96*bush.5$se.link)
bush.5$uci.link<-bush.5$fit.link+(1.96*bush.5$se.link)

bush.5$fit<-invlogit(bush.5$fit.link)
bush.5$se<-invlogit(bush.5$se.link)
bush.5$lci<-invlogit(bush.5$lci.link)
bush.5$uci<-invlogit(bush.5$uci.link)

##Thylogale stigmatica----

sti<-glm(ifelse(singledata[,which(colnames(singledata)=='stigmatica')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

stinull<-glm(ifelse(singledata[,which(colnames(singledata)=='stigmatica')]==0,0,1)~1,data=singledata,family=binomial )

stilrt<-anova(stinull,sti,test="LRT")
stip<-round(stilrt$`Pr(>Chi)`[2],3)

stitab<-list("FHC"=sti,"Null"=stinull)
aictab(stitab)

###Predictions

sti.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

sti.4<-predict(object=sti,newdata=sti.1,se.fit =T,type="link")

sti.5<-data.frame(sti.1,fit.link=sti.4$fit,se.link=sti.4$se.fit)
sti.5$lci.link<-sti.5$fit.link-(1.96*sti.5$se.link)
sti.5$uci.link<-sti.5$fit.link+(1.96*sti.5$se.link)

sti.5$fit<-invlogit(sti.5$fit.link)
sti.5$se<-invlogit(sti.5$se.link)
sti.5$lci<-invlogit(sti.5$lci.link)
sti.5$uci<-invlogit(sti.5$uci.link)

##Thylogale thetis----

the<-glm(ifelse(singledata[,which(colnames(singledata)=='thetis')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

thenull<-glm(ifelse(singledata[,which(colnames(singledata)=='thetis')]==0,0,1)~1,data=singledata,family=binomial )

thelrt<-anova(thenull,the,test="LRT")
thep<-round(thelrt$`Pr(>Chi)`[2],3)

thetab<-list("FHC"=the,"Null"=thenull)
aictab(thetab)

###Predictions

the.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

the.4<-predict(object=the,newdata=the.1,se.fit =T,type="link")

the.5<-data.frame(the.1,fit.link=the.4$fit,se.link=the.4$se.fit)
the.5$lci.link<-the.5$fit.link-(1.96*the.5$se.link)
the.5$uci.link<-the.5$fit.link+(1.96*the.5$se.link)

the.5$fit<-invlogit(the.5$fit.link)
the.5$se<-invlogit(the.5$se.link)
the.5$lci<-invlogit(the.5$lci.link)
the.5$uci<-invlogit(the.5$uci.link)

##Possum----

pos<-glm(ifelse(singledata[,which(colnames(singledata)=='possum')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

posnull<-glm(ifelse(singledata[,which(colnames(singledata)=='possum')]==0,0,1)~1,data=singledata,family=binomial )

poslrt<-anova(posnull,pos,test="LRT")
posp<-round(poslrt$`Pr(>Chi)`[2],3)

postab<-list("FHC"=pos,"Null"=posnull)
aictab(postab)

###Predictions

pos.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

pos.4<-predict(object=pos,newdata=pos.1,se.fit =T,type="link")

pos.5<-data.frame(pos.1,fit.link=pos.4$fit,se.link=pos.4$se.fit)
pos.5$lci.link<-pos.5$fit.link-(1.96*pos.5$se.link)
pos.5$uci.link<-pos.5$fit.link+(1.96*pos.5$se.link)

pos.5$fit<-invlogit(pos.5$fit.link)
pos.5$se<-invlogit(pos.5$se.link)
pos.5$lci<-invlogit(pos.5$lci.link)
pos.5$uci<-invlogit(pos.5$uci.link)

###Graphing----

dev.new(height=21,width=14,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,2),mgp=c(2.5,1,0)) 

plot(1:3,ante.6$fit,type="p", ylim=c(min(0,min(ante.6$lci)),max(1,max(ante.6$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n", main=expression(italic("Antechinus stuartii")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=ante.6$lci,x1=1:3, y1=ante.6$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",antep,sep=""))
mtext("(a)",3,0.4,F,0)

plot(1:3,bush.5$fit,type="p", ylim=c(min(0,min(bush.5$lci)),max(1,max(bush.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Rattus fuscipes")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=2:3, y0=bush.5$lci[2:3],x1=2:3, y1=bush.5$uci[2:3],angle=90,length=0.1, code=3, lwd=2)
text(3,0.1,labels=paste("p=",bushp,sep=""))
mtext("(b)",3,0.4,F,0)
points(1.2,0.92,pch="*",cex=2)

plot(1:3,the.5$fit,type="p", ylim=c(min(0,min(the.5$lci)),max(1,max(the.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale thetis")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=the.5$lci,x1=1:3, y1=the.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.95,labels=paste("p=",thep,sep=""))
mtext("(c)",3,0.4,F,0)

plot(1:3,sti.5$fit,type="p", ylim=c(min(0,min(sti.5$lci)),max(1,max(sti.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale stigmatica")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=sti.5$lci,x1=1:3, y1=sti.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",stip,sep=""))
mtext("(d)",3,0.4,F,0)

plot(1:3,pos.5$fit,type="p", ylim=c(min(0,min(pos.5$lci)),max(1,max(pos.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Trichosurus caninus")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=pos.5$lci,x1=1:3, y1=pos.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",posp,sep=""))
mtext("(e)",3,0.4,F,0)

#MOVEMENT----

site<-rawsite
site[,c(2,3,6,7,8,9)]<-NULL
names(site)<-c("Site","Elevation","Fire_habitat_catergory","Dis_to_road","Dis_to_path","Understory","Canopy_Cover","Midstory","Surrounding_Unburnt","Surrounding_Low","Surrounding_Moderate","Surrounding_High","Surrounding_Extreme","Surrounding_Burn", "Dis_to_rainforest_boundry", "Habitat")


table(mammals$Behaviour,mammals$Identification)

head(mammals,3);dim(mammals)
head(site,3);dim(site)

list(
  melomys = table(move.mel2$Fire_habitat_catergory[move.mel2$active==1])/table(move.mel2$Fire_habitat_catergory),
  rattus = table(move.fus2$Fire_habitat_catergory[move.fus2$active==1])/table(move.fus2$Fire_habitat_catergory),
  stigmatica = table(move.stig3$Fire_habitat_catergory[move.stig3$active==1])/table(move.stig3$Fire_habitat_catergory),
  thylogale = table(move.the2$Fire_habitat_catergory[move.the2$active==1])/table(move.the2$Fire_habitat_catergory),
  possum = table(move.can2$Fire_habitat_catergory[move.can2$active==1])/table(move.can2$Fire_habitat_catergory)
)



##melomys data----

move.mel<-mammals[c(which(mammals$Identification=="Melomys_cervinipes")),]

move.mel<-move.mel[c(which(move.mel$Behaviour==c("Foraging_Eating","No_")),which(move.mel$Behaviour==c("Walking_Running"))),]

move.mel1<-data.frame(
  inactive=ifelse(move.mel$Behaviour==c("Foraging_Eating"),1,0),
  active=ifelse(move.mel$Behaviour== c("Walking_Running"),1,0),
  Site=move.mel$Site)

move.mel2<-merge(x=move.mel1,y=site)

head(move.mel2,3);dim(move.mel2)

table(move.mel2$inactive==move.mel2$active)

str(move.mel2)
move.mel2$Fire_habitat_catergory<-factor(move.mel2$Fire_habitat_catergory,levels=c("Unburnt_Rainforest", "Burnt_Rainforest" ,  "Burnt_Sclerophyll" ))
unique(move.mel2$Fire_habitat_catergory)

table(move.mel2$Fire_habitat_catergory[move.mel2$active==1])/table(move.mel2$Fire_habitat_catergory)


###melomys modeling----

modm.null<-glmer(active~1+(1|Site),data=move.mel2,family=binomial)

modm.FHC<-glmer(active~Fire_habitat_catergory+(1|Site),data=move.mel2,family=binomial)

modm.1a<-glmer(active~Fire_habitat_catergory+Elevation+(1|Site),data=move.mel2,family=binomial)
modm.1b<-glmer(active~Fire_habitat_catergory*Elevation+(1|Site),data=move.mel2,family=binomial)
AICc(modm.1a);AICc(modm.1b)

modm.2a<-glmer(active~Fire_habitat_catergory+Dis_to_road+(1|Site),data=move.mel2,family=binomial)
modm.2b<-glmer(active~Fire_habitat_catergory*Dis_to_road+(1|Site),data=move.mel2,family=binomial)
AICc(modm.2a);AICc(modm.2b)

modm.3a<-glmer(active~Fire_habitat_catergory+Dis_to_path+(1|Site),data=move.mel2,family=binomial)
modm.3b<-glmer(active~Fire_habitat_catergory*Dis_to_path+(1|Site),data=move.mel2,family=binomial)
AICc(modm.3a);AICc(modm.3b)

modm.4a<-glmer(active~Fire_habitat_catergory+Dis_to_rainforest_boundry+(1|Site),data=move.mel2,family=binomial)
modm.4b<-glmer(active~Fire_habitat_catergory*Dis_to_rainforest_boundry+(1|Site),data=move.mel2,family=binomial)
AICc(modm.4a);AICc(modm.4b)

modm.5a<-glmer(active~Fire_habitat_catergory+Midstory+(1|Site),data=move.mel2,family=binomial)
modm.5b<-glmer(active~Fire_habitat_catergory*Midstory+(1|Site),data=move.mel2,family=binomial)
AICc(modm.5a);AICc(modm.5b)

modm.6a<-glmer(active~Fire_habitat_catergory+Understory+(1|Site),data=move.mel2,family=binomial)
modm.6b<-glmer(active~Fire_habitat_catergory*Understory+(1|Site),data=move.mel2,family=binomial)
AICc(modm.6a);AICc(modm.6b)

modm.list<-c('null'=modm.null,'FHC'=modm.FHC,'elevation'=modm.1a,"road"=modm.2a,'path'=modm.3a,'boundary'=modm.4b,'midstory'=modm.5a,'understory'=modm.6b)

aictab(modm.list)

anova(modm.FHC,modm.4b)

summary(modm.4b)

###melomys predictions----

modm.4b
class(modm.4b)
min(move.mel2$Dis_to_rainforest_boundry)
max(move.mel2$Dis_to_rainforest_boundry)
unique(move.mel2$Site)

move.melpred<-data.frame(Fire_habitat_catergory=factor(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),rep(50,3)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Dis_to_rainforest_boundry=seq(min(move.mel2$Dis_to_rainforest_boundry),max(move.mel2$Dis_to_rainforest_boundry),length.out=50))

summary(modm.4b)
anova(modm.4b,modm.null)


summary(model)$coefficients
summary(modm.4b)$coefficients

move.mel2$Site<-as.factor(move.mel2$Site)
summary(modm.4b)



move.melpred3<-predict(object=model,newdata=move.melpred,se.fit =TRUE,type="link")
head(move.melpred3)

move.melpred4<-data.frame(move.melpred,fit.link=move.melpred3$fit,se.link=move.melpred3$se.fit)
move.melpred4$lci.link<-move.melpred4$fit.link-(1.96*move.melpred4$se.link)
move.melpred4$uci.link<-move.melpred4$fit.link+(1.96*move.melpred4$se.link)

move.melpred4$fit<-invlogit(move.melpred4$fit.link)
move.melpred4$se<-invlogit(move.melpred4$se.link)
move.melpred4$lci<-invlogit(move.melpred4$lci.link)
move.melpred4$uci<-invlogit(move.melpred4$uci.link)

head(move.melpred4)

move.melpred5<-move.melpred4



move.melpred3<-predictSE(mod=modm.4b,newdata=move.melpred,se.fit =TRUE,type="link")
head(move.melpred3)

move.melpred4<-data.frame(move.melpred,fit.link=move.melpred3$fit,se.link=move.melpred3$se.fit)
move.melpred4$lci.link<-move.melpred4$fit.link-(1.96*move.melpred4$se.link)
move.melpred4$uci.link<-move.melpred4$fit.link+(1.96*move.melpred4$se.link)

move.melpred4$fit<-invlogit(move.melpred4$fit.link)
move.melpred4$se<-invlogit(move.melpred4$se.link)
move.melpred4$lci<-invlogit(move.melpred4$lci.link)
move.melpred4$uci<-invlogit(move.melpred4$uci.link)


head(move.melpred4[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest",])


###melomys graphing----

head(move.melpred4[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest",])
head(move.melpred4[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest",])
head(move.melpred4[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll",])

las=1,cex=2,pch=19,xaxt="n"
ylim=c(min(move.melpred4$lci),max(move.melpred4$uci))

move.melpred4[77:100,7:10]
?replace
move.melfit<-replace(x=move.melpred4$fit,list=77:100,values = NA)
move.melse<-replace(x=move.melpred4$se,list=77:100,values = NA)
move.mellci<-replace(x=move.melpred4$lci,list=77:100,values = NA)
move.meluci<-replace(x=move.melpred4$uci,list=77:100,values = NA)

move.melpred4$fit<-move.melfit
move.melpred4$se<-move.melse
move.melpred4$lci<-move.mellci
move.melpred4$uci<-move.meluci



dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],type="l",lwd=2,ylab=" ",xlab="",main="Unburnt Rainforest",cex.main=1.2,ylim=c(0,1),las=1)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
mtext(side=1,line=2,'Distance To Rainforest Boundry(m)',cex=0.9)
mtext(side=2,line=2.3,"Probability Of Movement",cex=0.9)
mtext("(a)",3,0.1,F,0)

plot(  x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,ylab=" ",xlab=" ",type="l",main="Burnt Rainforest",cex.main=1.2,xlim=c(min(move.melpred4$Dis_to_rainforest_boundry),max(move.melpred4$Dis_to_rainforest_boundry)),las=1)
lines(x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,lty=2)
lines(x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,lty=2)
mtext(side=1,line=2,at=160,'Distance To Rainforest Boundry(m)',cex=0.9)
mtext(side=2,line=2.3,"Probability Of Movement",cex=0.9)
mtext("(b)",3,0.1,F,0)

plot(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,ylab=" ",xlab=" ",type="l",main="Burnt Sclerophyll",cex.main=1.2,ylim=c(0,1),las=1)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
mtext(side=1,line=2,'Distance To Rainforest Boundry(m)',cex=0.9)
mtext(side=2,line=2.3,"Probability Of Movement",cex=0.9)
mtext("(c)",3,0.1,F,0)
mtext("(d)",3,0.1,F,550)

head(move.mel2)

tapply(move.mel2$Dis_to_rainforest_boundry,move.mel2$Fire_habitat_catergory,FUN=summary)

##rattus data----

table(mammals$Behaviour,mammals$Identification)


move.fus<-mammals[c(which(mammals$Identification=="Rattus_fuscipes")),]

dim(move.fus)
head(rawdata,30)

move.fus<-move.fus[c(
  which(move.fus$Behaviour==c("Foraging_Eating")),
  which(move.fus$Behaviour==c("Walking_Running")),
  which(move.fus$Behaviour==c("Grooming")),
  which(move.fus$Behaviour==c("Hopping_Jumping")),
  which(move.fus$Behaviour==c("Alert"))
),]

dim(move.fus)

move.fus1<-data.frame(
  inactive=ifelse(move.fus$Behaviour=="Foraging_Eating"|move.fus$Behaviour=="Grooming",1,0),
  active=ifelse(move.fus$Behaviour=="Walking_Running"| move.fus$Behaviour== "Hopping_Jumping"| move.fus$Behaviour== "Alert",1,0),
  Site=move.fus$Site)

dim(move.fus1)
sum(move.fus1$active)
sum(move.fus1$inactive)

move.fus2<-merge(x=move.fus1,y=site)

move.fus2$Fire_habitat_catergory<-factor(move.fus2$Fire_habitat_catergory,levels=c("Unburnt_Rainforest", "Burnt_Rainforest" ,  "Burnt_Sclerophyll" ))

head(move.fus2);dim(move.fus2)

trat<-table(move.fus2$active,move.fus2$Fire_habitat_catergory)
trat[2,]/sum(trat[1,],trat[2,])

table(move.fus2$Fire_habitat_catergory[move.fus2$active==1])/table(move.fus2$Fire_habitat_catergory)

###rattus modeling----

modr.null<-glmer(active~1+(1|Site),data=move.fus2,family=binomial)

modr.FHC<-glmer(active~Fire_habitat_catergory+(1|Site),data=move.fus2,family=binomial)
AICc(modr.null);AICc(modr.FHC)

modr.1a<-glmer(active~Fire_habitat_catergory+Elevation+(1|Site),data=move.fus2,family=binomial)
modr.1b<-glmer(active~Fire_habitat_catergory*Elevation+(1|Site),data=move.fus2,family=binomial)
AICc(modr.1a);AICc(modr.1b)

modr.2a<-glmer(active~Fire_habitat_catergory+Dis_to_road+(1|Site),data=move.fus2,family=binomial)
modr.2b<-glmer(active~Fire_habitat_catergory*Dis_to_road+(1|Site),data=move.fus2,family=binomial)
AICc(modr.2a);AICc(modr.2b)

modr.3a<-glmer(active~Fire_habitat_catergory+Dis_to_path+(1|Site),data=move.fus2,family=binomial)
modr.3b<-glmer(active~Fire_habitat_catergory*Dis_to_path+(1|Site),data=move.fus2,family=binomial)
AICc(modr.3a);AICc(modr.3b)

modr.4a<-glmer(active~Fire_habitat_catergory+Dis_to_rainforest_boundry+(1|Site),data=move.fus2,family=binomial)
modr.4b<-glmer(active~Fire_habitat_catergory*Dis_to_rainforest_boundry+(1|Site),data=move.fus2,family=binomial)
AICc(modr.4a);AICc(modr.4b)

modr.5a<-glmer(active~Fire_habitat_catergory+Midstory+(1|Site),data=move.fus2,family=binomial)
modr.5b<-glmer(active~Fire_habitat_catergory*Midstory+(1|Site),data=move.fus2,family=binomial)
AICc(modr.5a);AICc(modr.5b)

modr.6a<-glmer(active~Fire_habitat_catergory+Understory+(1|Site),data=move.fus2,family=binomial)
modr.6b<-glmer(active~Fire_habitat_catergory*Understory+(1|Site),data=move.fus2,family=binomial)
AICc(modr.6a);AICc(modr.6b)

modr.list<-c('null'=modr.null,'FHC'=modr.FHC,'elevation'=modr.1a,"road"=modr.2a,'path'=modr.3a,'boundary'=modr.4a,'midstory'=modr.5a,'understory'=modr.6a)

aictab(modr.list)

anova(modr.FHC,modr.6a)

summary(modr.6a)

###rattus predictions----

move.fus2$Understory<-mapping[move.fus2$Understory]

modr.6a

head(move.fus2,3);dim(move.fus2)

unique(move.fus2$Fire_habitat_catergory)

move.fuspredunder<-data.frame(Fire_habitat_catergory=factor(rep('Unburnt_Rainforest',4), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Understory=c(0,1,2,3))

table(move.fus2$Fire_habitat_catergory)

move.fuspredFHC<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Understory=c(2,2,2))

table(move.fus2$Fire_habitat_catergory,move.fus2$Understory)


move.fuspred5<-predictSE(mod=modr.6a,newdata=move.fuspredunder,se.fit =TRUE,type="link")

move.fuspred6<-data.frame(move.fuspredunder,fit.link=move.fuspred5$fit,se.link=move.fuspred5$se.fit)
move.fuspred6$lci.link<-move.fuspred6$fit.link-(1.96*move.fuspred6$se.link)
move.fuspred6$uci.link<-move.fuspred6$fit.link+(1.96*move.fuspred6$se.link)

move.fuspred6$fit<-invlogit(move.fuspred6$fit.link)
move.fuspred6$se<-invlogit(move.fuspred6$se.link)
move.fuspred6$lci<-invlogit(move.fuspred6$lci.link)
move.fuspred6$uci<-invlogit(move.fuspred6$uci.link)

move.fuspred7<-predictSE(mod=modr.6a,newdata=move.fuspredFHC,se.fit =TRUE,type="link")

move.fuspred8<-data.frame(move.fuspredFHC,fit.link=move.fuspred7$fit,se.link=move.fuspred7$se.fit)
move.fuspred8$lci.link<-move.fuspred8$fit.link-(1.96*move.fuspred8$se.link)
move.fuspred8$uci.link<-move.fuspred8$fit.link+(1.96*move.fuspred8$se.link)

move.fuspred8$fit<-invlogit(move.fuspred8$fit.link)
move.fuspred8$se<-invlogit(move.fuspred8$se.link)
move.fuspred8$lci<-invlogit(move.fuspred8$lci.link)
move.fuspred8$uci<-invlogit(move.fuspred8$uci.link)


###rattus graphing----


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(1:4,move.fuspred6$fit,type="p", ylim=c(min(0,min(move.fuspred6$lci)),max(1,max(move.fuspred6$uci))),xlab="Understory",ylab="Probability of Movement",las=1,cex=2,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=c(0,1,2,3))
arrows(x0=1:4, y0=move.fuspred6$lci,x1=1:4, y1=move.fuspred6$uci,angle=90,length=0.1, code=3, lwd=2)
mtext("(a)",3,0.2,F,0.5)

plot(1:3,move.fuspred8$fit,type="p", ylim=c(min(0,min(move.fuspred8$lci)),max(1,max(move.fuspred8$uci))),xlab="Fire Habitat Catergory",ylab="Probability of Movement",las=1,cex=2,pch=19,xaxt="n")
axis(side=1,at=1:3,labels=c("UBR","BR","BS"))
arrows(x0=1:3, y0=move.fuspred8$lci,x1=1:3, y1=move.fuspred8$uci,angle=90,length=0.1, code=3, lwd=2)
mtext("(b)",3,0.2,F,0.5)

##stigmatica data----

table(mammals$Behaviour,mammals$Identification)


move.stig<-mammals[c(which(mammals$Identification=="Thylogale_stigmatica")),]

dim(move.stig)

unique(move.stig$Behaviour)

move.stig<-move.stig[c(
  which(move.stig$Behaviour==c("Foraging_Eating")),
  which(move.stig$Behaviour==c("Hopping_Jumping")),
  which(move.stig$Behaviour==c("Alert"))
),]

dim(move.stig)

move.stig1<-data.frame(
  inactive=ifelse(move.stig$Behaviour=="Foraging_Eating",1,0),
  active=ifelse(move.stig$Behaviour== "Hopping_Jumping"| move.stig$Behaviour== "Alert",1,0),
  Site=move.stig$Site)


dim(move.stig1)
sum(move.stig1$active)
sum(move.stig1$inactive)

move.stig2<-merge(x=move.stig1,y=site)

head(move.stig2,3);dim(move.stig2)

move.stig2$Fire_habitat_catergory<-factor(move.stig2$Fire_habitat_catergory,levels=c("Unburnt_Rainforest", "Burnt_Rainforest","Burnt_Sclerophyll" ))

tstig<-table(move.stig2$active,move.stig2$Fire_habitat_catergory)
tstig[2,]/sum(tstig[1,],tstig[2,])

tstig2<-table(move.stig2$inactive,move.stig2$Fire_habitat_catergory)

move.stig3<- move.stig2

table(move.stig3$Fire_habitat_catergory[move.stig3$active==1])/table(move.stig3$Fire_habitat_catergory)

move.stig2[which(move.stig2$Fire_habitat_catergory==c("Burnt_Sclerophyll")),]


move.stig2<-move.stig2[c(
  which(move.stig2$Fire_habitat_catergory==c("Unburnt_Rainforest")),
  which(move.stig2$Fire_habitat_catergory==c("Burnt_Rainforest"))
),]
dim(move.stig2)
str(move.stig2)

levels.FHC<-list("Unburnt_Rainforest","Burnt_Rainforest")

move.stig2$Fire_habitat_catergory<-factor(move.stig2$Fire_habitat_catergory,levels=levels.FHC)

str(move.stig2)

###stigmatica modeling----

mods.null<-glmer(active~1+(1|Site),data=move.stig2,family=binomial)

mods.FHC<-glmer(active~Fire_habitat_catergory+(1|Site),data=move.stig2,family=binomial)

AICc(mods.null);AICc(mods.FHC)
summary(mods.FHC)
anova(mods.FHC)

mods.1a<-glmer(active~Fire_habitat_catergory+Elevation+(1|Site),data=move.stig2,family=binomial)
mods.1b<-glmer(active~Fire_habitat_catergory*Elevation+(1|Site),data=move.stig2,family=binomial)
AICc(mods.1a);AICc(mods.1b)

mods.2a<-glmer(active~Fire_habitat_catergory+Dis_to_road+(1|Site),data=move.stig2,family=binomial)
mods.2b<-glmer(active~Fire_habitat_catergory*Dis_to_road+(1|Site),data=move.stig2,family=binomial)
AICc(mods.2a);AICc(mods.2b)

mods.3a<-glmer(active~Fire_habitat_catergory+Dis_to_path+(1|Site),data=move.stig2,family=binomial)
mods.3b<-glmer(active~Fire_habitat_catergory*Dis_to_path+(1|Site),data=move.stig2,family=binomial)
AICc(mods.3a);AICc(mods.3b)

mods.4a<-glmer(active~Fire_habitat_catergory+Dis_to_rainforest_boundry+(1|Site),data=move.stig2,family=binomial)
mods.4b<-glmer(active~Fire_habitat_catergory*Dis_to_rainforest_boundry+(1|Site),data=move.stig2,family=binomial)
AICc(mods.4a);AICc(mods.4b)

mods.5a<-glmer(active~Fire_habitat_catergory+Midstory+(1|Site),data=move.stig2,family=binomial)
mods.5b<-glmer(active~Fire_habitat_catergory*Midstory+(1|Site),data=move.stig2,family=binomial)
AICc(mods.5a);AICc(mods.5b)

mods.6a<-glmer(active~Fire_habitat_catergory+Understory+(1|Site),data=move.stig2,family=binomial)
mods.6b<-glmer(active~Fire_habitat_catergory*Understory+(1|Site),data=move.stig2,family=binomial)
AICc(mods.6a);AICc(mods.6b)

mods.list<-c('null'=mods.null,'FHC'=mods.FHC,'elevation'=mods.1a,"road"=mods.2a,'path'=mods.3b,'boundary'=mods.4b,'midstory'=mods.5a,'understory'=mods.6b)

aictab(mods.list)

anova(mods.FHC,mods.1a)

summary(mods.1a)

head(move.stig2)

###stigmatica predictions----

mods.1a

move.stigpredelv<-data.frame(Fire_habitat_catergory=factor(rep('Unburnt_Rainforest',50), levels=c('Unburnt_Rainforest','Burnt_Rainforest')),Elevation=seq(min(move.stig2$Elevation),max(move.stig2$Elevation),length.out=50))

move.stigpredFHC<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest','Burnt_Rainforest'), levels=c('Unburnt_Rainforest','Burnt_Rainforest')),Elevation=mean(move.stig2$Elevation))

head(move.stigpred)
table(move.stig2$Fire_habitat_catergory)


move.stigpred5<-predictSE(mod=mods.1a,newdata=move.stigpredelv,se.fit =TRUE,type="link")
head(move.stigpred3)
warnings()

move.stigpred6<-data.frame(move.stigpredelv,fit.link=move.stigpred5$fit,se.link=move.stigpred5$se.fit)
move.stigpred6$lci.link<-move.stigpred6$fit.link-(1.96*move.stigpred6$se.link)
move.stigpred6$uci.link<-move.stigpred6$fit.link+(1.96*move.stigpred6$se.link)

move.stigpred6$fit<-invlogit(move.stigpred6$fit.link)
move.stigpred6$se<-invlogit(move.stigpred6$se.link)
move.stigpred6$lci<-invlogit(move.stigpred6$lci.link)
move.stigpred6$uci<-invlogit(move.stigpred6$uci.link)

head(move.stigpred6)

move.stigpred7<-predictSE(mod=mods.1a,newdata=move.stigpredFHC,se.fit =TRUE,type="link")
head(move.stigpred7)

move.stigpred8<-data.frame(move.stigpredFHC,fit.link=move.stigpred7$fit,se.link=move.stigpred7$se.fit)
move.stigpred8$lci.link<-move.stigpred8$fit.link-(1.96*move.stigpred8$se.link)
move.stigpred8$uci.link<-move.stigpred8$fit.link+(1.96*move.stigpred8$se.link)

move.stigpred8$fit<-invlogit(move.stigpred8$fit.link)
move.stigpred8$se<-invlogit(move.stigpred8$se.link)
move.stigpred8$lci<-invlogit(move.stigpred8$lci.link)
move.stigpred8$uci<-invlogit(move.stigpred8$uci.link)

head(move.stigpred8)



###stigmatica graphing----

head(move.stigpred4[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest",])
head(move.stigpred4[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest",])
head(move.stigpred4[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll",])

#Pub quality graph----
dev.new(height=4,width=8,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(1,2))

plot(x=move.stigpred6$Elevation,y=move.stigpred6$fit,type="l",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="",las=1)
lines(x=move.stigpred6$Elevation,y=move.stigpred6$lci,lwd=2,lty=2)
lines(x=move.stigpred6$Elevation,y=move.stigpred6$uci,lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=1.2)
mtext(side=2,line=2.3,"Probability Of Movement",cex=1.2)
mtext("(a)",3,0.2,F,765)

plot(1:2,move.stigpred8$fit,type="p", ylim=c(min(0,min(move.stigpred8$lci)),max(1,max(move.stigpred8$uci))),xlim=c(min(0.5),max(2.5)),las=1,cex=2,pch=19,xaxt="n",xlab = "",ylab = "")
axis(side=1,at=1:2,labels=c('UBR','BR'))
arrows(x0=1:2, y0=move.stigpred8$lci,x1=1:2, y1=move.stigpred8$uci,angle=90,length=0.1, code=3, lwd=2)
mtext(side=2,line=2.3,"Probability Of Movement",cex=1.2)
mtext(side=1,line=2,"Fire Habitat Category",cex=1.2)
mtext("(b)",side = 3,line = 0.2,outer = F,at = 0.1)

##thylogale data----

table(mammals$Behaviour,mammals$Identification)


move.the<-mammals[c(which(mammals$Identification=="Thylogale_thetis")),]

dim(move.the)

unique(move.the$Behaviour)

move.the<-move.the[c(
  which(move.the$Behaviour==c("Foraging_Eating")),
  which(move.the$Behaviour==c("Hopping_Jumping")),
  which(move.the$Behaviour==c("Walking_Running"))
),]

dim(move.the)

move.the1<-data.frame(
  inactive=ifelse(move.the$Behaviour=="Foraging_Eating",1,0),
  active=ifelse(move.the$Behaviour== "Hopping_Jumping"| move.the$Behaviour== "Walking_Running",1,0),
  Site=move.the$Site)


dim(move.the1)
sum(move.the1$active)
sum(move.the1$inactive)

move.the2<-merge(x=move.the1,y=site)

head(move.the2);dim(move.the2)

move.the2$Fire_habitat_catergory<-factor(move.the2$Fire_habitat_catergory,levels=c("Unburnt_Rainforest", "Burnt_Rainforest" ,  "Burnt_Sclerophyll" ))


tthe<-table(move.the2$active,move.the2$Fire_habitat_catergory)
tthe[2,]/sum(tthe[1,],tthe[2,])

table(move.the2$Fire_habitat_catergory[move.the2$active==1])/table(move.the2$Fire_habitat_catergory)


###thylogale modeling----

modt.null<-glmer(active~1+(1|Site),data=move.the2,family=binomial)

modt.FHC<-glmer(active~Fire_habitat_catergory+(1|Site),data=move.the2,family=binomial)

modt.1a<-glmer(active~Fire_habitat_catergory+Elevation+(1|Site),data=move.the2,family=binomial)
modt.1b<-glmer(active~Fire_habitat_catergory*Elevation+(1|Site),data=move.the2,family=binomial)
AICc(modt.1a);AICc(modt.1b)

modt.2a<-glmer(active~Fire_habitat_catergory+Dis_to_road+(1|Site),data=move.the2,family=binomial)
modt.2b<-glmer(active~Fire_habitat_catergory*Dis_to_road+(1|Site),data=move.the2,family=binomial)
AICc(modt.2a);AICc(modt.2b)

modt.3a<-glmer(active~Fire_habitat_catergory+Dis_to_path+(1|Site),data=move.the2,family=binomial)
modt.3b<-glmer(active~Fire_habitat_catergory*Dis_to_path+(1|Site),data=move.the2,family=binomial)
AICc(modt.3a);AICc(modt.3b)

modt.4a<-glmer(active~Fire_habitat_catergory+Dis_to_rainforest_boundry+(1|Site),data=move.the2,family=binomial)
modt.4b<-glmer(active~Fire_habitat_catergory*Dis_to_rainforest_boundry+(1|Site),data=move.the2,family=binomial)
AICc(modt.4a);AICc(modt.4b)

modt.5a<-glmer(active~Fire_habitat_catergory+Midstory+(1|Site),data=move.the2,family=binomial)
modt.5b<-glmer(active~Fire_habitat_catergory*Midstory+(1|Site),data=move.the2,family=binomial)
AICc(modt.5a);AICc(modt.5b)

modt.6a<-glmer(active~Fire_habitat_catergory+Understory+(1|Site),data=move.the2,family=binomial)
modt.6b<-glmer(active~Fire_habitat_catergory*Understory+(1|Site),data=move.the2,family=binomial)
AICc(modt.6a);AICc(modt.6b)

modt.list<-c('null'=modt.null,'FHC'=modt.FHC,'elevation'=modt.1a,"road"=modt.2a,'path'=modt.3a,'boundary'=modt.4a,'midstory'=modt.5a,'understory'=modt.6a)

aictab(modt.list)


##possum data----

table(mammals$Behaviour,mammals$Identification)


move.can<-mammals[c(which(mammals$Identification=="Trichosurus_caninus")),]

dim(move.can)

unique(move.can$Behaviour)

move.can<-move.can[c(
  which(move.can$Behaviour==c("Foraging_Eating")),
  which(move.can$Behaviour==c("Hopping_Jumping")),
  which(move.can$Behaviour==c("Walking_Running"))
),]

dim(move.can)

move.can1<-data.frame(
  inactive=ifelse(move.can$Behaviour=="Foraging_Eating",1,0),
  active=ifelse(move.can$Behaviour== "Hopping_Jumping"| move.can$Behaviour== "Walking_Running",1,0),
  Site=move.can$Site)


dim(move.can1)
sum(move.can1$active)
sum(move.can1$inactive)

move.can2<-merge(x=move.can1,y=site)

head(move.can2);dim(move.can2)

move.can2$Fire_habitat_catergory<-factor(move.can2$Fire_habitat_catergory,levels=c("Unburnt_Rainforest", "Burnt_Rainforest" ,  "Burnt_Sclerophyll" ))

table(move.can2$Fire_habitat_catergory[move.can2$active==1])/table(move.can2$Fire_habitat_catergory)

###possum modeling----

modp.null<-glmer(active~1+(1|Site),data=move.can2,family=binomial)

modp.FHC<-glmer(active~Fire_habitat_catergory+(1|Site),data=move.can2,family=binomial)

modp.1a<-glmer(active~Fire_habitat_catergory+Elevation+(1|Site),data=move.can2,family=binomial)
modp.1b<-glmer(active~Fire_habitat_catergory*Elevation+(1|Site),data=move.can2,family=binomial)
AICc(modp.1a);AICc(modp.1b)

modp.2a<-glmer(active~Fire_habitat_catergory+Dis_to_road+(1|Site),data=move.can2,family=binomial)
modp.2b<-glmer(active~Fire_habitat_catergory*Dis_to_road+(1|Site),data=move.can2,family=binomial)
AICc(modp.2a);AICc(modp.2b)

modp.3a<-glmer(active~Fire_habitat_catergory+Dis_to_path+(1|Site),data=move.can2,family=binomial)
modp.3b<-glmer(active~Fire_habitat_catergory*Dis_to_path+(1|Site),data=move.can2,family=binomial)
AICc(modp.3a);AICc(modp.3b)

modp.4a<-glmer(active~Fire_habitat_catergory+Dis_to_rainforest_boundry+(1|Site),data=move.can2,family=binomial)
modp.4b<-glmer(active~Fire_habitat_catergory*Dis_to_rainforest_boundry+(1|Site),data=move.can2,family=binomial)
AICc(modp.4a);AICc(modp.4b)

modp.5a<-glmer(active~Fire_habitat_catergory+Midstory+(1|Site),data=move.can2,family=binomial)
modp.5b<-glmer(active~Fire_habitat_catergory*Midstory+(1|Site),data=move.can2,family=binomial)
AICc(modp.5a);AICc(modp.5b)

modp.6a<-glmer(active~Fire_habitat_catergory+Understory+(1|Site),data=move.can2,family=binomial)
modp.6b<-glmer(active~Fire_habitat_catergory*Understory+(1|Site),data=move.can2,family=binomial)
AICc(modp.6a);AICc(modp.6b)

modp.list<-c('null'=modp.null,'FHC'=modp.FHC,'elevation'=modp.1a,"road"=modp.2a,'midstory'=modp.5a,'understory'=modp.6a)

aictab(modp.list)


#END OF SCRIPT----