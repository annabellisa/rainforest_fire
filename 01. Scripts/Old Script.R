#begining

#packages to load in----
library("dplyr")

library("corrplot")

library("ggpubr")

library("AICcmodavg")

library("ggiraphExtra")

library("sjPlot")

library("patchwork")

library("MASS")

install.packages("MASS")


citation("AICcmodavg")

#raw image data----
raw1<-read.csv("Image_Final.CSV")
head(raw1,2); dim(raw1)

#no false triggers----
nofalse<- filter(raw1, Pres_Abs==1)


#one line per animal Id----
totalanim<-nofalse %>% distinct(Animal_ID, .keep_all = TRUE)
head(totalanim,2);dim(totalanim)


#species list----
specieslist<-totalanim %>% distinct(Identification, .keep_all=TRUE) %>% 
  select(Habitat_Fire,Identification,Taxinomic_Rank)

head(specieslist,2);dim(specieslist)

specieslist$Identification
table(totalanim$Identification)



#exporting to excel then reimport new species data----

install.packages("writexl")

library(writexl)

write_xlsx(specieslist, "~/2021/Honours project/Data/Analysis/specieslist.xlsx")

bird_mammal<-read.csv("birdvsmammal.csv",header=T, check.names=FALSE)
colnames(bird_mammal)[1]<-c("Habitat_Fire") 
bird_mammal2<-bird_mammal %>% select(-Habitat_Fire)

#combining total animm with bird vs mammal colomn ----

raw_class<-totalanim %>% full_join(bird_mammal2)

#bringing in site data ----

rawsite<-read.csv("Site_Final.csv",header=T, check.names=FALSE)
colnames(rawsite)[1]<-c("Site")

#joining animal data and site data----

alldata<-raw_class %>% full_join(rawsite)
head(alldata,2); dim(alldata)

#culling down all data table----

culleddata<-alldata %>% select(-Latitude,-Longitude,-Vertical_Cam,-Horizontal_Cam,-Confidence,-Photo_Series,-Trigger_Group,-Trigger,-File,-Cam_ID, -Pres_Abs,-No_Animals,-Animal_Image,-Fire_habitat_catergory,-Height_CM,-To_Slope_CM, -Image_ID, -Cam_Type, -Orientation)


#sorting into just mammals ----
DFall<-culleddata %>% 
  select(Identification,Habitat_Fire,Class)

mammals<-DFall %>% distinct(Identification, Habitat_Fire, .keep_all = TRUE)%>% filter(Class=="Mammal")

head(mammals,2); dim(mammals)


#ggpubr----

install.packages("ggpubr")
library(ggpubr)

#quick note ----
#apperently this function(group_by(...) %>% count()) tallies up the observations based on a grouping variable

library(tidyr)


#with Annabel creating data to model with----

head(mammalsalldata,2); dim(mammalsalldata)
MAD<-data.frame(mammalsalldata)

head(MAD,2); dim(MAD)
unique(MAD$genus)

MAD$genus<-NA
us_ind<-unlist(gregexpr("_",MAD$Identification))
head(us_ind,10)
MAD$genus[which(us_ind<0)]<-MAD$Identification [which(us_ind<0)]
MAD$genus[which(us_ind>0)]<-substr(MAD$Identification,1,us_ind-1)[which(us_ind>1)]

MAD$species<-NA
MAD$species[which(us_ind>0)]<-substr(MAD$Identification,us_ind+1,nchar(MAD$Identification))[which(us_ind>1)]

MADS<-MAD[which(MAD$Taxinomic_Rank=="Species"),]
head(MADS,2);dim(MADS)

?which

#MAD= taxinomic richness and MADS = species richness

head(rawsite,2);dim(rawsite)
head(MAD,2); dim(MAD)

aggregate(genus~Site,data=MAD,FUN = function(x)length(unique(x)))

SP_sum<-data.frame(Site=aggregate(genus~Site,data=MAD,FUN = function(x)length(unique(x)))$Site,genus_rich=aggregate(genus~Site,data=MAD,FUN = function(x)length(unique(x)))$genus)
head(SP_sum,2);dim(SP_sum)

speciesdata<-aggregate(species~Site,data=MAD,FUN = function(x)length(unique(x)))
head(speciesdata,4);dim(speciesdata)

SP_sum<-merge(SP_sum,speciesdata,by="Site",all.x=T,all.y=F)
SP_sum$species[which(is.na(SP_sum$species))]<-0


modelsite<-merge(rawsite,SP_sum,by="Site",all.x=T,all.y=F) %>% 
  arrange(desc(Fire_habitat_catergory))

head(modelsite,3);dim(modelsite)



#Box Plots----


dev.new(height=6,width=6,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
ggboxplot(modelsite,x= "Fire_habitat_catergory" ,y="species", xlab="Habitat Fire Category",ylab="Species Richness")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  scale_x_discrete(limits = c("Unburnt_Rainforest", "Burnt_Sclerophyll", "Burnt_Rainforest"),
                   labels = c("Unburnt Rainforest", "Burnt Sclerophyll", 
                              "Burnt Rainforest"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))


#TOTAL count of individual species table----
table(mammalsalldata$Identification)

Species<-c('Aepyprymnus_rufescens','Mus_musculus','Perameles_nasuta','Phascolarctos_cinereus','Pseudomys_gracilicaudatus','Rattus_fuscipes','Rattus_lutreolus','Rattus_norvegicus','Rattus_rattus','Sminthopsis_murina','Sus_scrofa','Thylogale_stigmatica','Thylogale_thetis','Trichosurus_caninus')

Total<-c(1,10,11,1,82,161,17,59,47,2,4,256,79,39)

SpeciesTotal<-data.frame(Species,Total) #highest numbers for bush rat and red legged pademelon

capture.output(SpeciesTotal,file="Total species caught table")


#Step 1 Model selection (Best Burn Catergory)----

summary(modnull<-glm.nb(species~1,data=modelsite))

summary(modhab<-glm.nb(species~Fire_habitat_catergory,data=modelsite))

summary(modunburn<-glm.nb(species~Surrounding_Unburnt,data=modelsite))

summary(modlow<-glm.nb(species~Surrounding_Low,data=modelsite))

summary(modmoder<-glm.nb(species~Surrounding_Moderate,data=modelsite))

summary(modhigh<-glm.nb(species~Surrounding_High,data=modelsite))

summary(modextr<-glm.nb(species~Surrounding_Extreme,data=modelsite))

summary(modtotal<-glm.nb(species~Surrounding_Burn,data=modelsite))


modburnlist<-list("null"=modnull,"habitat"=modhab,"unburnt"=modunburn,"low"=modlow,"moderate"=modmoder,"high"=modhigh,"extreme"=modextr,"total"=modtotal)

aictab(modburnlist)

capture.output(aictab(modburnlist),file="AICC Table - Burn")


#Corelation analysis----

modelsite$habitat<-substr(modelsite$Fire_habitat_catergory,unlist(gregexpr("_",modelsite$Fire_habitat_catergory))+1,nchar(as.character(modelsite$Fire_habitat_catergory)))

cordata<-modelsite %>% select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Underbrush,Midstory,Canopy_Cover,Surrounding_Burn)

cordata<-modelsite[,c("Elevation","Dis_to_road","Dis_to_path", "Dis_to_rainforest_boundry","Underbrush","Midstory","Canopy_Cover","Surrounding_Burn")]

head(modelsite,3);dim(modelsite)

str(cordata)
mapping<-c("Open" = 0, "Sparse" = 1, "Moderate" = 2,"Dense" = 3)
cordata$Underbrush <- mapping[cordata$Underbrush]
cordata$Midstory <- mapping[cordata$Midstory]
cordata$Canopy_Cover <- mapping[cordata$Canopy_Cover]

cr<-cor(cordata, method="pearson")
head(cr,3);dim(cr)
summary(cr)
colnames(cr)<-c("Elevation" ,"Distance To Road" ,"Distance To Path", "Dist. To Rainforest","Understory","Midstory","Canopy Cover", "Surrounding Burn")
rownames(cr)<-c("Elevation" ,"Distance To Road" ,"Distance To Path", "Dist. To Rainforest","Understory","Midstory","Canopy Cover", "Surrounding Burn")

corrplot(cr,method="color")
corrplot(cr,method="number")
corrplot(cr,method="pie")

corrplot(cr,method="color",order="AOE",tl.col="Black")

?corrplot

corrplot(cr,method="color",  
         type="upper", order="hclust", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE)

#Step 2 (addative vs interaction models)----

modelsite$Underbrush <- mapping[modelsite$Underbrush]
modelsite$Midstory <- mapping[modelsite$Midstory]
modelsite$Canopy_Cover <- mapping[modelsite$Canopy_Cover]

modelsite<-modelsite2
summary(mod1a<-glm.nb(species~Surrounding_Burn*habitat+Elevation,data=modelsite))



summary(mod1a<-glm.nb(species~Surrounding_Burn+Elevation,data=modelsite))
summary(mod1b<-glm.nb(species~Surrounding_Burn*Elevation,data=modelsite))
AICc(mod1a);AICc(mod1b)

summary(mod2a<-glm.nb(species~Surrounding_Burn+Dis_to_road,data=modelsite))
summary(mod2b<-glm.nb(species~Surrounding_Burn*Dis_to_road,data=modelsite))
AICc(mod2a);AICc(mod2b)

summary(mod3a<-glm.nb(species~Surrounding_Burn+Dis_to_path,data=modelsite))
summary(mod3b<-glm.nb(species~Surrounding_Burn*Dis_to_path,data=modelsite))
AICc(mod3a);AICc(mod3b)

summary(mod4a<-glm.nb(species~Surrounding_Burn+Dis_to_rainforest_boundry,data=modelsite))
summary(mod4b<-glm.nb(species~Surrounding_Burn*Dis_to_rainforest_boundry,data=modelsite))
AICc(mod4a);AICc(mod4b)

summary(mod5a<-glm.nb(species~Surrounding_Burn+Underbrush,data=modelsite))
summary(mod5b<-glm.nb(species~Surrounding_Burn*Underbrush,data=modelsite))
AICc(mod5a);AICc(mod5b)

summary(mod6a<-glm.nb(species~Surrounding_Burn+Midstory,data=modelsite))
summary(mod6b<-glm.nb(species~Surrounding_Burn*Midstory,data=modelsite))
AICc(mod6a);AICc(mod6b)

summary(mod7a<-glm.nb(species~Surrounding_Burn+habitat,data=modelsite))
summary(mod7b<-glm.nb(species~Surrounding_Burn*habitat,data=modelsite))
AICc(mod7a);AICc(mod7b)

summary(mod8a<-glm.nb(species~Surrounding_Burn+Canopy_Cover,data=modelsite))
summary(mod8b<-glm.nb(species~Surrounding_Burn*Canopy_Cover,data=modelsite))
AICc(mod8a);AICc(mod8b)



#Step 3 (Final comparison of models)----

modfinal<-list("Null Model"=modnull,"Total Surrounding Burn"=modtotal,"Total Burn + Elevation"=mod1a,"Total Burn + Distance To Road"=mod2a,"Total Burn x Distance To Path"=mod3b,"Total Burn x Distance To Rainforest Boundary"=mod4b,"Total Burn + Level Of Underbrush"=mod5a,"Total Burn x Level Of Midstory"=mod6b,"Total Burn + Habitat Type"=mod7a, "Total Burn + Levle Of Canopy Cover"=mod8a)

modaictable<-aictab(modfinal)

capture.output(aictab(modfinal),file="AICC Table")


#Step 4 (plot model)----
totalburn<-data.frame(Surrounding_Burn=seq(0,100,by=1))

head(totalburn,3);dim(totalburn);tail(totalburn,3)

pr2<-predict(object = modtotal,newdata = totalburn,se.fit=T,type="response")

summary(modelsite$Surrounding_Burn)

totalburn1<-data.frame(totalburn,fit=pr2$fit,se=pr2$se.fit)
totalburn1$lci<-totalburn1$fit-(1.96*totalburn1$se)
head(totalburn1,3);dim(totalburn1);tail(totalburn1,3)
totalburn1$uci<-totalburn1$fit+(1.96*totalburn1$se)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(totalburn1$Surrounding_Burn,totalburn1$fit,type="l",ylim=c(min(0),max(c(totalburn1$uci,modelsite$species))),xlab="Proportion Of Surrounding Area Burnt",ylab="Species Richness",las=1)+ 
  lines(totalburn1$Surrounding_Burn,totalburn1$lci,lty=2)+
  lines(totalburn1$Surrounding_Burn,totalburn1$uci,lty=2)
points(modelsite$Surrounding_Burn,modelsite$species,pch=20)

head(modelsite,3);dim(modelsite)
?plot

#path----

summary(mod3b)

hist(modelsite$species,breaks = 10)


nd1<-data.frame(Dis_to_path=c(rep(as.numeric(summary(modelsite$Dis_to_path)[2]),101),rep(as.numeric(summary(modelsite$Dis_to_path)[5]),101)),Surrounding_Burn=seq(0,100,by=1))

head(nd1,3);dim(nd1);tail(nd1,3)

pr1<-predict(object = mod3b,newdata = nd1,se.fit=T, type="response")

nd2<-data.frame(nd1,fit=pr1$fit,se=pr1$se.fit)
nd2$lci<-nd2$fit-(1.96*nd2$se)
head(nd2,3);dim(nd2);tail(nd2,3)
nd2$uci<-nd2$fit+(1.96*nd2$se)

str(nd2)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))

plot(nd2$Surrounding_Burn[nd2$Dis_to_path==17.265],nd2$fit[nd2$Dis_to_path==17.265],type="l",ylim=c(min(nd2$lci),max(nd2$uci)),xlab="Proportion Of surrounding area burnt",ylab="Predicted species richness",las=1)+
  lines(nd2$Surrounding_Burn[nd2$Dis_to_path==17.265],nd2$lci[nd2$Dis_to_path==17.265],lty=2)+
  lines(nd2$Surrounding_Burn[nd2$Dis_to_path==17.265],nd2$uci[nd2$Dis_to_path==17.265],lty=2)

lines(nd2$Surrounding_Burn[102:202],nd2$fit[102:202],col="red")
lines(nd2$Surrounding_Burn[102:202],nd2$lci[102:202],col="red", lty=2)
lines(nd2$Surrounding_Burn[102:202],nd2$uci[102:202],col="red",lty=2)
legend("topright",legend=c("1st Q","3rd Q"),col=c("black","red"),lty=1)


range(modelsite$Surrounding_Burn)
range(modelsite$Dis_to_path)
summary(modelsite$Dis_to_path)

nd2$Dis_to_path==6.95

#midstory----

midstory<-data.frame(Midstory=c(rep(0,101),rep(1,101),rep(2,101),rep(3,101)),Surrounding_Burn=seq(0,100,by=1))

head(midstory,3);dim(midstory);tail(midstory,3)

pr3<-predict(object = mod6b,newdata = midstory,se.fit=T)

midstory1<-data.frame(midstory,fit=pr3$fit,se=pr3$se.fit)
midstory1$lci<-midstory1$fit-(1.96*midstory1$se)
head(midstory1,3);dim(midstory1);tail(midstory1,3)
midstory1$uci<-midstory1$fit+(1.96*midstory1$se)



plot(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$fit[midstory1$Midstory==0],type="l",ylim=c(min(midstory1$lci),max(midstory1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Predicted Species Richness",las=1)

lines(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$lci[midstory1$Midstory==0],lty=2)
lines(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$uci[midstory1$Midstory==0],lty=2)

lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$fit[midstory1$Midstory==3],col="red")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$lci[midstory1$Midstory==3],lty=2,col="red")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$uci[midstory1$Midstory==3],lty=2,col="red")

legend("bottomleft",legend=c("Open","Dense"),col=c("black","red"),lty=1)


plot(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$fit[midstory1$Midstory==0],type="l",ylim=c(min(midstory1$lci),max(midstory1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Predicted Species Richness",las=1)

lines(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$lci[midstory1$Midstory==0],lty=2)
lines(midstory1$Surrounding_Burn[midstory1$Midstory==0],midstory1$uci[midstory1$Midstory==0],lty=2)

lines(midstory1$Surrounding_Burn[midstory1$Midstory==1],midstory1$fit[midstory1$Midstory==1],col="blue")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==1],midstory1$lci[midstory1$Midstory==1],lty=2,col="blue")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==1],midstory1$uci[midstory1$Midstory==1],lty=2,col="blue")

lines(midstory1$Surrounding_Burn[midstory1$Midstory==2],midstory1$fit[midstory1$Midstory==2],col="purple")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==2],midstory1$lci[midstory1$Midstory==2],lty=2,col="purple")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==2],midstory1$uci[midstory1$Midstory==2],lty=2,col="purple")

lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$fit[midstory1$Midstory==3],col="red")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$lci[midstory1$Midstory==3],lty=2,col="red")
lines(midstory1$Surrounding_Burn[midstory1$Midstory==3],midstory1$uci[midstory1$Midstory==3],lty=2,col="red")

legend("bottomleft",legend=c("Open","Sparse","Moderate","Dense"),col=c("black","blue","purple","red"),lty=1)


#PCA----

modelspecies1<-modelspecies
summary(modelsite$Surrounding_Burn)


head(modelspecies,3);dim(modelspecies)


?princomp

pca1<-princomp(~.,data=modelspecies1[,which(colnames(modelspecies1)=="A_ruf"):ncol(modelspecies1),],cor=T)

pca2<-princomp(~Surrounding_Burn,data=modelspecies1[,c(which(colnames(modelspecies1)=="Surrounding_Burn"),which(colnames(modelspecies1)=="A_ruf"):ncol(modelspecies1)),],cor=T)

plot(pca1,cex.axis=0.6,las=1)
summary(pca1)
str(pca1)
str(summary(pca1))
summary(pca1)$importance[2,]
PoV<-pca1$sdev^2/sum(pca1$sdev^2)
PoV1<-as.data.frame(PoV)
PoV1$Comp<-c(1:14)
structure(PoV)
head(pca1);dim(pca1)
pca1$loadings

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
barplot(PoV1$PoV,ylab="Propotion Varience Explained",ylim = c(0.00,0.30),names.arg=PoV1$Comp)

dev.new(height=9,width=9,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,3,3),mfrow=c(2,2))
biplot(pca1, xlab="Component 1",ylab="Component 2",ylim=c(-0.6,0.4),xlim=c(-0.2,0.7))

plot(modelspecies$pca.comp1,modelspecies$pca.comp2,col=as.factor(modelspecies$Fire_habitat_catergory),pch=20, xlab="Component 1",ylab="Component 2")

plot(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$fit[road1$Dis_to_road==72.58],type="l",ylim=c(min(road1$lci),max(road1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Principal Component 1",las=1)+
  lines(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$lci[road1$Dis_to_road==72.58],lty=2)+
  lines(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$uci[road1$Dis_to_road==72.58],lty=2)+
  lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$fit[road1$Dis_to_road==1183.36],col="red")+
  lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$lci[road1$Dis_to_road==1183.36],lty=2,col="red")+
  lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$uci[road1$Dis_to_road==1183.36],lty=2,col="red")+
  legend("bottomright",legend=c("Min","Max"),col=c("black","red"),lty=1)




?biplot



head(modelspecies,3);dim(modelspecies)
modelspecies$pca.comp1<-pca1$scores[,1]
modelspecies$pca.comp2<-pca1$scores[,2]

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(modelspecies$pca.comp1,modelspecies$pca.comp2,col=as.factor(modelspecies$Fire_habitat_catergory),pch=20, xlab="Component 1",ylab="Component 2")

legend("bottomright",legend=levels(as.factor(modelspecies$Fire_habitat_catergory)),col=1:3,pch=20)

str(modelspecies)


#comp1 models----

modelspecies$Underbrush <- mapping[modelspecies$Underbrush]
modelspecies$Midstory <- mapping[modelspecies$Midstory]
modelspecies$Canopy_Cover <- mapping[modelspecies$Canopy_Cover]

summary(comp1null<-lm(pca.comp1~1,data=modelspecies))

summary(comp1a<-lm(pca.comp1~Surrounding_Burn+Elevation,data=modelspecies))
summary(comp1b<-lm(pca.comp1~Surrounding_Burn*Elevation,data=modelspecies))
AICc(comp1a);AICc(comp1b)

summary(comp1c<-lm(pca.comp1~Surrounding_Burn+Dis_to_road,data=modelspecies))
summary(comp1d<-lm(pca.comp1~Surrounding_Burn*Dis_to_road,data=modelspecies))
AICc(comp1c);AICc(comp1d)

summary(comp1e<-lm(pca.comp1~Surrounding_Burn+Dis_to_path,data=modelspecies))
summary(comp1f<-lm(pca.comp1~Surrounding_Burn*Dis_to_path,data=modelspecies))
AICc(comp1e);AICc(comp1f)

summary(comp1g<-lm(pca.comp1~Surrounding_Burn+Dis_to_rainforest_boundry,data=modelspecies))
summary(comp1h<-lm(pca.comp1~Surrounding_Burn*Dis_to_rainforest_boundry,data=modelspecies))
AICc(comp1g);AICc(comp1h)

summary(comp1i<-lm(pca.comp1~Surrounding_Burn+Underbrush,data=modelspecies))
summary(comp1j<-lm(pca.comp1~Surrounding_Burn*Underbrush,data=modelspecies))
AICc(comp1i);AICc(comp1j)

summary(comp1k<-lm(pca.comp1~Surrounding_Burn+Midstory,data=modelspecies))
summary(comp1l<-lm(pca.comp1~Surrounding_Burn*Midstory,data=modelspecies))
AICc(comp1k);AICc(comp1l)

summary(comp1m<-lm(pca.comp1~Surrounding_Burn+habitat,data=modelspecies))
summary(comp1n<-lm(pca.comp1~Surrounding_Burn*habitat,data=modelspecies))
AICc(comp1m);AICc(comp1n)

summary(comp1o<-lm(pca.comp1~Surrounding_Burn+Canopy_Cover,data=modelspecies))
summary(comp1p<-lm(pca.comp1~Surrounding_Burn*Canopy_Cover,data=modelspecies))
AICc(comp1o);AICc(comp1p)

trial<-lm(pca.comp1~Surrounding_Burn,data=modelspecies)


comp1final<-list("Null Model"=comp1null,"Total Burn + Elevation"=comp1a,"Total Burn X Distance To Road"=comp1d,"Total Burn x Distance To Path"=comp1f,"Total Burn + Distance To Rainforest Boundary"=comp1g,"Total Burn + Level Of Underbrush"=comp1i,"Total Burn + Level Of Midstory"=comp1k,"Total Burn + Habitat Type"=comp1m, "Total Burn + Levle Of Canopy Cover"=comp1o,"Total Burn"=trial)

comp1aictable<-aictab(comp1final)

capture.output(aictab(comp1final),file="AHHAHAHA")



summary(modelspecies$Dis_to_road)

road<-data.frame(Dis_to_road=c(rep(as.numeric(summary(modelsite$Dis_to_road)[1]),101),rep(as.numeric(summary(modelsite$Dis_to_road)[6]),101)),Surrounding_Burn=seq(0,100,by=1))

head(road,3);dim(road);tail(road,3)

pr4<-predict(object = comp1d,newdata = road,se.fit=T)

road1<-data.frame(road,fit=pr4$fit,se=pr4$se.fit)
road1$lci<-road1$fit-(1.96*road1$se)
head(road1,3);dim(road1);tail(road1,3)
road1$uci<-road1$fit+(1.96*road1$se)



plot(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$fit[road1$Dis_to_road==72.58],type="l",ylim=c(min(road1$lci),max(road1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Principal Component 1",las=1)

lines(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$lci[road1$Dis_to_road==72.58],lty=2)
lines(road1$Surrounding_Burn[road1$Dis_to_road==72.58],road1$uci[road1$Dis_to_road==72.58],lty=2)

lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$fit[road1$Dis_to_road==1183.36],col="red")
lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$lci[road1$Dis_to_road==1183.36],lty=2,col="red")
lines(road1$Surrounding_Burn[road1$Dis_to_road==1183.36],road1$uci[road1$Dis_to_road==1183.36],lty=2,col="red")

legend("bottomright",legend=c("Min","Max"),col=c("black","red"),lty=1)




totalb<-data.frame(Surrounding_Burn=seq(0,100,by=1))

head(totalb,3);dim(totalb);tail(totalb,3)

pr5<-predict(object = trial,newdata = totalb,se.fit=T)

totalb1<-data.frame(totalb,fit=pr5$fit,se=pr5$se.fit)
totalb1$lci<-totalb1$fit-(1.96*totalb1$se)
head(totalb1,3);dim(totalb1);tail(totalb1,3)
totalb1$uci<-totalb1$fit+(1.96*totalb1$se)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))

plot(totalb1$Surrounding_Burn,totalb1$fit,type="l",ylim=c(min(totalb1$lci),max(totalb1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Principal Component 1",las=1)

lines(totalb1$Surrounding_Burn,totalb1$lci,lty=2)
lines(totalb1$Surrounding_Burn,totalb1$uci,lty=2)




#comp2 models----

summary(comp2null<-lm(pca.comp2~1,data=modelspecies))

summary(comp2burn<-lm(pca.comp2~Surrounding_Burn,data=modelspecies))

summary(comp2a<-lm(pca.comp2~Surrounding_Burn+Elevation,data=modelspecies))
summary(comp2b<-lm(pca.comp2~Surrounding_Burn*Elevation,data=modelspecies))
AICc(comp2a);AICc(comp2b)

summary(comp2c<-lm(pca.comp2~Surrounding_Burn+Dis_to_road,data=modelspecies))
summary(comp2d<-lm(pca.comp2~Surrounding_Burn*Dis_to_road,data=modelspecies))
AICc(comp2c);AICc(comp2d)

summary(comp2e<-lm(pca.comp2~Surrounding_Burn+Dis_to_path,data=modelspecies))
summary(comp2f<-lm(pca.comp2~Surrounding_Burn*Dis_to_path,data=modelspecies))
AICc(comp2e);AICc(comp2f)

summary(comp2g<-lm(pca.comp2~Surrounding_Burn+Dis_to_rainforest_boundry,data=modelspecies))
summary(comp2h<-lm(pca.comp2~Surrounding_Burn*Dis_to_rainforest_boundry,data=modelspecies))
AICc(comp2g);AICc(comp2h)

summary(comp2i<-lm(pca.comp2~Surrounding_Burn+Underbrush,data=modelspecies))
summary(comp2j<-lm(pca.comp2~Surrounding_Burn*Underbrush,data=modelspecies))
AICc(comp2i);AICc(comp2j)

summary(comp2k<-lm(pca.comp2~Surrounding_Burn+Midstory,data=modelspecies))
summary(comp2l<-lm(pca.comp2~Surrounding_Burn*Midstory,data=modelspecies))
AICc(comp2k);AICc(comp2l)

summary(comp2m<-lm(pca.comp2~Surrounding_Burn+habitat,data=modelspecies))
summary(comp2n<-lm(pca.comp2~Surrounding_Burn*habitat,data=modelspecies))
AICc(comp2m);AICc(comp2n)

summary(comp2o<-lm(pca.comp2~Surrounding_Burn+Canopy_Cover,data=modelspecies))
summary(comp2p<-lm(pca.comp2~Surrounding_Burn*Canopy_Cover,data=modelspecies))
AICc(comp2o);AICc(comp2p)


comp2final<-list("Null Model"=comp2null,"Total Burn + Elevation"=comp2a,"Total Burn + Distance To Road"=comp2c,"Total Burn + Distance To Path"=comp2e,"Total Burn + Distance To Rainforest Boundary"=comp2g,"Total Burn + Level Of Underbrush"=comp2i,"Total Burn + Level Of Midstory"=comp2k,"Total Burn + Habitat Type"=comp2m, "Total Burn + Levle Of Canopy Cover"=comp2o, "Total Burn"=comp2burn)

comp2aictable<-aictab(comp2final)

capture.output(aictab(comp2final),file="AAAAHHHH2")


#total burn predict/plot
burn<-data.frame(Surrounding_Burn=seq(0,100,by=1))

head(burn,3);dim(burn);tail(burn,3)

pr6<-predict(object = comp2burn,newdata = burn,se.fit=T)

burn1<-data.frame(burn,fit=pr6$fit,se=pr6$se.fit)
burn1$lci<-burn1$fit-(1.96*burn1$se)
head(burn1,3);dim(burn1);tail(burn1,3)
burn1$uci<-burn1$fit+(1.96*burn1$se)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))

plot(burn1$Surrounding_Burn,burn1$fit,type="l",ylim=c(min(burn1$lci),max(burn1$uci)),xlab="Proportion Of Surrounding Area Burnt",ylab="Principal Component 2",las=1)

lines(burn1$Surrounding_Burn,burn1$lci,lty=2)
lines(burn1$Surrounding_Burn,burn1$uci,lty=2)



#species specific added to data---- 

Species0<-filter(mammalsalldata,Taxinomic_Rank=="Species")
head(Species0,4); dim(Species0)
table(Species0$Identification)


Species0$A_ruf<-NA
Species0$Mus_mus<-NA
Species0$P_nas<-NA
Species0$Ph_cin<-NA
Species0$Pseudo<-NA
Species0$Rat_fus<-NA
Species0$Rat_lut<-NA
Species0$Rat_nor<-NA
Species0$Rat_rat<-NA
Species0$S_mur<-NA
Species0$S_scr<-NA
Species0$Th_stig<-NA
Species0$Th_th<-NA
Species0$Tr_can<-NA
Species0$present<-1
Species0$absent<-0

Aruf<-unlist(gregexpr("Aepyprymnus_rufescens",Species0$Identification))
Musmus<-unlist(gregexpr("Mus_musculus",Species0$Identification))
Pnas<-unlist(gregexpr("Perameles_nasuta",Species0$Identification))
Phcin<-unlist(gregexpr("Phascolarctos_cinereus",Species0$Identification))
Pseudo<-unlist(gregexpr("Pseudomys_gracilicaudatus",Species0$Identification))
Ratfus<-unlist(gregexpr("Rattus_fuscipes",Species0$Identification))
Ratlut<-unlist(gregexpr("Rattus_lutreolus",Species0$Identification))
Ratnor<-unlist(gregexpr("Rattus_norvegicus",Species0$Identification))
Ratrat<-unlist(gregexpr("Rattus_rattus",Species0$Identification))
Smur<-unlist(gregexpr("Sminthopsis_murina",Species0$Identification))
Sscr<-unlist(gregexpr("Sus_scrofa",Species0$Identification))
Thstig<-unlist(gregexpr("Thylogale_stigmatica",Species0$Identification))
Thth<-unlist(gregexpr("Thylogale_thetis",Species0$Identification))
Trcan<-unlist(gregexpr("Trichosurus_caninus",Species0$Identification))

Species0$A_ruf[which(Aruf>0)]<-Species0$present[which(Aruf>0)]
Species0$A_ruf[which(Aruf<0)]<-Species0$absent[which(Aruf<0)]

Species0$Mus_mus[which(Musmus>0)]<-Species0$present[which(Musmus>0)]
Species0$Mus_mus[which(Musmus<0)]<-Species0$absent[which(Musmus<0)]

Species0$Ph_cin[which(Phcin>0)]<-Species0$present[which(Phcin>0)]
Species0$Ph_cin[which(Phcin<0)]<-Species0$absent[which(Phcin<0)]

Species0$P_nas[which(Pnas>0)]<-Species0$present[which(Pnas>0)]
Species0$P_nas[which(Pnas<0)]<-Species0$absent[which(Pnas<0)]

Species0$Pseudo[which(Pseudo>0)]<-Species0$present[which(Pseudo>0)]
Species0$Pseudo[which(Pseudo<0)]<-Species0$absent[which(Pseudo<0)]

Species0$Rat_fus[which(Ratfus>0)]<-Species0$present[which(Ratfus>0)]
Species0$Rat_fus[which(Ratfus<0)]<-Species0$absent[which(Ratfus<0)]

Species0$Rat_lut[which(Ratlut>0)]<-Species0$present[which(Ratlut>0)]
Species0$Rat_lut[which(Ratlut<0)]<-Species0$absent[which(Ratlut<0)]

Species0$Rat_nor[which(Ratnor>0)]<-Species0$present[which(Ratnor>0)]
Species0$Rat_nor[which(Ratnor<0)]<-Species0$absent[which(Ratnor<0)]

Species0$Rat_rat[which(Ratrat>0)]<-Species0$present[which(Ratrat>0)]
Species0$Rat_rat[which(Ratrat<0)]<-Species0$absent[which(Ratrat<0)]

Species0$S_mur[which(Smur>0)]<-Species0$present[which(Smur>0)]
Species0$S_mur[which(Smur<0)]<-Species0$absent[which(Smur<0)]

Species0$S_scr[which(Sscr>0)]<-Species0$present[which(Sscr>0)]
Species0$S_scr[which(Sscr<0)]<-Species0$absent[which(Sscr<0)]

Species0$Th_stig[which(Thstig>0)]<-Species0$present [which(Thstig>0)]
Species0$Th_stig[which(Thstig<0)]<-Species0$absent [which(Thstig<0)]
Species0$Th_stig<-as.numeric(Species0$Th_stig)

Species0$Th_th[which(Thth>0)]<-Species0$present[which(Thth>0)]
Species0$Th_th[which(Thth<0)]<-Species0$absent[which(Thth<0)]

Species0$Tr_can[which(Trcan>0)]<-Species0$present[which(Trcan>0)]
Species0$Tr_can[which(Trcan<0)]<-Species0$absent[which(Trcan<0)]

str(Species0)

head(Species0,4);dim(Species0)

as.data.frame(Species0)

species2<-Species0 %>% select(Site,A_ruf,Mus_mus,P_nas,Ph_cin,Pseudo,Rat_fus,Rat_lut,Rat_nor,Rat_rat,S_mur,S_scr,Th_stig,Th_th,Tr_can)


species3<-aggregate(.~Site,species2,sum)


rawsite$habitat<-substr(rawsite$Fire_habitat_catergory,unlist(gregexpr("_",rawsite$Fire_habitat_catergory))+1,nchar(as.character(rawsite$Fire_habitat_catergory)))


modelspecies<-merge(rawsite,species3,by="Site",all.x =T,All.y=F) %>% 
  arrange(desc(Fire_habitat_catergory))

head(modelspecies,3);dim(modelspecies)

mapping<-c("Open" = 0, "Sparse" = 1, "Moderate" = 2,"Dense" = 3)
modelspecies$Underbrush <- mapping[modelspecies$Underbrush]
modelspecies$Midstory <- mapping[modelspecies$Midstory]
modelspecies$Canopy_Cover <- mapping[modelspecies$Canopy_Cover]

modelspecies[is.na(modelspecies)] = 0
#this makes the vertical cam hight and number into 0 to since they where NA




#pres abs----
presabs<-modelspecies


presabs$Mus_mus[presabs$Mus_mus>1]<-1
presabs$P_nas[presabs$P_nas>1]<-1
presabs$Pseudo[presabs$Pseudo>1]<-1
presabs$Rat_fus[presabs$Rat_fus>1]<-1
presabs$Rat_lut[presabs$Rat_lut>1]<-1
presabs$Rat_nor[presabs$Rat_nor>1]<-1
presabs$Rat_rat[presabs$Rat_rat>1]<-1
presabs$S_mur[presabs$S_mur>1]<-1
presabs$S_scr[presabs$S_scr>1]<-1
presabs$Th_stig[presabs$Th_stig>1]<-1
presabs$Th_th[presabs$Th_th>1]<-1
presabs$Tr_can[presabs$Tr_can>1]<-1

presabs$habitat<-substr(presabs$Fire_habitat_catergory,unlist(gregexpr("_",presabs$Fire_habitat_catergory))+1,nchar(as.character(presabs$Fire_habitat_catergory)))



#Bottom----



install.packages("corrplot")
install.packages("dplyr")
library ("corrplot")
library ("dplyr")


cordata<-select(modelsite,Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Underbrush,Midstory,Canopy_Cover)

str(cordata)

cr<-cor(cordata, method="pearson")

corrplot(cr,method="color")
corrplot(cr,method="number")


