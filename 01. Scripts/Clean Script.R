#Packages----

library("dplyr")
library("ggpubr")
library("patchwork")
library("AICcmodavg")
library("MASS")
library("corrplot")
library('lmerTest')

#Sorting/filtering Data----
#using dplyr package

rawdata<-read.csv("Raw.Image.Data.CSV")
names(rawdata)[1]<-"Image_ID"
head(rawdata,2); dim(rawdata)

nofalse<-filter(rawdata, Pres_Abs==1)
head(nofalse,2);dim(nofalse)

uniqueanimID<-nofalse %>% distinct(Animal_ID,.keep_all = TRUE)
head(uniqueanimID,2);dim(uniqueanimID)

onlyspecies<-filter(uniqueanimID,Taxinomic_Rank=="Species")
table(onlyspecies$Identification)
head(onlyspecies,3);dim(onlyspecies)

specieslist<-onlyspecies %>% distinct(Identification, .keep_all=TRUE) %>% 
  select(Identification)
head(specieslist,2);dim(specieslist)

mammals<-filter(onlyspecies,Bird_Mammal=="Mammal")
head(mammals,3);dim(mammals)
table(mammals$Identification)

#adding site data
rawsite<-read.csv("Site.Data.csv")
names(rawsite)[1]<-"Site"
head(rawsite,3);dim(rawsite)

image_site<-merge(rawsite,mammals,by="Site",all.x=T,all.y=F)
head(image_site,3);dim(image_site)

#Creating A Species Richness Data frame (modelrich)----

MAD<-data.frame(mammals)
head(MAD,2); dim(MAD)

us_ind<-unlist(gregexpr("_",MAD$Identification))
head(us_ind,10)

MAD$species<-NA
MAD$species[which(us_ind>0)]<-substr(MAD$Identification,us_ind+1,nchar(MAD$Identification))[which(us_ind>1)]

head(rawsite,2);dim(rawsite)
head(MAD,2); dim(MAD)

speciesdata<-aggregate(species~Site,data=MAD,FUN = function(x)length(unique(x)))
head(speciesdata,4);dim(speciesdata)

modelrich<-merge(rawsite,speciesdata,by="Site",all.x=T,all.y=F) %>% 
  arrange(desc(Fire_habitat_catergory))

head(modelrich,3);dim(modelrich)
str(modelrich)

#Richness at each site and orientation

culled<-image_site %>% select(-Latitude,-Longitude,-Vertical_Cam,-Horizontal_Cam,To_Slope_CM,-Height_CM,-Confidence,-Photo_Series,-Trigger_Group,-Trigger,-File,-Cam_ID, -Pres_Abs,-No_Animals,-Animal_Image,-Image_ID, -Cam_Type,-Old_ID,-Old_Rank,-Date,-Time,-Temp,-Taxinomic_Rank,-Bird_Mammal) %>% 
  na.omit()
head(culled,4);dim(culled)

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

richness3<-filter(culled, Orientation=="Horizontal")

#exploring total number of all animals caught at sites using orientation----
Abundances<-culled %>% group_by(Site,Orientation) %>% summarise(abundance=n())
summary(Abundances$abundance)

ggboxplot(Abundances,x= "Orientation" ,y="abundance",ylab="Total Species Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

abundances1<-filter(Abundances,Orientation=="Horizontal")
summary(abundances1$abundance)
abundances2<-filter(Abundances,Orientation=="Vertical")
summary(abundances2$abundance)

Abund<-culled %>% group_by(Site,Orientation,Identification) %>% summarise(abun=n())
head(Abund,3);dim(Abund)

df1<-Abund %>% filter(Identification=="Thylogale_stigmatica")

df2<-Abund %>% filter(Identification=="Thylogale_thetis")
head(df2,3);dim(df2)

df3<-Abund %>% filter(Identification=="Melomys_cervinipes")%>% 
  arrange(Orientation)
head(df3,3);dim(df3)

df4<-Abund %>% filter(Identification=="Sus_scrofa")
head(df4,3);dim(df4)

df5<-Abund %>% filter(Identification=="Perameles_nasuta") %>% 
  arrange(Orientation)
head(df5,3);dim(df5)

df6<-Abund %>% filter(Identification=="Rattus_fuscipes") %>% 
  arrange(Orientation)
head(df6,3);dim(df6)

df7<-Abund %>% filter(Identification=="Trichosurus_caninus")
head(df7,3);dim(df7)

df8<- Abund %>% filter(Identification=="Antechinus_stuartii")
head(df8,3);dim(df8)

df9<- Abund %>% filter(Identification=="Phascolarctos_cinereus")
head(df9,3);dim(df9)


#created graphic using patchwork package

(p1 | plot_spacer()|p2)/  
  (p3 | p4 | p5) /
  (p11 | p6 | p8) /
  (p9 | p10 | p7)

summary(lm(richness1$Identification~richness1$Orientation)

p1<-ggboxplot(data=richness1,x="Orientation" ,y="Identification",ylab="Species Richness",width = 0.5)
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=14))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(Abundances$abundance~Abundances$Orientation))

p2<-ggboxplot(Abundances,x= "Orientation" ,y="abundance",ylab="Species Abundance",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df1$abun~df1$Orientation))

p3<-ggboxplot(df1,x= "Orientation" ,y="abun",ylab="Red-Legged Pademelon",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

#no p value as only one orientation

p4<-ggboxplot(df2,x= "Orientation" ,y="abun",ylab="Red-Necked Pademelon",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df3$abun~df3$Orientation))

p5<-ggboxplot(df3,x= "Orientation" ,y="abun",ylab="Fawn-footed Melomys",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

#no p value as only one orientation

p6<-ggboxplot(df4,x= "Orientation" ,y="abun",ylab="Feral Pig",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df5$abun~df5$Orientation))

p7<-ggboxplot(df5,x= "Orientation" ,y="abun",ylab="Long Nosed Bandicoot",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df6$abun~df6$Orientation))

p8<-ggboxplot(df6,x= "Orientation" ,y="abun",ylab="Bush Rat",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df7$abun~df7$Orientation))

p9<-ggboxplot(df7,x= "Orientation" ,y="abun",ylab="Short Eared Possum",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

summary(lm(df8$abun~df8$Orientation))

p10<-ggboxplot(df8,x= "Orientation" ,y="abun",ylab="Brown Antechinus",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

#no p value as only one orientation

p11<-ggboxplot(df9,x= "Orientation" ,y="abun",ylab="Koala",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

#Preparing for species richness modeling----

mapping<-c("Open" = 0, "Sparse" = 1, "Moderate" = 2,"Dense" = 3)
modelrich$Underbrush <- mapping[modelrich$Underbrush]
modelrich$Midstory <- mapping[modelrich$Midstory]
modelrich$Canopy_Cover <- mapping[modelrich$Canopy_Cover]
names(modelrich)[12]<-"Understory"

modelrich<-modelrich %>%  select(-Latitude,-Longitude,-Vertical_Cam,-Height_CM,-Horizontal_Cam,-To_Slope_CM)

mapping2<-c("Rainforest"="R","Sclerophyll"="S")
modelrich$Habitat<-mapping2[modelrich$Habitat]
modelrich$Habitat<-factor(modelrich$Habitat)
modelrich$Fire_habitat_catergory<-factor(modelrich$Fire_habitat_catergory)

modelrich[is.na(modelrich)]<-0
?is.na

#Species Richness Modelling = burn(+ OR *)habitat----

summary(modnull1<-glm.nb(species~1+Habitat,data=modelrich))
summary(modnull2<-glm.nb(species~1*Habitat,data=modelrich))
AICc(modnull1);AICc(modnull2)

summary(modunburnt1<-glm.nb(species~Surrounding_Unburnt+Habitat,data=modelrich))
summary(modunburnt2<-glm.nb(species~Surrounding_Unburnt*Habitat,data=modelrich))
AICc(modunburnt1);AICc(modunburnt2)

summary(modlow1<-glm.nb(species~Surrounding_Low+Habitat,data=modelrich))
summary(modlow2<-glm.nb(species~Surrounding_Low*Habitat,data=modelrich))
AICc(modlow1);AICc(modlow2)

summary(modmoderate1<-glm.nb(species~Surrounding_Moderate+Habitat,data=modelrich))
summary(modmoderate2<-glm.nb(species~Surrounding_Moderate*Habitat,data=modelrich))
AICc(modmoderate1);AICc(modmoderate2)

summary(modhigh1<-glm.nb(species~Surrounding_High+Habitat,data=modelrich))
summary(modhigh2<-glm.nb(species~Surrounding_High*Habitat,data=modelrich))
AICc(modhigh1);AICc(modhigh2)

summary(modextreme1<-glm.nb(species~Surrounding_Extreme+Habitat,data=modelrich))
summary(modextreme2<-glm.nb(species~Surrounding_Extreme*Habitat,data=modelrich))
AICc(modextreme1);AICc(modextreme2)

summary(modtotal1<-glm.nb(species~Surrounding_Burn+Habitat,data=modelrich))
summary(modtotal2<-glm.nb(species~Surrounding_Burn*Habitat,data=modelrich))
AICc(modtotal1);AICc(modtotal2)

burnlist<-list("Null"=modnull2,"Unburnt"=modunburnt1,"Low"=modlow1,"Moderate"=modmoderate1,"High"=modhigh1,"Extreme"=modextreme1,"Total Burn"=modtotal1)

aictab(burnlist)

#surrounding burn chosen to continue as its more applicable for replication then high intensity burns

#Correlation analysis --> including surrounding burn----

cordata1<-modelrich %>% dplyr::select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Understory,Midstory,Canopy_Cover,Habitat,Surrounding_Burn)

mapping3<-c("R"=0,"S"=1)
cordata1$Habitat<-mapping3[cordata1$Habitat]

cr1<-cor(cordata1, method="pearson")

corrplot(cr1,method="color",  
         type="upper", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE)

#Species richness modelling with co-varients----

plot(modelrich$Habitat,modelrich$Surrounding_Burn)

cor.test(as.numeric(modelrich$Fire_habitat_catergory),modelrich$Surrounding_Burn)
as.numeric(modelrich$Fire_habitat_catergory)
plot(modelrich$Fire_habitat_catergory,modelrich$Surrounding_Burn)


summary(mod1a<-glm.nb(species~Surrounding_Burn+Habitat+Elevation,data=modelrich))
summary(mod1b<-glm.nb(species~Surrounding_Burn+Habitat*Elevation,data=modelrich))
summary(mod1c<-glm.nb(species~Habitat+Surrounding_Burn*Elevation,data=modelrich))
AICc(mod1a);AICc(mod1b);AICc(mod1c)

summary(mod2a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_path,data=modelrich))
summary(mod2b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_path,data=modelrich))
summary(mod2c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_path,data=modelrich))
AICc(mod2a);AICc(mod2b);AICc(mod2c)

summary(mod3a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_rainforest_boundry,data=modelrich))
summary(mod3b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_rainforest_boundry,data=modelrich))
summary(mod3c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_rainforest_boundry,data=modelrich))
AICc(mod3a);AICc(mod3b);AICc(mod3c)

summary(mod4a<-glm.nb(species~Surrounding_Burn+Habitat+Understory,data=modelrich))
summary(mod4b<-glm.nb(species~Surrounding_Burn+Habitat*Understory,data=modelrich))
summary(mod4c<-glm.nb(species~Habitat+Surrounding_Burn*Understory,data=modelrich))
AICc(mod4a);AICc(mod4b);AICc(mod4c)

summary(mod5a<-glm.nb(species~Surrounding_Burn+Habitat+Midstory,data=modelrich))
summary(mod5b<-glm.nb(species~Surrounding_Burn+Habitat*Midstory,data=modelrich))
summary(mod5c<-glm.nb(species~Habitat+Surrounding_Burn*Midstory,data=modelrich))
AICc(mod5a);AICc(mod5b);AICc(mod5c)

summary(mod6a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_road,data=modelrich))
summary(mod6b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_road,data=modelrich))
summary(mod6c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_road,data=modelrich))
AICc(mod6a);AICc(mod6b);AICc(mod6c)

summary(modtotal3<-glm.nb(species~Surrounding_Burn+Habitat,data=modelrich))
summary(modnull3<-glm.nb(species~1*Habitat,data=modelrich))

mods<-list("Null"=modnull3,"total burn + Habitat"=modtotal3,"Elevation"=mod1a,"Dis to path"=mod2a,"Dis to Rainforest boundary"=mod3a,"Understory"=mod4a,"Midstory"=mod5a,"Dis to road"=mod6a)

aictab(mods)


#modelling without site 15 to see if the lack of a vertical camera on that site is skewing the data

no15<-modelrich[-c(19), ]

summary(modnull4<-glm.nb(species~1+Habitat,data=no15))
summary(modnull5<-glm.nb(species~1*Habitat,data=no15))
AICc(modnull4);AICc(modnull5)

summary(modtotal4<-glm.nb(species~Surrounding_Burn+Habitat,data=no15))
summary(modtotal5<-glm.nb(species~Surrounding_Burn*Habitat,data=no15))
AICc(modtotal4);AICc(modtotal5)

summary(mod7a<-glm.nb(species~Surrounding_Burn+Habitat+Elevation,data=no15))
summary(mod7b<-glm.nb(species~Surrounding_Burn+Habitat*Elevation,data=no15))
summary(mod7c<-glm.nb(species~Habitat+Surrounding_Burn*Elevation,data=no15))
AICc(mod7a);AICc(mod7b);AICc(mod7c)

summary(mod8a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_path,data=no15))
summary(mod8b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_path,data=no15))
summary(mod8c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_path,data=no15))
AICc(mod8a);AICc(mod8b);AICc(mod8c)

summary(mod9a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_rainforest_boundry,data=no15))
summary(mod9b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_rainforest_boundry,data=no15))
summary(mod9c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_rainforest_boundry,data=no15))
AICc(mod9a);AICc(mod9b);AICc(mod9c)

summary(mod10a<-glm.nb(species~Surrounding_Burn+Habitat+Understory,data=no15))
summary(mod10b<-glm.nb(species~Surrounding_Burn+Habitat*Understory,data=no15))
summary(mod10c<-glm.nb(species~Habitat+Surrounding_Burn*Understory,data=no15))
AICc(mod10a);AICc(mod10b);AICc(mod10c)

summary(mod11a<-glm.nb(species~Surrounding_Burn+Habitat+Midstory,data=no15))
summary(mod11b<-glm.nb(species~Surrounding_Burn+Habitat*Midstory,data=no15))
summary(mod11c<-glm.nb(species~Habitat+Surrounding_Burn*Midstory,data=no15))
AICc(mod11a);AICc(mod11b);AICc(mod11c)

summary(mod12a<-glm.nb(species~Surrounding_Burn+Habitat+Dis_to_road,data=no15))
summary(mod12b<-glm.nb(species~Surrounding_Burn+Habitat*Dis_to_road,data=no15))
summary(mod12c<-glm.nb(species~Habitat+Surrounding_Burn*Dis_to_road,data=no15))
AICc(mod12a);AICc(mod12b);AICc(mod12c)

nosite15<-list("Null"=modnull4,"total burn + Habitat"=modtotal3,"Elevation"=mod7a,"Dis to path"=mod8a,"Dis to Rainforest boundary"=mod9a,"Understory"=mod10a,"Midstory"=mod11a,"Dis to road"=mod12a)

aictab(nosite15)

#species richness modelling with fire habitat category----

cordata2<-modelrich %>% dplyr::select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Understory,Midstory,Canopy_Cover,Fire_habitat_catergory)
?select

mapping3<-c("R"=0,"S"=1)
cordata2$Fire_habitat_catergory<-mapping4[cordata2$Fire_habitat_catergory]

cr2<-cor(cordata2, method="pearson")

corrplot(cr2,method="color",  
         type="upper", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE)

summary(modnull1<-glm.nb(species~1,data=modelrich))
AICc(modnull1)

summary(mod0<-glm.nb(species~Fire_habitat_catergory,data=modelrich))

summary(mod1a<-glm.nb(species~Fire_habitat_catergory+Elevation,data=modelrich))
summary(mod1b<-glm.nb(species~Fire_habitat_catergory*Elevation,data=modelrich))
AICc(mod1a);AICc(mod1b)

summary(mod2a<-glm.nb(species~Fire_habitat_catergory+Dis_to_path,data=modelrich))
summary(mod2b<-glm.nb(species~Fire_habitat_catergory*Dis_to_path,data=modelrich))
AICc(mod2a);AICc(mod2b)

summary(mod3a<-glm.nb(species~Fire_habitat_catergory+Dis_to_rainforest_boundry,data=modelrich))
summary(mod3b<-glm.nb(species~Fire_habitat_catergory*Dis_to_rainforest_boundry,data=modelrich))
AICc(mod3a);AICc(mod3b)

summary(mod4a<-glm.nb(species~Fire_habitat_catergory+Understory,data=modelrich))
summary(mod4b<-glm.nb(species~Fire_habitat_catergory*Understory,data=modelrich))
AICc(mod4a);AICc(mod4b)

summary(mod5a<-glm.nb(species~Fire_habitat_catergory+Midstory,data=modelrich))
summary(mod5b<-glm.nb(species~Fire_habitat_catergory*Midstory,data=modelrich))
AICc(mod5a);AICc(mod5b)

summary(mod6a<-glm.nb(species~Fire_habitat_catergory+Dis_to_road,data=modelrich))
summary(mod6b<-glm.nb(species~Fire_habitat_catergory*Dis_to_road,data=modelrich))
AICc(mod6a);AICc(mod6b)

models<-list("Null"=modnull1,"Fire Habitat Catergory"=mod0,"Elevation"=mod1a,"Dis to path"=mod2a,"Dis to Rainforest boundary"=mod3a,"Understory"=mod4a,"Midstory"=mod5a,"Dis to road"=mod6a)

aictab(models)

lrt<-anova(modnull1,mod0,test="LRT")
pvalue<-round(lrt[2,8],3)

#Species richness model predictions----

predictions<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

pr<-predict(object=mod0,newdata=predictions,se.fit =T,type="response")

firehabitat<-data.frame(predictions,fit=pr$fit,se=pr$se.fit)
firehabitat$lci<-firehabitat$fit-(1.96*firehabitat$se)
firehabitat$uci<-firehabitat$fit+(1.96*firehabitat$se)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:3,firehabitat$fit,type="p", ylim=c(min(0),max(firehabitat$uci)),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Species Richness",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:3,labels=c("UB R.forest","B R.forest","B Sclero"))
arrows(x0=c(1,2,3), y0=firehabitat$lci,x1=c(1,2,3), y1=firehabitat$uci,angle=90,length=0.1, code=3, col="black", lwd=2)


#up to pca in main script

