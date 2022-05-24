#PACKAGES----

#testing git


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





install.packages("glmmADMB")

install.packages("R2admb")

install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")


#SORTING/FILTERING DATA----

rawdata<-read.csv("CSV_DATA.csv")
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
rawsite<-read.csv("Site.csv")
names(rawsite)[1]<-"Site"
head(rawsite,3);dim(rawsite)

image_site<-merge(rawsite,mammals,by="Site",all.x=T,all.y=F)
head(image_site,3);dim(image_site)


##Creating A Species Richness Data frame----

MAD<-data.frame(mammals)
head(MAD,2); dim(MAD)

us_ind<-unlist(gregexpr("_",MAD$Identification))
head(us_ind,10)

MAD$species<-NA
MAD$species[which(us_ind>0)]<-substr(MAD$Identification,us_ind+1,nchar(MAD$Identification))[which(us_ind>1)]

?gregexpr

head(rawsite,2);dim(rawsite)
head(MAD,2); dim(MAD)

speciesdata<-aggregate(species~Site,data=MAD,FUN = function(x)length(unique(x)))
head(speciesdata,4);dim(speciesdata)

modelrich<-merge(rawsite,speciesdata,by="Site",all.x=T,all.y=F) %>% 
  arrange(desc(Fire_habitat_catergory))

head(modelrich,3);dim(modelrich)
str(modelrich)

##summary stats ----

culled<-image_site %>% select(-Latitude,-Longitude,-Vertical_Cam,-Horizontal_Cam,To_Slope_CM,-Height_CM,-Confidence,-Photo_Series,-Trigger_Group,-Trigger,-File,-Cam_ID, -Pres_Abs,-No_Animals,-Animal_Image,-Image_ID, -Cam_Type,-Old_ID,-Old_Rank,-Date,-Time,-Temp,-Taxinomic_Rank,-Bird_Mammal) %>% 
  na.omit()
head(culled,4);dim(culled)

#richness at each site and orientation

richness<-data.frame(culled)

richness1<-aggregate(data=richness,
                     Identification~Site+Orientation,
                     function(x)length(unique(x)))
?aggregate

ggboxplot(richness1,x= "Orientation" ,y="Identification",ylab="Species Richness")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

summary(richness1$Identification)

richness2<-filter(richness1,Orientation=="Vertical")
summary(richness2$Identification)


#horizontal species list

list<-filter(culled, Orientation=="Vertical")
table(list$Identification)


##abundance of all animals caught at sites----

Abundances<-culled %>% group_by(Site,Orientation) %>% summarise(abundance=n())
summary(Abundances$abundance)
head(Abundances,3);dim(Abundances)

ggboxplot(Abundances,x= "Orientation" ,y="abundance",ylab="Total Species Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

abundances1<-filter(Abundances,Orientation=="Horizontal")
summary(abundances1$abundance)
head(abundances1,3);dim(abundances1)

abundances2<-filter(Abundances,Orientation=="Vertical")
summary(abundances2$abundance)
head(abundances2,3);dim(abundances2)

Abund<-culled %>% group_by(Site,Orientation,Identification) %>% summarise(abun=n())
head(Abund,3);dim(Abund)

#Species specific abundance boxplots

df1<-Abund %>% filter(Identification=="Thylogale_stigmatica")
head(df1,3);dim(df1)

ggboxplot(df1,x= "Orientation" ,y="abun",ylab="Red Legged Pademelon Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df2<-Abund %>% filter(Identification=="Thylogale_thetis")
head(df2,3);dim(df2)

ggboxplot(df2,x= "Orientation" ,y="abun",ylab="Red Necked Pademelon Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df3<-Abund %>% filter(Identification=="Melomys_cervinipes")%>% 
  arrange(Orientation)
head(df3,3);dim(df3)

ggboxplot(df3,x= "Orientation" ,y="abun",ylab="Fawn-footed Melomys Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df4<-Abund %>% filter(Identification=="Sus_scrofa")
head(df4,3);dim(df4)

ggboxplot(df4,x= "Orientation" ,y="abun",ylab="Feral Pig Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))


df5<-Abund %>% filter(Identification=="Perameles_nasuta") %>% 
  arrange(Orientation)
head(df5,3);dim(df5)

ggboxplot(df5,x= "Orientation" ,y="abun",ylab="Long Nosed Bandicoot Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df6<-Abund %>% filter(Identification=="Rattus_fuscipes") %>% 
  arrange(Orientation)
head(df6,3);dim(df6)

ggboxplot(df6,x= "Orientation" ,y="abun",ylab="Bush Rat Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df7<-Abund %>% filter(Identification=="Trichosurus_caninus")
head(df7,3);dim(df7)

ggboxplot(df7,x= "Orientation" ,y="abun",ylab="Short Eared Possum Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df8<- Abund %>% filter(Identification=="Antechinus_stuartii")
head(df8,3);dim(df8)

ggboxplot(df8,x= "Orientation" ,y="abun",ylab="Brown Antechinus Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

df9<- Abund %>% filter(Identification=="Phascolarctos_cinereus")
head(df9,3);dim(df9)

ggboxplot(df9,x= "Orientation" ,y="abun",ylab="Koala Abundance")+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))

##creating graphic----

dev.new(height=27,width=27,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,0,0))

(p1 | plot_spacer()|p2)/  
(p3 | p4 | p5) /
  (p11 | p6 | p8) /
  (p9 | p10 | p7)
  
?ggpar 
length(unique(richness1$Site))
length(unique(modelrich$Site))
range(modelrich$species,na.rm=T)
head(modelrich,3)
head(richness1);dim(richness1)

#boxplot(richness1$Orientation,richness1$Identification)

dev.new(height=22,width=27,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(4,3),mar=c(4,4,0,0))

t.test(richness1$Orientation,richness1$Identification)

?t.test
str(richness1)

summary(lm(richness1$Identification~richness1$Orientation))

p1<-ggboxplot(richness1,x= "Orientation" ,y="Identification",ylab="Species Richness",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=14))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(Abundances,3)

ab2<-as.data.frame(Abundances)
head(ab2,4);dim(ab2)

summary(lm(ab2$abundance~ab2$Orientation))

p2<-ggboxplot(Abundances,x= "Orientation" ,y="abundance",ylab="Species Abundance",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df1,4)
df1<-as.data.frame(df1)

summary(lm(df1$abun~df1$Orientation))

p3<-ggboxplot(df1,x= "Orientation" ,y="abun",ylab="Red-Legged Pademelon",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

p4<-ggboxplot(df2,x= "Orientation" ,y="abun",ylab="Red-Necked Pademelon",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df3,4)
df3<-as.data.frame(df3)

summary(lm(df3$abun~df3$Orientation))

p5<-ggboxplot(df3,x= "Orientation" ,y="abun",ylab="Fawn-footed Melomys",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

p6<-ggboxplot(df4,x= "Orientation" ,y="abun",ylab="Feral Pig",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df5,4)
df5<-as.data.frame(df5)

summary(lm(df5$abun~df5$Orientation))

p7<-ggboxplot(df5,x= "Orientation" ,y="abun",ylab="Long Nosed Bandicoot",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df6,4)
df6<-as.data.frame(df6)

summary(lm(df6$abun~df6$Orientation))

p8<-ggboxplot(df6,x= "Orientation" ,y="abun",ylab="Bush Rat",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df7,4)
df7<-as.data.frame(df7)

summary(lm(df7$abun~df7$Orientation))

p9<-ggboxplot(df7,x= "Orientation" ,y="abun",ylab="Short Eared Possum",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

head(df8,4)
df8<-as.data.frame(df8)

summary(lm(df8$abun~df8$Orientation))

p10<-ggboxplot(df8,x= "Orientation" ,y="abun",ylab="Brown Antechinus",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

p11<-ggboxplot(df9,x= "Orientation" ,y="abun",ylab="Koala",width = 0.5)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 15))

#SPECIES RICHNESS----

mapping<-c("Open" = 0, "Sparse" = 1, "Moderate" = 2,"Dense" = 3)
modelrich$Underbrush <- mapping[modelrich$Underbrush]
modelrich$Midstory <- mapping[modelrich$Midstory]
modelrich$Canopy_Cover <- mapping[modelrich$Canopy_Cover]
names(modelrich)[12]<-"Understory"

str(modelrich)
modelrich<-modelrich %>%  select(-Latitude,-Longitude,-Vertical_Cam,-Height_CM,-Horizontal_Cam,-To_Slope_CM)


mapping2<-c("Rainforest"="R","Sclerophyll"="S")
modelrich$Habitat<-mapping2[modelrich$Habitat]
modelrich$Habitat<-factor(modelrich$Habitat)

modelrich$Fire_habitat_catergory<-factor(modelrich$Fire_habitat_catergory)

modelrich[is.na(modelrich)]<-0
?is.na

##modeling --> burn+habitat+covarients----

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

head(modelrich,3);dim(modelrich)

##correlation----


head(modelrich,3);dim(modelrich)

cordata1<-modelrich %>% dplyr::select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Understory,Midstory,Canopy_Cover,Habitat,Surrounding_Burn)
?select

mapping3<-c("R"=0,"S"=1)
cordata1$Habitat<-mapping3[cordata1$Habitat]
head(cordata1,3);dim(cordata)


cr1<-cor(cordata1, method="pearson")

corrplot(cr1,method="color",  
         type="upper", 
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE)
head(modelrich,3);dim(modelrich)
modelrich$Fire_habitat_catergory<-factor(modelrich$Fire_habitat_catergory,levels=c("Unburnt_Rainforest","Burnt_Rainforest","Burnt_Sclerophyll"))

cordata2<-modelrich %>% dplyr::select(Habitat,Surrounding_Burn)

cordata2$Habitat<-mapping3[cordata2$Habitat]
head(cordata2,3);dim(cordata2)


cr2<-cor(cordata2, method="pearson")

corrplot(cr2,method="color",  
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         sig.level = 0.01, insig = "blank", 
         diag=FALSE)

#Modeling with co-varients

head(modelrich,3);dim(modelrich)
plot(modelrich$Habitat,modelrich$Surrounding_Burn)
str(modelrich)
cor.test(as.numeric(modelrich$Fire_habitat_catergory),modelrich$Surrounding_Burn)
as.numeric(modelrich$Fire_habitat_catergory)
plot(modelrich$Fire_habitat_catergory,modelrich$Surrounding_Burn)


summary(mod1a<-glm.nb(species~Surrounding_Burn+Habitat+Elevation,data=modelrich))
summary(mod1b<-glm.nb(species~Surrounding_Burn+Habitat*Elevation,data=modelrich))
summary(mod1c<-glm.nb(species~Habitat+Surrounding_Burn*Elevation,data=modelrich))
AICc(mod1a);AICc(mod1b);AICc(mod1c)
summary(mod1a<-glm.nb(species~Surrounding_Burn+Habitat+Elevation,data=modelrich))
#Habitat*Elevation = Habitat+Elevation+Habitat:Elevation

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

models<-list("Null"=modnull2,"total burn + Habitat"=modtotal1,"Elevation"=mod1a,"Dis to path"=mod2a,"Dis to Rainforest boundary"=mod3a,"Understory"=mod4a,"Midstory"=mod5a,"Dis to road"=mod6a)

aictab(models)

#modelling without site 15

no15<-modelrich[-c(19), ]

dim(modelrich)
dim(no15)
head(modelrich,3);dim(modelrich)

head(no15,3);dim(no15)

summary(modnull3<-glm.nb(species~1+Habitat,data=no15))
summary(modnull4<-glm.nb(species~1*Habitat,data=no15))
AICc(modnull3);AICc(modnull4)

summary(modtotal3<-glm.nb(species~Surrounding_Burn+Habitat,data=no15))
summary(modtotal4<-glm.nb(species~Surrounding_Burn*Habitat,data=no15))
AICc(modtotal3);AICc(modtotal4)

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
aictab(models)


##modelling with fire habitat category----

cordata1<-modelrich %>% dplyr::select(Elevation,Dis_to_road,Dis_to_path, Dis_to_rainforest_boundry,Understory,Midstory,Canopy_Cover,Fire_habitat_catergory)
?select

mapping3<-c("R"=0,"S"=1)
cordata1$Fire_habitat_catergory<-mapping4[cordata1$Fire_habitat_catergory]
head(cordata1,3);dim(cordata1)


cr1<-cor(cordata1, method="pearson")

corrplot(cr1,method="color",  
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

summary(mod0)$coefficients
aicc<-list("FHC"=mod0,"Null"=modnull1)
aictab(aicc)



##predictions----


predictions<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

pr<-predict(object=mod0,newdata=predictions,se.fit =T,type="response")

firehabitat<-data.frame(predictions,fit=pr$fit,se=pr$se.fit)
firehabitat$lci<-firehabitat$fit-(1.96*firehabitat$se)
firehabitat$uci<-firehabitat$fit+(1.96*firehabitat$se)

head(firehabitat,3);dim(firehabitat)

#mapping4<-c("Unburnt_Rainforest"=1,"Burnt_Rainforest"=2,"Burnt_Sclerophyll"=3)

#firehabitat$Fire_habitat_catergory<-mapping4[firehabitat$Fire_habitat_catergory]


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
  par(mar=c(4,4,1,1))
  plot(1:3,firehabitat$fit,type="p", ylim=c(min(0),max(firehabitat$uci)),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Species Richness",las=1,cex=3,pch=19,xaxt="n")
  axis(side=1,at=1:3,labels=c("UB R.forest","B R.forest","B Sclero"))
    arrows(x0=c(1,2,3), y0=firehabitat$lci,x1=c(1,2,3), y1=firehabitat$uci,angle=90,length=0.1, code=3, lwd=2)
str(firehabitat)
head(modelrich,3);dim(modelrich)
str(modelrich)

?arrows


#PCA----

head(modelrich,3);dim(modelrich)
head(MAD,3);dim(MAD)
table(list$Identification)

pcadata<-modelrich
xx<-MAD

count<-MAD %>% group_by(Site) %>% count(Identification)
count<-as.data.frame(count)

head(xx,5);dim(xx)

xx$antechinus<-ifelse(xx$Identification=="Antechinus_stuartii",1,0)
xx$melomys<-ifelse(xx$Identification=="Melomys_cervinipes",1,0)
xx$bandicoot<-ifelse(xx$Identification=="Perameles_nasuta",1,0)
xx$koala<-ifelse(xx$Identification=="Phascolarctos_cinereus",1,0)
xx$bushrat<-ifelse(xx$Identification=="Rattus_fuscipes",1,0)
xx$pig<-ifelse(xx$Identification=="Sus_scrofa",1,0)
xx$stigmatica<-ifelse(xx$Identification=="Thylogale_stigmatica",1,0)
xx$thetis<-ifelse(xx$Identification=="Thylogale_thetis",1,0)
xx$possum<-ifelse(xx$Identification=="Trichosurus_caninus",1,0)


antechinus<-aggregate(xx$antechinus,by=list(Site=xx$Site),FUN=sum)
melomys<-aggregate(xx$melomys,by=list(Site=xx$Site),FUN=sum)
bandicoot<-aggregate(xx$bandicoot,by=list(Site=xx$Site),FUN=sum)
koala<-aggregate(xx$koala,by=list(Site=xx$Site),FUN=sum)
bushrat<-aggregate(xx$bushrat,by=list(Site=xx$Site),FUN=sum)
pig<-aggregate(xx$pig,by=list(Site=xx$Site),FUN=sum)
stigmatica<-aggregate(xx$stigmatica,by=list(Site=xx$Site),FUN=sum)
thetis<-aggregate(xx$thetis,by=list(Site=xx$Site),FUN=sum)
possum<-aggregate(xx$possum,by=list(Site=xx$Site),FUN=sum)

head(pcadata,3);dim(pcadata)

pcadata<-full_join(pcadata,antechinus)
colnames(pcadata)[18]<-c("antechinus")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,melomys)
colnames(pcadata)[19]<-c("melomys")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,bandicoot)
colnames(pcadata)[20]<-c("bandicoot")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,koala)
colnames(pcadata)[21]<-c("koala")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,bushrat)
colnames(pcadata)[22]<-c("bushrat")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,pig)
colnames(pcadata)[23]<-c("pig")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,stigmatica)
colnames(pcadata)[24]<-c("stigmatica")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,thetis)
colnames(pcadata)[25]<-c("thetis")
pcadata[is.na(pcadata)] = 0

pcadata<-full_join(pcadata,possum)
colnames(pcadata)[26]<-c("possum")
pcadata[is.na(pcadata)] = 0


pca1<-princomp(~.,data=pcadata[,which(colnames(pcadata)=="antechinus"):which(colnames(pcadata)=="possum"),],cor=T)

pcadata2<-pcadata
pcadata2<-select(pcadata2,-pca.comp1,-pca.comp2)
pcadata2<- rename(pcadata2,A.sti=antechinus,M.cer=melomys,P.nas=bandicoot,P.cin=koala,R.fus=bushrat,S.scr=pig,T.sti=stigmatica,T.the=thetis,T.can=possum)

head(pcadata2,4)


pca2<-prcomp(~.,data=pcadata2[,which(colnames(pcadata2)=="A.sti"):which(colnames(pcadata2)=="T.can"),])

summary(pca1)
head(pca2,4)
str(pca1)
summary(pca2)

#this is how you get proportion of varience explained for princomp 
PoV<-pca1$sdev^2/sum(pca1$sdev^2)
PoV1<-as.data.frame(PoV)
PoV1$Comp<-c(1:9)
#this is how you get proportion of varience explained for prcomp 
PoV2<-summary(pca2)$importance[2,]
str(pca2)
pca2$x

pcadata$pca.comp1<-pca2$x[,1]
pcadata$pca.comp2<-pca2$x[,2]

pcadata2$pca.comp1<-pca2$x[,1]
pcadata2$pca.comp2<-pca2$x[,2]


summary(pca2)
names(pca2$center)


dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0))
plot(x=1:length(PoV2),y=PoV2,ylab="Propotion Varience Explained",xlab="Components",type="p")
lines(x=1:length(PoV2),y=PoV2)
mtext("(a)",3,0.7,F,0)


biplot(pca2,xlab="Component 1",ylab="Component 2",col=c("grey40","black"),var.axes=T,arrow.len=0.1,ylim=c(-0.3,0.6))
       
mtext("(b)",3,0.7,F,adj = 0,at=-200)

ylim=c(-0.3,0.6),xlim=c(-0.7,0.4)

str(pca2)
pcadata

pca2$x
dimnames(pca2$x)[[1]]
dimnames(pca2$rotation)[[1]]
dimnames(pca2$rotation)[[1]][5]<-'R.fuscipes'
dimnames(pca2$rotation)[[1]][7]<-'T.stigmatica'
dimnames(pca2$rotation)[[1]][1]<-"A.stu"
str(pca2)
?biplot

pca2
summary(pca2)
str(pca2)

dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0))

plot(x=1:length(PoV2),y=PoV2,ylab="Propotion Varience Explained",xlab="Components",type="p")
lines(x=1:length(PoV2),y=PoV2)
mtext("(a)",3,0.7,F,0)

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
  


 pc2rotaion<- pca2$rotation[,1:2]

shapes<-c(15,17,16)
shapes<-shapes[as.factor(pcadata$Fire_habitat_catergory)]
col.1<-c("grey20","grey40","grey60")
col.1<-col.1[as.factor(pcadata$Fire_habitat_catergory)]

ggplot_pca(pca2)

?ggplot



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

head(pcadata,4)
anova(comp2l)
comp2final<-list("Null"=comp2null,"FHC + Elevation"=comp2a,"FHC + Distance To Road"=comp2c,"FHC + Distance To Path"=comp2e,"FHC + Distance To Rainforest Boundary"=comp2g,"FHC + Level Of Underbrush"=comp2i,"FHC X Level Of Midstory"=comp2l,"FHC"=comp2fhc)

aictab(comp2final)

lrt1<-anova(comp2null,comp2l,test="LRT")
pvalue1<-round(lrt[2,8],3)

summary(comp2l)$coefficients
aicc1<-list("FHC"=comp2l,"Null"=comp2null)
aictab(aicc1)

###component 2 predictions----

predictions1<-data.frame(Fire_habitat_catergory=factor(c(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),4)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Midstory=c(0,0,0,1,1,1,2,2,2,3,3,3))


pr1<-predict(object=comp2l,newdata=predictions1,se.fit =T,type="response")

midstory<-data.frame(predictions1,fit=pr1$fit,se=pr1$se.fit)
midstory$lci<-midstory$fit-(1.96*midstory$se)
midstory$uci<-midstory$fit+(1.96*midstory$se)

head(midstory,3);dim(midstory)

###component 2 graphing----

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
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

legend("bottomright",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=c("grey20","grey40","grey60"),pch=c(15,17,16),pt.cex=2,title="Fire Habitat Catergory")

?legend
?arrows




?axis
str(midstory)
midstory<-arrange(midstory,Fire_habitat_catergory,Midstory)
?arrange
head(pcadata,3);dim(pcadata)


#SPECIES PRESENCE/ABSENCE DATA----

length(which(pcadata$bandicoot==0))/length(pcadata$bandicoot)
table(pcadata$bandicoot==0)[2]/sum(table(pcadata$bandicoot==0))

propzero<-data.frame(species=names(apply(pcadata[,which(colnames(pcadata)=="antechinus"):which(colnames(pcadata)=="possum")],2,FUN = function(x)table(x==0)[2]/sum(table(x==0)))),propzero=apply(pcadata[,which(colnames(pcadata)=="antechinus"):which(colnames(pcadata)=="possum")],2,FUN = function(x)table(x==0)[2]/sum(table(x==0))))
#only use species with prop zero less then 80

propzero$abundance<-apply(pcadata[,which(colnames(pcadata)=="antechinus"):which(colnames(pcadata)=="possum")],2,FUN = function(x)sum(x))
propzero$binomial<-ifelse(propzero$propzero<0.9,1,0)
binomsp<-propzero$species[propzero$binomial==1]


head(pcadata,3);dim(pcadata)
singledata<-pcadata
head(singledata,3);dim(singledata)

singledata<-select(singledata,-koala,-pig,-pca.comp1,-pca.comp2)


##Loop modelling----

AIC.out<-list()
best.out<-list()

for(i in 1:length(binomsp)){

sp.thisrun<-binomsp[i]
sp.datathisrun<-ifelse(singledata[,sp.thisrun]==0,0,1)

null<-glm(sp.datathisrun~1,data=singledata,family=binomial )
modnull.type<-c("null")

mod0<-glm(sp.datathisrun~Fire_habitat_catergory,data=singledata,family=binomial )
mod0.type<-c("--")

mod1a<-glm(sp.datathisrun~Fire_habitat_catergory+Elevation,data=singledata,family=binomial )
mod1b<-glm(sp.datathisrun~Fire_habitat_catergory*Elevation,data=singledata,family=binomial )
mod1.mods<-data.frame(mods=c("mod1a","mod1b"),AICc=c(AICc(mod1a),AICc(mod1b)))
mod1<-mod1.mods$mods[which(mod1.mods$AICc==min(mod1.mods$AICc))]
mod1<-get(mod1)
mod1.type<-ifelse(which(mod1.mods$AICc==min(mod1.mods$AICc))==1,"add","int")

mod2a<-glm(sp.datathisrun~Fire_habitat_catergory+Dis_to_path,data=singledata,family=binomial )
mod2b<-glm(sp.datathisrun~Fire_habitat_catergory*Dis_to_path,data=singledata,family=binomial )
mod2.mods<-data.frame(mods=c("mod2a","mod2b"),AICc=c(AICc(mod2a),AICc(mod2b)))
mod2<-mod1.mods$mods[which(mod2.mods$AICc==min(mod2.mods$AICc))]
mod2<-get(mod2)
mod2.type<-ifelse(which(mod2.mods$AICc==min(mod2.mods$AICc))==1,"add","int")

mod3a<-glm(sp.datathisrun~Fire_habitat_catergory+Dis_to_rainforest_boundry,data=singledata,family=binomial )
mod3b<-glm(sp.datathisrun~Fire_habitat_catergory*Dis_to_rainforest_boundry,data=singledata,family=binomial )
mod3.mods<-data.frame(mods=c("mod3a","mod3b"),AICc=c(AICc(mod3a),AICc(mod3b)))
mod3<-mod3.mods$mods[which(mod3.mods$AICc==min(mod3.mods$AICc))]
mod3<-get(mod3)
mod3.type<-ifelse(which(mod3.mods$AICc==min(mod3.mods$AICc))==1,"add","int")

mod4a<-glm(sp.datathisrun~Fire_habitat_catergory+Understory,data=singledata,family=binomial )
mod4b<-glm(sp.datathisrun~Fire_habitat_catergory*Understory,data=singledata,family=binomial )
mod4.mods<-data.frame(mods=c("mod4a","mod4b"),AICc=c(AICc(mod4a),AICc(mod4b)))
mod4<-mod4.mods$mods[which(mod4.mods$AICc==min(mod4.mods$AICc))]
mod4<-get(mod4)
mod4.type<-ifelse(which(mod4.mods$AICc==min(mod4.mods$AICc))==1,"add","int")

mod5a<-glm(sp.datathisrun~Fire_habitat_catergory+Midstory,data=singledata,family=binomial )
mod5b<-glm(sp.datathisrun~Fire_habitat_catergory*Midstory,data=singledata,family=binomial )
mod5.mods<-data.frame(mods=c("mod5a","mod5b"),AICc=c(AICc(mod5a),AICc(mod5b)))
mod5<-mod5.mods$mods[which(mod5.mods$AICc==min(mod5.mods$AICc))]
mod5<-get(mod5)
mod5.type<-ifelse(which(mod5.mods$AICc==min(mod5.mods$AICc))==1,"add","int")

mod6a<-glm(sp.datathisrun~Fire_habitat_catergory+Dis_to_road,data=singledata,family=binomial )
mod6b<-glm(sp.datathisrun~Fire_habitat_catergory*Dis_to_road,data=singledata,family=binomial )
mod6.mods<-data.frame(mods=c("mod6a","mod6b"),AICc=c(AICc(mod6a),AICc(mod6b)))
mod6<-mod6.mods$mods[which(mod6.mods$AICc==min(mod6.mods$AICc))]
mod6<-get(mod6)
mod6.type<-ifelse(which(mod6.mods$AICc==min(mod6.mods$AICc))==1,"add","int")

mod.sum<-data.frame(species=sp.thisrun,mods=c('mod0',"mod1","mod2","mod3","mod4","mod5","mod6","null"),type=c(mod0.type,mod1.type,mod2.type,mod3.type,mod4.type,mod5.type,mod6.type,modnull.type),AICc=c(AICc(mod0), AICc(mod1),AICc(mod2),AICc(mod3),AICc(mod4),AICc(mod5),AICc(mod6),AICc(null)))
mod.sum<-mod.sum[order(mod.sum$AICc),]

AIC.out[[i]]<-mod.sum
best.out[[i]]<-get(mod.sum$mods[mod.sum$AICc==min(mod.sum$AICc)])



} #close for


summary(best.out[[1]])

AIC.res<-do.call(rbind,AIC.out)
#stigmatica top model is null but mod 4 (understory) is within 2 AICc's of this

best.res<-do.call(rbind,best.out)


summary(best.out[[i]])



##pred loop----



?do.call
?get
best.out[[2]]

best.out[[i]]$call
best.out[[1]]$terms

str(best.res[[i,24]])

head(best.res,7);dim(best.res)
str(best.res)

str(best.out[[1]])
str(best.res)

best.out[[1]]$data
best.out[[1]]$model[[2]]
best.out[[1]]$model[[3]]





binomsp
best.out
AIC.out

best.out.1<-best.out[-c(4,5,6)]
binomsp1<-list("antechinus","melomys","bandicoot","possum")

pred<-list()

for(i in 1:length(binomsp1)){


col1<-c(list(colnames(best.out.1[[i]]$model))[[1]][2])
col2<-c(list(colnames(best.out.1[[i]]$model))[[1]][3])
col.names<-list(col1,col2)

pr<-data.frame(fact1=factor((c(rep(levels(best.out.1[[i]]$model[[2]]),4)))),fact2=c(rep(0,3),rep(1,3),rep(2,3),rep(3,3)))

colnames(pr)<-col.names

pr1<-predict(object=best.out.1[[i]],newdata=pr,se.fit =T,type="response")

pr2<-data.frame(pr,fit=pr1$fit,se=pr1$se.fit)
pr2$lci<-pr2$fit-(1.96*pr2$se)
pr2$uci<-pr2$fit+(1.96*pr2$se)


pred[[i]]<-pr2

} #close

##plotting preds from loop----


###antecinus----

pred[[1]]
summary(best.out[[1]])


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:4,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(pred[[1]]$lci),max(pred[[1]]$uci)),xlim=c(min(0.5),max(4.5)),xlab="Understory",ylab="Pres/Abs of Antechinus stuartii",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=0:3,xlab="Understory")
lines(1:4,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Unburnt_Rainforest"])
arrows(x0=1:4, y0=pred[[1]]$lci[pred[[1]]$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=pred[[1]]$uci[pred[[1]]$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="black", lwd=2)

points((1:4)-0.2,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=19, col="2")
lines((1:4)-0.2,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Burnt_Rainforest"],col="2")
arrows(x0=(1:4)-0.2, y0=pred[[1]]$lci[pred[[1]]$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=pred[[1]]$uci[pred[[1]]$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="2", lwd=2)

points((1:4)+0.2,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=19,col="3")
lines((1:4)+0.2,pred[[1]]$fit[pred[[1]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='3')
arrows((1:4)+0.2, y0=pred[[1]]$lci[pred[[1]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=pred[[1]]$uci[pred[[1]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="3", lwd=2)

legend("topleft",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=1:3,pch=19,pt.cex=2,title="Fire Habitat Catergory")


###melomys----
pred[[2]]
str(summary(best.out.1[[2]]))
model.matrix(best.out.1[[2]])
AIC.out

head(singledata)
table(singledata$Fire_habitat_catergory,singledata$melomys>0)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:4,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(pred[[2]]$lci),max(pred[[2]]$uci)),xlim=c(min(0.5),max(4.5)),xlab="Midstory",ylab="Pres/Abs of Melomys crevinipes",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=0:3,xlab="Midstory")
lines(1:4,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Unburnt_Rainforest"])
arrows(x0=1:4, y0=pred[[2]]$lci[pred[[2]]$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=pred[[2]]$uci[pred[[2]]$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="black", lwd=2)
 
points((1:4)-0.2,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=19, col="2")
lines((1:4)-0.2,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Burnt_Rainforest"],col="2")
arrows(x0=(1:4)-0.2, y0=pred[[2]]$lci[pred[[2]]$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=pred[[2]]$uci[pred[[2]]$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="2", lwd=2)

points((1:4)+0.2,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=19,col="3")
lines((1:4)+0.2,pred[[2]]$fit[pred[[2]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='3')
arrows((1:4)+0.2, y0=pred[[2]]$lci[pred[[2]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=pred[[2]]$uci[pred[[2]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="3", lwd=2)

legend("bottomright",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=1:3,pch=19,pt.cex=2,title="Fire Habitat Catergory")

###bandicoot ----

pred[[3]]

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:4,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(pred[[3]]$lci),max(pred[[3]]$uci)),xlim=c(min(0.5),max(4.5)),xlab="Understory",ylab="Pres/abs of Perameles nasuta",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=0:3,xlab="Understory")
lines(1:4,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Unburnt_Rainforest"])
arrows(x0=1:4, y0=pred[[3]]$lci[pred[[3]]$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=pred[[3]]$uci[pred[[3]]$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="black", lwd=2)

points((1:4)-0.2,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=19, col="2")
lines((1:4)-0.2,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Burnt_Rainforest"],col="2")
arrows(x0=(1:4)-0.2, y0=pred[[3]]$lci[pred[[3]]$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=pred[[3]]$uci[pred[[3]]$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="2", lwd=2)

points((1:4)+0.2,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=19,col="3")
lines((1:4)+0.2,pred[[3]]$fit[pred[[3]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='3')
arrows((1:4)+0.2, y0=pred[[3]]$lci[pred[[3]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=pred[[3]]$uci[pred[[3]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="3", lwd=2)

legend("bottomright",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=1:3,pch=19,pt.cex=2,title="Fire Habitat Catergory")

###possum----

pred[[4]]

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:4,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(pred[[4]]$lci),max(pred[[4]]$uci)),xlim=c(min(0.5),max(4.5)),xlab="Midstory",ylab="Pres/Abs of Trichosurus caninus",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=0:3,xlab="Midstory")
lines(1:4,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Unburnt_Rainforest"])
arrows(x0=1:4, y0=pred[[4]]$lci[pred[[4]]$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=pred[[4]]$uci[pred[[4]]$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="black", lwd=2)

points((1:4)-0.2,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=19, col="2")
lines((1:4)-0.2,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Burnt_Rainforest"],col="2")
arrows(x0=(1:4)-0.2, y0=pred[[4]]$lci[pred[[4]]$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=pred[[4]]$uci[pred[[4]]$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="2", lwd=2)

points((1:4)+0.2,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=19,col="3")
lines((1:4)+0.2,pred[[4]]$fit[pred[[4]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='3')
arrows((1:4)+0.2, y0=pred[[4]]$lci[pred[[4]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=pred[[4]]$uci[pred[[4]]$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="3", lwd=2)

legend("topleft",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=1:3,pch=19,pt.cex=2,title="Fire Habitat Catergory")



###stigmataica----
AIC.out


#had to rerun loop with i=5 to get this
stig<-get(mod.sum$mods[5])

stig1<-data.frame(Fire_habitat_catergory=factor(c(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),4)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Understory=c(0,0,0,1,1,1,2,2,2,3,3,3))


stig2<-predict(object=stig,newdata=stig1,se.fit =T,type="response")

stig3<-data.frame(stig1,fit=stig2$fit,se=stig2$se.fit)
stig3$lci<-stig3$fit-(1.96*stig3$se)
stig3$uci<-stig3$fit+(1.96*stig3$se)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:4,stig3$fit[stig3$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p", ylim=c(min(stig3$lci),max(stig3$uci)),xlim=c(min(0.5),max(4.5)),xlab="Understory",ylab="Pres/Abs of Thylogale stigmatica",las=1,cex=3,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=0:3,xlab="Understory")
lines(1:4,stig3$fit[stig3$Fire_habitat_catergory=="Unburnt_Rainforest"])
arrows(x0=1:4, y0=stig3$lci[stig3$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=stig3$uci[stig3$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, col="black", lwd=2)

points((1:4)-0.2,stig3$fit[stig3$Fire_habitat_catergory=="Burnt_Rainforest"],cex=3,pch=19, col="2")
lines((1:4)-0.2,stig3$fit[stig3$Fire_habitat_catergory=="Burnt_Rainforest"],col="2")
arrows(x0=(1:4)-0.2, y0=stig3$lci[stig3$Fire_habitat_catergory=="Burnt_Rainforest"],x1=(1:4)-0.2, y1=stig3$uci[stig3$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, col="2", lwd=2)

points((1:4)+0.2,stig3$fit[stig3$Fire_habitat_catergory=="Burnt_Sclerophyll"],cex=3,pch=19,col="3")
lines((1:4)+0.2,stig3$fit[stig3$Fire_habitat_catergory=="Burnt_Sclerophyll"],col='3')
arrows((1:4)+0.2, y0=stig3$lci[stig3$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=(1:4)+0.2, y1=stig3$uci[stig3$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, col="3", lwd=2)

legend("bottomleft",legend=c("Unburnt Rainforest","Burnt Rainforest","Burnt Sclerophyll"),col=1:3,pch=19,pt.cex=2,title="Fire Habitat Catergory")


##Plotting simpler models----

binomsp

head(singledata,6)



###Antechinus----

ante<-glm(ifelse(singledata[,which(colnames(singledata)=='antechinus')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

antenull<-glm(ifelse(singledata[,which(colnames(singledata)=='antechinus')]==0,0,1)~1,data=singledata,family=binomial )
antelrt<-anova(antenull,ante,test="LRT")
antep<-round(antelrt$`Pr(>Chi)`[2],3)

summary(ante)$coefficients
antetab<-list("FHC"=ante,"Null"=antenull)
aictab(antetab)

head(singledata,3);dim(singledata)

ante.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

#new predictions

class(ante)
family(ante)

head(ante.1)

ante.4<-predict(object=ante,newdata=ante.1,se.fit =T,type="link")

ante.6<-data.frame(ante.1,fit.link=ante.4$fit,se.link=ante.4$se.fit)
ante.6$lci.link<-ante.6$fit.link-(1.96*ante.6$se.link)
ante.6$uci.link<-ante.6$fit.link+(1.96*ante.6$se.link)

ante.6$fit<-invlogit(ante.6$fit.link)
ante.6$se<-invlogit(ante.6$se.link)
ante.6$lci<-invlogit(ante.6$lci.link)
ante.6$uci<-invlogit(ante.6$uci.link)

install.packages("arm")
library("arm")

?invlogit


#new graph

anova(ante)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,ante.6$fit,type="p", ylim=c(min(0,min(ante.6$lci)),max(1,max(ante.6$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n", main=expression(italic("Antechinus stuartii")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=ante.6$lci,x1=1:3, y1=ante.6$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",antep,sep=""))
mtext("(a)",3,0.4,F,0)


#old predictions and graph

ante.2<-predict(object=ante,newdata=ante.1,se.fit =T,type="response")

ante.3<-data.frame(ante.1,fit=ante.2$fit,se=ante.2$se.fit)
ante.3$lci<-ante.3$fit-(1.96*ante.3$se)
ante.3$uci<-ante.3$fit+(1.96*ante.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,ante.3$fit,type="p", ylim=c(min(0,min(ante.3$lci)),max(1,max(ante.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n", main=expression(italic("Antechinus stuartii")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=ante.3$lci,x1=1:3, y1=ante.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",antep,sep=""))
mtext("(a)",3,0.4,F,0)


###Melomys----

mel<-glm(ifelse(singledata[,which(colnames(singledata)=='melomys')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

melnull<-glm(ifelse(singledata[,which(colnames(singledata)=='melomys')]==0,0,1)~1,data=singledata,family=binomial )
mellrt<-anova(melnull,mel,test="LRT")
melp<-round(mellrt$`Pr(>Chi)`[2],3)

summary(mel)$coefficients
meltab<-list("FHC"=mel,"Null"=melnull)
aictab(meltab)

head(singledata)
table(singledata$Fire_habitat_catergory,singledata$melomys>0)


mel.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))


#new predictions

family(mel)

mel.4<-predict(object=mel,newdata=mel.1,se.fit =T,type="link")

mel.5<-data.frame(mel.1,fit.link=mel.4$fit,se.link=mel.4$se.fit)
mel.5$lci.link<-mel.5$fit.link-(1.96*mel.5$se.link)
mel.5$uci.link<-mel.5$fit.link+(1.96*mel.5$se.link)

mel.5$fit<-invlogit(mel.5$fit.link)
mel.5$se<-invlogit(mel.5$se.link)
mel.5$lci<-invlogit(mel.5$lci.link)
mel.5$uci<-invlogit(mel.5$uci.link)


#new graph

anova(mel)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,mel.5$fit,type="p", ylim=c(min(0,min(mel.5$lci)),max(1,max(mel.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Melomys crevinipes")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=2:3, y0=mel.5$lci[2:3],x1=2:3, y1=mel.5$uci[2:3],angle=90,length=0.1, code=3, lwd=2)
text(3,1,labels=paste("p=",melp,sep=""))
mtext("(b)",3,0.4,F,0)
points(1.2,0.92,pch="*",cex=2)


#old predictions and graph


mel.2<-predict(object=mel,newdata=mel.1,se.fit =T,type="response")

mel.3<-data.frame(mel.1,fit=mel.2$fit,se=mel.2$se.fit)
mel.3$lci<-mel.3$fit-(1.96*mel.3$se)
mel.3$uci<-mel.3$fit+(1.96*mel.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,mel.3$fit,type="p", ylim=c(min(0,min(mel.3$lci)),max(1,max(mel.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Melomys crevinipes")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=mel.3$lci,x1=1:3, y1=mel.3$uci,angle=90,length=0.1, code=3, lwd=2)
  text(3,1,labels=paste("p=",melp,sep=""))
  mtext("(b)",3,0.4,F,0)


###Bandicoot----

band<-glm(ifelse(singledata[,which(colnames(singledata)=='bandicoot')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial )

bandnull<-glm(ifelse(singledata[,which(colnames(singledata)=='bandicoot')]==0,0,1)~1,data=singledata,family=binomial )
bandlrt<-anova(bandnull,band,test="LRT")
bandp<-round(bandlrt$`Pr(>Chi)`[2],3)


summary(band)$coefficients
bandtab<-list("FHC"=band,"Null"=bandnull)
aictab(bandtab)

band.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

#new predictions

family(band)

band.4<-predict(object=band,newdata=band.1,se.fit =T,type="link")

band.5<-data.frame(band.1,fit.link=band.4$fit,se.link=band.4$se.fit)
band.5$lci.link<-band.5$fit.link-(1.96*band.5$se.link)
band.5$uci.link<-band.5$fit.link+(1.96*band.5$se.link)

band.5$fit<-invlogit(band.5$fit.link)
band.5$se<-invlogit(band.5$se.link)
band.5$lci<-invlogit(band.5$lci.link)
band.5$uci<-invlogit(band.5$uci.link)


#new graph

anova(mel)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,band.5$fit,type="p", ylim=c(min(0,min(band.5$lci)),max(1,max(band.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Perameles nasuta")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1, y0=band.5$lci[1],x1=1, y1=band.5$uci[1],angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",bandp,sep=""))
mtext("(c)",3,0.4,F,0)
points(2.2,0.1,pch="*",cex=2)
points(3.2,0.1,pch="*",cex=2)


#old graph and predictions

band.2<-predict(object=band,newdata=band.1,se.fit =T,type="response")

band.3<-data.frame(band.1,fit=band.2$fit,se=band.2$se.fit)
band.3$lci<-band.3$fit-(1.96*band.3$se)
band.3$uci<-band.3$fit+(1.96*band.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,band.3$fit,type="p", ylim=c(min(0,min(band.3$lci)),max(1,max(band.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Perameles nasuta")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=band.3$lci,x1=1:3, y1=band.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",bandp,sep=""))
mtext("(c)",3,0.4,F,0)


###Bushrat----

bush<-glm(ifelse(singledata[,which(colnames(singledata)=='bushrat')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

bushnull<-glm(ifelse(singledata[,which(colnames(singledata)=='bushrat')]==0,0,1)~1,data=singledata,family=binomial )
bushlrt<-anova(bushnull,bush,test="LRT")
bushp<-round(bushlrt$`Pr(>Chi)`[2],3)

summary(bush)$coefficients
bushtab<-list("FHC"=bush,"Null"=bushnull)
aictab(bushtab)


bush.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

#new predictions

family(bush)

bush.4<-predict(object=bush,newdata=bush.1,se.fit =T,type="link")

bush.5<-data.frame(bush.1,fit.link=bush.4$fit,se.link=bush.4$se.fit)
bush.5$lci.link<-bush.5$fit.link-(1.96*bush.5$se.link)
bush.5$uci.link<-bush.5$fit.link+(1.96*bush.5$se.link)

bush.5$fit<-invlogit(bush.5$fit.link)
bush.5$se<-invlogit(bush.5$se.link)
bush.5$lci<-invlogit(bush.5$lci.link)
bush.5$uci<-invlogit(bush.5$uci.link)


#new graph

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,bush.5$fit,type="p", ylim=c(min(0,min(bush.5$lci)),max(1,max(bush.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Rattus fuscipes")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=2:3, y0=bush.5$lci[2:3],x1=2:3, y1=bush.5$uci[2:3],angle=90,length=0.1, code=3, lwd=2)
text(3,0.1,labels=paste("p=",bushp,sep=""))
mtext("(d)",3,0.4,F,0)
points(1.2,0.92,pch="*",cex=2)




#old graph and predicitons

bush.2<-predict(object=bush,newdata=bush.1,se.fit =T,type="response")

bush.3<-data.frame(bush.1,fit=bush.2$fit,se=bush.2$se.fit)
bush.3$lci<-bush.3$fit-(1.96*bush.3$se)
bush.3$uci<-bush.3$fit+(1.96*bush.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,bush.3$fit,type="p", ylim=c(min(0,min(bush.3$lci)),max(1,max(bush.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Rattus fuscipes")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=bush.3$lci,x1=1:3, y1=bush.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.1,labels=paste("p=",bushp,sep=""))
mtext("(d)",3,0.4,F,0)

###Stigmatica----

sti<-glm(ifelse(singledata[,which(colnames(singledata)=='stigmatica')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

stinull<-glm(ifelse(singledata[,which(colnames(singledata)=='stigmatica')]==0,0,1)~1,data=singledata,family=binomial )
stilrt<-anova(stinull,sti,test="LRT")
stip<-round(stilrt$`Pr(>Chi)`[2],3)

summary(sti)$coefficients
stitab<-list("FHC"=sti,"Null"=stinull)
aictab(stitab)

sti.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))



#new predictions

family(sti)

sti.4<-predict(object=sti,newdata=sti.1,se.fit =T,type="link")

sti.5<-data.frame(sti.1,fit.link=sti.4$fit,se.link=sti.4$se.fit)
sti.5$lci.link<-sti.5$fit.link-(1.96*sti.5$se.link)
sti.5$uci.link<-sti.5$fit.link+(1.96*sti.5$se.link)

sti.5$fit<-invlogit(sti.5$fit.link)
sti.5$se<-invlogit(sti.5$se.link)
sti.5$lci<-invlogit(sti.5$lci.link)
sti.5$uci<-invlogit(sti.5$uci.link)


#new graph

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,sti.5$fit,type="p", ylim=c(min(0,min(sti.5$lci)),max(1,max(sti.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale stigmatica")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=sti.5$lci,x1=1:3, y1=sti.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",stip,sep=""))
mtext("(e)",3,0.4,F,0)

#old graph and predictions

sti.2<-predict(object=sti,newdata=sti.1,se.fit =T,type="response")

sti.3<-data.frame(sti.1,fit=sti.2$fit,se=sti.2$se.fit)
sti.3$lci<-sti.3$fit-(1.96*sti.3$se)
sti.3$uci<-sti.3$fit+(1.96*sti.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,sti.3$fit,type="p", ylim=c(min(0,min(sti.3$lci)),max(1,max(sti.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale stigmatica")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=sti.3$lci,x1=1:3, y1=sti.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",stip,sep=""))
mtext("(e)",3,0.4,F,0)

###Thetis----

the<-glm(ifelse(singledata[,which(colnames(singledata)=='thetis')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

thenull<-glm(ifelse(singledata[,which(colnames(singledata)=='thetis')]==0,0,1)~1,data=singledata,family=binomial )
thelrt<-anova(thenull,the,test="LRT")
thep<-round(thelrt$`Pr(>Chi)`[2],3)

summary(the)$coefficients
thetab<-list("FHC"=the,"Null"=thenull)
aictab(thetab)

the.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))


#new predictions

family(the)

the.4<-predict(object=the,newdata=the.1,se.fit =T,type="link")

the.5<-data.frame(the.1,fit.link=the.4$fit,se.link=the.4$se.fit)
the.5$lci.link<-the.5$fit.link-(1.96*the.5$se.link)
the.5$uci.link<-the.5$fit.link+(1.96*the.5$se.link)

the.5$fit<-invlogit(the.5$fit.link)
the.5$se<-invlogit(the.5$se.link)
the.5$lci<-invlogit(the.5$lci.link)
the.5$uci<-invlogit(the.5$uci.link)


#new graph

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(1:3,the.5$fit,type="p", ylim=c(min(0,min(the.5$lci)),max(1,max(the.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale thetis")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=the.5$lci,x1=1:3, y1=the.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.95,labels=paste("p=",thep,sep=""))
mtext("(f)",3,0.4,F,0)



#old graph and predicitons

the.2<-predict(object=the,newdata=the.1,se.fit =T,type="response")

the.3<-data.frame(the.1,fit=the.2$fit,se=the.2$se.fit)
the.3$lci<-the.3$fit-(1.96*the.3$se)
the.3$uci<-the.3$fit+(1.96*the.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,the.3$fit,type="p", ylim=c(min(0,min(the.3$lci)),max(1,max(the.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Thylogale thetis")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=the.3$lci,x1=1:3, y1=the.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.95,labels=paste("p=",thep,sep=""))
mtext("(f)",3,0.4,F,0)

###Possum----


pos<-glm(ifelse(singledata[,which(colnames(singledata)=='possum')]==0,0,1)~Fire_habitat_catergory,data=singledata,family=binomial)

posnull<-glm(ifelse(singledata[,which(colnames(singledata)=='possum')]==0,0,1)~1,data=singledata,family=binomial )
poslrt<-anova(posnull,pos,test="LRT")
posp<-round(poslrt$`Pr(>Chi)`[2],3)

summary(pos)$coefficients
postab<-list("FHC"=pos,"Null"=posnull)
aictab(postab)

pos.1<-data.frame(Fire_habitat_catergory=factor(c('Unburnt_Rainforest',"Burnt_Rainforest","Burnt_Sclerophyll"),levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')))

#new predictions

family(pos)

pos.4<-predict(object=pos,newdata=pos.1,se.fit =T,type="link")

pos.5<-data.frame(pos.1,fit.link=pos.4$fit,se.link=pos.4$se.fit)
pos.5$lci.link<-pos.5$fit.link-(1.96*pos.5$se.link)
pos.5$uci.link<-pos.5$fit.link+(1.96*pos.5$se.link)

pos.5$fit<-invlogit(pos.5$fit.link)
pos.5$se<-invlogit(pos.5$se.link)
pos.5$lci<-invlogit(pos.5$lci.link)
pos.5$uci<-invlogit(pos.5$uci.link)


#new graph

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,pos.5$fit,type="p", ylim=c(min(0,min(pos.5$lci)),max(1,max(pos.5$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Trichosurus caninus")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=pos.5$lci,x1=1:3, y1=pos.5$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",posp,sep=""))
mtext("(g)",3,0.4,F,0)



#old graph and predicitons

pos.2<-predict(object=pos,newdata=pos.1,se.fit =T,type="response")

pos.3<-data.frame(pos.1,fit=pos.2$fit,se=pos.2$se.fit)
pos.3$lci<-pos.3$fit-(1.96*pos.3$se)
pos.3$uci<-pos.3$fit+(1.96*pos.3$se)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1)) 
plot(1:3,pos.3$fit,type="p", ylim=c(min(0,min(pos.3$lci)),max(1,max(pos.3$uci))),xlim=c(min(0.5),max(3.5)),xlab="Fire Habitat Catergory",ylab="Probability of Occurence",las=1,cex=2,pch=19,xaxt="n",main=expression(italic("Trichosurus caninus")),font.main=1,cex.lab=1.2)
axis(side=1,at=1:3,labels=c('UBR','BR','BS'))
arrows(x0=1:3, y0=pos.3$lci,x1=1:3, y1=pos.3$uci,angle=90,length=0.1, code=3, lwd=2)
text(3,0.9,labels=paste("p=",posp,sep=""))
mtext("(g)",3,0.4,F,0)



dev.new(height=12,width=12,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,3,3),mfrow=c(3,3))


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

table(mammals$Behaviour,mammals$Identification)

move.mel<-mammals[c(which(mammals$Identification=="Melomys_cervinipes")),]

dim(move.mel)

unique(move.mel$Behaviour)

move.mel<-move.mel[c(which(move.mel$Behaviour==c("Foraging_Eating")),which(move.mel$Behaviour==c("Walking_Running"))),]

dim(move.mel)

move.mel1<-data.frame(
  inactive=ifelse(move.mel$Behaviour==c("Foraging_Eating"),1,0),
  active=ifelse(move.mel$Behaviour== c("Walking_Running"),1,0),
  Site=move.mel$Site)

dim(move.mel1)
sum(move.mel1$active)
sum(move.mel1$inactive)


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


head(move.mel2);dim(move.mel2)
summary(modm.4b)

###melomys predictions----

modm.4b
min(move.mel2$Dis_to_rainforest_boundry)
max(move.mel2$Dis_to_rainforest_boundry)
unique(move.mel2$Site)

head(move.mel2,3);dim(move.mel2)


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

plot(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],type="l",lwd=2,ylab=" ",xlab="",main="Unburnt Rainforest",cex.main=1.1,ylim=c(0,1))
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(  x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,ylab=" ",xlab=" ",type="l",main="Burnt Rainforest",cex.main=1.1,xlim=c(min(move.melpred4$Dis_to_rainforest_boundry),max(move.melpred4$Dis_to_rainforest_boundry)))
lines(x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,lty=2)
lines(x=subset(move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),y=subset(move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Burnt_Rainforest"],move.melpred4$Dis_to_rainforest<213),lwd=2,lty=2)
mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$fit[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,ylab=" ",xlab=" ",type="l",main="Burnt Sclerophyll",cex.main=1.1,ylim=c(0,1))
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$lci[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
lines(x=move.melpred4$Dis_to_rainforest_boundry[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.melpred4$uci[move.melpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

head(move.mel2)

tapply(move.mel2$Dis_to_rainforest_boundry,move.mel2$Fire_habitat_catergory,FUN=summary)


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))
hist(move.mel2$Dis_to_rainforest_boundry[move.mel2$Fire_habitat_catergory=="Unburnt_Rainforest"],xlab="Dis to Rainforest Boundary",main="Unburnt Rainforest")
hist(move.mel2$Dis_to_rainforest_boundry[move.mel2$Fire_habitat_catergory=="Burnt_Rainforest"],xlab="Dis to Rainforest Boundary",main="Burnt Rainforest")
hist(move.mel2$Dis_to_rainforest_boundry[move.mel2$Fire_habitat_catergory=="Burnt_Sclerophyll"],xlab="Dis to Rainforest Boundary",main="Burnt Sclerophyll")



###melomys graphing  and predictions old----

move.melpred<-data.frame(Fire_habitat_catergory=factor(c(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),each=6192)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Dis_to_rainforest_boundry=rep(c(0:386),time=48),Site=rep(c(5,6,7,9,10,11,12,13,14,15,16,17,21,22,23,24),each=387,times=48))


head(move.melpred,6);dim(move.melpred)
tail(move.melpred,6);dim(move.melpred)

move.melpred1<-predictSE(mod=modm.4b,newdata=move.melpred,se.fit =TRUE,type="response")
head(move.melpred2)

str(move.melpred1)

move.melpred2<-data.frame(move.melpred,fit=move.melpred1$fit,se=move.melpred1$se.fit)
move.melpred2$lci<-move.melpred2$fit-(1.96*move.melpred2$se)
move.melpred2$uci<-move.melpred2$fit+(1.96*move.melpred2$se)


std.error(move.melpred1)



head(move.melpred2[101:150,])

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))

plot(x=move.melpred2[1:50,2],y=move.melpred2[1:50,3],type="l",ylim=c(0:1),col="grey20",lwd=2,ylab="Probability Of Active Behaviour",xlab="Distance To Rainforest Boundry (m)")
  lines(x=move.melpred2[1:50,2],y=move.melpred2[1:50,5],lwd=2,col="grey20",lty=2)
  lines(x=move.melpred2[1:50,2],y=move.melpred2[1:50,6],lwd=2,col="grey20",lty=2)
  
  lines(x=move.melpred2[51:100,2],y=move.melpred2[51:100,3],lwd=2,col="grey40")
  lines(x=move.melpred2[51:100,2],y=move.melpred2[51:100,5],lwd=2,col="grey40",lty=2)
  lines(x=move.melpred2[51:100,2],y=move.melpred2[51:100,6],lwd=2,col="grey40",lty=2)
  
  lines(x=move.melpred2[101:150,2],y=move.melpred2[101:150,3],lwd=2,col="grey60")
  lines(x=move.melpred2[101:150,2],y=move.melpred2[101:150,5],lwd=2,col="grey60",lty=2)
  lines(x=move.melpred2[101:150,2],y=move.melpred2[101:150,6],lwd=2,col="grey60",lty=2)
  
  legend(x=250,y=0.8,legend=c("Unburnt Rainforest","Burnt Rainforest", "Burnt Sclerophyll"),bty="n",text.col = c("grey20","grey40","grey60"))


?plot

  las=1,cex=2,pch=19,xaxt="n"
  
dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))
  
plot(x=move.melpred2[1:50,2],y=move.melpred2[1:50,3],type="l",ylim=c(0,1),lwd=2,ylab=" ",xlab="",main="Unburnt Rainforest",cex.main=1.1)
  lines(x=move.melpred2[1:50,2],y=move.melpred2[1:50,5],lwd=2,lty=2)
  lines(x=move.melpred2[1:50,2],y=move.melpred2[1:50,6],lwd=2,lty=2)
  mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
  mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)
  
plot(x=move.melpred2[51:100,2],y=move.melpred2[51:100,3],lwd=2,ylim=c(0,1),ylab=" ",xlab=" ",type="l",main="Burnt Rainforest",cex.main=1.1)
  lines(x=move.melpred2[51:100,2],y=move.melpred2[51:100,5],lwd=2,lty=2)
  lines(x=move.melpred2[51:100,2],y=move.melpred2[51:100,6],lwd=2,lty=2)
  mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
  mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)
  
plot(x=move.melpred2[101:150,2],y=move.melpred2[101:150,3],lwd=2,ylim=c(0,1.5),ylab=" ",xlab=" ",type="l",main="Burnt Sclerophyll",cex.main=1.1)
  lines(x=move.melpred2[101:150,2],y=move.melpred2[101:150,5],lwd=2,lty=2)
  lines(x=move.melpred2[101:150,2],y=move.melpred2[101:150,6],lwd=2,lty=2)
  mtext(side=1,line=2,'Distance To Rainforest Boundry (m)',cex=0.9)
  mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)
  
  
  
##rattus data----

table(mammals$Behaviour,mammals$Identification)


move.fus<-mammals[c(which(mammals$Identification=="Rattus_fuscipes")),]

dim(move.fus)

unique(move.fus$Behaviour)

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

#predicitons for an addative model

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



#interactive predictions ie don't use


move.fuspred<-data.frame(Fire_habitat_catergory=factor(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),4), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Understory=c(0,0,0,1,1,1,2,2,2,3,3,3))


head(move.fuspred)
str(move.fus2)
str(move.fuspred)
summary(modr.6a)

move.fuspred3<-predictSE(mod=modr.6a,newdata=move.fuspred,se.fit =TRUE,type="link")
head(move.fuspred3)

move.fuspred4<-data.frame(move.fuspred,fit.link=move.fuspred3$fit,se.link=move.fuspred3$se.fit)
move.fuspred4$lci.link<-move.fuspred4$fit.link-(1.96*move.fuspred4$se.link)
move.fuspred4$uci.link<-move.fuspred4$fit.link+(1.96*move.fuspred4$se.link)

move.fuspred4$fit<-invlogit(move.fuspred4$fit.link)
move.fuspred4$se<-invlogit(move.fuspred4$se.link)
move.fuspred4$lci<-invlogit(move.fuspred4$lci.link)
move.fuspred4$uci<-invlogit(move.fuspred4$uci.link)

head(move.stigpred4)

###rattus graphing----


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(1:4,move.fuspred6$fit,type="p", ylim=c(min(0,min(move.fuspred6$lci)),max(1,max(move.fuspred6$uci))),xlab="Understory",ylab="Probability of Active Behaviour",las=1,cex=2,pch=19,xaxt="n")
axis(side=1,at=1:4,labels=c(0,1,2,3))
arrows(x0=1:4, y0=move.fuspred6$lci,x1=1:4, y1=move.fuspred6$uci,angle=90,length=0.1, code=3, lwd=2)

plot(1:3,move.fuspred8$fit,type="p", ylim=c(min(0,min(move.fuspred8$lci)),max(1,max(move.fuspred8$uci))),xlab="Fire Habitat Catergory",ylab="Probability of Active Behaviour",las=1,cex=2,pch=19,xaxt="n")
axis(side=1,at=1:3,labels=c("UBR","BR","BS"))
arrows(x0=1:3, y0=move.fuspred8$lci,x1=1:3, y1=move.fuspred8$uci,angle=90,length=0.1, code=3, lwd=2)


#interactive ie wrong kind of pred graph

head(move.fuspred4[move.fuspred4$Fire_habitat_catergory=="Unburnt_Rainforest",])
head(move.fuspred4[move.fuspred4$Fire_habitat_catergory=="Burnt_Rainforest",])
head(move.fuspred4[move.fuspred4$Fire_habitat_catergory=="Burnt_Sclerophyll",])

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(x=1:4,y=move.fuspred4$fit[move.fuspred4$Fire_habitat_catergory=="Unburnt_Rainforest"],type="p",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="Unburnt Rainforest",main.cex=1,pch=19,cex=2,xaxt="n")
arrows(x0=1:4, y0=move.fuspred4$lci[move.fuspred4$Fire_habitat_catergory=="Unburnt_Rainforest"],x1=1:4, y1=move.fuspred4$uci[move.fuspred4$Fire_habitat_catergory=="Unburnt_Rainforest"],angle=90,length=0.1, code=3, lwd=2)
axis(side=1,at=1:4,labels=c('1','2','3','4'))
mtext(side=1,line=2,'Understory',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=1:4,y=move.fuspred4$fit[move.fuspred4$Fire_habitat_catergory=="Burnt_Rainforest"],type="p",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="Burnt Rainforest",main.cex=1,pch=19,cex=2,xaxt="n")
arrows(x0=1:4, y0=move.fuspred4$lci[move.fuspred4$Fire_habitat_catergory=="Burnt_Rainforest"],x1=1:4, y1=move.fuspred4$uci[move.fuspred4$Fire_habitat_catergory=="Burnt_Rainforest"],angle=90,length=0.1, code=3, lwd=2)
axis(side=1,at=1:4,labels=c('1','2','3','4'))
mtext(side=1,line=2,'Understory',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=1:4,y=move.fuspred4$fit[move.fuspred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],type="p",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="Burnt Sclerophyll",main.cex=1,pch=19,cex=2,xaxt="n")
arrows(x0=1:4, y0=move.fuspred4$lci[move.fuspred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],x1=1:4, y1=move.fuspred4$uci[move.fuspred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],angle=90,length=0.1, code=3, lwd=2)
axis(side=1,at=1:4,labels=c('1','2','3','4'))
mtext(side=1,line=2,'Understory',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)


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

#interaction DO NOT USE
move.stigpred<-data.frame(Fire_habitat_catergory=factor(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),rep(50,3)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Elevation=seq(min(move.stig2$Elevation),max(move.stig2$Elevation),length.out=50))

move.stigpred3<-predictSE(mod=mods.1a,newdata=move.stigpred,se.fit =TRUE,type="link")
head(move.stigpred3)

move.stigpred4<-data.frame(move.stigpred,fit.link=move.stigpred3$fit,se.link=move.stigpred3$se.fit)
move.stigpred4$lci.link<-move.stigpred4$fit.link-(1.96*move.stigpred4$se.link)
move.stigpred4$uci.link<-move.stigpred4$fit.link+(1.96*move.stigpred4$se.link)

move.stigpred4$fit<-invlogit(move.stigpred4$fit.link)
move.stigpred4$se<-invlogit(move.stigpred4$se.link)
move.stigpred4$lci<-invlogit(move.stigpred4$lci.link)
move.stigpred4$uci<-invlogit(move.stigpred4$uci.link)



#additive model predictions

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


dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(x=move.stigpred6$Elevation,y=move.stigpred6$fit,type="l",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="")
lines(x=move.stigpred6$Elevation,y=move.stigpred6$lci,lwd=2,lty=2)
lines(x=move.stigpred6$Elevation,y=move.stigpred6$uci,lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(1:2,move.stigpred8$fit,type="p", ylim=c(min(0,min(move.stigpred8$lci)),max(1,max(move.stigpred8$uci))),xlim=c(min(0.5),max(2.5)),xlab="Fire Habitat Catergory",ylab="Probability of Active Behaviour",las=1,cex=2,pch=19,xaxt="n")
axis(side=1,at=1:2,labels=c('UB Rainforest','B Rainforest'))
arrows(x0=1:2, y0=move.stigpred8$lci,x1=1:2, y1=move.stigpred8$uci,angle=90,length=0.1, code=3, lwd=2)














plot(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.stigpred4$fit[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],type="l",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="Unburnt Rainforest",main.cex=1)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.stigpred4$lci[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],y=move.stigpred4$uci[move.stigpred4$Fire_habitat_catergory=="Unburnt_Rainforest"],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],y=move.stigpred4$fit[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],lwd=2,ylim=c(0,1),ylab=" ",xlab=" ",type="l",main="Burnt Rainforest",main.cex=1)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],y=move.stigpred4$lci[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],lwd=2,lty=2)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],y=move.stigpred4$uci[move.stigpred4$Fire_habitat_catergory=="Burnt_Rainforest"],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.stigpred4$fit[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,ylim=c(0,1),ylab=" ",xlab=" ",type="l",main="Burnt Sclerophyll",main.cex=1)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.stigpred4$lci[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
lines(x=move.stigpred4$Elevation[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],y=move.stigpred4$uci[move.stigpred4$Fire_habitat_catergory=="Burnt_Sclerophyll"],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))
hist(move.stig2$Elevation[move.stig2$Fire_habitat_catergory=="Unburnt_Rainforest"],xlab="Dis to Rainforest Boundary",main="Unburnt Rainforest")
hist(move.stig2$Elevation[move.stig2$Fire_habitat_catergory=="Burnt_Rainforest"],xlab="Dis to Rainforest Boundary",main="Burnt Rainforest")
hist(move.stig2$Elevation[move.stig2$Fire_habitat_catergory=="Burnt_Sclerophyll"],xlab="Dis to Rainforest Boundary",main="Burnt Sclerophyll")


###stigmatica graphing  and predictions old----



mods.1a
min(move.stig2$Elevation)
max(move.stig2$Elevation)


head(move.stig2,3);dim(move.stig2)

move.stigpred<-data.frame(Fire_habitat_catergory=factor(rep(c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll'),rep(50,3)), levels=c('Unburnt_Rainforest','Burnt_Rainforest','Burnt_Sclerophyll')),Elevation=seq(min(move.stig2$Elevation),max(move.stig2$Elevation),length.out=50))

head(move.stigpred,6);dim(move.stigpred)


move.stigpred1<-predictSE(mod=mods.1a,newdata=move.stigpred,se.fit =TRUE,type="response")
head(move.stigpred1)
str(move.stigpred1)

move.stigpred2<-data.frame(move.stigpred,fit=move.stigpred1$fit,se=move.stigpred1$se.fit)
move.stigpred2$lci<-move.stigpred2$fit-(1.96*move.stigpred2$se)
move.stigpred2$uci<-move.stigpred2$fit+(1.96*move.stigpred2$se)

head(move.stigpred2)


head(move.stigpred2[101:150,])
dim(move.stigpred2)

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1))

plot(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,3],type="l",ylim=c(0:1),col="grey20",lwd=2,ylab="Probability Of Active Behaviour",xlab="Elevation (m)")
lines(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,5],lwd=2,col="grey20",lty=2)
lines(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,6],lwd=2,col="grey20",lty=2)

lines(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,3],lwd=2,col="grey40")
lines(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,5],lwd=2,col="grey40",lty=2)
lines(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,6],lwd=2,col="grey40",lty=2)

lines(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,3],lwd=2,col="grey60")
lines(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,5],lwd=2,col="grey60",lty=2)
lines(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,6],lwd=2,col="grey60",lty=2)

legend(x=700,y=1,legend=c("Unburnt Rainforest","Burnt Rainforest", "Burnt Sclerophyll"),bty="n",text.col = c("grey20","grey40","grey60"))



?plot

dev.new(height=7,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mfrow=c(2,2))

plot(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,3],type="l",ylim=c(0,1),lwd=2,ylab=" ",xlab=" ",main="Unburnt Rainforest",main.cex=1)
lines(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,5],lwd=2,lty=2)
lines(x=move.stigpred2[1:50,2],y=move.stigpred2[1:50,6],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,3],lwd=2,ylim=c(0,1),ylab=" ",xlab=" ",type="l",main="Burnt Rainforest",main.cex=1)
lines(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,5],lwd=2,lty=2)
lines(x=move.stigpred2[51:100,2],y=move.stigpred2[51:100,6],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)

plot(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,3],lwd=2,ylim=c(0,1.5),ylab=" ",xlab=" ",type="l",main="Burnt Sclerophyll",main.cex=1)
lines(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,5],lwd=2,lty=2)
lines(x=move.stigpred2[101:150,2],y=move.stigpred2[101:150,6],lwd=2,lty=2)
mtext(side=1,line=2,'Elevation (m)',cex=0.9)
mtext(side=2,line=2,"Probability Of Active Behaviour",cex=0.9)




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








#bottom----























#Time active preliminary----

table(startsWith(mammals$Time,'0')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'1')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'2')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'3')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'4')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'5')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'6')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'7')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'8')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'9')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'10')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'11')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'12')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'13')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'14')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'15')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'16')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'17')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'18')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'19')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'20')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'21')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'22')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'23')==TRUE,mammals$Identification)
table(startsWith(mammals$Time,'24')==TRUE,mammals$Identification)

time.ante<-data.frame(time=c(0:24),activity.level=c(0,15,6,1,4,7,1,0,0,0,0,0,0,0,0,1,1,1,4,7,2,2,1,1,0))
time.mel<-data.frame(time=c(0:24),activity.level=c(3,52,58,9,12,8,0,0,0,0,1,0,0,0,0,0,0,4,23,12,6,13,18,7,0))
time.ban<-data.frame(time=c(0:24),activity.level=c(2,3,3,1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0))
time.bush<-data.frame(time=c(0:24),activity.level=c(18,128,146,19,20,18,3,1,0,1,0,1,0,0,0,0,0,2,62,31,40,35,23,28,0))
time.sti<-data.frame(time=c(0:24),activity.level=c(7,103,46,13,11,14,32,14,7,7,4,6,10,7,3,10,14,27,7,5,5,14,7,11,0))
time.the<-data.frame(time=c(0:24),activity.level=c(5,22,12,1,8,9,7,4,7,5,0,1,3,2,0,1,3,7,2,2,1,3,2,1,0))
time.pos<-data.frame(time=c(0:24),activity.level=c(4,13,15,5,2,1,0,0,0,0,0,0,0,0,0,0,0,4,5,4,4,2,1,1,0))



unique(mammals$Identification)

dev.new(height=12,width=12,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,3,3),mfrow=c(3,3))

plot(time.ante,type="l",main="Antechinus")
rect(-0.5,-0.4,6,15.5,col="grey90",border="grey90")
rect(18,-0.4,24.7,15.5,col="grey90",border="grey90")
lines(time.ante)

plot(time.mel,type="l", main="Melomys")
rect(-0.7,-1.5,6,59.5,col="grey90",border="grey90")
rect(18,-1.5,24.7,59.5,col="grey90",border="grey90")
lines(time.mel)

plot(time.ban,type="l",main="Bandicoot")
rect(-0.5,-0.09,6,3.1,col="grey90",border="grey90")
rect(18,-0.09,24.7,3.1,col="grey90",border="grey90")
lines(time.ban)

plot(time.bush,type="l",main="Bushrat")
rect(-0.5,-4,6,150.5,col="grey90",border="grey90")
rect(18,-4,24.7,150.5,col="grey90",border="grey90")
lines(time.bush)

plot(time.sti,type="l",main="Stigmatica")
rect(-0.7,-3.1,6,106,col="grey90",border="grey90")
rect(18,-3,24.7,106,col="grey90",border="grey90")
lines(time.sti)

plot(time.the,type="l",main="thylogale")
rect(-0.5,-0.7,6,22.5,col="grey90",border="grey90")
rect(18,-0.5,24.7,22.5,col="grey90",border="grey90")
lines(time.the)

plot(time.pos,type="l",main="possum")
rect(-0.5,-0.4,6,15.5,col="grey90",border="grey90")
rect(18,-0.4,24.7,15.5,col="grey90",border="grey90")
lines(time.pos)



