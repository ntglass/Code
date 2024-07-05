#ThreshEx 1 - Nick G
#Welcome! This is my code for ThreshEx 1, see ReadMe for more information
#Citation: Havstad, K. and B. Bestelmeyer. 2019. Perennial grass recovery following livestock overgazing and shrub removal: an experiment at the Jornada Experimental Range (Jornada Basin LTER), 1996-2016 ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/cc0b40f7a34e52b612f384ff0f246e06 (Accessed 2024-07-01).

############################################# PT. 1: OBTAIN DATA
library("tidyverse")
rm(list = ls()) #clear directory
# Package ID: knb-lter-jrn.210461001.17 Cataloging System:https://pasta.edirepository.org.
# Data set title: Perennial grass recovery following livestock overgazing and shrub removal: an experiment at the Jornada Experimental Range (Jornada Basin LTER), 1996-2016.
# Data set creator:  Kris Havstad - USDA ARS Jornada Experimental Range (JER) 
# Data set creator:  Brandon Bestelmeyer - USDA ARS Jornada Experimental Range (JER) 
# Contact:  Data Manager -  Jornada Basin LTER  - datamanager.jrn.lter@gmail.com 
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@Virginia.edu 

infile1 <- trimws("https://pasta.lternet.edu/package/data/eml/knb-lter-jrn/210461001/17/12137791f6cade12e21e48dfcb8c55d9") 
infile1 <-sub("^https","http",infile1)
# This creates a tibble named: dt1 
dt1 <-read_delim(infile1  
                 ,delim=","   
                 ,skip=1 
                 , col_names=c( 
                   "year",   
                   "date",   
                   "exclosure",   
                   "block",   
                   "shrub_treatment",   
                   "grazing_treatment",   
                   "line",   
                   "point",   
                   "layer",   
                   "USDA_code",   
                   "live_dead",   
                   "scientific_name",   
                   "common_name",   
                   "habit",   
                   "form"   ), 
                 col_types=list(
                   col_character(),  
                   col_date("%Y-%m-%d"),   
                   col_character(),  
                   col_character(),  
                   col_character(),  
                   col_character(), 
                   col_number() , 
                   col_number() ,  
                   col_character(),  
                   col_character(),  
                   col_character(),  
                   col_character(),  
                   col_character(),  
                   col_character(),  
                   col_character()), 
                 na=c(" ",".","NA","")  )

# Convert Missing Values to NA for individual vectors 
dt1$date <- ifelse((trimws(as.character(dt1$date))==trimws(".")),NA,dt1$date)               
suppressWarnings(dt1$date <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$date))==as.character(as.numeric("."))),NA,dt1$date))
dt1$live_dead <- ifelse((trimws(as.character(dt1$live_dead))==trimws(".")),NA,dt1$live_dead)               
suppressWarnings(dt1$live_dead <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$live_dead))==as.character(as.numeric("."))),NA,dt1$live_dead))
dt1$scientific_name <- ifelse((trimws(as.character(dt1$scientific_name))==trimws(".")),NA,dt1$scientific_name)               
suppressWarnings(dt1$scientific_name <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$scientific_name))==as.character(as.numeric("."))),NA,dt1$scientific_name))
dt1$habit <- ifelse((trimws(as.character(dt1$habit))==trimws(".")),NA,dt1$habit)               
suppressWarnings(dt1$habit <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$habit))==as.character(as.numeric("."))),NA,dt1$habit))
dt1$form <- ifelse((trimws(as.character(dt1$form))==trimws(".")),NA,dt1$form)               
suppressWarnings(dt1$form <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(dt1$form))==as.character(as.numeric("."))),NA,dt1$form))

# Observed issues when reading the data. An empty list is good!
problems(dt1) 
# Here is the structure of the input data tibble: 
glimpse(dt1) 
# And some statistical summaries of the data 
summary(dt1) 
# Get more details on character variables
summary(as.factor(dt1$exclosure)) 
summary(as.factor(dt1$block)) 
summary(as.factor(dt1$shrub_treatment)) 
summary(as.factor(dt1$grazing_treatment)) 
summary(as.factor(dt1$layer)) 
summary(as.factor(dt1$USDA_code)) 
summary(as.factor(dt1$live_dead)) 
summary(as.factor(dt1$scientific_name)) 
summary(as.factor(dt1$common_name)) 
summary(as.factor(dt1$habit)) 
summary(as.factor(dt1$form))

############################################# PT. 1: CREATE IDEAL DATASET
#Calculating percents
dt1.counts<-dt1%>% group_by(year,exclosure) %>% select(year,exclosure,line,point) %>%unique()%>%summarise(count=length(point))
dat<-merge(dt1,dt1.counts,by=c('year','exclosure')) #working dataset
dat<-dat%>%group_by(year,exclosure,USDA_code,live_dead)%>%summarise(n=length(USDA_code),
                                                          Count=mean(count))%>%mutate(per=(n/Count)*100) #as calculated in Bestelmeyer 2013

#Removing dead biomass
summary(as.factor(dat$live_dead))
dat<-dat%>%filter(live_dead!='dead')%>%select(-live_dead)

#Add metadata
meta<-dt1%>%select(block,exclosure,shrub_treatment,grazing_treatment)%>%arrange(exclosure)%>%unique() #metadata
dat<-merge(meta,dat,by='exclosure')

#Prepare additional variables for working dataset DF:
DF<-dat%>%dplyr::select(exclosure,year,block,shrub_treatment,grazing_treatment)%>%arrange(exclosure,year)%>%unique()

#diversity measures: richness and shannon
library(vegan)
spe<-dat%>%dplyr::select(-block,-shrub_treatment,-grazing_treatment,-Count,-per)%>%spread(USDA_code,n) #creating species dataframe
spe<-spe%>%dplyr::select(-exclosure,-year)
spe[is.na(spe)]<-0 #no counts set to zero
columns_to_remove <- c('2LC', '2LTR', '2RF', 'S', 'D', '2LTRWL') #removing non-plant data
spe<- spe %>%select(-one_of(columns_to_remove)) #remove columns if they are there
DF$rich<-specnumber(spe) #species richness
DF$shan<-diversity(spe,index='shannon') #shannon diversity

#create explanatory variables: pre and post-treatment diversity
shannon.96<-DF%>%filter(year==1996)%>%dplyr::select(exclosure,rich,shan)
names(shannon.96)<-c('exclosure','rich.96','shan.96')
shannon.02<-DF%>%filter(year==2002)%>%dplyr::select(exclosure,rich,shan)
names(shannon.02)<-c('exclosure','rich.02','shan.02')
DF<-merge(DF,shannon.96,by='exclosure')
DF<-merge(DF,shannon.02,by='exclosure')

#create response variable: recovery rate = BOER (Black Grama) cover in relation to control (no grazing, shrub intact)
o.boer<-dat%>%filter(shrub_treatment=='intact',grazing_treatment=='none',USDA_code=='BOER4')%>%dplyr::select(block,year,per)
names(o.boer)[3]<-'boer.o'
boer<-filter(dat,USDA_code=='BOER4') 
boer<-merge(boer,o.boer,by=c('block','year'))
boer<-mutate(boer,rr=per/boer.o) ##rr is recovery proportion
for.DF<-boer%>% dplyr::select(exclosure,year,per,rr) 
names(for.DF)[3]<-'per.boer'
DF<-merge(DF,for.DF,by=c('exclosure','year'))

#PRGL (Honey Mesquite) foliar cover
prgl<-filter(dat,USDA_code=='PRGL2') 
names(prgl)[9]<-'per.prgl'
for.DF<-prgl%>% dplyr::select(exclosure,year,per.prgl) 
DF<-merge(DF,for.DF,by=c('exclosure','year'))

#create explanatory variables: pre and post-treatment BOER and PRGL cover
boer.96<-boer%>%filter(year==1996)%>%dplyr::select(exclosure,per)
names(boer.96)[2]<-'boer.96'
boer.02<-boer%>%filter(year==2002)%>%dplyr::select(exclosure,per)
names(boer.02)[2]<-'boer.02'
prgl.96<-dat%>%filter(year==1996,USDA_code=='PRGL2')%>%dplyr::select(exclosure,per)
names(prgl.96)[2]<-'prgl.96'
prgl.02<-dat%>%filter(year==2002,USDA_code=='PRGL2')%>%dplyr::select(exclosure,per)
names(prgl.02)[2]<-'prgl.02'
DF<-merge(DF,boer.96,by='exclosure')
DF<-merge(DF,boer.02,by='exclosure')
DF<-merge(DF,prgl.96,by='exclosure')
DF<-merge(DF,prgl.02,by='exclosure')

#create grazing severity = percent reduction in BOER cover to initial cover
DF<-mutate(DF,gs=((boer.02-boer.96)/boer.96)*-1)
DF[DF<0]<-0 #for grazing severity values, controls gained biomass

#Calculating plant net primary productivity (NPP) and rain use efficiency (RUE)
DF<-DF%>%mutate(boer.npp=per.boer*264.56,
                prgl.npp=per.prgl*184.61) #equations from Reichmann et al 2011
df.96<-filter(DF,year==1996)
df.02<-filter(DF,year==2002)
df.09<-filter(DF,year==2009)
df.16<-filter(DF,year==2016)
df.96$ppt<-178.3 #growing season (May thru Sept) ppt from IBP well
df.02$ppt<-118.6
df.09$ppt<-126.2
df.16$ppt<-178.4 #growing season ppt from IBP-G well
DF<-bind_rows(df.96,df.02,df.09,df.16)
DF<-DF%>%mutate(rue.boer=boer.npp/ppt,
                rue.prgl=prgl.npp/ppt)

#create explanatory variable: years since disturbance
DF$year<-as.numeric(as.character(DF$year))
DF<-mutate(DF,yr2=year-2002)

#REMOVING 1996 DATA
DF.96<-DF
DF<-filter(DF,year!=1996)

#creating a few more alternative datasets
DF.exclosure<-DF%>%group_by(exclosure,shrub_treatment,grazing_treatment,gs,prgl.02,rich.02,shan.02)%>%summarise(r=mean(rr),prgl=mean(per.prgl),boer=mean(per.boer),rich=mean(rich),shan=mean(shan),rue.b=mean(rue.boer),
                                                                                                       rue.p=mean(rue.prgl)) #for testing effects of grazing severity and 2002 levels of PRGL cover, richness and Shannon diversity
DF.grazed<-filter(DF,grazing_treatment!='none') #For recovery analysis
DF.exclosure.grazed<-filter(DF.exclosure,grazing_treatment!='none')
DF.initial<-filter(DF,year==2002)
DF.initial.grazed<-filter(DF.initial,grazing_treatment!='none')

############################################# PT. 2: ANALYZE DATA
library(ggfortify)
library(lme4)
library(lmerTest)

#library(ggplot2)
#library(drc)

str(DF) #the dataframe contains observations of BOER and PRGL foliar cover (per.boer and per.prgl) by exclosure for 2002-2023 in 7 year intervals. RR is the quotient of the BOER foliar cover of an exclosure divided by the BOER foliar cover of the control within the corresponding block, where the control is the shrub-intact and no-grazing paddock. Diversity metrics include rich (species richness) and shan (Shannon diversity) for each exclosure. This response variables change through time. In addition, constant variables are present: the BOER and PRGL foliar coverage and the species richness and Shannon diversity in 1996 and 2002 are given for each exclosure. GS is also constant and represents the percent change in BOER cover from 1196 to 2002. Finally, sat.MM is constant for each combination of the grazing and shrub treatment and represents the BOER saturation level (expressed in % foliar coverage) as predicted by Michaelis-Menten models. Time is expressed in 4 variables: year = sampling year, yr3 = years since treatment began (1996), yr2 = years since treatment ceased (2002), and yrlog is the log of yr2.
test<-DF%>%filter(year==2002,block==1)
print(test)

#analyze initial treatment effects
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)} #a function for testing dispersion in poisson models

#grazing severity
hist(DF.initial.grazed$gs)
shapiro.test(DF.initial.grazed$gs) #yes
boxplot(gs~grazing_treatment,data=DF.initial.grazed) #no effect
boxplot(gs~shrub_treatment,data=DF.initial.grazed) #no effect
hist(DF.initial$rue.prgl)
shapiro.test(DF.initial$rue.prgl) #nope
boxplot(rue.prgl~grazing_treatment,data=DF.initial) #no effect
boxplot(rue.prgl~shrub_treatment,data=DF.initial) #effect
#models
mod<-lm(gs~grazing_treatment+shrub_treatment,data=DF.initial.grazed)
shapiro.test(resid(mod))
autoplot(mod) #sure
summary(mod) 
anova(mod) #No effect of treatments on grazing severity
#This means that ecological resistance of Black Grama was not impacted by treatments

#shrubs
hist(DF.initial$prgl.02)
shapiro.test(DF.initial$prgl.02) #nope
boxplot(prgl.02~grazing_treatment,data=DF.initial) #no effect
#models
mod<-lm(prgl.02~grazing_treatment,data=DF.initial)
autoplot(mod) 
shapiro.test(resid(mod)) #barely
summary(mod)#.9
DF.initial.shrub<-filter(DF.initial,shrub_treatment=='intact')
DF.initial.noshrub<-filter(DF.initial,shrub_treatment=='removed')
mean((DF.initial.noshrub$per.prgl-DF.initial.shrub$per.prgl)/DF.initial.shrub$per.prgl)
#shrub treatment successfully removed 90% of PRGL biomass

#richness
hist(DF.initial$rich)
shapiro.test(DF.initial$rich) #yep
boxplot(rich~grazing_treatment,data=DF.initial) #no effect
boxplot(rich~shrub_treatment,data=DF.initial) #effect
#models
mod<-lm(rich~grazing_treatment+shrub_treatment,data=DF.initial)
autoplot(mod) 
summary(mod)
anova(mod)

#shannon
hist(DF.initial$shan)
shapiro.test(DF.initial$shan) #yep
boxplot(shan~grazing_treatment,data=DF.initial) #effect
boxplot(shan~shrub_treatment,data=DF.initial) #no effect
#models
mod<-lm(shan~grazing_treatment+shrub_treatment,data=DF.initial)
autoplot(mod) 
shapiro.test(resid(mod)) 
summary(mod)
#Shrub treatment effected number of species but not Shannon diversity.
#Grazing treatment effected Shannon diversity but not number of species.

#BOER recovery response to treatments
boxplot(rr~grazing_treatment,data=DF) #effect
boxplot(rr~shrub_treatment,data=DF) #effect
hist(DF.grazed$rr)
shapiro.test(DF.grazed$rr) #normal 
mod<-lmer(rr~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF.grazed)
shapiro.test(resid(mod)) #valid
drop1(mod)
mod2<-lmer(rr~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF.grazed)
shapiro.test(resid(mod2)) #valid
plot(mod2)
#models - quadratic
DF.grazed<-mutate(DF.grazed,yrsq=yr2^2)
quad<-lmer(rr~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF.grazed)
shapiro.test(resid(quad)) #valid
drop1(quad)
quad2<-lmer(rr~grazing_treatment*yrsq+grazing_treatment*yr2+shrub_treatment*yrsq+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF.grazed)
shapiro.test(resid(quad2)) #valid
plot(quad2)

#check model with visualization
library(gridExtra) 
DF.grazed$rr.preds.linear<-predict(mod2,DF.grazed)
DF.grazed$rr.preds.quad<-predict(quad2,DF.grazed)
DF.grazed2<-DF.grazed%>%group_by(grazing_treatment,shrub_treatment,yr2)%>%summarise(r=mean(rr),p.boer=mean(per.boer),rr.pred.lin=mean(rr.preds.linear),
                                                                                        rr.pred.quad=mean(rr.preds.quad))
a<-ggplot(DF.grazed2,aes(x=yr2,y=r))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=rr.pred.lin,color=grazing_treatment,linetype=shrub_treatment))+
  geom_hline(yintercept=1,size=2)+
  theme_classic()
b<-ggplot(DF.grazed2,aes(x=yr2,y=r))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=rr.pred.quad,color=grazing_treatment,linetype=shrub_treatment))+
  geom_hline(yintercept=1,size=2)+
  theme_classic()
grid.arrange(a,b,nrow=2) #visualizations suggest linear model is best

#get effect sizes of treatments
print(summary(mod2),correlation=TRUE)
modg<-lmer(rr~grazing_treatment*yr2+(1|block),data=DF.grazed)
mods<-lmer(rr~shrub_treatment*yr2+(1|block),data=DF.grazed)
modni<-lmer(rr~grazing_treatment*yr2+shrub_treatment*yr2+(1|block),data=DF.grazed)
anova(mod2,modg) #shrub treatment important p = .006
anova(mod2,mods) #grazing treatment important p = .001
anova(mod2,modni) #no treatment interaction p = .6
s <- as.data.frame(t(coef(mod)$block[1,]))
s2  <- as.data.frame(confint(mod)[3:10,])
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#BOER recovery response to grazing severity, post-treatment PRGL, post-treatment diversity
cor(DF.grazed$shan.02,DF.grazed$rich.02) #0.7, chance of variance inflation is high
perm <- adonis2(DF.grazed$rr ~ gs*prgl.02*shan.02*rich.02*yr2, data=DF.grazed, permutations = 999, method="euclidean") #start with PERMANOVA
perm #important variables: gs, prgl.02, shan.02, yr2
mod<-lm(r~gs+prgl.02+shan.02,data=DF.exclosure.grazed) 
autoplot(mod) #maybe violates assumptions
shapiro.test(resid(mod)) #valid
pois<-glm(r~gs+prgl.02+shan.02,family=poisson,data=DF.exclosure.grazed)
overdisp_fun(pois) #underdispersed
quas<-glm(r~gs+prgl.02+shan.02,family=quasipoisson,data=DF.exclosure.grazed)
plot(quasi)
summary(mod) #gs and prgl.02 predict, Shannon doesnt
summary(quas) #quasipoisson agrees
mod2<-lm(r~gs*prgl.02,data=DF.exclosure.grazed) #now there's room for interaction term
shapiro.test(resid(mod2)) #valid
autoplot(mod2) #also iffy
quas2<-glm(r~gs*prgl.02,family=quasipoisson,data=DF.exclosure.grazed)
summary(mod2) #no interaction
summary(quas2) #interaction, slightly more error
plot(quasi2)

#get effect sizes: lets do the quasipoisson to be safe
summary(quas2)
s <- as.data.frame(coef(mod2))
s2  <- as.data.frame(confint(mod2))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#BOER recovery response to PRGL
plot(rue.boer~rue.prgl,data=DF)
plot(rue.b~rue.p,data=DF.exclosure.grazed)
mod<-lm(rue.b~rue.p,data=DF.exclosure.grazed)
shapiro.test(resid(mod)) #valid
autoplot(mod) #good
#effect sizes
summary(mod)
s <- as.data.frame(coef(mod))
s2  <- as.data.frame(confint(mod))
s <- bind_cols(s, s2)

#BOER RUE response to treatments
boxplot(rue.boer~grazing_treatment,data=DF) #effect
boxplot(rue.boer~shrub_treatment,data=DF) #no effect
hist(DF$rue.boer)
shapiro.test(DF$rue.boer) #not normal
mod<-lmer(rue.boer~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF)
shapiro.test(resid(mod)) #valid
drop1(mod)
mod2<-lmer(rue.boer~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(mod2)) #valid
plot(mod2) #looks decent
#models - quadratic
DF<-mutate(DF,yrsq=yr2^2)
quad<-lmer(rue.boer~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF)
shapiro.test(resid(quad)) #valid
drop1(quad)
quad2<-lmer(rue.boer~grazing_treatment*yrsq+grazing_treatment*yr2+shrub_treatment*yrsq+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(quad2)) #valid
plot(quad2)

#check model with visualization
DF$rue.boer.preds.linear<-predict(mod2,DF)
DF$rue.boer.preds.quad<-predict(quad2,DF)
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2)%>%summarise(rue=mean(rue.boer),p.boer=mean(per.boer),rr.pred.lin=mean(rue.boer.preds.linear),
                                                                                        rr.pred.quad=mean(rue.boer.preds.quad))
a<-ggplot(DF2,aes(x=yr2,y=rue))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=rr.pred.lin,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
b<-ggplot(DF2,aes(x=yr2,y=rue))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=rr.pred.quad,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
grid.arrange(a,b,nrow=2) #visualizations suggest linear model is best

#get effect sizes of treatments
print(summary(mod2),correlation=TRUE)
modg<-lmer(rue.boer~grazing_treatment*yr2+(1|block),data=DF)
mods<-lmer(rue.boer~shrub_treatment*yr2+(1|block),data=DF)
modni<-lmer(rue.boer~grazing_treatment*yr2+shrub_treatment*yr2+(1|block),data=DF)
anova(mod2,modg) #P>.001
anova(mod2,mods)  #P>.001
anova(mod2,modni) #P=.26
s <- as.data.frame(t(coef(mod2)$block[1,]))
s2  <- as.data.frame(confint(mod2)[3:12,])
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#BOER RUE response to grazing severity, post-treatment PRGL, post-treatment diversity
perm <- adonis2(DF$rue.boer ~ gs*prgl.02*shan.02*rich.02*yr2, data=DF, permutations = 999, method="euclidean") #start with PERMANOVA
perm #important variables: gs, prgl.02, shan.02, rich.02
mod<-lm(rue.b~gs+prgl.02+shan.02+rich.02,data=DF.exclosure) 
autoplot(mod) #Normal Q-Q plot looks decent
shapiro.test(resid(mod)) #valid
drop1(mod) #prgl not significant
mod<-lm(rue.b~gs+shan.02+rich.02,data=DF.exclosure) 
autoplot(mod) 
shapiro.test(resid(mod)) #valid
summary(mod) #all significant

#get effect sizes 
summary(mod)
s <- as.data.frame(coef(mod))
s2  <- as.data.frame(confint(mod))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#PRGL productivity response to treatments
boxplot(per.prgl~grazing_treatment,data=DF) #no effect
boxplot(per.prgl~shrub_treatment,data=DF) #effect
hist(DF$per.prgl)
shapiro.test(DF$per.prgl) #non-normal proportion data
mod<-lmer(per.prgl~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF) #singular
shapiro.test(resid(mod)) #NOT valid (non-linear correct way to go)
#models - logarithmic 
emod<-lmer(log(per.prgl)~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF) #singular
emod2<-lmer(log(per.prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF) #still singular
shapiro.test(resid(emod2))
plot(emod2)

#models- quadratic
quad<-lmer(log(per.prgl)~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF) #singular
quad2<-lmer(log(per.prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment*yrsq+shrub_treatment*yrsq+grazing_treatment:shrub_treatment+(1|block),data=DF) #still singular
shapiro.test(resid(quad2)) #valid

#because there's a lot of singularity happening, I will simplify the dataset
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2,yrsq)%>%summarise(r=mean(rr),p.boer=mean(per.boer),prgl=mean(per.prgl))
mod<-lm(prgl~grazing_treatment*shrub_treatment*yr2,data=DF2)
autoplot(mod) #iffy
shapiro.test(resid(mod))
drop1(mod,test='Chisq')
mod2<-lm(prgl~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
shapiro.test(resid(mod2)) #valid
AIC(mod,mod2) #mod2 has lower AIC
#model - trying exponential on simplified data
emod<-lm(log(prgl)~grazing_treatment*shrub_treatment*yr2,data=DF2)
autoplot(emod) #not great
shapiro.test(resid(emod)) #not valid
emod2<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
shapiro.test(resid(emod2)) #valid
autoplot(emod2) #much better!
#model - trying quadratic on simplified data
quad2<-lm(prgl~grazing_treatment*yr2+grazing_treatment*yrsq+shrub_treatment*yr2+shrub_treatment*yrsq+grazing_treatment:shrub_treatment,data=DF2)
autoplot(quad2) #Q-Q plot isn't as good as emod2
shapiro.test(resid(quad2)) #valid

#check model with visualization
DF2$preds.e.prgl<-exp(predict(emod2,DF2))
DF2$preds.linear.prgl<-predict(mod2,DF2)
e<-ggplot(DF2,aes(x=yr2,y=prgl))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=preds.linear.prgl,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
f<-ggplot(DF2,aes(x=yr2,y=prgl))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=preds.e.prgl,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
grid.arrange(e,f,nrow=2) #Q-Q plots, and visualize suggest logarithmic model is best

#get effect sizes of treatments 
summary(emod2)
emodg<-lm(log(prgl)~grazing_treatment*yr2,data=DF2)
emods<-lm(log(prgl)~shrub_treatment*yr2,data=DF2)
emodni<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2,data=DF2)
anova(emod2,emodg) #P>.001 for shrub trt
anova(emod2,emods) #P=.2 for grazing trt
anova(emod2,emodni) #P=.08 for grazing*shrub interaction
#because interaction term is marginally significant, I will include it in some more anovas
emodg2<-lm(log(prgl)~grazing_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
emods2<-lm(log(prgl)~shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
anova(emod2,emodg2) #still P>.001 for shrub trt
anova(emod2,emods2) #P=.4 for grazing trt when interaction still included
s <- as.data.frame(coef(emod2))
s2  <- as.data.frame(confint(emod2))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#looking further into shrub treatment effects over time
emod<-lm(log(prgl)~shrub_treatment:yr2,data=DF2)
shapiro.test(resid(emod)) #valid
anova(emod) #P > .01, rate of growth differs between shrub treatments

#PRGL productivity response to grazing severity, post-treatment PRGL, post-treatment diversity
perm <- adonis2(DF$per.prgl ~ gs*prgl.02*shan.02*rich.02*yr2, data=DF, permutations = 999, method="euclidean") #start with PERMANOVA
perm #important variables: prgl.02, gs
mod<-lm(prgl~gs*prgl.02,data=DF.exclosure) 
autoplot(mod) #Normal Q-Q plot could be better
shapiro.test(resid(mod)) #barely valid
summary(mod) #only prgl.02 predicts
mod2<-lm(prgl~prgl.02,data=DF.exclosure)
shapiro.test(resid(mod2)) #hmmm
quasi<-glm(prgl~prgl.02,family=quasipoisson,data=DF.exclosure)
summary(quasi)
plot(quasi)

#get effect sizes  
s <- as.data.frame(coef(quasi))
s2  <- as.data.frame(confint(quasi))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#PRGL RUE response to treatments
boxplot(rue.prgl~grazing_treatment,data=DF) #no effect
boxplot(rue.prgl~shrub_treatment,data=DF) #effect
hist(DF$rue.prgl)
shapiro.test(DF$rue.prgl) #non-normal
mod<-lmer(rue.prgl~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF)
shapiro.test(resid(mod)) #valid
plot(mod) #nope, inconstant variance
#models - logarithmic 
emod<-lmer(log(rue.prgl)~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF) #singular
drop1(emod)
emod2<-lmer(log(rue.prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF) #still singular
shapiro.test(resid(emod2)) #valid
plot(emod2)
#models- quadratic
quad<-lmer(log(rue.prgl)~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF) #singular
drop1(quad)
quad2<-lmer(log(rue.prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment*yrsq+shrub_treatment*yrsq+grazing_treatment:shrub_treatment+(1|block),data=DF) #still singular
shapiro.test(resid(quad2)) #valid
plot(quad2)
AIC(quad,quad2) #quad2 has lower AIC

#because there's a lot of singularity happening, I will simplify the dataset
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2,yrsq)%>%summarise(r=mean(rr),p.boer=mean(per.boer),prgl=mean(rue.prgl))
mod<-lm(prgl~grazing_treatment*shrub_treatment*yr2,data=DF2)
autoplot(mod) #Q-Q not great
drop1(mod,test='Chisq')
mod2<-lm(prgl~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
shapiro.test(resid(mod2)) #valid
autoplot(mod2) #better
#model - trying exponential on simplified data
emod<-lm(log(prgl)~grazing_treatment*shrub_treatment*yr2,data=DF2)
drop1(emod,test='Chisq')
emod2<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
shapiro.test(resid(emod2)) #valid
autoplot(emod2)

#check model with visualization
DF2$preds.e.prgl<-exp(predict(emod2,DF2))
DF2$preds.linear.prgl<-predict(mod2,DF2)
e<-ggplot(DF2,aes(x=yr2,y=prgl))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=preds.linear.prgl,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
f<-ggplot(DF2,aes(x=yr2,y=prgl))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=preds.e.prgl,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
grid.arrange(e,f,nrow=2) #similar to productivity, we'll go with logarithmic again

#get effect sizes of treatments 
summary(emod2)
emodg<-lm(log(prgl)~grazing_treatment*yr2,data=DF2)
emods<-lm(log(prgl)~shrub_treatment*yr2,data=DF2)
emodni<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2,data=DF2)
anova(emod2,emodg) #P>.001 for shrub trt
anova(emod2,emods) #P=.13 for grazing trt
anova(emod2,emodni) #P=.05 for grazing*shrub interaction
#because interaction term is marginally significant, I will include it in some more anovas
emodg2<-lm(log(prgl)~grazing_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
emods2<-lm(log(prgl)~shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
anova(emod2,emodg2) #still P>.001 for shrub trt
anova(emod2,emods2) #P=.4 for grazing trt when interaction still included
DF2.grazed<-filter(DF2,grazing_treatment!='none') #For getting difference between winter and summer grazing
emod2<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2.grazed)
emods<-lm(log(prgl)~shrub_treatment*yr2,data=DF2.grazed)
anova(emod2,emods) #grazing unimportant for winter vs summer at P = .5
s <- as.data.frame(coef(emod2))
s2  <- as.data.frame(confint(emod2))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#PRGL RUE response to grazing severity, post-treatment PRGL, post-treatment diversity
perm <- adonis2(DF$rue.prgl ~ gs*prgl.02*shan.02*rich.02*yr2, data=DF, permutations = 999, method="euclidean") #start with PERMANOVA
perm #important variables: gs, prgl.02
mod<-lm(rue.p~gs*prgl.02,data=DF.exclosure) 
autoplot(mod) #Normal Q-Q plot could be better
shapiro.test(resid(mod)) #not valid
quasi<-glm(rue.p~gs*prgl.02,family=quasipoisson,data=DF.exclosure) 
summary(quasi) #gs does not predict
quasi<-glm(rue.p~prgl.02,family=quasipoisson,data=DF.exclosure) 
plot(quasi)

#get effect sizes  
summary(quasi)
s <- as.data.frame(coef(quasi))
s2  <- as.data.frame(confint(quasi))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

##Analyze treatment effects on diversity
hist(DF$rich)
shapiro.test(DF$rich) #normal
hist(DF$shan)
shapiro.test(DF$shan) #normal
boxplot(rich~grazing_treatment,data=DF) #no trend
boxplot(rich~shrub_treatment,data=DF) #no trend
boxplot(shan~grazing_treatment,data=DF) #greater in grazed plots
boxplot(shan~shrub_treatment,data=DF) #no trend

#shannon diversity
mod<-lmer(shan~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF)
shapiro.test(resid(mod)) #valid
drop1(mod)
mod2<-lmer(shan~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(mod2)) #valid
plot(mod2)
#model - exponential
emod<-lmer(log(shan)~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF)
shapiro.test(resid(emod)) #valid
drop1(emod)
emod2<-lmer(log(shan)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(emod2)) #valid
plot(emod2)
#model - quadratic
quad<-lmer(shan~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF)
shapiro.test(resid(quad)) #valid
drop1(quad)
quad2<-lmer(shan~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment*yrsq+shrub_treatment*yrsq+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(quad2)) #valid
plot(quad2)

#check model with visualization
DF$preds.shan<-predict(mod2,DF)
DF$preds.shan.e<-exp(predict(emod2,DF))
DF$preds.shan.quad<-predict(quad2,DF)
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2,yrsq)%>%summarise(pred.shan=mean(preds.shan),pred.shan.e=mean(preds.shan.e),pred.shan.quad=mean(preds.shan.quad),s=mean(shan))
o<-ggplot(DF2,aes(x=yr2,y=s))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.shan,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
p<-ggplot(DF2,aes(x=yr2,y=s))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.shan.e,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
q<-ggplot(DF2,aes(x=yr2,y=s))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.shan.quad,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
grid.arrange(o,p,q,nrow=3) #similar, use linear

#check effect sizes of treatments
summary(mod2)
mods<-lmer(shan~shrub_treatment*yr2+(1|block),data=DF)
modg<-lmer(shan~grazing_treatment*yr2+(1|block),data=DF)
modni<-lmer(shan~grazing_treatment*yr2+shrub_treatment*yr2+(1|block),data=DF)
anova(mod2,modg) #P=.3 for shrub trt
anova(mod2,mods) #P>.001 for grazing trt
anova(mod2,modni) #P=.7 for grazing*shrub interaction
mod2<-lmer(shan~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF.grazed)
mods<-lmer(shan~shrub_treatment*yr2+(1|block),data=DF.grazed)
anova(mod2,mods) #grazing season unimportant at P = .1
s <- as.data.frame(t(coef(mod2)$block[1,]))
s2  <- as.data.frame(confint(mod2)[3:12,])
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#Shannon response to grazing severity, post-treatment PRGL, post-treatment diversity
perm<-adonis2(DF$shan~gs*yr2*rich.02*shan.02*prgl.02,data=DF,permutations = 999, method="euclidean")
perm #important variables: gs, shan.02
mod<-lm(shan~gs*shan.02,data=DF.exclosure)
autoplot(mod) #Normal Q-Q plot decent
shapiro.test(resid(mod)) #valid
summary(mod) #only shan important
mod<-lm(shan~shan.02,data=DF.exclosure)
shapiro.test(resid(mod))

#check effect sizes 
summary(mod)
s <- as.data.frame(coef(mod))
s2  <- as.data.frame(confint(mod))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#species richness
mod<-lmer(rich~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF) #singular
shapiro.test(resid(mod)) #valid
drop1(mod)
mod2<-lmer(rich~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF) #singular
shapiro.test(resid(mod)) #valid
plot(mod2)
#model - exponential
emod<-lmer(log(rich)~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF) #singular
shapiro.test(resid(emod)) #valid
drop1(emod)
emod2<-lmer(log(rich)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF) #singular
shapiro.test(resid(emod2)) #valid
plot(emod2)
#model - quadratic
quad<-lmer(rich~grazing_treatment*shrub_treatment*yr2+grazing_treatment*shrub_treatment*yrsq+(1|block),data=DF)
shapiro.test(resid(quad)) #valid
drop1(quad)
quad2<-lmer(rich~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment*yrsq+shrub_treatment*yrsq+grazing_treatment:shrub_treatment+(1|block),data=DF)
shapiro.test(resid(quad2)) #valid
plot(quad2)

#check model with visualization
DF$preds.rich.lin<-predict(mod2,DF)
DF$preds.rich.e<-exp(predict(emod2,DF))
DF$preds.rich.quad2<-predict(quad2,DF)
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2,yrsq)%>%summarise(pred.rich=mean(preds.rich.lin),pred.rich.e=mean(preds.rich.e),pred.rich.q=mean(preds.rich.quad2),r=mean(rich))
w<-ggplot(DF2,aes(x=yr2,y=r))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.rich,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
x<-ggplot(DF2,aes(x=yr2,y=r))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.rich.e,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
y<-ggplot(DF2,aes(x=yr2,y=r))+
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=2)+
  geom_line(aes(y=pred.rich.q,color=grazing_treatment,linetype=shrub_treatment))+
  theme_classic()
grid.arrange(w,x,y,nrow=3) #similar, go with linear

#check effect sizes of treatments
print(summary(mod2),correlation=TRUE)
mods<-lmer(rich~shrub_treatment*yr2+(1|block),data=DF)
modg<-lmer(rich~grazing_treatment*yr2+(1|block),data=DF)
modni<-lmer(rich~grazing_treatment*yr2+shrub_treatment*yr2+(1|block),data=DF)
anova(mod2,modg) #P>.001 for shrub trt
anova(mod2,mods) #P=.16 for grazing trt
anova(mod2,modni) #P=.07 for grazing*shrub interaction
#because interaction term is marginally significant, I will include it in some more anovas
modg2<-lmer(rich~grazing_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
mods2<-lmer(rich~shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
anova(mod2,modg2) #still P>.001 for shrub trt
anova(mod2,mods2) #P=.2 for grazing trt when interaction still included
s <- as.data.frame(t(coef(mod2)$block[1,]))
s2  <- as.data.frame(confint(mod2)[3:12,])
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

#richness response to grazing severity, post-treatment PRGL, post-treatment diversity 
perm<-adonis2(DF$rich~gs*yr2*rich.02*shan.02*prgl.02,data=DF,permutations = 999, method="euclidean")
perm #only richness is important
mod<-lm(rich~rich.02,data=DF.exclosure)
shapiro.test(resid(mod))

#check effect sizes 
summary(mod)
s <- as.data.frame(coef(mod))
s2  <- as.data.frame(confint(mod))
s <- bind_cols(s, s2) #effect sizes and 95% confidence intervals

############################################# PT. 3: FIGURES
library(Manu)

se <- function(x, na.rm = TRUE) {sd(x, na.rm = TRUE)/sqrt(length(x))} #a function for computing standard errors
print_pal(get_pal('Hihi'))
get_pal('Hihi') #color palattes for figures

#creating modeling dataset
formods<-data.frame(yr2=rep(seq(0,300,.1),18),
                    yrlog=log(rep(seq(0,300,.1),18)),
                    grazing_treatment=rep(c(rep('none',3001),rep('summer',3001),rep('winter',3001)),6),
                    shrub_treatment=rep(c(rep('removed',3001*3),rep('intact',3001*3)),6),
                    block=c(rep(1,36012),rep(2,36012),rep(3,36012)))
formods<-mutate(formods,yrsq=yr2^2)
formods.grazed<-filter(formods,grazing_treatment!='none')

#Figure 1: additional variables by treatment (grazing severity, diversity)
forgg<-DF.exclosure%>%group_by(grazing_treatment,shrub_treatment)%>%summarise(GS=mean(gs),GS.se=se(gs),
                                                                              prgl=mean(prgl.02),prgl.se=se(prgl.02),
                                                                              rich=mean(rich.02),rich.se=se(rich.02),
                                                                              shan=mean(shan.02),shan.se=se(shan.02))
forgg.grazed<-filter(forgg,grazing_treatment!='none')
dodge <- position_dodge(width=0.9)

gg1 <- ggplot(forgg.grazed, aes(x = grazing_treatment, y = GS, fill = shrub_treatment)) +   
  geom_bar(aes(fill = shrub_treatment), color='black', position = "dodge", stat="identity") +
  labs(x = "Grazing treatment", y = "Grazing Severity") +
  geom_errorbar(aes(ymin = GS-GS.se, ymax=GS+GS.se), width=0.25,
                position = dodge) +
  theme_classic() +
  scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by=0.1)) +
  scale_fill_manual(name = "Shrub Treatment",
                    labels = c("intact", "removed"),
                    values = c("dark grey", "white")) +
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.8, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 1, label = "A", size = 14)
gg1

gg2 <- ggplot(forgg, aes(x = grazing_treatment, y = prgl, fill = shrub_treatment)) +   
  geom_bar(aes(fill = shrub_treatment), color='black', position = "dodge", stat="identity") +
  labs(x = "Grazing treatment", y = "Honey Mesquite Cover '02 (%)") +
  geom_errorbar(aes(ymin = prgl-prgl.se, ymax=prgl+prgl.se), width=0.25,
                position = dodge) +
  theme_classic() +
  scale_y_continuous(limits=c(0,5.17), breaks = seq(0, 5, by=0.5)) +
  scale_fill_manual(name = "Shrub Treatment",
                    labels = c("intact", "removed"),
                    values = c("dark grey", "white")) +
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = 'none',
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 5.17, label = "B", size = 14)
gg2

gg3 <- ggplot(forgg, aes(x = grazing_treatment, y = rich, fill = shrub_treatment)) +   
  geom_bar(aes(fill = shrub_treatment), color='black', position = "dodge", stat="identity") +
  labs(x = "Grazing treatment", y = "Species richness '02") +
  geom_errorbar(aes(ymin = rich-rich.se, ymax=rich+rich.se), width=0.25,
                position = dodge) +
  theme_classic() +
  scale_y_continuous(limits=c(0,24), breaks = seq(0, 24, by=3)) +
  scale_fill_manual(name = "Shrub Treatment",
                    labels = c("intact", "removed"),
                    values = c("dark grey", "white")) +
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = 'none',
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 24, label = "C", size = 14)
gg3

gg4 <- ggplot(forgg, aes(x = grazing_treatment, y = shan, fill = shrub_treatment)) +   
  geom_bar(aes(fill = shrub_treatment), color='black', position = "dodge", stat="identity") +
  labs(x = "Grazing treatment", y = "Shannon diversity '02") +
  geom_errorbar(aes(ymin = shan-shan.se, ymax=shan+shan.se), width=0.25,
                position = dodge) +
  theme_classic() +
  scale_y_continuous(limits=c(0,2), breaks = seq(0, 2, by=0.2)) +
  scale_fill_manual(name = "Shrub Treatment",
                    labels = c("intact", "removed"),
                    values = c("dark grey", "white")) +
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = 'none',
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 2, label = "D", size = 14)
gg4

grid.arrange(gg1,gg2,gg3,gg4,nrow=2) #export 1200 x 1005

#Resilience figures (BOER recovery, PRGL, diversity)
boer.mod<-lmer(rue.boer~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
DF2<-DF%>%group_by(grazing_treatment,shrub_treatment,yr2,yrsq,year)%>%summarise(r=mean(rr),rse=se(rr),p.boer=mean(per.boer),prgl=mean(rue.prgl),prglse=se(rue.prgl),
                                                                           s=mean(shan),sse=se(shan),Rich=mean(rich),Richse=se(rich),GS=mean(gs))
prgl.mod<-lm(log(prgl)~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment,data=DF2)
recovery.mod<-lmer(rr~grazing_treatment*shrub_treatment*yr2+(1|block),data=DF.grazed)
shan.mod<-lmer(shan~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF)
rich.mod<-lmer(rich~grazing_treatment*yr2+shrub_treatment*yr2+grazing_treatment:shrub_treatment+(1|block),data=DF) #singular

formods$b.preds<-predict(boer.mod,formods,type='response')
formods$p.preds<-exp(predict(prgl.mod,formods,type='response'))
formods.grazed$r.preds<-predict(recovery.mod,formods.grazed,type='response')
formods$s.preds<-predict(shan.mod,formods,type='response')
formods$r.preds<-predict(rich.mod,formods,type='response')
formods2.grazed<-formods.grazed%>%group_by(grazing_treatment,shrub_treatment,yr2)%>%summarise(rp=mean(r.preds))
formods2<-formods%>%group_by(grazing_treatment,shrub_treatment,yr2)%>%summarise(bp=mean(b.preds),pp=mean(p.preds),
                                                                                sp=mean(s.preds),rp=mean(r.preds))

DF2.grazed<-filter(DF2,grazing_treatment!='none')
gb<- ggplot(DF2.grazed,aes(x=yr2,y=r)) +
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=5) +
  geom_errorbar(aes(ymin = r-rse, ymax=r+rse,color=grazing_treatment), size=1, width=.25)+
  geom_line(data=formods2.grazed,aes(y=rp,color=grazing_treatment,linetype=shrub_treatment),size=2) +
  geom_hline(yintercept=1,color='red',size=3)+
  theme_classic() +
  scale_color_manual(name='Grazing Treatment',
                     values=c("#F9E211", "#A8ACAD"))+
  scale_shape_manual(name='Shrub Treatment',
                     values=c('circle', 'triangle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  scale_x_continuous(limits=c(-0.25,30), breaks=seq(0,30,5))+
  scale_y_continuous(limits=c(0,1.4), breaks=seq(0,1.4,.2))+
  labs(x="Years Since Disturbance", y="Black Grama Recovery")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 12))
gb #export in 1000 x 709

gp<- ggplot(DF2,aes(x=yr2)) +
  geom_point(aes(y=prgl,color=grazing_treatment,shape=shrub_treatment),size=5) +
  geom_errorbar(aes(ymin = prgl-prglse, ymax=prgl+prglse,color=grazing_treatment), size=1, width=.25)+
  geom_line(data=formods2,aes(y=pp,color=grazing_treatment,linetype=shrub_treatment),size=2) +
  theme_classic() +
  scale_color_manual(name='Grazing Treatment',
                     values=c("dark green","#F9E211", "#A8ACAD"))+
  scale_shape_manual(name='Shrub Treatment',
                     values=c('circle', 'triangle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  scale_x_continuous(limits=c(-0.25,21.5), breaks=seq(0,30,5))+
  scale_y_continuous(limits=c(0,20), breaks=seq(0,20,4))+
  labs(x="Years Since Disturbance", y="Honey Mesquite RUE (g/mm)")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))
gp #export in 1000 x 709

gs<- ggplot(DF2,aes(x=yr2,y=s)) +
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=5) +
  geom_errorbar(aes(ymin =s-sse, ymax=s+sse,color=grazing_treatment), size=1, width=.25)+
  geom_line(data=formods2,aes(y=sp,color=grazing_treatment,linetype=shrub_treatment),size=2) +
  theme_classic() +
  scale_color_manual(name='Grazing Treatment',
                     values=c("dark green","#F9E211", "#A8ACAD"))+
  scale_shape_manual(name='Shrub Treatment',
                     values=c('circle', 'triangle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  scale_x_continuous(limits=c(-0.25,30), breaks=seq(0,30,5))+
  scale_y_continuous(limits=c(0,3), breaks=seq(0,3,.5))+
  labs(x="Years Since Disturbance", y="Shannon Diversity Index")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = 'none',
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 3, label = "B", size = 14)
gs

gr<- ggplot(DF2,aes(x=yr2,y=Rich)) +
  geom_point(aes(color=grazing_treatment,shape=shrub_treatment),size=5) +
  geom_errorbar(aes(ymin =Rich-Richse, ymax=Rich+Richse,color=grazing_treatment), size=1, width=.25)+
  geom_line(data=formods2,aes(y=rp,color=grazing_treatment,linetype=shrub_treatment),size=2) +
  theme_classic() +
  scale_color_manual(name='Grazing Treatment',
                     values=c("dark green","#F9E211", "#A8ACAD"))+
  scale_shape_manual(name='Shrub Treatment',
                     values=c('circle', 'triangle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  scale_x_continuous(limits=c(-0.25,30), breaks=seq(0,30,5))+
  scale_y_continuous(limits=c(0,40), breaks=seq(0,40,5))+
  labs(x="Years Since Disturbance", y="Species Richness")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.7, 0.2),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 0.6, y = 40, label = "A", size = 14)
gr

grid.arrange(gr,gs,nrow=1) #export in 1000 x 709

#More figures (inspired by Bestelmeyer et al 2013)
DF.decline2<-DF%>%mutate(bchange=per.boer-boer.96,
                         pchange=per.prgl-prgl.96)%>%group_by(grazing_treatment,shrub_treatment,year)%>%
  summarise(bc=mean(bchange),bcse=se(bchange),
            pc=mean(pchange),pcse=se(pchange))
DF2<-DF2%>%unite('trt',grazing_treatment:shrub_treatment,sep=' ',remove=FALSE)
DF3<-merge(DF.decline2,DF2,by=c('grazing_treatment','shrub_treatment','year'))%>%dplyr::select(grazing_treatment,shrub_treatment,year,bc,bcse,
                                                                                               pc,pcse,trt)
forDF3<-data.frame(trt=c('none intact', 'none removed', 'summer intact', 'summer removed', 'winter intact', 'winter removed'),
                   year=rep(1996,6),
                   grazing_treatment=c('none', 'none', 'summer', 'summer', 'winter', 'winter'),
                   shrub_treatment=rep(c('intact','removed'),3),
                   bc=rep(0,6),bcse=rep(0,6),pc=rep(0,6),pcse=rep(0,6))
DF3<-bind_rows(DF3,forDF3)

gb2<- ggplot(DF3,aes(x=year,y=bc)) +
  geom_point(aes(color=trt,shape=grazing_treatment),size=5) +
  geom_errorbar(aes(ymin = bc-bcse, ymax=bc+bcse,color=trt), size=1, width=.25)+
  geom_line(aes(color=trt,linetype=shrub_treatment),size=2) +
  geom_hline(yintercept = 0,color='dark green',linetype='dotted')+
  theme_classic() +
  scale_color_manual(name='Treatment',
                     values=c("light green","dark green", "orange",'brown','light blue','blue'))+
  scale_shape_manual(name='Grazing Treatment',
                     values=c('square', 'triangle','circle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  scale_y_continuous(limits=c(-10,60), breaks=seq(-10,60,10))+
  labs(x="Year", y="Mean % foliar cover of Black Grama (relative to 1996)")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))
gb2 #export in 1000 x 837

gp2<- ggplot(DF3,aes(x=year,y=pc)) +
  geom_point(aes(color=trt,shape=grazing_treatment),size=5) +
  geom_errorbar(aes(ymin = pc-pcse, ymax=pc+pcse,color=trt), size=1, width=.25)+
  geom_line(aes(color=trt,linetype=shrub_treatment),size=2) +
  geom_hline(yintercept = 0,color='dark green',linetype='dotted')+
  theme_classic() +
  scale_color_manual(name='Treatment',
                     values=c("light green","dark green", "orange",'brown','light blue','blue'))+
  scale_shape_manual(name='Grazing Treatment',
                     values=c('square', 'triangle','circle'))+
  scale_linetype_manual(name='Shrub Treatment',
                        values=c('solid','dashed'))+
  #scale_x_continuous(limits=c(-0.25,21.5), breaks=seq(0,30,5))+
  scale_y_continuous(limits=c(-10,60), breaks=seq(-10,60,10))+
  labs(x="Year", y="Mean % foliar cover of Mesquite (relative to 1996)")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))
gp2 #export in 1000 x 837

#Initial cover plots
DF.decline16<-filter(DF,year==2016)%>%mutate(bchange=per.boer-boer.02)%>%filter(grazing_treatment!='none')
plot(bchange~boer.02,data=DF.decline16)
mod2<-lm(bchange~boer.02,data=DF.decline16)
autoplot(mod2)
shapiro.test(resid(mod2)) #everything looks good
DF.decline16$preds<-predict(mod2,DF.decline16,type='response')
summary(mod2) #relationship at P<.001, R2=.8

g16<- ggplot(DF.decline16,aes(x=boer.02,y=bchange)) +
  geom_point(size=5) +
  geom_line(aes(y=preds),size=2) +
  theme_classic() +
  scale_x_continuous(limits=c(0,16), breaks=seq(0,16,2))+
  scale_y_continuous(limits=c(0,45), breaks=seq(0,45,5))+
  labs(x="% foliar cover of Black Grama after petrurbation", y="Change in % cover of Black Grama during recovery")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 12, y = 14, label = "2016", size = 12)+
  annotate("text", x = 12, y = 8, label = "R-squared = .8", size = 8)
g16 #export 1000 x 837

DF.decline09<-filter(DF,year==2009)%>%mutate(bchange=per.boer-boer.02)%>%filter(grazing_treatment!='none')
plot(bchange~boer.02,data=DF.decline09)
mod3<-lm(bchange~boer.02,data=DF.decline09)
autoplot(mod3)
shapiro.test(resid(mod3)) #everything looks good
DF.decline09$preds<-predict(mod3,DF.decline09,type='response')
summary(mod3) #relationship at P=.003, R2=.56

g09<- ggplot(DF.decline09,aes(x=boer.02,y=bchange)) +
  geom_point(size=5) +
  geom_line(aes(y=preds),size=2) +
  theme_classic() +
  scale_x_continuous(limits=c(0,16), breaks=seq(0,16,2))+
  scale_y_continuous(limits=c(0,45), breaks=seq(0,45,5))+
  labs(x="% foliar cover of Black Grama after petrurbation", y="Change in % cover of Black Grama during recovery")+ 
  theme(axis.text.x=element_text(angle=0, vjust=0.65)) +
  theme(
    legend.position = c(.2, 0.85),
    legend.box.just = "right") +
  theme(axis.text=element_text(colour = "black", size = 20)) +
  theme(axis.title=element_text(colour = "black", size=24)) +
  theme(
    legend.title = element_text(colour = "black", size = 15),
    legend.text = element_text(colour = "black", size = 15))+
  annotate("text", x = 12, y = 40, label = "2009", size = 12)+
  annotate("text", x = 12, y = 34, label = "R-squared = .56", size = 8)
g09 #export 1000 x 837

#End

