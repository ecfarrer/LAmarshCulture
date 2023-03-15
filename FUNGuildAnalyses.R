#FUNGuild analyses

head(dat7)
dat7$plantpathogenpercent<-dat7$plantpathogen/dat7$totalabundance*100
dat7$plantpathogentaxapercent<-dat7$plantpathogentaxa/dat7$SR*100
dat7$plantpathogenbroadpercent<-dat7$plantpathogenbroad/dat7$totalabundance*100
dat7$plantpathogenbroadprobablehppercent<-dat7$plantpathogenbroadprobablehp/dat7$totalabundance*100
dat7$symbiotrophpercent<-dat7$symbiotroph/dat7$totalabundance*100
dat7$plantpathogenbroadprobablehptaxapercent<-dat7$plantpathogenbroadprobablehptaxa/dat7$SR*100
dat7$symbiotrophtaxapercent<-dat7$symbiotrophtaxa/dat7$SR*100
dat7$symbiotrophprobablehptaxapercent<-dat7$symbiotrophprobablehp/dat7$SR*100

options(contrasts=c("contr.helmert","contr.poly"))


pathmean<-dat7%>%
  group_by(HostPlant,Site)%>%
  summarise(mean=mean(plantpathogentaxapercent),se=std.error(plantpathogentaxapercent),count=n())

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/richness.pdf",width=6.6,height=3.5)
ggplot(pathmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="Path") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) +
  facet_wrap(vars(Site),strip.position = "bottom")
#dev.off()

pathmean<-dat7%>%
  group_by(HostPlant)%>%
  summarise(mean=mean(plantpathogentaxapercent),se=std.error(plantpathogentaxapercent),count=n())

ggplot(pathmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="Path") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) 


m1<-lme(plantpathogentaxapercent~HostPlant+Site,random=~1|Year,data=dat7,na.action=na.omit)
anova(m1,type="marginal")
hist(resid(m1))
boxplot(resid(m1)~dat7$HostPlantSite)

m2<-lme(plantpathogentaxapercent~HostPlant*Site,random=~1|Year,data=phragspartina,na.action=na.omit)
anova(m2,type="marginal")




##### Pie charts by species #####
#this didn't seem very useful so I discontinued
#get the full list of taxa from each species, and make a pie chart
otulist<-dat7%>%
  group_by(HostPlant)%>%
  summarise(across(OTU0:OTU9,sum))

ind<-which(guilds2$guild=="Plant Pathogen")
pathos<-guilds2[ind,]$OTU

Junroe<-colnames(otulist)[which(otulist[1,]>0)]
Junroe2<-length(which(pathos%in%Junroe))
Junroe3<-length(Junroe)
pie(c(Junroe2/Junroe3,1-Junroe2/Junroe3),col=c("black","green"))

Phraus<-colnames(otulist)[which(otulist[2,]>0)]
Phraus2<-length(which(pathos%in%Phraus))
Phraus3<-length(Phraus)
pie(c(Phraus2/Phraus3,1-Phraus2/Phraus3),col=c("black","green"))

