#Rarefactions


#By Host Plant
rar1<-dat6%>%
  dplyr::select(HostPlant,OTU0:OTU9)%>%
  group_by(HostPlant)%>%
  summarise(across(OTU0:OTU9,sum))
rar1<-data.frame(rar1)

rownames(rar1)<-rar1$HostPlant
rar1[,1]<-NULL
out<-rarecurve(rar1,tidy=T)
colnames(out)[1]<-"HostPlant"

#note the fig looks a little funny with too large text
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/rarefactionbyhostplant.pdf",width=3.5,height=1.75)
ggplot(out,aes(x=Sample,y=Species,color=HostPlant))+
  theme_classic()+
  labs(x="Number of isolates",y="Accumulated taxa")+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),strip.placement = "outside",panel.spacing=unit(0,"cm"))+
  geom_line(linewidth=1)
dev.off()


#By Site
rar2<-dat6%>%
  dplyr::select(Site,OTU0:OTU9)%>%
  group_by(Site)%>%
  summarise(across(OTU0:OTU9,sum))
rar2<-data.frame(rar2)

rownames(rar2)<-rar2$Site
rar2[,1]<-NULL
out2<-rarecurve(rar2,tidy=T)
colnames(out2)[1]<-"Site"

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/rarefactionbysite.pdf",width=3.5,height=1.75)
ggplot(out2,aes(x=Sample,y=Species,color=Site))+
  theme_classic()+
  labs(x="Number of isolates",y="Accumulated taxa")+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),strip.placement = "outside",panel.spacing=unit(0,"cm"))+
  geom_line(linewidth=1)
dev.off()


#### most abundant taxa ####
genusspecies
dat6

abunotus<-dat6%>%
  dplyr::select(OTU0:OTU9)
abunotus2<-colSums(abunotus)
abunotus3<-data.frame(otu=colnames(abunotus),abun=abunotus2)
abunotus4<-abunotus3%>%
  plyr::join(genusspecies)%>%
  arrange(desc(abun))
sort(abunotus4$genusspecies)
