#Ordination 


# Because some species were only collected in some years, I can't put year in TC models (phrag only from 2018, patens only 2017), CF (patens only in 2017), LU (patens only 2017)


###### Year - does year matter ######
#upshot - yes year sometimes matters, esp for phrag, not as much for spartina alterniflora

sample_data(datp)$Year<-factor(sample_data(datp)$Year)

#have multiple years: phragCF, phragLU, SaCF, SaLU
datpyear<-datp%>%
  subset_samples(Site!="Turtle Cove")%>%
  subset_samples(HostPlant!="Juncus roemerianus")%>%
  subset_samples(HostPlant!="Spartina patens")%>%
  subset_samples(HostPlant!="Sagittaria lancifolia")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracyear<-subset_dist(datpyear, unifracp)
sample_data(datpyear)$Year<-factor(sample_data(datpyear)$Year)
sample_data(datpyear)$HostPlantYear<-paste(sample_data(datpyear)$HostPlant,sample_data(datpyear)$Year)
sample_data(datpyear)$HostPlantSiteYear<-paste(sample_data(datpyear)$HostPlant,sample_data(datpyear)$Site,sample_data(datpyear)$Year)

mynmdsyearj <- ordinate(datpyear, "CAP",distance(datpyear, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Site+Year))#
anova(mynmdsyearj,by="margin",permutations = how(nperm=9999))#blocks=sample_data(datpTC)$HostPlant,

mynmdsyearu <- ordinate(datpyear, "CAP",distance=unifracyear,formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsyearu,by="margin",permutations = how(nperm=99999))

plot_ordination(datpyear, mynmdsyearj, type="samples", color="HostPlantSiteYear",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlantSiteYear),level=.95)


#subsetting
datpyear<-datp%>%
  subset_samples(Site=="CERF"&HostPlant=="Spartina alterniflora")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracyear<-subset_dist(datpyear, unifracp)

mynmdsyearj <- ordinate(datpyear, "CAP",distance(datpyear, method = "jaccard", binary = TRUE),formula=as.formula(~Year))#
anova(mynmdsyearj,by="margin",permutations = how(nperm=9999))

mynmdsyearu <- ordinate(datpyear, "CAP",distance=unifracyear,formula=as.formula(~Year))
anova(mynmdsyearu,by="margin",permutations = how(nperm=9999))

plot_ordination(datpyear, mynmdsyearj, type="samples", color="Year",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Year),level=.95)




#calculating MDP (mean pairwise distance) matrix for use in ordination


MDPdist<-comdist(data.frame(t(otu_table(datp))), cophenetic(phy_tree(datp)), abundance.weighted=F)



##### By Site #####
#Does host plant affect composition within each site?
#blocks=sample_data(datpTC)$HostPlant

#TURTLE COVE

datpTC<-datp%>%
  subset_samples(Site=="Turtle Cove")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracTC<-subset_dist(datpTC, unifracp)
MDPdistTC<-subset_dist(datpTC, MDPdist)

mynmdsTCj <- ordinate(datpTC, "CAP",distance(datpTC, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsTCj,by="margin",permutations = how(nperm=9999))

mynmdsTCu <- ordinate(datpTC, "CAP",distance=unifracTC,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsTCu,by="terms",permutations = how(nperm=9999))

mynmdsTCm <- ordinate(datpTC, "CAP",distance=MDPdistTC,formula=as.formula(~HostPlant))
anova(mynmdsTCm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpTC, mynmdsTCm, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)


#CERF

datpCF<-datp%>%
  subset_samples(Site=="CERF")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracCF<-subset_dist(datpCF, unifracp)
MDPdistTC<-subset_dist(datpCF, MDPdist)

mynmdsCFj <- ordinate(datpCF, "CAP",distance(datpCF, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFj,by="terms",permutations = how(nperm=9999))

mynmdsCFu <- ordinate(datpCF, "CAP",distance=unifracCF,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFu,by="terms",permutations = how(nperm=9999))#

mynmdsCFm <- ordinate(datpCF, "CAP",distance=MDPdistTC,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpCF, mynmdsCFj, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)


#LUMCON

datpLU<-datp%>%
  subset_samples(Site=="LUMCON")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracLU<-subset_dist(datpLU, unifracp)
MDPdistLU<-subset_dist(datpLU, MDPdist)

mynmdsLUj <- ordinate(datpLU, "CAP",distance(datpLU, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUj,by="terms",permutations = how(nperm=9999))

mynmdsLUu <- ordinate(datpLU, "CAP",distance=unifracLU,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUu,by="terms",permutations = how(nperm=9999))

mynmdsLUm <- ordinate(datpLU, "CAP",distance=MDPdistLU,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpLU, mynmdsLUj, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)





##### By Host #####
#Does site affect composition within each site?


#Phragmites 

datpPa<-datp%>%
  subset_samples(HostPlant=="Phragmites australis")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracPa<-subset_dist(datpPa, unifracp)
MDPdistPa<-subset_dist(datpPa, MDPdist)

mynmdsPaj <- ordinate(datpPa, "CAP",distance(datpPa, method = "jaccard", binary = TRUE),formula=as.formula(~Site+Condition(Year)))
anova(mynmdsPaj,by="terms",permutations = how(nperm=9999))

mynmdsPau <- ordinate(datpPa, "CAP",distance=unifracPa,formula=as.formula(~Site+Condition(Year)))#
anova(mynmdsPau,by="terms",permutations = how(nperm=9999))

mynmdsPam <- ordinate(datpPa, "CAP",distance=MDPdistPa,formula=as.formula(~Site+Condition(Year)))#+Condition(Year)
anova(mynmdsPam,by="terms",permutations = how(nperm=9999))

plot_ordination(datpPa, mynmdsPaj, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)


#Spartina alterniflora 

datpSa<-datp%>%
  subset_samples(HostPlant=="Spartina alterniflora")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracSa<-subset_dist(datpSa, unifracp)
MDPdistSa<-subset_dist(datpSa, MDPdist)

mynmdsSaj <- ordinate(datpSa, "CAP",distance(datpSa, method = "jaccard", binary = TRUE),formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSaj,by="margin",permutations = how(nperm=9999))

mynmdsSau <- ordinate(datpSa, "CAP",distance=unifracSa,formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSau,by="margin",permutations = how(nperm=9999))

mynmdsSam <- ordinate(datpSa, "CAP",distance=MDPdistSa,formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSam,by="margin",permutations = how(nperm=9999))

plot_ordination(datpSa, mynmdsSam, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)


#Spartina patens 
#only collected in 2017

datpSp<-datp%>%
  subset_samples(HostPlant=="Spartina patens")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracSp<-subset_dist(datpSp, unifracp)
MDPdistSp<-subset_dist(datpSp, MDPdist)

mynmdsSpj <- ordinate(datpSp, "CAP",distance(datpSp, method = "jaccard", binary = TRUE),formula=as.formula(~Site))
anova(mynmdsSpj,by="terms",permutations = how(nperm=9999))

mynmdsSpu <- ordinate(datpSp, "CAP",distance=unifracSp,formula=as.formula(~Site))
anova(mynmdsSpu,by="terms",permutations = how(nperm=9999))

mynmdsSpm <- ordinate(datpSp, "CAP",distance=MDPdistSp,formula=as.formula(~Site))
anova(mynmdsSpm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpSp, mynmdsSpm, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)




