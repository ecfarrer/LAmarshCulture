##### Loading packages #####

library(tidyverse)
library(picante)
library(plotrix)
library(phyloseq)
library(QsRutils)

#install.packages("remotes")
#remotes::install_github("jfq3/QsRutils") #distance matrix subsetting function

library(remotes)
library(ape)
library(phytools)


save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture/workspace1.Rdata")  # 

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture/workspace1.Rdata")  # 


##### Reading in data #####

##### T-BAS data, 99% similarity ####

#OTUold is the otu name that was input for the isolate name in T-BAS, every isolate has a unique OTUold. The numbers of OTUold are from Nelle's original T-BAS run and then we added A, B C, etc to them to make them unique.
#Query.sequence is also unique to each isolate. it is the OTUold plus the species name that Nelle's original T-BAS run came up with
OTUreport<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned3/tbas21_archiveCERVNAZT_/assignments_report_addvoucherCERVNAZT.csv",stringsAsFactors = T)
cbind(OTUreport$otu,OTUreport$Query.sequence)
OTUreport2<-OTUreport%>%
  dplyr::select(Query.sequence,Phylum,Taxon.assignment,Genusspecies,otu,OTU_CERVNAZT:Juncus_roemerianus_CERVNAZT,Trophic.Mode,Guild)%>%
  rename(OTUold=OTU_CERVNAZT,HostPlant=HostPlant_CERVNAZT,Site=Site_CERVNAZT,TurtleCove=TurtleCove_CERVNAZT,LUMCON=LUMCON_CERVNAZT,CERF=CERF_CERVNAZT,Spartina_patens=Spartina_patens_CERVNAZT,Spartina_alterniflora=Spartina_alterniflora_CERVNAZT,Phragmites_australis=Phragmites_australis_CERVNAZT,Sagittaria_lancifolia=Sagittaria_lancifolia_CERVNAZT,Juncus_roemerianus=Juncus_roemerianus_CERVNAZT,OTU=otu)%>%  mutate(Site = factor(Site, levels = c("Turtle Cove", "CERF", "LUMCON")))

head(OTUreport2)


##### Culture data ####

#Culture IDs are unique to the sequence
#Otus 18 and 46 did not have an ITS region, remove
culturefile<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Copy of Mastersheet - July 29, 2020, 11_55 AM.csv",stringsAsFactors = T)
culturefile2<-culturefile%>%
  filter(is.na(CultureID)==F,OTU!="OTU18",OTU!="OTU46")%>%
  dplyr::select(-Sequence)%>%
  unite(Query.sequence,OTU,TBAS_Species,remove=F)%>%
  dplyr::select(Query.sequence,Year,CultureID,PlantIndividual,TBAS_Species)%>%
  separate(PlantIndividual,c("Plant","Individual"),sep=" ")%>%
  unite(PlantIndividual,Plant,Individual,remove=T)
head(culturefile2)

length(culturefile2$CultureID)

dat<-full_join(OTUreport2,culturefile2,by="Query.sequence")
head(dat)
dim(dat)

#Collapse identical OTUs into one row. Looking at abundances of OTUs in the whole dataset. there are many many 1's (singletons). also this is needed to create the community dataset below
dat2<-dat%>%
  group_by(HostPlant,Site,Year,PlantIndividual,OTU,Genusspecies)%>%
  summarise(abundance=n())
as.data.frame(dat2)
data.frame(dat2$HostPlant,dat2$PlantIndividual,dat2$Year,dat2$OTU,dat2$abundance)

#make wide
dat3<-dat2%>%
  ungroup()%>%
  unite(PlantIndividualYear,PlantIndividual,Year,remove=F)%>%
  unite(HostPlantSite,HostPlant,Site,remove=F)%>%
  dplyr::select(-Genusspecies)%>%
  spread(OTU,abundance,fill=0)
dat.comm<-data.frame(dat3[,7:66])
row.names(dat.comm)<-dat3$PlantIndividualYear


##### Phylogenetic tree #####
#note - the uncleaned tree in this directory contains the genus species names attached to the OTU numbers
tree<-read.newick("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned3/Farrertreecleaned.nwk")

plot(tree)





##### Faith's Phylogenetic distance #####
PD<-pd(dat.comm,tree)

dat4<-data.frame(dat3[,1:6],PD,dat3[,7:66])





##### MPD by plantindividual #####
phydist <- cophenetic(tree)
ses.mpd.result.notweighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999) #takes 5 min with 999
ses.mpd.result.notweighted
ses.mpd.result.notweighted$PlantIndividualYear<-rownames(ses.mpd.result.notweighted)
ses.mpd.result.notweighted1<-ses.mpd.result.notweighted%>%
  select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.notweighted=mpd.obs.z)

ses.mpd.result.weighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=TRUE, runs=999) #takes 5 min with 999
ses.mpd.result.weighted
ses.mpd.result.weighted$PlantIndividualYear<-rownames(ses.mpd.result.weighted)
ses.mpd.result.weighted1<-ses.mpd.result.weighted%>%
  select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.weighted=mpd.obs.z)

dat5<-dat4%>%
  full_join(ses.mpd.result.notweighted1)%>%
  full_join(ses.mpd.result.weighted1)
  
dat6<-data.frame(dat5[,1:8],dat5[,69:70],dat5[,9:68])
head(dat6)








##### phyloseq object ####
otus<-dat6[,11:70]
otus2<-t(otus)
sampleotus<-dat6[,c(1:10)]
taxonomyotus<-as.matrix(data.frame(Kingdom=row.names(otus2),Phylum=row.names(otus2),Class=row.names(otus2),Order=row.names(otus2),Class=row.names(otus2),Family=row.names(otus2),Genus=row.names(otus2),Species=row.names(otus2)))
rownames(taxonomyotus)<-row.names(otus2)

datp <- merge_phyloseq(otu_table(otus2,taxa_are_rows = T), tax_table(taxonomyotus), sample_data(sampleotus),tree)

#calculate unifrac distances
unifracp<-unifrac(otus,tree)



#### Summary of files #####
datp #phyloseq object
head(dat6) #big data frame, wide data format, dat3 plus PD and MPD data
dat2 #long dataformat 



##### git hub token stuff #####
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
 #first when it asks to enter password or token I put my computer password
 #then do gitcreds_set() again and select 2, then paste my token
