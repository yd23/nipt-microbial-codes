####################### ####################### ####################### #######################
# Load Count Matrix
library("data.table")
tdb = fread("CountMatrix.gz")
nLeadCol <- 2
annot = read.table("SpeciesInfo.txt", header = T, sep='\t', stringsAsFactors = F)
save(tdb, nLeadCol, annot, file=file.path('./','niptmicrobial.Rda'),version = 2)


####################### ####################### ####################### #######################
# calculating RPKM and output
# util functions
writeGzFile <- function(df, gzfile){
  gzcon <- gzfile(gzfile, "w")
  write.table(df, gzcon, row.names=F, col.names=T, sep='\t', quote=F)
  close(gzcon)
}
maxCore <- 16
nLeadCol <- 2
rpkmTdb<-do.call(cbind, lapply((nLeadCol+1):ncol(tdb), function(i) as.numeric(tdb[,i])/tdb$libSize*1E6/annot$genomeSize[i-nLeadCol]*1000))
tdbFreq<-apply(tdb[,-(1:nLeadCol)],2, function(x) sum(x>0)/length(x))
rpkmTdbMax<-apply(rpkmTdb, 2, function(x) max(x, na.rm=T))
rpkmTdbMax[is.infinite(rpkmTdbMax)]<-0
rpkmTdbRef<-apply(rpkmTdb, 2, function(x) {
  tmm<-AnnoroadPD:::trim(x, 0.01)
  v1=quantile(tmm,probs=0.95, na.rm=T)
  ifelse(v1 <= 0, quantile(tmm, 0.99, na.rm=T), v1)
})
rpkmTdbRef[is.infinite(rpkmTdbRef)]<-0

save(annot, tdb, nLeadCol,
     rpkmTdb, tdbFreq, rpkmTdbMax, rpkmTdbRef,
     file=file.path('./','niptmicrobial.Rda'))

writeGzFile(data.frame(lib=tdb[,1], rpkmTdb, stringsAsFactors = F), file.path('niptmicrobial.all.RPKMmat.tsv.gz'))
writeGzFile(data.frame(annot,rpkmTdbRef, rpkmTdbMax, tdbFreq), file.path('niptmicrobial.rpkmRef.popFreq.tsv.gz'))


####################### ####################### ####################### #######################
# Subjects and general information (Fig1C/D)
library(dplyr)
require(reshape2)
library(ggplot2)
library(RColorBrewer)
library(egg)
require(grid)

load('niptmicrobial.Rda')
tdb = as.data.frame.array(tdb)
#all <- colSums(tdb[, (nLeadCol+1):ncol(tdb)])
#No0Sp <- names(sort(all,decreasing = T)[sort(all,decreasing = T)>0])
#tdb = as.data.frame.array(tdb)[No0Sp]


TopFreqSp <- names(sort(tdbFreq,decreasing = T))
topn <- 15
skipn <- 6 # skip highest 6 BACT manually
TopBACT<-TopFreqSp[grepl('BACT_.*',TopFreqSp)][(skipn+1):(skipn+topn)]
TopEUKY<-TopFreqSp[grepl('EUKY_.*',TopFreqSp)][1:topn]
TopVIR<-TopFreqSp[grepl('virus',TopFreqSp)][1:topn]
TopOther<- TopFreqSp[!(grepl('virus',TopFreqSp)|grepl('EUKY_.*',TopFreqSp)|grepl('BACT_.*',TopFreqSp))][1:15]

BactSp <- sapply(strsplit(annot[match(TopBACT,annot$head),3],split = ' '),function(x) paste(x[1:2],collapse = ' '))
EukySp <- sapply(strsplit(annot[match(TopEUKY,annot$head),3],split = ' '),function(x) paste(x[1:2],collapse = ' '))
OtherSp <- sapply(strsplit(annot[match(TopOther,annot$head),3],split = ' '),function(x) paste(x[1:2],collapse = ' '))

rpkmTDB <- data.frame(rpkmTdb)
#rpkmTDB <- data.frame(rpkmTDB)
colnames(rpkmTDB) <-names(tdbFreq)
TopBactRPKM <- rpkmTDB[TopBACT]
names(TopBactRPKM) = BactSp
TopEukyRPKM <- rpkmTDB[TopEUKY]
names(TopEukyRPKM) <- EukySp
TopVirRPKM <- rpkmTDB[TopVIR]
TopOtherRPKM <- rpkmTDB[TopOther]
names(TopOtherRPKM) <- OtherSp

TopBactAb <- tdb[TopBACT]
names(TopBactAb) = BactSp
TopEukyAb<- tdb[TopEUKY]
names(TopEukyAb) <- EukySp
TopVirAb <- tdb[TopVIR]
TopOtherAb <- tdb[TopOther]
names(TopOtherAb) <- OtherSp

save(tdbFreq,
     TopBACT, TopEUKY, TopVIR, TopOther,
     TopBactRPKM,TopEukyRPKM,TopVirRPKM,TopOtherRPKM,
     TopBactAb,TopVirAb,TopEukyAb,TopOtherAb,
     BactSp, EukySp, OtherSp,
     file = 'TopN.Rda')

load('TopN.Rda')

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Mycol = rev(getPalette(15))

# RPKM data
BactPdata <-melt(TopBactRPKM,measure.vars  = rev(colnames(TopBactRPKM)))
EukyPdata <-melt(TopEukyRPKM,measure.vars  = rev(colnames(TopEukyRPKM)))
VirPdata <-melt(TopVirRPKM,measure.vars  = rev(colnames(TopVirRPKM)))
OtherPdata <-melt(TopOtherRPKM,measure.vars  = rev(colnames(TopOtherRPKM)))

Prpkm <- function(Pdata){
  p <- ggplot(data=Pdata,aes(x=variable,y=log10(value)))+
    labs(x = NULL,y='log10(RPKM)')+
    scale_fill_manual(values = Mycol)+
    coord_flip()+
    geom_boxplot(aes(fill=variable))+
    guides(fill=FALSE)+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=15,color="black"),
          axis.title.x = element_text(size=15,color="black"))
}
BactPrpkm <- Prpkm(BactPdata)
EukyPrpkm<-Prpkm(EukyPdata)
VirPrpkm<-Prpkm(VirPdata)
OtherPrpkm<-Prpkm(OtherPdata)

# Abundance data
BactPdata <-melt(TopBactAb,measure.vars  = rev(colnames(TopBactAb)))
BactPdata <- BactPdata[BactPdata$value>0,]
EukyPdata <-melt(TopEukyAb,measure.vars  = rev(colnames(TopEukyAb)))
EukyPdata <- EukyPdata[EukyPdata$value>0,]
VirPdata <-melt(TopVirAb,measure.vars  = rev(colnames(TopVirAb)))
VirPdata <- VirPdata[VirPdata$value>0,]
OtherPdata <-melt(TopOtherAb,measure.vars  = rev(colnames(TopOtherAb)))
OtherPdata <- OtherPdata[OtherPdata$value>0,]

Pab <- function(Pdata){
  p<-ggplot(data=Pdata,aes(x=variable,y=value))+
    labs(x = NULL,y='Abundance')+
    coord_flip()+
    #ylim(10,10000)+
    #scale_y_log10(limits = c(10,5000))+
    guides(fill=FALSE)+
    #geom_violin(aes(fill=variable))+
    scale_color_manual(values = Mycol)+
    geom_jitter(shape=20,stat = 'identity',aes(color=variable),show.legend = FALSE)+
    #scale_fill_manual(values = Mycol)+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=15,color="black"),
          axis.title.x = element_text(size=15,color="black"),
          legend.position="none")
}
Pab_B <- function(Pdata){
  p<-ggplot(data=Pdata,aes(x=variable,y=value))+
    labs(x = NULL,y='Abundance')+
    coord_flip()+
    #ylim(10,10000)+
    scale_y_log10(limits = c(10,5000))+
    guides(fill=FALSE)+
    #geom_violin(aes(fill=variable))+
    scale_color_manual(values = Mycol)+
    geom_jitter(shape=20,stat = 'identity',aes(color=variable),show.legend = FALSE)+
    #scale_fill_manual(values = Mycol)+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=15,color="black"),
          axis.title.x = element_text(size=15,color="black"),
          legend.position="none")
}
Pab_V <- function(Pdata){
  p<-ggplot(data=Pdata,aes(x=variable,y=value))+
    labs(x = NULL,y='Abundance')+
    coord_flip()+
    #ylim(10,10000)+
    scale_y_log10(limits = c(1,500))+
    guides(fill=FALSE)+
    #geom_violin(aes(fill=variable))+
    scale_color_manual(values = Mycol)+
    geom_jitter(shape=20,stat = 'identity',aes(color=variable),show.legend = FALSE)+
    #scale_fill_manual(values = Mycol)+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=15,color="black"),
          axis.title.x = element_text(size=15,color="black"),
          legend.position="none")
}
BactPab<-Pab_B(BactPdata)
#EukyPab<-Pab(EukyPdata)
VirPab<-Pab_V(VirPdata)
#OtherPab<-Pab(OtherPdata)

# Freq data
BactFreq <- data.frame(Species=BactSp, Freq=tdbFreq[TopBACT],stringsAsFactors = F)
EukyFreq <- data.frame(Species=EukySp,Freq=tdbFreq[TopEUKY],stringsAsFactors = F)
VirFreq <- data.frame(Species=TopVIR,Freq=tdbFreq[TopVIR],stringsAsFactors = F)
OtherFreq <- data.frame(Species=OtherSp,Freq=tdbFreq[TopOther],stringsAsFactors = F)

Pfreq <- function(Pdata){
  P = ggplot(data=Pdata, mapping=aes(x=reorder(Species,Freq), y=Freq*100, fill=reorder(Species,Freq,decreasing = F)))+
    labs(x = NULL,y = 'prevalence (%)')+
    guides(fill=FALSE)+
    geom_bar(stat= 'identity')+
    scale_fill_manual(values = Mycol)+
    coord_flip()+
    theme(axis.text.y = element_text(size=15,color="black"),
          axis.text.x = element_text(size=15,color="black"),
          axis.title.x = element_text(size=15,color="black")
    )
}
BactPfreq <- Pfreq(BactFreq)
EukyPfreq = Pfreq(EukyFreq)
VirPfreq = Pfreq(VirFreq)
OtherPfreq = Pfreq(OtherFreq)

# Draw
Mycol_bak=rev(c('red3', 'deeppink', 'orchid4', 'mediumslateblue', 
                "royalblue1", 'lightseagreen','green4',
                'paleturquoise3','yellow2', 'springgreen2', 'yellowgreen', 
                'darkgoldenrod1', 'darksalmon', 'tomato','orangered'))
Mycol_using = c("#999999", "#CE8BAE", "#EB7AA9", "#BD6253", "#BF862B", "#F2E631", "#FFC81D", "#FF7F00", "#C4625D"
                ,"#8D5B96", "#629363", "#46A169", "#3A85A8", "#815375", "#E41A1C")

BactPfreq <- set_panel_size(BactPfreq, width = unit(10, "cm"), height = unit(10, "cm"))
BactPrpkm <- set_panel_size(BactPrpkm, width = unit(10, "cm"), height = unit(10, "cm"))
BactPab <- set_panel_size(BactPab, width = unit(10, "cm"), height = unit(10, "cm"))

VirPfreq <- set_panel_size(VirPfreq, width = unit(10, "cm"), height = unit(10, "cm"))
VirPrpkm <- set_panel_size(VirPrpkm, width = unit(10, "cm"), height = unit(10, "cm"))
VirPab <- set_panel_size(VirPab, width = unit(10, "cm"), height = unit(10, "cm"))

fig <- grid.arrange(BactPfreq, BactPrpkm,BactPab,VirPfreq, VirPrpkm,VirPab, ncol=3)
ggsave(file='F1CD.generalInformation.pdf',plot=fig,width = 25,height = 20)


####################### ####################### ####################### #######################
# Simulation of potential saturation effect (Fig 1A)
load('niptmicrobial.Rda')
maxCore <- 1
Ns<- seq(1, 2001, length.out = 201)
#Ns<- seq(1, 101, length.out = 11)
n <- 100
set.seed(1000)
libi<-seq_len(nrow(tdb))
library(parallel)
simfig2<-mclapply(Ns, function(N) {
  message(N)
  mergs<-do.call(rbind, lapply(1:n, function(i) {
    ii<-sample(libi, size = N, replace = T)
    apply(tdb[ii, -1],2, function(x) sum(as.integer(x), na.rm=T)) 
  }))
  apply(mergs, 1, function(x) sum(x>0, na.rm=T) )
}, mc.cores=maxCore)

pdf(file.path('F1A.detectionRateSimulation.pdf'), width = 8, height=5)
boxplot(simfig2, xlab='Number of individuals sampled', ylab='Number of distinct species detected', names=Ns, cex=0.5)
dev.off()


####################### ####################### ####################### #######################
# Geographic-associated
library('openxlsx')
library('vegan')
library('data.table')
library("ggplot2")
library("plyr")
library("ggthemes")
library("rgdal")
library("dplyr")
library(RColorBrewer)

load('niptmicrobial.Rda')
all <- colSums(tdb[,(nLeadCol+1):ncol(tdb)])
# 按reads多少进行sort微生物,只保留有检出的微生物物种。1933个微生物被检出过。
Spid <- names(sort(all,decreasing = T)[sort(all,decreasing = T)>0]) 
tdb <- as.data.frame(tdb)[c('Province',Spid)]
tdb[,2:ncol(tdb)][tdb[,2:ncol(tdb)]>0] <-1 #用样本检出数作为物种的检出数，不是reads数。############ Notice!

MyRed = '#E41A1C'
MyBlue = '#3A85A8'
MyYellow = "#F2E631"

china_map<-readOGR("lib/CHN_adm/bou2_4p.shp",stringsAsFactors=FALSE)
mydata<-china_map@data["NAME"]
mydata$id<-0:924
mydata[mydata$id==898,"NAME"]<-"澳门特别行政区"
mapdata<-fortify(china_map)
mapdata$id<-as.numeric(mapdata$id)
mapdata<-merge(mapdata,mydata,all.x=TRUE)
names(mapdata)[names(mapdata) == 'NAME'] <- 'region'

Province <- unique(mydata$NAME)
t <- data.frame()
for (x in Province){
  t <- rbind(t,apply(tdb[,2:ncol(tdb)][which(tdb[,1]==x),],2,sum))
}
regionSamplesCount<-as.data.frame(table(tdb$Province))
names(t) <- names(tdb[,2:ncol(tdb)])
SampleNum<-regionSamplesCount[match(Province,regionSamplesCount$Var1),2]
countTab <- rbind(SampleNum,t(t))
colnames(countTab) <- Province

skipn = 6
ShannonIndex <- apply(t[(skipn+1):ncol(t)], 1, diversity) # 去掉了Top6的物种
SIData <- data.frame(Province,SampleNum,ShannonIndex)
#SIData[SIData==0] <- NA
SIMapData <- data.frame(mapdata,ShannonIndex=SIData[match(mapdata$region,SIData[,1]),][,3],row.names=NULL)
SIMapData$ShannonIndex <- as.numeric(as.character(SIMapData$ShannonIndex))


captiallist <- 'capital.xlsx'
capName <- read.xlsx(captiallist)
capName[,1] <- paste0(capName[,1],'省')
capName[capName[,1]=='新疆省',1] <- '新疆维吾尔自治区'
capName[capName[,1]=='西藏省',1] <- '西藏自治区'
capName[capName[,1]=='内蒙古省',1] <- '内蒙古自治区'
capName[capName[,1]=='广西省',1] <- '广西壮族自治区'
capName[capName[,1]=='宁夏省',1] <- '宁夏回族自治区'
capName[capName[,1]=='香港省',1] <- '香港特别行政区'
capName[capName[,1]=='澳门省',1] <- '澳门特别行政区'
capName[capName[,1]=='天津省',1] <- '天津市'
capName[capName[,1]=='北京省',1] <- '北京市'
capName[capName[,1]=='上海省',1] <- '上海市'
capName[capName[,1]=='重庆省',1] <- '重庆市'

capName$group <- SIMapData[match(capName$Province,SIMapData$region),7]
num <- SIData[match(capName$Province,SIData$Province),2]
si <- round(SIData[match(capName$Province,SIData$Province),3], digits = 1)
SampleMapData <- data.frame(capName,num=num,knum=round(num/1000,digits = 1),si=si,stringsAsFactors = F)

MapData <- cbind(SIMapData,SampleMapData[match(SIMapData$region,SampleMapData$Province),5:7])
SampleMapData[is.na(SampleMapData)] <- 0

MyCol = c('grey88', MyYellow, MyRed)
getPalette = colorRampPalette(MyCol)
pdf("F2.SampleNum.pdf")
ggplot(MapData, aes(x = long, y = lat, group=group, fill=knum)) +
  geom_polygon(colour="LightSlateGray") +
  scale_fill_gradient(low=MyYellow, high=MyRed, na.value = 'gray88') +
  coord_map("polyconic")+
  guides(fill=guide_legend(title="Sample number",title.position='bottom',reverse = T))+
  geom_text(aes(label=knum),size =3,data = SampleMapData,group=Province)+
  theme_void()+
  theme(legend.position = c(0.9,0.5))
dev.off()
pdf("F2.ShannonIndex.pdf")
ggplot(MapData, aes(x = long, y = lat, group=group, fill=factor(si))) +
  geom_polygon(colour="LightSlateGray") +
  #scale_fill_gradient(low = 'red', high = 'blue', na.value = 'gray88',limits=c(4.8,5.6)) +
  scale_fill_manual(values = getPalette(8))+
  coord_map("polyconic")+
  guides(fill=guide_legend(title="Diversity index",title.position='bottom',reverse = T))+
  geom_text(aes(label=si),size =3,data = SampleMapData, group=Province)+
  theme_void()+theme(legend.position = c(0.9,0.5))
dev.off()


####################### ####################### ####################### #######################
# Heat Map Fig3/Fig4
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(pheatmap)

load('niptmicrobial.Rda')
#annot$description[which(annot$head=='BACT_5005')] <- 'Actinomyces sp.'

tranName <- function(head){
  sp <- annot$description[match(head,annot$head)] %>%
    strsplit(.,split = ' ')%>%
    sapply(., function(x) {
      if(length(x)>=3){
        paste(x[1:3],collapse = ' ')
      }else{
        paste(x[1:2],collapse = ' ')
      }
    })
  return(sp)
}
province <- unique(na.omit(as.character(tdb$Province)))%>%
  .[.!='全国省']%>%
  .[.!='黑吉辽省']%>%
  .[.!='青岛省']
tranProvince <- function(province){
  province2 <- substr(province,1,2)
  province2[which(province2=='黑龙')]<-'黑龙江'
  cp =c('江苏','广东','上海','福建','北京',
        '山东','河北','陕西','湖北','重庆',
        '四川','甘肃','浙江','安徽','河南',
        '辽宁','新疆','江西','山西','天津',
        '广西','吉林','云南','湖南','海南',
        '贵州','黑龙江','宁夏'
  )
  ep = c('JS','GD','SH','FJ','BJ',
         'SD','HEB','SXI','HUB','CQ',
         'SC','GS','ZJ','AH','HEN',
         'LN','XJ','JX','SX','TJ',
         'GX','JL','YN','HUN','HN',
         'GZ','HLJ','NX')
  Province <- data.frame(cp,ep)
  province2 <- Province$ep[match(province2,Province$cp)]
  return(province2)
}

# Heat Map
MyRed = '#E41A1C'
MyBlue = '#3A85A8'
MyYellow = "#F2E631"

MyHeatColor <- c(MyBlue, "lightsteelblue2","lightsteelblue1", "white", "yellow", 'orange', MyRed)

# Top BACT Heat Map
topN = 50
skipN = 6
n = skipN+topN-1
topNFreq <- lapply(province, function(x){
  message(x)
  apply(tdb[which(tdb$Province == x),match(names(sort(tdbFreq,decreasing = T))[(skipN+1):n],colnames(tdb))], 2, function(x){
    sum(x>0)/length(x)
  })
})%>%
  do.call(rbind,.)
rownames(topNFreq) <- tranProvince(province)
colnames(topNFreq) <- tranName(colnames(topNFreq))
heatData <- t(topNFreq)
hp1 <- pheatmap(heatData,
                color=colorRampPalette(MyHeatColor)(20))
pdf('F3.top50BactHeatMap.pdf')
hp1
dev.off()

# Pearson Matrix
bacts <- names(sort(tdbFreq,decreasing = T)) %>%
  .[grepl('BACT.*',.)]%>%
  .[(skipN+1):n+1]

bactFreq <- lapply(province, function(x){
  message(x)
  apply(tdb[which(tdb$Province == x),match(bacts,colnames(tdb))], 2, function(x){
    sum(x>0)/length(x)
  })
})%>%
  do.call(rbind,.)
rownames(bactFreq) <- tranProvince(province)
BactsProvincePrevelance <- t(bactFreq)
PearsonBact <- cor(BactsProvincePrevelance)
# write.table(PearsonBact,file = "BactsProvincePrevelance_Pearson")
hp2 <- pheatmap(
  PearsonBact,
  color=colorRampPalette(MyHeatColor)(20),
  display_numbers = FALSE
)
pdf('F3.Top50BactEachProvincePearson.pdf', width = 8,height = 8)
hp2
dev.off()

# Virus Heat Map
virus <- names(sort(tdbFreq,decreasing = T)) %>%
  .[grepl('virus.*',.)|grepl('Virus.*',.)]%>%
  .[1:20]

virFreq <- lapply(province, function(x){
  message(x)
  apply(tdb[which(tdb$Province == x),match(virus,colnames(tdb))], 2, function(x){
    sum(x>0)/length(x)
  })
})%>%
  do.call(rbind,.)
rownames(virFreq) <- tranProvince(province)
heatData <- t(virFreq)
hp2<-pheatmap(heatData,
              color=colorRampPalette(MyHeatColor)(20))

pdf('top20virus_heatMap.pdf',
    width = 10,height = 5
)
hp2
dev.off()

VirusProvincePrevelance <- t(virFreq)
PearsonVirus <- cor(VirusProvincePrevelance)
#write.table(PearsonVirus,file = "Virus(TopN)ProvincePrevelance_Pearson")
hp1 <- pheatmap(PearsonVirus,
                color=colorRampPalette(MyHeatColor)(20))
pdf('F4.VirusEachProvincePearson.pdf',width = 8,height = 8)
hp1
dev.off()


####################### ####################### ####################### #######################
library('data.table')
library('openxlsx')
library("ggplot2")
library("plyr")
library("ggthemes")
library("rgdal")
library("dplyr")
load("niptmicrobial.Rda")


HBV <- tdb$`Hepatitis B virus`
HSV <- tdb$`Human herpesvirus 1` +tdb$`Human herpesvirus 2`
Toxo <- tdb$EUKY_26
CMV <- tdb$`Human herpesvirus 5`
HPV <- c('Human papillomavirus type 4',
         'Human papillomavirus type 63',
         'Human papillomavirus RTRX7',
         'Human papillomavirus type 115',
         'Human papillomavirus - 18',
         'Human papillomavirus type 5',
         'Human papillomavirus type 9',
         'Human papillomavirus type 104',
         'Human papillomavirus type 24',
         'Human papillomavirus type 49',
         'Human papillomavirus type 98',
         'Human papillomavirus SIBX-3a',
         'Human papillomavirus type 48',
         'Human papillomavirus type 105',
         'Human papillomavirus type 71',
         'Human papillomavirus - 1',
         'Human papillomavirus type 121',
         'Human papillomavirus type 10',
         'Human papillomavirus type 96',
         'Human papillomavirus type 116',
         'Human papillomavirus type 88',
         'Human papillomavirus type 60',
         'Human papillomavirus type 112',
         'Human papillomavirus - 2',
         'Human papillomavirus 54',
         'Human papillomavirus FA75/KI88-03',
         'Human papillomavirus type 100',
         'Human papillomavirus type 101',
         'Human papillomavirus type 103',
         'Human papillomavirus type 108',
         'Human papillomavirus type 109',
         'Human papillomavirus type 113',
         'Human papillomavirus type 114',
         'Human papillomavirus type 16',
         'Human papillomavirus type 26',
         'Human papillomavirus type 32',
         'Human papillomavirus type 34',
         'Human papillomavirus type 41',
         'Human papillomavirus type 50',
         'Human papillomavirus type 53',
         'Human papillomavirus type 61',
         'Human papillomavirus type 6b',
         'Human papillomavirus type 7',
         'Human papillomavirus type 90',
         'Human papillomavirus type 92',
         'Human papillomavirus type 99'
)
HPV = apply(as.data.frame(tdb)[HPV], 1, sum)
EBV <- tdb$`Human herpesvirus 4 type 1` # 15 samples
VIR <- cbind(Province=tdb$Province,HBV=HBV,HSV=HSV,Toxo=Toxo,CMV=CMV,HPV=HPV,EBV=EBV)
VIR[,2:ncol(VIR)][VIR[,2:ncol(VIR)]>0] <- 1

Prov_spnum <- as.data.frame(table(as.data.frame(VIR)$Province))
names(Prov_spnum) <- c('Province','Count')

VIR <- as.data.table(VIR)
Virus = VIR[,.(HBV=sum(as.integer(HBV)),HPV=sum(as.integer(HPV)),EBV=sum(as.integer(EBV))),by=Province]
Virus = cbind(Prov_spnum,Virus[match(Prov_spnum$Province,Virus$Province),][,-1])
ToRCH = VIR[,.(Toxo=sum(as.integer(Toxo)),CMV=sum(as.integer(CMV)),HSV=sum(as.integer(HSV))),by=Province]
ToRCH = cbind(Prov_spnum,ToRCH[match(Prov_spnum$Province,ToRCH$Province),][,-1])
VirusFreq <- do.call(cbind,lapply(3:ncol(Virus), function(i) Virus[,i]/Virus[2]))
names(VirusFreq) <- names(Virus[,3:ncol(Virus)])
VirusFreq <- cbind(Virus[,1:2],VirusFreq)
TorchFreq <- do.call(cbind,lapply(3:ncol(ToRCH), function(i) ToRCH[,i]/ToRCH[2]))
names(TorchFreq) <- names(ToRCH[,3:ncol(Virus)])
TorchFreq <- cbind(ToRCH[,1:2],TorchFreq)

china_map<-readOGR("lib/CHN_adm/bou2_4p.shp",stringsAsFactors=FALSE)
mydata<-china_map@data["NAME"]
mydata$id<-0:924
mydata[mydata$id==898,"NAME"]<-"澳门特别行政区"
mapdata<-fortify(china_map)
mapdata$id<-as.numeric(mapdata$id)
mapdata<-merge(mapdata,mydata,all.x=TRUE)
mapdata<- mapdata %>%rename(region=NAME)

VirusMapData <- cbind(mapdata,VirusFreq[match(mapdata$region,VirusFreq$Province),][3:5],row.names=NULL)
VirusMapData[is.na(VirusMapData)] <- 0
library("tidyr")
VirusFactor <- tidyr::gather(VirusMapData[,9:11],Species,Freq)
MapData2 <- rbind(VirusMapData[,1:8],VirusMapData[,1:8],VirusMapData[,1:8])
VirusMapData2 <- cbind(MapData2,VirusFactor)

TorchMapData <- cbind(mapdata,TorchFreq[match(mapdata$region,TorchFreq$Province),][3:5],row.names=NULL)
TorchMapData[is.na(TorchMapData)] <- 0
TorchFactor <- tidyr::gather(TorchMapData[,9:11],Species,Freq)
#TorchMapData2 <- cbind(MapData2,TorchFactor)
############## Split Virus
PV <- function(Pdata, f){
  p<-ggplot(Pdata, aes(x = long, y = lat,group=group,fill=f)) +
    geom_polygon(colour="LightSlateGray") +
    #facet_wrap(~Species)+
    scale_fill_gradient(low="white",high="red3") +
    coord_map("polyconic")+
    theme_void()%+replace% theme(legend.text.align=1)
}
p1 <- PV(VirusMapData, HBV)

library(RColorBrewer)
MyRed = '#E41A1C'
MyBlue = '#3A85A8'
MyYellow = "#F2E631"
MyOrange = "#FF7F00"
MyCol = c('white',MyYellow,MyOrange, MyRed)

p1<-ggplot(VirusMapData, aes(x = long, y = lat,group=group,fill=HBV)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  #scale_color_gradient(MyCol(50))+
  scale_fill_gradientn(colors = MyCol, 
                       limits=c(0,0.036), 
                       na.value = 'gray88') +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)
p2<-ggplot(VirusMapData, aes(x = long, y = lat,group=group,fill=HPV)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  #scale_fill_gradient(low="white",high="red3") +
  scale_fill_gradientn(colors = MyCol, 
                       limits=c(0,0.0064), 
                       na.value = 'gray88') +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)
p3<-ggplot(VirusMapData, aes(x = long, y = lat,group=group,fill=EBV)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  scale_fill_gradientn(colors = MyCol) +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)
p4<-ggplot(TorchMapData, aes(x = long, y = lat,group=group,fill=Toxo)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  scale_fill_gradientn(colors = MyCol) +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)
p5<-ggplot(TorchMapData, aes(x = long, y = lat,group=group,fill=CMV)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  scale_fill_gradientn(colors = MyCol) +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)
p6<-ggplot(TorchMapData, aes(x = long, y = lat,group=group,fill=HSV)) +
  geom_polygon(colour="LightSlateGray") +
  #facet_wrap(~Species)+
  scale_fill_gradientn(colors = MyCol) +
  coord_map("polyconic")+
  theme_void()%+replace% theme(legend.text.align=1)

library("gridExtra")
pdf('F4.VirusEachProvince.pdf', width = 9, height = 6)
grid.arrange(p1,p2,p3,
             p4,p5,p6,
             ncol=3,nrow=2)
dev.off()


####################### ####################### ####################### #######################
# Differential abundance analysis
library("limma")
library("edgeR")
library(ggplot2)
require("ggrepel")

load("niptmicrobial.Rda")
count <- as.data.frame(tdb)[3:ncol(tdb)]
count <- as.data.frame(lapply(count,as.integer))
count <- count[apply(count,1,sum)>0,]
count <- t(count)

funDE <- function(group){
  #group <- factor( count['EUKY_26',]>0)
  y <- DGEList(counts=count, group=group)
  keep <- rowSums(cpm(y)>1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  summary(decideTestsDGE(qlf))
  #plotSmear(qlf)
  res<-topTags(qlf,n = length(Spid))$table
  res$significant[res$logFC > 1.5] <- 'UP'
  res$significant[res$logFC < -1.5] <- 'Down'
  res$significant[res$FDR > 0.05] <- 'Unchanged'
  res$significant[is.na(res$significant)] <- 'Unchanged'
  return(res)
  
}
Toxo_edgeR <- funDE(factor(count['EUKY_26',]>0))
CMV_edgeR <- funDE(factor(count['Human herpesvirus 5',]>0))
HSV_edgeR <- funDE(factor(count['Human herpesvirus 1',]>0|count['Human herpesvirus 2',]>0))

HBV_edgeR <- funDE(factor(count['Hepatitis B virus',]>0))
EBV_edgeR <- funDE(factor(count['Human herpesvirus 4 type 1',]>0))
HPV <- c('Human papillomavirus type 4',
         'Human papillomavirus type 63',
         'Human papillomavirus RTRX7',
         'Human papillomavirus type 115',
         'Human papillomavirus - 18',
         'Human papillomavirus type 5',
         'Human papillomavirus type 9',
         'Human papillomavirus type 104',
         'Human papillomavirus type 24',
         'Human papillomavirus type 49',
         'Human papillomavirus type 98',
         'Human papillomavirus SIBX-3a',
         'Human papillomavirus type 48',
         'Human papillomavirus type 105',
         'Human papillomavirus type 71',
         'Human papillomavirus - 1',
         'Human papillomavirus type 121',
         'Human papillomavirus type 10',
         'Human papillomavirus type 96',
         'Human papillomavirus type 116',
         'Human papillomavirus type 88',
         'Human papillomavirus type 60',
         'Human papillomavirus type 112'
)
HPV_edgeR <- funDE(factor(apply(count[HPV,], 2, sum) > 0))

# write.table(Toxo_edgeR,file = 'Toxo_edgeR',quote = F,sep = '\t',row.names = T)
# write.table(CMV_edgeR,file = 'CMV_edgeR',quote = F,sep = '\t',row.names = T)
# write.table(HSV_edgeR,file = 'HSV_edgeR',quote = F,sep = '\t',row.names = T)
# write.table(HBV_edgeR,file = 'HBV_edgeR',quote = F,sep = '\t',row.names = T)
# write.table(EBV_edgeR,file = 'EBV_edgeR',quote = F,sep = '\t',row.names = T)
# write.table(HPV_edgeR,file = 'HPV_edgeR',quote = F,sep = '\t',row.names = T)
save(Toxo_edgeR,CMV_edgeR,HSV_edgeR,HBV_edgeR,EBV_edgeR,HPV_edgeR,annot,
     file=file.path('./','ToRCH-Virus_edgeR.Rda'))

# vocalno fig
load('ToRCH-Virus_edgeR.Rda')
MyRed = '#E41A1C'
MyBlue = '#3A85A8'
options(scipen=200)
Toxo_edgeR$id <-annot$description[match(row.names(Toxo_edgeR),annot$head)]
Toxo_edgeR$id <-  paste(strsplit2(Toxo_edgeR$id,split = ' ')[,1],strsplit2(Toxo_edgeR$id,split = ' ')[,2],sep = ' ')
Toxo_edgeR['BACT_302','id'] <- 'TM7b'
Toxo_edgeR['BACT_301','id'] <- 'TM7a'
#Toxo_edgeR$id[which(Toxo_edgeR$significant == 'No')] <- NA

pToxo<-ggplot(
  Toxo_edgeR[2:nrow(Toxo_edgeR),],
  aes(logFC, -1*log10(PValue))) +
  geom_point(aes(color = significant,
                 size = abs(logFC))) +
  scale_size_continuous(range = c(0,4))+
  scale_y_log10()+
  #geom_point()+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'Toxo')+
  scale_color_manual(values =c(MyBlue, MyRed))+
  geom_text_repel(aes(x = logFC-0.5,label = id),fontface = "italic",size=3.6,color="black",
                  data = subset(Toxo_edgeR[2:nrow(Toxo_edgeR),],significant != 'Unchanged' & FDR ==0))+
  #scale_size_continuous(range = c(0.5,4))+
  xlim(-2,3)+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2,2),linetype=3)+
  theme_bw()+theme(panel.grid=element_blank()) +
  geom_text_repel(
    data = subset(Toxo_edgeR[2:nrow(Toxo_edgeR),], significant != 'Unchanged' & FDR !=0),
    aes(label = id,x = logFC-0.5),
    fontface = "italic",size=3.6,color="black",
    box.padding = unit(0.2, "lines"),
    point.padding = unit(0.2, "lines"))
pToxo

CMV_edgeR$significant[which(CMV_edgeR$logFC<= -2 & CMV_edgeR$FDR< 0.05)] <- 'Down'
CMV_edgeR$significant[which(CMV_edgeR$logFC>= 2 & CMV_edgeR$FDR< 0.05)] <- 'UP'
CMV_edgeR$significant[which(CMV_edgeR$logFC< 2 & CMV_edgeR$logFC> -2)] <- 'Unchanged'
CMV_edgeR$id <-annot$description[match(row.names(CMV_edgeR),annot$head)]
CMV_edgeR$id <-  paste(strsplit2(CMV_edgeR$id,split = ' ')[,1],strsplit2(CMV_edgeR$id,split = ' ')[,2],sep = ' ')
pCMV<-ggplot(
  CMV_edgeR[2:nrow(CMV_edgeR),],
  aes(logFC, -1*log10(PValue))) +
  geom_point(aes(color = significant, size = abs(logFC))) +
  scale_size_continuous(range = c(0,4))+
  xlim(-2,3)+
  scale_y_log10()+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'CMV')+
  scale_color_manual(values =c(MyBlue, MyRed))+
  geom_text(aes(x = logFC,label = id),fontface = "italic",
            data = subset(CMV_edgeR[2:nrow(CMV_edgeR),], significant != 'Unchanged' & FDR !=0))+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2,2),linetype=3)+
  theme_bw()+theme(panel.grid=element_blank()) 
# geom_text_repel(
#   data = subset(CMV_edgeR, significant != 'Unchanged'),
#   aes(label = id,x = logFC),
#   fontface = "italic",
#   box.padding = unit(0.35, "lines"),
#   point.padding = unit(0.3, "lines")
pCMV

HSV_edgeR$significant[which(HSV_edgeR$logFC<= -2 & HSV_edgeR$FDR< 0.05)] <- 'Down'
HSV_edgeR$significant[which(HSV_edgeR$logFC>= 2 & HSV_edgeR$FDR< 0.05)] <- 'UP'
HSV_edgeR$significant[which(HSV_edgeR$logFC< 2 & HSV_edgeR$logFC> -2)] <- 'Unchanged'
HSV_edgeR$id <-annot$description[match(row.names(HSV_edgeR),annot$head)]
HSV_edgeR$id <-  paste(strsplit2(HSV_edgeR$id,split = ' ')[,1],strsplit2(HSV_edgeR$id,split = ' ')[,2],sep = ' ')

pHSV<-ggplot(
  HSV_edgeR[c(2:3,5:nrow(HSV_edgeR)),],
  aes(x = logFC, y=-1*log10(PValue))) +
  geom_point(aes(color = significant, size= abs(logFC)))+
  scale_size_continuous(range = c(0,4))+
  xlim(-2,3)+
  scale_y_log10()+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'HSV')+
  scale_color_manual(values =c(MyBlue, MyRed))+
  #geom_text(aes(x = logFC + 3,label = id),fontface = "italic",
  #data = subset(HSV_edgeR,significant != 'Unchanged' & FDR !=0))+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2,2),linetype=3)+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_text_repel(
    data = subset(HSV_edgeR, significant != 'Unchanged'),
    aes(label = id,x = logFC),
    fontface = "italic",
    #box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
pHSV

require(gridExtra)
pdf('F5.ToRCHVolcano.pdf',height = 14)
grid.arrange(pToxo,pCMV,pHSV,nrow=3)
dev.off()

HBV_edgeR$significant[which(HBV_edgeR$logFC<= -2 & HBV_edgeR$FDR< 0.05)] <- 'Down'
HBV_edgeR$significant[which(HBV_edgeR$logFC>= 2 & HBV_edgeR$FDR< 0.05)] <- 'UP'
HBV_edgeR$significant[which(HBV_edgeR$logFC< 2 & HBV_edgeR$logFC> -2)] <- 'Unchanged'
HBV_edgeR$id <-annot$description[match(row.names(HBV_edgeR),annot$head)]
HBV_edgeR$id <-  paste(strsplit2(HBV_edgeR$id,split = ' ')[,1],strsplit2(HBV_edgeR$id,split = ' ')[,2],sep = ' ')
HBV_edgeR <- HBV_edgeR[2:nrow(HBV_edgeR),]
pHBV<-ggplot(HBV_edgeR,
             aes(logFC,-1*log10(PValue))) +
  geom_point(aes(color =significant, size= abs(logFC)))+
  scale_size_continuous(range = c(0,4))+
  xlim(-2,2)+ #ylim(0,20)+
  scale_y_log10()+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'HBV')+
  scale_color_manual(values =c(MyBlue, "red3"))+
  #geom_text(aes(x = logFC + 3,label = id),fontface = "italic",
  #data = subset(HSV_edgeR,significant != 'Unchanged' & FDR !=0))+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2,2),linetype=3) +
  theme_bw()+theme(panel.grid=element_blank())+
  geom_text_repel(
    data = subset(HBV_edgeR, significant != 'Unchanged'),
    aes(label = id,x = logFC),
    fontface = "italic",
    #box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
pHBV

HPV_edgeR$significant[which(HPV_edgeR$logFC<= -2.5 & HPV_edgeR$FDR< 0.05)] <- 'Down'
HPV_edgeR$significant[which(HPV_edgeR$logFC>= 2.5 & HPV_edgeR$FDR< 0.05)] <- 'UP'
HPV_edgeR$significant[which(HPV_edgeR$logFC< 2.5 & HPV_edgeR$logFC> -2.5)] <- 'Unchanged'
HPV_edgeR$id <-annot$description[match(row.names(HPV_edgeR),annot$head)]
HPV_edgeR$id <-  paste(strsplit2(HPV_edgeR$id,split = ' ')[,1],strsplit2(HPV_edgeR$id,split = ' ')[,2],sep = ' ')
HPV_edgeR <- HPV_edgeR[c(1:2,4:7,9:23,25:39,41:nrow(HPV_edgeR)),]
pHPV<-ggplot(HPV_edgeR,aes(logFC,-1*log10(PValue))) +
  geom_point(aes(color =significant, size= abs(logFC)))+
  scale_size_continuous(range = c(0,4))+
  scale_y_log10()+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'HPV')+
  scale_color_manual(values =c(MyBlue, MyRed))+
  #geom_text(aes(x = logFC + 3,label = id),fontface = "italic",
  #data = subset(HSV_edgeR,significant != 'Unchanged' & FDR !=0))+
  xlim(-2,5)+ #ylim(0,20)+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2.5,2.5),linetype=3)+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_text_repel(
    data = subset(HPV_edgeR, significant != 'Unchanged'),
    aes(label = id,x = logFC),
    fontface = "italic",
    #box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
pHPV

EBV_edgeR$significant[which(EBV_edgeR$logFC<= -2 & EBV_edgeR$FDR< 0.05)] <- 'Down'
EBV_edgeR$significant[which(EBV_edgeR$logFC>= 2 & EBV_edgeR$FDR< 0.05)] <- 'UP'
EBV_edgeR$significant[which(EBV_edgeR$logFC< 2 & EBV_edgeR$logFC> -2)] <- 'Unchanged'
EBV_edgeR$id <-annot$description[match(row.names(EBV_edgeR),annot$head)]
EBV_edgeR$id <-  paste(strsplit2(EBV_edgeR$id,split = ' ')[,1],strsplit2(EBV_edgeR$id,split = ' ')[,2],sep = ' ')
EBV_edgeR <- EBV_edgeR[2:nrow(EBV_edgeR),]
EBV_edgeR['Human herpesvirus 6B','id'] <-'Human herpesvirus 6B'
pEBV<-ggplot(EBV_edgeR,aes(logFC,-1*log10(PValue))) +
  geom_point(aes(color =significant, size= abs(logFC)))+
  scale_size_continuous(range = c(0, 4))+
  labs(x ='log2(Fold Change)',y='-log10(PValue)',title = 'EBV')+
  #scale_y_log10()+
  scale_color_manual(values =c(MyBlue, MyRed))+
  #geom_text(aes(x = logFC + 3,label = id),fontface = "italic",
  #data = subset(HSV_edgeR,significant != 'Unchanged' & FDR !=0))+
  xlim(-2,5)+ #ylim(0,20)+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-2,2),linetype=3)+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_text_repel(
    data = subset(EBV_edgeR, significant != 'Unchanged'),
    aes(label = id,x = logFC),
    fontface = "italic",
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
pEBV

pdf('F5.VirusVolcano.pdf',height = 14)
grid.arrange(pHBV,pHPV,pEBV,nrow=3)
dev.off()


####################### ####################### ####################### #######################
library(dplyr)
library(tidyr)
library(data.table)
library(WGCNA)

load("niptmicrobial.Rda")

# Merge Bact species by max RPKM of each Genus
rpkmTdb = as.data.frame(rpkmTdb)
colnames(rpkmTdb) <- annot$head
bact_annot = annot[grep('BACT', annot$head),]
g = sapply(strsplit(as.character(bact_annot$description),' '), "[", 1)
bact_annot$g = g
merged_bact_rpkm = data.frame(row.names = rownames(rpkmTdb))
for(g in unique(sort(bact_annot$g))){
  message(g)
  g_bact = bact_annot$head[which(bact_annot$g == g)]
  t = as.data.frame(apply(as.data.frame(rpkmTdb[,g_bact]),1,max))
  colnames(t) = g
  merged_bact_rpkm = cbind(merged_bact_rpkm,t)
}

vir_rpkm = rpkmTdb[,grep('BACT|EUKY|ARCH', colnames(rpkmTdb), invert = TRUE)]
merged_rpkm_by_max = cbind(merged_bact_rpkm, vir_rpkm)
bactORvir = rbind(data.frame(head=colnames(merged_bact_rpkm),bact=1),data.frame(head=colnames(vir_rpkm),bact=0))
save(merged_rpkm_by_max,bactORvir,file='merged_by_max.Rda')

# WGCNA
edata <- merged_rpkm_by_max
# check missing and bad sample
gsg<-goodSamplesGenes(edata, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(edata)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(edata)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  edata <-edata[gsg$goodSamples, gsg$goodGenes]
}
#  ..Excluding 3242 genes from the calculation due to too many missing samples or zero variance.

library(flashClust)
nGenes<-ncol(edata)
nSamples<-nrow(edata)
## 1st run
# Choose a set of soft-thresholding powers
powers<-c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft<-pickSoftThreshold(edata, powerVector = powers, verbose = 5)

pdf('network.topology.pdf', paper='special', width=9, height=5)#, horizontal=F)
par(mfrow = c(1,2))
cex1<-0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# minimal power with R^2 > 90%
## for gene optimal power is 12
softPower <- 2

adjacency <- adjacency(edata, power = softPower);
# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Call the hierarchical clustering function
geneTree <- flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf('Gene.cluster.TOM-based.pdf', paper='special', width=12, height=9)#, horizontal=F)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
## for gene the min module size is 30
minModuleSize<-30

# Module identification using dynamic tree cut:
dynamicMods<-cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf('Gene.dendrogram.module.colors.pdf', paper='special', width=8, height=6)#, horizontal=F)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList <- moduleEigengenes(edata, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('Cluster.module.eigengenes.pdf', paper='special', width=7, height=6)#, horizontal=F)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(edata, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
pdf('Gene.dendrogram.module.colors.with.merge.pdf', paper='special', width=12, height=9)#, horizontal=F)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

# save so far
## for gene
save(edata, MEs, moduleLabels, dynamicColors, moduleColors, geneTree, file = file.path("niptmetagenWGCNA.RData"))

# names (colors) of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(edata, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
probes <- names(edata)	
#read in the trait data and match the samples for which they were measured to the expression samples.
geneInfo <- data.frame(probeset_id = probes, moduleColor = moduleColors)
# for external network visualization
#annot = annot[sapply(colnames(edata), function(x) match(x, annot$head)),]
#annot$description[is.na(annot$description)]=annot$head[is.na(annot$description)]
#colnames(edata)<-annot[gsg$goodGenes, 1]

# Select modules
modules = unique(moduleColors[!moduleColors %in% 'grey'])
# Select module probes
probes = names(edata)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = modProbes
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM, 
                               edgeFile = file.path(paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep="")), 
                               nodeFile = file.path(paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep="")), 
                               weighted = TRUE, threshold = 0.02,  altNodeNames= modProbes, nodeNames = colnames(edata)[!moduleColors %in% 'grey'], 
                               nodeAttr = data.frame(bactvir=as.integer(bactORvir$bact[match(colnames(edata), bactORvir$head)][!moduleColors %in% 'grey']), moduleColors[inModule]));

# if hub gene is defined by mod mem...
geneModuleMembership <- as.data.frame(cor(edata, MEs, use = "p"))
rownames(geneModuleMembership)<-probes
write.table(cbind(bactORvir$head[match(colnames(edata), bactORvir$head)], moduleColors, geneModuleMembership), file=file.path("geneModuleMembership.csv"), sep='\t', col.names=T, row.names=F, quote=F)

select_edges = read.csv('CytoscapeInput-edges-blue-turquoise.txt', stringsAsFactors=F, header=T,sep='\t')
vir_row = sapply(select_edges$toNode, function(x) !as.logical(bactORvir$bact[match(x,bactORvir$head)]))
write.table(select_edges[vir_row,], sep = '\t', file = 'to_vir_edges.txt',row.names =F, quote = F)
to_vir = select_edges[vir_row,]
bact_row = sapply(to_vir$fromNode, function(x) as.logical(bactORvir$bact[match(x,bactORvir$head)]))
bact_to_vir = to_vir[bact_row,]
write.table(bact_to_vir , sep = '\t', file = 'bact_to_vir_edges.txt',row.names =F, quote = F)

for(modcol in unique(moduleColors)){
  modname<-paste('MM', modcol, sep='')
  cat('processing module ', modcol, '\n')	
  # Select module probes
  inModule = (moduleColors==modcol);
  modProbes = probes[inModule];
  modN<-sum(inModule)
  
  # output hub gene for each mod
  # this is based on WGCNA's hub gene defination, and could match exactly with cor, but not always with TOM
  nTopH = 10
  IMConn = softConnectivity(edata[, modProbes]);
  hgtopidx<-which(rank(-IMConn)<=nTopH)
  write.table(cbind(modProbes[hgtopidx], bactORvir[match(modProbes[hgtopidx],bactORvir[,1]),]), file = file.path(paste0("wgcnaHubgene-", modcol, "-topH", nTopH, ".txt")), sep='\t', col.names=F, row.names=F, quote=F)
  
  #this is based on plos paper's own hubgene, high module membership
  # this is not always consitant with TOM or cor
  hgtopidx<-which(rank(-abs(geneModuleMembership[inModule,modname]))<=nTopH)
  write.table(annot[match(modProbes[hgtopidx],bactORvir[,1]),], file = file.path(paste0("plosHubgene-", modcol, "-topH", nTopH, ".txt")), sep='\t', col.names=F, row.names=F, quote=F)
}




####################### ####################### ####################### #######################
# sim for baseline prediction poison model lambda
# start resampling
# we generate n length of N idxs and sum them up repectively
# then for each column, we do a pois estimate, the maximum likelyhood estimator is just the mean.

library(parallel)
n <- 1000
N<- seq(1, 1001, length.out = 21) # each sample size N can be considered a function related to sequencing depth, as relatively to the depth of our nipt data, which is typically 0.1X of the human genome.
nrep<-1000
set.seed(1000)
libi<-seq_len(nrow(tdb))
tdbSizeFactors<-tdb$libSize/median(tdb$libSize)
tdbM<-as.matrix(tdb[, -(1:nLeadCol)])
sim<-parallel::mclapply(N, function(sN) {
  message(sN)
  tmp<-do.call(rbind, parallel::mclapply(seq_len(nrep), function(j) {
    message(j, ' j of sN ', sN)
    mergs<-do.call(rbind, lapply(1:n, function(i) {
      message(i, ' i of n ', n)
      ii<-sample(libi, size = sN, replace = T)
      colSums(tdb[ii, -(1:nLeadCol)]*tdbSizeFactors[ii])
    }))
    apply(mergs, 2, function(x) mean(x, na.rm=T))
  },mc.cores=maxCore))
  apply(tmp, 2, function(x) c( mean(x, na.rm=T), sd(x, na.rm=T) ))
}, mc.cores=maxCore)
###
simMean<-do.call(rbind, lapply(sim, function(x) x[1,]))
simSd<-do.call(rbind, lapply(sim, function(x) x[2,]))
simSd<-cbind(simN=N, simSd)
simMean<-cbind(simN=N, simMean)
writeGzFile(simMean,file.path('simPoisLambdaMean.1-1001.n1000.tsv.gz'))
writeGzFile(simSd,file.path('simPoisLambdaSd.1-1001.n1000.tsv.gz'))