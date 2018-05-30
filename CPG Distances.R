## CpG distances ##

setwd("~/Year 1 MRC PhD/Rotation 3 Gabrielle/Data")
load('DataWithMars.rData')

library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(methylKit)
library(karyoploteR)
#------------------------------------------------------------Indexing D1 D2-----------------------------------------------------------------------------------#
WT = DD2$WT
KO = DD2$KO
D.max = 20 #only interested in neighbours close by
idx.N1 = precede(POS,POS) #neighbours
idx.N2 = follow(POS,POS)
ii = which(!is.na(idx.N1)&!is.na(idx.N2)) #ignore NaN precede and folllow neghtbours 

D1 = distance(POS[ii],POS[idx.N1[ii]])
D2 = distance(POS[ii],POS[idx.N2[ii]])

ii.2 = which(abs(D1)<D.max & abs(D2)<D.max)
d1 = D1[ii.2]
d2 = D2[ii.2]

Triplets = data.frame(
  WT.p = WT[idx.N1[ii[ii.2]]],
  WT.0 = WT[ii[ii.2]], 
  WT.n = WT[idx.N2[ii[ii.2]]],
  KO.p = KO[idx.N1[ii[ii.2]]],
  KO.0 = KO[ii[ii.2]],
  KO.n = KO[idx.N2[ii[ii.2]]])

Triplets = cbind(Triplets,WT.nh = rowMeans(Triplets[,c('WT.p','WT.n')]),
                 KO.nh = rowMeans(Triplets[,c('KO.p','KO.n')]))

Triplets.t = Triplets
Triplets.t[Triplets < cuts[4]] = 'UM'
Triplets.t[Triplets > cuts[4]]='IM'
Triplets.t[Triplets > cuts[3]]='FM'   #creating a maxtrix 

colnames(Triplets.t) = NULL
TripType = rep(0,length(ii.2))

##Permutations of states of items 
# n^r = 3^3 #items can repeat states = 27 permuations
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "UM","UM"))] = 'uuu' #1
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "UM","IM"))] = 'uui' #2
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "UM","FM"))] = 'uuf' #3
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "IM","UM"))] = 'uiu' #4
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "IM","IM"))] = 'uii' #5
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "IM","FM"))] = 'uif' #6
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "FM","UM"))] = 'ufu' #7
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "FM","IM"))] = 'ufi' #8
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("UM", "FM","FM"))] = 'uff' #9
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "UM","UM"))] = 'iuu' #10
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "UM","IM"))] = 'iui' #11
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "UM","FM"))] = 'iuf' #12
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "IM","UM"))] = 'iiu' #13
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "IM","IM"))] = 'iii' #14
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "IM","FM"))] = 'iif' #15
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "FM","UM"))] = 'ifu' #16
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "FM","IM"))] = 'ifi' #17
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("IM", "FM","FM"))] = 'iff' #18
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "UM","UM"))] = 'fuu' #19
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "UM","IM"))] = 'fui' #20
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "UM","FM"))] = 'fuf' #21
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "IM","UM"))] = 'fiu' #22
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "IM","IM"))] = 'fii' #23
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "IM","FM"))] = 'fif' #24
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "FM","UM"))] = 'ffu' #25
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "FM","IM"))] = 'ffi' #26
TripType[apply(Triplets.t[,c(1:3)], 1, identical, c("FM", "FM","FM"))] = 'fff' #27

levels = levels=c('uuu','uui','uuf','uiu','uii',
                  'uif','ufu','ufi','uff','iuu','iui','iuf','iiu','iii',
                  'iif','ifu','ifi','iff','fuu','fui','fuf','fiu','fii', 'fif', 'ffu', 'ffi', 'fff')

colnames(Triplets.t)=colnames(Triplets)

DD3_Dist_WT = data.frame(WT.t=Triplets.t[,2],type=factor(TripType,levels), D1 = D1[ii.2], D2 = D2[ii.2])# all Triplet Distances WT

#-----------------------------------------------------------Plot Triplet Fequency Distances ------------------------------------------------------------------------------------------------

####-------WT 'fff'- FREQUENCY PLOT & Linear Reg Plot------######
selected <- c('fff')
WT_fff <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fff <- WT_fff[, !(names(WT_fff) %in% drops)]
WT_fff1 <- melt(WT_fff)
P1 <- ggplot(WT_fff1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FFF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,25000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
         panel.grid.major=element_blank(),
         panel.background = element_rect(fill = "white"),
         axis.text.y = element_text(face="bold", size=10),
         axis.text.x  = element_text(face = "bold", size=10),
         title = element_text(face = "bold", size = 14),
         plot.title = element_text(hjust = 0.5),
         axis.title.y = element_text(face="bold", size=14),
         axis.title.x = element_text(face = "bold", size=14),
         axis.line.x = element_line(color="black", size = 1),
         axis.line.y = element_line(color="black", size = 1),
         legend.text = element_text(face = "bold", size = 12))
print(P1)
#no relationship between D1 and D2! = same Distances 



####-------WT 'uuu'- FREQUENCY PLOT------######
selected <- c('uuu')
WT_uuu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uuu <- WT_uuu[, !(names(WT_uuu) %in% drops)]
WT_uuu1 <- melt(WT_uuu)
P2 <- ggplot(WT_uuu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UUU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,500000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P2)

####-------WT 'III'- FREQUENCY PLOT------######
selected <- c('iii')
WT_iii <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iii <- WT_iii[, !(names(WT_iii) %in% drops)]
WT_iii1 <- melt(WT_iii)
P3 <- ggplot(WT_iii1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT III Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,20000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P3)

####-------WT 'UUF'- FREQUENCY PLOT------######
selected <- c('uuf')
WT_uuf <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uuf <- WT_uuf[, !(names(WT_uuf) %in% drops)]
WT_uuf1 <- melt(WT_uuf)
P4 <- ggplot(WT_uuf1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UUF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,20), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P4)


####-------WT 'FUU'- FREQUENCY PLOT------######
selected <- c('fuu')
WT_fuu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fuu <- WT_fuu[, !(names(WT_fuu) %in% drops)]
WT_fuu1 <- melt(WT_fuu)
P5 <- ggplot(WT_fuu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FUU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,20), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P5)
####-------WT 'UII'- FREQUENCY PLOT------######
selected <- c('uii')
WT_uii <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uii <- WT_uii[, !(names(WT_uii) %in% drops)]
WT_uii1 <- melt(WT_uii)
P6 <- ggplot(WT_uii1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UII Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,1500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P6)


####-------WT 'IIU'- FREQUENCY PLOT------######
selected <- c('iiu')
WT_iiu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iiu <- WT_iuu[, !(names(WT_iiu) %in% drops)]
WT_iiu1 <- melt(WT_iiu)
P7 <- ggplot(WT_iiu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IIU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0, 1500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P7)

####-------WT 'UUI'- FREQUENCY PLOT------######
selected <- c('uui')
WT_uui <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uui <- WT_uui[,!(names(WT_uui) %in% drops)]
WT_uui1 <- melt(WT_uui)
P8 <- ggplot(WT_uui1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UUI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0, 1500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P8)


####-------WT 'IUU'- FREQUENCY PLOT------######
selected <- c('iuu')
WT_iuu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iuu <- WT_iuu[, !(names(WT_iuu) %in% drops)]
WT_iuu1 <- melt(WT_iuu)
P9 <- ggplot(WT_iuu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IUU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0, 1500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P9)

####-------WT 'UIU'- FREQUENCY PLOT------######
selected <- c('uiu')
WT_uiu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uiu <- WT_uiu[, !(names(WT_uiu) %in% drops)]
WT_uiu1 <- melt(WT_uiu)
P10 <- ggplot(WT_uiu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UIU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0, 750), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P10)

####-------WT 'IUI'- FREQUENCY PLOT------######
selected <- c('iui')
WT_iui <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iui <- WT_iui[, !(names(WT_iuu) %in% drops)]
WT_iui1 <- melt(WT_iui)
P11 <- ggplot(WT_iui1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IUI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0, 700), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P11)

####-------WT 'UIF'- FREQUENCY PLOT------######
selected <- c('uif')
WT_uif <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uif <- WT_uif[, !(names(WT_uif) %in% drops)]
WT_uif1 <- melt(WT_uif)
P12 <- ggplot(WT_uif1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UIF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,20), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P12)

####-------WT 'FIU'- FREQUENCY PLOT------######
selected <- c('fiu')
WT_fiu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fiu <- WT_fiu[, !(names(WT_fiu) %in% drops)]
WT_fiu1 <- melt(WT_fiu)
P13 <- ggplot(WT_fiu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FIU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,10), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P13)

####-------WT 'UFU'- FREQUENCY PLOT------######
selected <- c('ufu')
WT_ufu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ufu <- WT_ufu[, !(names(WT_ufu) %in% drops)]
WT_ufu1 <- melt(WT_ufu)
P14 <- ggplot(WT_ufu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UFU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,5), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P14)

####-------WT 'FUF'- FREQUENCY PLOT------######
selected <- c('fuf')
WT_fuf <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fuf <- WT_fuf[, !(names(WT_fuf) %in% drops)]
WT_fuf1 <- melt(WT_fuf)
P15 <- ggplot(WT_fuf1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FUF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,10), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P15)

####-------WT 'UFF'- FREQUENCY PLOT------######
selected <- c('uff')
WT_uff <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_uff <- WT_uff[, !(names(WT_uff) %in% drops)]
WT_uff1 <- melt(WT_uff)
P16 <- ggplot(WT_uff1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UFF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,100), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P16)

####-------WT 'FFU'- FREQUENCY PLOT------######
selected <- c('ffu')
WT_ffu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ffu <- WT_ffu[, !(names(WT_ffu) %in% drops)]
WT_ffu1 <- melt(WT_ffu)
P17 <- ggplot(WT_ffu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FFU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,100), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P17)

####-------WT 'UFI'- FREQUENCY PLOT------######
selected <- c('ufi')
WT_ufi <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ufi <- WT_ufi[, !(names(WT_ufi) %in% drops)]
WT_ufi1 <- melt(WT_ufi)
P18 <- ggplot(WT_ufi1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT UFI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,15), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P18)

####-------WT 'IFU'- FREQUENCY PLOT------######
selected <- c('ifu')
WT_ifu <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ifu <- WT_ifu[, !(names(WT_ifu) %in% drops)]
WT_ifu1 <- melt(WT_ifu)
P19 <- ggplot(WT_ifu1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IFU Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,20), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P19)

####-------WT 'IUF'- FREQUENCY PLOT------######
selected <- c('iuf')
WT_iuf <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iuf <- WT_iuf[, !(names(WT_iuf) %in% drops)]
WT_iuf1 <- melt(WT_iuf)
P20 <- ggplot(WT_iuf1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IUF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,4), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P20)

####-------WT 'FUI'- FREQUENCY PLOT------######
selected <- c('fui')
WT_fui <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fui <- WT_fui[, !(names(WT_fui) %in% drops)]
WT_fui1 <- melt(WT_fui)
P21 <- ggplot(WT_fui1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FUI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,5), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P21)

####-------WT 'IIF'- FREQUENCY PLOT------######
selected <- c('iif')
WT_iif <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iif <- WT_iif[, !(names(WT_iif) %in% drops)]
WT_iif1 <- melt(WT_iif)
P22 <- ggplot(WT_iif1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IIF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,3000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P22)

####-------WT 'FII'- FREQUENCY PLOT------######
selected <- c('fii')
WT_fii <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fii <- WT_fii[, !(names(WT_fii) %in% drops)]
WT_fii1 <- melt(WT_fii)
P23 <- ggplot(WT_fii1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FII Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,3000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P23)

####-------WT 'IFI'- FREQUENCY PLOT------######
selected <- c('ifi')
WT_ifi <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ifi <- WT_ifi[, !(names(WT_ifi) %in% drops)]
WT_ifi1 <- melt(WT_ifi)
P24 <- ggplot(WT_ifi1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IFI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,1500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P24)

####-------WT 'FIF'- FREQUENCY PLOT------######
selected <- c('fif')
WT_fif <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_fif <- WT_fif[, !(names(WT_fif) %in% drops)]
WT_fif1 <- melt(WT_fif)
P25 <- ggplot(WT_fif1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FIF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,2500), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P25)

####-------WT 'IFF'- FREQUENCY PLOT------######
selected <- c('iff')
WT_iff <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_iff <- WT_iff[, !(names(WT_iff) %in% drops)]
WT_iff1 <- melt(WT_iff)
P26 <- ggplot(WT_iff1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT IFF Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,4000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P26)

####-------WT 'FFI'- FREQUENCY PLOT------######
selected <- c('ffi')
WT_ffi <- DD3_Dist_WT[DD3_Dist_WT$type %in% selected,]
drops <- c("WT.t", "type")
WT_ffi <- WT_ffi[, !(names(WT_ffi) %in% drops)]
WT_ffi1 <- melt(WT_ffi)
P27 <- ggplot(WT_ffi1, aes(value, fill=variable)) + 
  geom_histogram(position="dodge", colour="white",binwidth=2) +
  labs(title = "WT FFI Triplet Distances",x = "Distance", y = "Frequency", fill = "") +
  coord_cartesian( xlim=c(0,20), expand = FALSE )+
  coord_cartesian( ylim=c(0,4000), expand = FALSE ) + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(face="bold", size=10),
        axis.text.x  = element_text(face = "bold", size=10),
        title = element_text(face = "bold", size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.text = element_text(face = "bold", size = 12))
print(P27)

######################################---------KAROTYPES Plots--------################################################################
#Interested in triplet pairs -EXHIBIT similiar relationships in D1 and D2, why?
#1. IIF & FII 4000>Granges
#2. IIU & UII
#3. FFI & IFF
#4. UUI & IUU
#5. FFU & UFF

#Distances suggest these triplets could be close to one another- based on histograms- 
# Aim - plot karotype plot for each pair types (1-5) to determine any distance patterns??

POS.0 <- POS[ii[ii.2]] #only middle CpGs nrow matches DD3_Dist_WT

####----- Pair 1 - IIF and FII ,rows same as DD3 Extract Granges
#iff
row_iif <- rownames(WT_iif)
class(row_iif) <- "numeric"
GR_WT_iif<- POS.0[row_iif]
seqlevels(GR_WT_iif)
#Keep norm Chrs
keep_seq = c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",     
             "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7",      
             "chr8","chr9","chrX","chrY")
GR_WT_iif <-keepSeqlevels(GR_WT_iif, keep_seq, pruning.mode = "coarse")

#fii
row_fii <- rownames(WT_fii)
class(row_fii) <- "numeric"
GR_WT_fii<- POS.0[row_fii]
seqlevels(GR_WT_fii)
GR_WT_fii <-keepSeqlevels(GR_WT_fii, keep_seq, pruning.mode = "coarse")

kp1I <- plotKaryotype(genome="mm10")
kpPlotDensity(kp1I, data = GR_WT_iif, col="#FF0000", r0=0, r1=1)
kp1I <-legend("right",legend=c("IIF"), fill=c("red"), par(xpd=FALSE)+
               theme(axis.text.y = element_text(face="bold", size=13),
                     legend.text = element_text(face = "bold", size = 12)))
print(kp1I)

kp1II <- plotKaryotype(genome="mm10")
kpPlotDensity(kp1II, data = GR_WT_fii, col="#0000FF", r0=0, r1=1)
kp1II <-legend("right",legend=c("FII"), fill=c("blue"), par(xpd=FALSE)+
               theme(axis.text.y = element_text(face="bold", size=13),
                     legend.text = element_text(face = "bold", size = 12)))
print(kp1II)

#TOO many CpGs so you cannot see anything clearly.....

####----- Pair 2 - IIU and UII ,rows same as DD3 Extract Granges
#IIU
row_iiu <- rownames(WT_iiu)
class(row_iiu) <- "numeric"
GR_WT_iiu<- POS.0[row_iiu]
seqlevels(GR_WT_iiu)
#Keep norm Chrs
keep_seq = c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",     
             "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7",      
             "chr8","chr9","chrX","chrY")
GR_WT_iiu <-keepSeqlevels(GR_WT_iiu, keep_seq, pruning.mode = "coarse")

#UII
row_uii <- rownames(WT_uii)
class(row_uii) <- "numeric"
GR_WT_uii<- POS.0[row_uii]
seqlevels(GR_WT_uii)
GR_WT_uii <-keepSeqlevels(GR_WT_uii, keep_seq, pruning.mode = "coarse")

#Do separate plots as superimpose one another 
kp2I <- plotKaryotype(genome="mm10")
kpPlotDensity(kp2I, data = GR_WT_iiu, col="#FF0000", r0=0, r1=1)
kp2I <-legend("right",legend=c("IIU"), fill=c("red"), par(xpd=FALSE)+
                theme(axis.text.y = element_text(face="bold", size=13),
                      legend.text = element_text(face = "bold", size = 12)))
print(kp2I)

kp2II <- plotKaryotype(genome="mm10")
kpPlotDensity(kp2II, data = GR_WT_uii, col="#0000FF", r0=0, r1=1)
kp2II <-legend("right",legend=c("UII"), fill=c("blue"), par(xpd=FALSE)+
                 theme(axis.text.y = element_text(face="bold", size=13),
                       legend.text = element_text(face = "bold", size = 12)))
print(kp2II)

####----- Pair 3 - FFI & IFF,rows same as DD3 Extract Granges
#FFI
row_ffi <- rownames(WT_ffi)
class(row_ffi) <- "numeric"
GR_WT_ffi<- POS.0[row_ffi]
seqlevels(GR_WT_ffi)
#Keep norm Chrs
keep_seq = c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",     
             "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7",      
             "chr8","chr9","chrX","chrY")
GR_WT_ffi <-keepSeqlevels(GR_WT_ffi, keep_seq, pruning.mode = "coarse")

#IIF
row_iif <- rownames(WT_iif)
class(row_iif) <- "numeric"
GR_WT_iif<- POS.0[row_iif]
seqlevels(GR_WT_iif)
GR_WT_iif <-keepSeqlevels(GR_WT_iif, keep_seq, pruning.mode = "coarse")

#Do separate plots as superimpose one another 
kp3I <- plotKaryotype(genome="mm10")
kpPlotDensity(kp3I, data = GR_WT_ffi, col="#FF0000", r0=0, r1=1)
kp3I <-legend("right",legend=c("FFI"), fill=c("red"), par(xpd=FALSE)+
                theme(axis.text.y = element_text(face="bold", size=13),
                      legend.text = element_text(face = "bold", size = 12)))
print(kp3I)

kp3II <- plotKaryotype(genome="mm10")
kpPlotDensity(kp3II, data = GR_WT_uii, col="#0000FF", r0=0, r1=1)
kp3II <-legend("right",legend=c("IFF"), fill=c("blue"), par(xpd=FALSE)+
                 theme(axis.text.y = element_text(face="bold", size=13),
                       legend.text = element_text(face = "bold", size = 12)))
print(kp3II)

####----- Pair 4 - UUI & IUU,rows same as DD3 Extract Granges
#UUI
row_uui <- rownames(WT_uui)
class(row_uui) <- "numeric"
GR_WT_uui<- POS.0[row_uui]
seqlevels(GR_WT_uui)
#Keep norm Chrs
keep_seq = c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",     
             "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7",      
             "chr8","chr9","chrX","chrY")
GR_WT_uui <-keepSeqlevels(GR_WT_uui, keep_seq, pruning.mode = "coarse")

#IUU
row_iuu <- rownames(WT_iuu)
class(row_iuu) <- "numeric"
GR_WT_iuu<- POS.0[row_iuu]
seqlevels(GR_WT_iuu)
GR_WT_iuu <-keepSeqlevels(GR_WT_iuu, keep_seq, pruning.mode = "coarse")

#Do separate plots as superimpose one another 
kp4I <- plotKaryotype(genome="mm10")
kpPlotDensity(kp4I, data = GR_WT_uui, col="#FF0000", r0=0, r1=1)
kp4I <-legend("right",legend=c("UUI"), fill=c("red"), par(xpd=FALSE)+
                theme(axis.text.y = element_text(face="bold", size=13),
                      legend.text = element_text(face = "bold", size = 12)))
print(kp4I)

kp4II <- plotKaryotype(genome="mm10")
kpPlotDensity(kp4II, data = GR_WT_iiu, col="#0000FF", r0=0, r1=1)
kp4II <-legend("right",legend=c("IIU"), fill=c("blue"), par(xpd=FALSE)+
                 theme(axis.text.y = element_text(face="bold", size=13),
                       legend.text = element_text(face = "bold", size = 12)))
print(kp4II)
####----- Pair 5 - FFU and UUF ,rows same as DD3 Extract Granges
#FFU
row_ffu <- rownames(WT_ffu)
class(row_ffu) <- "numeric"
GR_WT_ffu<- POS.0[row_ffu]
seqlevels(GR_WT_ffu)#check
GR_WT_ffu <-keepSeqlevels(GR_WT_ffu, keep_seq, pruning.mode = "coarse")

#UUF
row_uff <- rownames(WT_uff)
class(row_uff) <- "numeric"
GR_WT_uff<- POS.0[row_uff]
seqlevels(GR_WT_uff)
GR_WT_uff <-keepSeqlevels(GR_WT_uff, keep_seq, pruning.mode = "coarse")

kp5 <- plotKaryotype(genome="mm10")
kpPlotDensity(kp5, data = GR_WT_ffu, col="#FF0000", r0=0, r1=1)
kpPlotDensity(kp5, data = GR_WT_uff, col="#0000FF", r0=0, r1=1)
print(kp5)
kpPlotRegions(kp5, data = GR_WT_ffu, r=0, r1=1)
kpPlotRegions(kp5, data = GR_WT_uff, r=0, r1=1)
kp5 <-legend("right",legend=c("FFU", "UUF"), fill=c("red", "blue"), par(xpd=FALSE)+
theme(axis.text.y = element_text(face="bold", size=13),
             legend.text = element_text(face = "bold", size = 12)))
print(kp5)
