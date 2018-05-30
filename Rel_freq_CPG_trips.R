###Relative Frequency of CpG island Triplets#####
# Interested in probabilities of triplet types in WT vs KO 

setwd("~/Year 1 MRC PhD/Rotation 3 Gabrielle/Data")
load('DataWithMars.rData')
#run whole script = source("Rel_freq_CPG_trips.R", echo = TRUE)
library(plyr)
library(ggplot2)
library(scales)
library(karyoploteR)

WT = DD2$WT
KO = DD2$KO


D.max = 20 

idx.N1 = precede(POS,POS) 
idx.N2 = follow(POS,POS)

ii = which(!is.na(idx.N1)&!is.na(idx.N2)) 

D1 = distance(POS[ii],POS[idx.N1[ii]])
D2 = distance(POS[ii],POS[idx.N2[ii]])

ii.2 = which(abs(D1)<D.max & abs(D2)<D.max)
d1 = D1[ii.2]
d2 = D2[ii.2]

#-----------------------------------------------------------Extract WT Triplets+ WT Triplet Frequency---------------------------------------------------------------------------------

WT_DF = data.frame(
  WT.p = WT[idx.N1[ii[ii.2]]],
  WT.0 = WT[ii[ii.2]], 
  WT.n = WT[idx.N2[ii[ii.2]]])
  
WT_DF.t = WT_DF
WT_DF.t[WT_DF < cuts[4]] ='UM'
WT_DF.t[WT_DF > cuts[4]]='IM'
WT_DF.t[WT_DF > cuts[3]]='FM'  

colnames(WT_DF.t) = NULL
TripWT = rep(0,length(ii.2))

##Permutations of states of items 
# n^r = 3^3 #items can repeat states = 27 permuations
#ONLY APPLIES TO WT
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "UM","UM"))] = 'uuu' #1
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "UM","IM"))] = 'uui' #2
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "UM","FM"))] = 'uuf' #3
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "IM","UM"))] = 'uiu' #4
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "IM","IM"))] = 'uii' #5
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "IM","FM"))] = 'uif' #6
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "FM","UM"))] = 'ufu' #7
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "FM","IM"))] = 'ufi' #8
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("UM", "FM","FM"))] = 'uff' #9
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "UM","UM"))] = 'iuu' #10
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "UM","IM"))] = 'iui' #11
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "UM","FM"))] = 'iuf' #12
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "IM","UM"))] = 'iiu' #13
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "IM","IM"))] = 'iii' #14
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "IM","FM"))] = 'iif' #15
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "FM","UM"))] = 'ifu' #16
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "FM","IM"))] = 'ifi' #17
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("IM", "FM","FM"))] = 'iff' #18
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "UM","UM"))] = 'fuu' #19
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "UM","IM"))] = 'fui' #20
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "UM","FM"))] = 'fuf' #21
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "IM","UM"))] = 'fiu' #22
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "IM","IM"))] = 'fii' #23
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "IM","FM"))] = 'fif' #24
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "FM","UM"))] = 'ffu' #25
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "FM","IM"))] = 'ffi' #26
TripWT[apply(WT_DF.t[,c(1:3)], 1, identical, c("FM", "FM","FM"))] = 'fff' #27

colnames(WT_DF.t)=colnames(WT_DF)


lv_WT = levels=c('uuu','uui','uuf','uiu','uii',
                  'uif','ufu','ufi','uff','iuu','iui','iuf','iiu','iii',
                  'iif','ifu','ifi','iff','fuu','fui','fuf','fiu','fii', 'fif', 'ffu', 'ffi', 'fff')

WT_trip = data.frame(WT=WT_DF[,'WT.0'],WT.t=WT_DF.t[,2], type=factor(TripWT,lv_WT))

#Triplet Frequency for WT (Relative Frequency = frequency X/Total Frequency)


WT_sum <- nrow(WT_trip)
WT_T <- data.frame(table(WT_trip$type))
options(scipen = 999)

#WT Triplet Frequencies
WT_freq<- data.frame(WT_T, Rel_Freq = WT_T[,2]/WT_sum)            
CPG_per<- WT_freq[,3]*100
lg<-log10(CPG_per)
WT_freq <- cbind(WT_freq, CPG_per, lg)
#WT_freq =  relative frequency of triplets for all possible CpG triplets

#----------------------------------------------- Extract KO Triplets + KO Triplet Frequency--------------------------------------------------------------------------------------------#
#Need to extract KO Triplets separate --->>>Doesn't work with same dataframe...

KO_DF = data.frame(KO.p = KO[idx.N1[ii[ii.2]]],
                    KO.0 = KO[ii[ii.2]],
                    KO.n = KO[idx.N2[ii[ii.2]]])

KO_DF.t = KO_DF
KO_DF.t[KO_DF < cuts[4]] = 'UM'
KO_DF.t[KO_DF > cuts[4]]='IM'
KO_DF.t[KO_DF > cuts[3]]='FM'  

colnames(KO_DF.t) = NULL
TripKO = rep(0,length(ii.2))


TripKO = rep(0,length(ii.2))
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "UM","UM"))] = 'uuu' #1
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "UM","IM"))] = 'uui' #2
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "UM","FM"))] = 'uuf' #3
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "IM","UM"))] = 'uiu' #4
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "IM","IM"))] = 'uii' #5
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "IM","FM"))] = 'uif' #6
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "FM","UM"))] = 'ufu' #7
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "FM","IM"))] = 'ufi' #8
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("UM", "FM","FM"))] = 'uff' #9
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "UM","UM"))] = 'iuu' #10
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "UM","IM"))] = 'iui' #11
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "UM","FM"))] = 'iuf' #12
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "IM","UM"))] = 'iiu' #13
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "IM","IM"))] = 'iii' #14
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "IM","FM"))] = 'iif' #15
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "FM","UM"))] = 'ifu' #16
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "FM","IM"))] = 'ifi' #17
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("IM", "FM","FM"))] = 'iff' #18
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "UM","UM"))] = 'fuu' #19
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "UM","IM"))] = 'fui' #20
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "UM","FM"))] = 'fuf' #21
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "IM","UM"))] = 'fiu' #22
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "IM","IM"))] = 'fii' #23
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "IM","FM"))] = 'fif' #24
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "FM","UM"))] = 'ffu' #25
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "FM","IM"))] = 'ffi' #26
TripKO[apply(KO_DF.t[,c(1:3)], 1, identical, c("FM", "FM","FM"))] = 'fff' #27
table(TripKO)#check KO levels DIFFERENT FROM WT


lv_KO = levels=c('fff', 'ffi', 'fif', 'fii', 'fiu', 'fui', 
                 'iff', 'ifi', 'ifu', 'iif', 'iii', 'iiu', 
                 'iuf', 'iui', 'iuu', 'uif', 'uii', 'uiu', 'uui', 'uuu')


KO_trip = data.frame(KO= KO_DF[,'KO.0'],K.t= KO_DF.t[,2],type=factor(TripKO,lv_KO))
KO_sum <- nrow(KO_trip)
KO_T <- data.frame(table(KO_trip$type))
options(scipen = 999)

#KO Triplet Frequencies
KO_freq<- data.frame(KO_T, Rel_Freq = KO_T[,2]/KO_sum)            
CPG_per2<- KO_freq[,3]*100
lg2<-log10(CPG_per2)
KO_freq<- cbind(KO_freq, CPG_per2, lg2)
#KO_freq =  relative frequency of triplets for all possible CpG triplets

##---------------------------------------PLOT Different Triplet Frequencies for KO and WT-------------------------------------------------------------------##

#WT_freq and KO_freq for %-long way-----if triplets absent = 0, EXTRACT log vals for KO and WT 
#Long data frame as WT_freq & KO_freq are different lengths

Trips <- c('uuu','uui','uuf','uiu','uii',
           'uif','ufu','ufi','uff','iuu','iui','iuf','iiu','iii',
           'iif','ifu','ifi','iff','fuu','fui','fuf','fiu','fii', 'fif', 'ffu', 'ffi', 'fff')

WT <- c(rbind(WT_freq$Freq[1],
                WT_freq$Freq[2],
                WT_freq$Freq[3],
                WT_freq$Freq[4],
                WT_freq$Freq[5],
                WT_freq$Freq[6],
                WT_freq$Freq[7],
                WT_freq$Freq[8],
                WT_freq$Freq[9],
                WT_freq$Freq[10],
                WT_freq$Freq[11],
                WT_freq$Freq[12],
                WT_freq$Freq[13],
                WT_freq$Freq[14],
                WT_freq$Freq[15],
                WT_freq$Freq[16],
                WT_freq$Freq[17],
                WT_freq$Freq[18],
                WT_freq$Freq[19],
                WT_freq$Freq[20],
                WT_freq$Freq[21],
                WT_freq$Freq[22],
                WT_freq$Freq[23],
                WT_freq$Freq[24],
                WT_freq$Freq[25],
                WT_freq$Freq[26],
                WT_freq$Freq[27]))
        
KO <- c(rbind(KO_freq$Freq[20], 
              KO_freq$Freq[19],
               0,
               KO_freq$Freq[18],
              KO_freq$Freq[17],
              KO_freq$Freq[16],
              0,
              0,
              0,
              KO_freq$Freq[15],
              KO_freq$Freq[14],
              KO_freq$Freq[13],
              KO_freq$Freq[12],
              KO_freq$Freq[11],
              KO_freq$Freq[10],
              KO_freq$Freq[9],
              KO_freq$Freq[8],
              KO_freq$Freq[7],
              0,
               KO_freq$Freq[6],
              0,
              KO_freq$Freq[5],
              KO_freq$Freq[4],
              KO_freq$Freq[3],
              0,
              KO_freq$Freq[2],
              KO_freq$Freq[1]))

data <- data.frame(Trips = Trips, WT = WT, KO = KO)
data1 <- data.frame(Trips = Trips, WT = WT, KO = KO)
Per_diff <- (data$WT-data$KO)/data$WT *100
data1 <- cbind(data, Per_diff)#contains all numbers
drops <- c("WT", "KO")
data1<- data1[, !(names(data1) %in% drops)]

data1$colour <- ifelse(data1$Per_diff < 0, "negative","positive")
trip_plot <- ggplot(data1, aes(Trips, Per_diff)) +
  geom_bar(stat = "identity",position="identity", aes(fill= colour))+ 
  scale_fill_manual(values=c(positive="firebrick1",negative="steelblue"))+
  labs(title = "WT vs KO CpG Coverage",x = "Triplet Types", y = "Percentage Difference(%)", fill = "")+
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
        axis.line.y = element_line(color="black", size = 1))
trip_plot +theme(legend.position="none")

###Check Numbers of pairs
#Rearrange data----
Three <- data %>% filter(Trips %in% c("fff", "uuu", "iii")) #pairs1
UUF_FFU <- data %>% filter(Trips %in% c("uuf", "ffu"))      #pairs2
UII_IIU <- data %>% filter(Trips %in% c("uii", "iiu"))      #pairs3 #SIMILIAR FOR WT&KO
UUI_IUU <- data %>% filter(Trips %in% c("uui", "iuu"))      #pairs4 #SIMILIAR FOR WT&KO
UIU_IUI <- data %>% filter(Trips %in% c("uiu", "iui"))      #pairs5 #SIMILIAR FOR WT&KO
UIF_FIU <- data %>% filter(Trips %in% c("uif", "fiu"))      #pairs6 #SIMILIAR FOR WT&KO
UFU_FUF <- data %>% filter(Trips %in% c("ufu", "fuf"))      #pairs7
UFF_FFU <- data %>% filter(Trips %in% c("uff", "ffu"))      #pairs8  #SIMILIAR FOR WT&KO
UFI_IFU <- data %>% filter(Trips %in% c("ufi", "ifu"))      #pairs9  #SIMILIAR FOR WT&KO
IUF_FUI <- data %>% filter(Trips %in% c("iuf", "fui"))      #pairs10 #SIMILIAR FOR WT&KO
IFF_FII <- data %>% filter(Trips %in% c("iff", "fii"))      #pairs11
IFI_FIF <- data %>% filter(Trips %in% c("ifi", "fif"))      #pairs12
IFF_FFI <- data %>% filter(Trips %in% c("iff", "ffi"))      #pairs13 #SIMILIAR FOR WT&KO



