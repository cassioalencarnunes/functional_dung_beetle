#### Script to run simulations of species loss of Undisturbed Forests ####

## This script is adapted from scripts created by Leitão et al. (2016).

## The first part of the script runs the simulations of species loss at the regional scale
## The second part of the script runs the simulations of local species loss

# LEITÃO, R.P., J. ZUANON, S. VILLÉGER, S.E. WILLIAMS, C. BARALOTO, C. FORTUNEL, F.P. MENDONÇA,
# and D. MOUILLOT. 2016. Rare species contribute disproportionately to the functional structure
# of species assemblages. Proc. R. Soc. B Biol. Sci. 283: 20160084. 
# Available at: https://doi.org/10.1098/rspb.2016.0084



### List of packages used
library(geometry)


## Loading data

Local.pa <- read.csv("matrix_species.csv", h=T, row.names = 1)
RI <- read.csv("rarity_all_species.csv", h=T)
Axes <- read.csv("PCoA_axes.csv", h=T, row.names = 1)

RI_UF<-RI[RI$sp %in% row.names(Axes),]
Axes<-Axes[order(row.names(Axes)),]

##############################################################################################################################

#### Regional scale ####

## Functional Richness (FRic)

# Input: PCoA axes + RI

ERS<-cbind(Axes,RI_UF$RI)
ERS<-as.data.frame(ERS)
names(ERS)<-c("Axis.1","Axis.2","Axis.3","RI_UF")

# Number of species (-4 = 3 axes + 1)
NoS<-nrow(Local)-4

# Scenario losing rarest species first

Simulat_rare<-ERS[order(ERS[,"RI_UF"]),]
tr_rare<-Simulat_rare[,1:NoD]
FRic_Simulat_rare=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  tr_rare=tr_rare[-1,]
  FRic_Simulat_rare[i]<-convhulln(tr_rare,"FA")$vol
}

# Scenario losing most common species first

Simulat_comm<-ERS[order(ERS[,"RI_UF"],decreasing=T),]
tr_comm<-Simulat_comm[,1:NoD]
FRic_Simulat_comm=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  tr_comm=tr_comm[-1,]
  FRic_Simulat_comm[i]<-convhulln(tr_comm,"FA")$vol
}

# Scenario losing species randomly (null model)

dataperm=ERS
randFRic=matrix(NA,nrow=NoS,ncol=1000)  #ncol= n randomizations
for (j in 1:1000){
  dataperm[,"RI_UF"]=sample(dataperm[,"RI_UF"])
  Null_Simulat<-dataperm[order(dataperm[,"RI_UF"]),]
  Null_tr<-Null_Simulat[,1:NoD]
  for(i in 1:NoS){
    Null_tr=Null_tr[-1,]
    randFRic[i,j]<-convhulln(Null_tr,"FA")$vol
    }
}

## Calculating median and quantiles from the null models
Null_Fric=matrix(NA,nrow=NoS,ncol=3)
for(i in 1:NoS)
{
  Null_Fric[i,1]=median(randFRic[i,])
  Null_Fric[i,2]=quantile(randFRic[i,],probs=0.025)
  Null_Fric[i,3]=quantile(randFRic[i,],probs=0.975)
}

# Plot regional species loss scenarios - FRic
Spp_erosion<-c(1:NoS)
plot(Spp_erosion,Null_Fric[,1]/max(Null_Fric[,1]),type="l",xlim=c(0,NoS),xlab="Regional Species
     loss (%)",ylab="FRic",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n")
axis(side=1,at=c(0,NoS/4,NoS/2,NoS/1.333,NoS),labels=c(0,25,50,75,100),line=F,tick=-0.3,cex.axis=1.2,
     mgp=c(3,1,0))
polygon(c(1:NoS,rev(1:NoS)),c(Null_Fric[,2]/max(Null_Fric[,2]),rev(Null_Fric[,3]/max(Null_Fric))),col="grey88",border=F)     
points(Null_Fric[,1]/max(Null_Fric[,1]),col="gray",cex=1.6,type="l",lty=1,lwd=4)
points(FRic_Simulat_rare/max(FRic_Simulat_rare),col="black",cex=1.6,type="l",lwd=4)
points(FRic_Simulat_comm/max(FRic_Simulat_comm),col="black",cex=1.6,type="l",lwd=4,lty=3)

##############################################################################################################################

## Functional Specialization (FSpe)

# Computing functional specialization for each species (speS)
O<-apply(Axes,2,mean)
speS<-apply(Axes,1,function(x){(sum((x-O)^2))^0.5})
speS<-speS/max(speS)
speS<-as.data.frame(speS)

ERS<-cbind(RI_UF$RI,speS)
names(ERS)<-c("RI_UF","speS")

# Number of species
NoS<-nrow(Local)-1

# Scenario losing rarest species first

Simulat_rare<-ERS[order(ERS[,"RI_UF"]),]
FSpe_Simulat_rare=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  Simulat_rare=Simulat_rare[-1,]
  FSpe_Simulat_rare[i]<-mean(Simulat_rare[,"speS"])
}

# Scenario losing most common species first

Simulat_comm<-ERS[order(ERS[,"RI_UF"],decreasing=T),]
FSpe_Simulat_comm=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  Simulat_comm=Simulat_comm[-1,]
  FSpe_Simulat_comm[i]<-mean(Simulat_comm[,"speS"])
}

# Scenario losing species randomly (null model)

dataperm=ERS
randSPE=matrix(NA,nrow=NoS,ncol=1000)  #ncol= n randomizations
for (j in 1:1000){
  dataperm[,1]=sample(dataperm[,1])
  Null_Simulat_rare<-dataperm[order(dataperm[,"RI_UF"]),]
  for(i in 1:NoS){
    Null_Simulat_rare=Null_Simulat_rare[-1,]
    randSPE[i,j]<-mean(Null_Simulat_rare[,"speS"])
  }
}

## Calculating median and quantiles from the null models
Null_FSpe=matrix(NA,nrow=NoS,ncol=3)
for(i in 1:NoS){
  Null_FSpe[i,1]=median(randSPE[i,])
  Null_FSpe[i,2]=quantile(randSPE[i,],probs=0.025)
  Null_FSpe[i,3]=quantile(randSPE[i,],probs=0.975)
}

# Plot regional species loss scenarios - FSpe
Spp_erosion<-c(1:NoS)
plot(Spp_erosion,Null_FSpe[,1],type="l",xlim=c(0,NoS),xlab="Regional Species
     loss (%)",ylab="FSpe",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n",
     ylim=c(0,1))
axis(side=1,at=c(0,NoS/4,NoS/2,NoS/1.333,NoS),labels=c(0,25,50,75,100),line=F,tick=-0.3,cex.axis=1.2,
     mgp=c(3,1,0))
polygon(c(1:NoS,rev(1:NoS)),c(Null_FSpe[,2],rev(Null_FSpe[,3])),col="grey88",border=F)     
points(Null_FSpe[,1],col="gray",cex=1.6,type="l",lty=1,lwd=4)
points(FSpe_Simulat_rare,col="black",cex=1.6,type="l",lwd=4)
points(FSpe_Simulat_comm,col="black",cex=1.6,type="l",lwd=4,lty=3)

##############################################################################################################################

## Functional Originality (FOri)

# Computing functional originality for each species (oriS)
dist_F<-as.matrix(dist(Axes, method="euclidean")) ; dist_F[which(dist_F==0)]<-NA
oriS<-apply(dist_F, 1, min, na.rm=T)
oriS<-oriS/max(oriS)
oriS<-as.data.frame(oriS)

ERS<-cbind(RI_UF$RI,oriS)
names(ERS)<-c("RI_UF","oriS")

# Scenario losing rarest species first

Simulat_rare<-ERS[order(ERS[,"RI_UF"]),]
FOri_Simulat_rare=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  Simulat_rare=Simulat_rare[-1,]
  FOri_Simulat_rare[i]<-mean(Simulat_rare[,"oriS"])
}

#Scenario losing most common species first

Simulat_comm<-ERS[order(ERS[,"RI_UF"],decreasing=T),]
FOri_Simulat_comm=matrix(NA,nrow=NoS,ncol=1)
for(i in 1:NoS){
  Simulat_comm=Simulat_comm[-1,]
  FOri_Simulat_comm[i]<-mean(Simulat_comm[,"oriS"])
}

#Scenario losing species randomly (null model)

dataperm=ERS
randORI=matrix(NA,nrow=NoS,ncol=1000)  #ncol= n randomizations
for (j in 1:1000){
  dataperm[,1]=sample(dataperm[,1])
  Null_Simulat_rare<-dataperm[order(dataperm[,"RI_UF"]),]
  for(i in 1:NoS)
  {
    Null_Simulat_rare=Null_Simulat_rare[-1,]
    randORI[i,j]<-mean(Null_Simulat_rare[,"oriS"])
  }
}

## Calculating median and quantiles from the null models

Null_FOri=matrix(NA,nrow=NoS,ncol=3)
for(i in 1:NoS){
  Null_FOri[i,1]=median(randORI[i,])
  Null_FOri[i,2]=quantile(randORI[i,],probs=0.025)
  Null_FOri[i,3]=quantile(randORI[i,],probs=0.975)
}

# Plot regional species loss scenarios - FSpe
Spp_erosion<-c(1:NoS)
plot(Spp_erosion,Null_FOri[,1],type="l",xlim=c(0,NoS),xlab="Regional Species
     loss (%)",ylab="FOri",font.lab=2,cex.lab=1.4,cex.axis=1.2,col="gray",pch=16,cex=1.5,xaxt="n",
     ylim=c(0,1))
axis(side=1,at=c(0,NoS/4,NoS/2,NoS/1.333,NoS),labels=c(0,25,50,75,100),line=F,tick=-0.3,cex.axis=1.2,
     mgp=c(3,1,0))
polygon(c(1:NoS,rev(1:NoS)),c(Null_FOri[,2],rev(Null_FOri[,3])),col="grey88",border=F)     
points(Null_FOri[,1],col="gray",cex=1.6,type="l",lty=1,lwd=4)
points(FOri_Simulat_rare,col="black",cex=1.6,type="l",lwd=4)
points(FOri_Simulat_comm,col="black",cex=1.6,type="l",lwd=4,lty=3)

##############################################################################################################################

##### Local scale #####

std.error<-function(x){
  sd(x)/sqrt(length(x))}

# Loading the functions to calculate functional indices and simulate species loss

source("indices_functions.R")


### Functional Richness (FRic) ###

# Input: Local-assemblage matrix + RI + PCoA axes

Local.pa<-Local.pa[order(row.names(Local.pa)),]

fullMat<-cbind(Local.pa,RI_UF$RI,Axes)
fullMat<-as.data.frame(fullMat)
names_mat<-names(fullMat)
names_mat[13]<-"RI"
names(fullMat)<-names_mat


# Scenarios losing most rare and most common species

## Rare

FRic_Local_RARE<-c(mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="rare"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="rare"),na.rm=T))

FRic_Local_RARE_SE<-c(std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="rare")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="rare")))


## Common

FRic_Local_COMM<-c(mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="common"),na.rm=T),
                   mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="common"),na.rm=T))

FRic_Local_COMM_SE<-c(std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="common")),
                      std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="common")))

# Scenario losing species randomly

FRic_Local_mean<-matrix(NA,nrow=7,ncol=1000) # ncol = n randomizations

for(i in 1:1000){
  FRic_Local_mean[,i]<-c(mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="random"),na.rm=T),
                         mean(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="random"),na.rm=T))
}


FRic_Local_se<-matrix(NA,nrow=7,ncol=1000) # ncol = n randomizations
for(i in 1:1000){
  FRic_Local_se[,i]<-c(std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.1,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.2,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.3,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.4,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.5,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.6,rule="random")),
                       std.error(calcFRic(fullMat=fullMat,nsite=ncol(fullMat)-NoD-1,prob=0.7,rule="random")))
}

FRic_Local_NULL_mean<-rowMeans(FRic_Local_mean)
FRic_Local_NULL_SE<-rowMeans(FRic_Local_se)

# Plot local species loss scenarios - FRic

erosion1<-seq(0,70,10)
plot(erosion1,c(1,FRic_Local_NULL_mean),type="l",ylim=c(0,1),xlim=c(0,70),col="gray",pch=16,cex=1.5,
     xlab="Local species loss (%)",ylab="FRic",font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion1[-c(1,11)],FRic_Local_NULL_mean,type="p",col="black",pch=21,bg="gray",cex=1)
arrows(x0 = erosion1, x1=erosion1, y0=c(0,FRic_Local_NULL_mean-FRic_Local_NULL_SE), y1=c(0,FRic_Local_NULL_mean
                                                                                         +FRic_Local_NULL_SE),code=0,angle=90,length = 0.05, col="gray")
points(erosion1,c(1,FRic_Local_RARE),type="l",ylim=c(0,1),col="black")
points(erosion1[-c(1,11)],FRic_Local_RARE,type="p",col="black",pch=21,bg="black",cex=1)
points(erosion1,c(1,FRic_Local_COMM),type="l",ylim=c(0,1),col="black",lty=2)
points(erosion1[-c(1,11)],FRic_Local_COMM,type="p",col="black",pch=21,bg="white",cex=1)
arrows(x0 = erosion1, x1=erosion1, y0=c(0,FRic_Local_RARE-FRic_Local_RARE_SE), y1=c(0,FRic_Local_RARE
                                                                                    +FRic_Local_RARE_SE),code=0,angle=90,length = 0.05)
arrows(x0 = erosion1, x1=erosion1, y0=c(0,FRic_Local_COMM-FRic_Local_COMM_SE), y1=c(0,FRic_Local_COMM
                                                                                    +FRic_Local_COMM_SE),code=0,angle=90,length = 0.05)
legend(10,0.2,legend=c("common","random","rare"),col=c("black","gray","black"),
       cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,
       bty="n")

##############################################################################################################################

### Functional Specialization (FSpe) ###

# Computing functional specialization for each species (speS)
O<-apply(Axes,2,mean)
speS<-apply(Axes,1,function(x){(sum((x-O)^2))^0.5})
speS<-speS/max(speS)
speS<-as.data.frame(speS)

# Input: Local-assemblage matrix + RI + speS
fullMat<-cbind(Local.pa,RI_UF$RI,speS)
names_mat<-names(fullMat)
names_mat[13]<-"RI"
names(fullMat)<-names_mat

# Scenario losing most rare and most common species

## Rare

FSpe_Local_RARE<-c(mean(calcSPE_total(fullMat=fullMat,nsite=12),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="rare"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="rare"),na.rm=T))

FSpe_Local_RARE_SE<-c(std.error(calcSPE_total(fullMat=fullMat,nsite=12)),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="rare")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="rare")))

## Common

FSpe_Local_COMM<-c(mean(calcSPE_total(fullMat=fullMat,nsite=12),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="common"),na.rm=T),
                   mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="common"),na.rm=T))

FSpe_Local_COMM_SE<-c(std.error(calcSPE_total(fullMat=fullMat,nsite=12)),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="common")),
                      std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="common")))

# Scenario losing species randomly

FSpe_Local_mean<-matrix(NA,nrow=10,ncol=1000) # ncol = n randomizations

for(i in 1:1000){
  FSpe_Local_mean[,i]<-c(mean(calcSPE_total(fullMat=fullMat,nsite=12),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="random"),na.rm=T),
                         mean(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="random"),na.rm=T))
}

FSpe_Local_se<-matrix(NA,nrow=10,ncol=1000) # ncol = n randomizations

for(i in 1:1000){
  FSpe_Local_se[,i]<-c(std.error(calcSPE_total(fullMat=fullMat,nsite=12)),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.1,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.2,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.3,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.4,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.5,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.6,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.7,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.8,rule="random")),
                       std.error(calcSPE(fullMat=fullMat,nsite=12,prob=0.9,rule="random")))
}


FSpe_Local_NULL_mean<-rowMeans(FSpe_Local_mean)
FSpe_Local_NULL_SE<-rowMeans(FSpe_Local_se)

#Plot local species loss scenarios - FSpe
erosion<-seq(0,90,10)
plot(erosion,FSpe_Local_NULL_mean,type="l",ylim=c(0.20,0.52),
     xlim=c(0,100),col="gray",pch=16,cex=1.5,xlab="Local species loss (%)",ylab="FSpe",
     font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion[-1],FSpe_Local_NULL_mean[-1],type="p",col="gray",pch=21,bg="gray",cex=1)
arrows(x0 = erosion, x1=erosion, y0=FSpe_Local_NULL_mean-FSpe_Local_NULL_SE, y1=FSpe_Local_NULL_mean
       +FSpe_Local_NULL_SE,code=0,angle=90,length = 0.05, col="gray")
points(erosion,FSpe_Local_RARE,type="l",ylim=c(0,1),col="black")
points(erosion[-1],FSpe_Local_RARE[-1],type="p",col="black",pch=21,bg="black",cex=1)
points(erosion,FSpe_Local_COMM,type="l",ylim=c(0,1),col="black",lty=2)
points(erosion[-1],FSpe_Local_COMM[-1],type="p",col="black",pch=21,bg="white",cex=1)
arrows(x0 = erosion, x1=erosion, y0=FSpe_Local_RARE-FSpe_Local_RARE_SE, y1=FSpe_Local_RARE
       +FSpe_Local_RARE_SE,code=0,angle=90,length = 0.05)
arrows(x0 = erosion, x1=erosion, y0=FSpe_Local_COMM-FSpe_Local_COMM_SE, y1=FSpe_Local_COMM
       +FSpe_Local_COMM_SE,code=0,angle=90,length = 0.05)
legend(0,0.28,legend=c("common","random","rare"),col=c("black","gray","black"),
       cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,bty="n")

##############################################################################################################################

## Functional Originality (FOri)

#Computing functional originality for each species (speS)
dist_F<-as.matrix(dist(Axes, method="euclidean")) ; dist_F[which(dist_F==0)]<-NA
oriS<-apply(dist_F, 1, min, na.rm=T)
oriS<-oriS/max(oriS)
oriS<-as.data.frame(oriS)

#Input: Local-assemblage matrix + RI + oriS
fullMat<-cbind(Local.pa,RI_UF$RI,oriS)
names_mat<-names(fullMat)
names_mat[13]<-"RI"
names(fullMat)<-names_mat

# Scenario losing most rare and most common species

## Rare

FOri_Local_RARE<-c(mean(calcORI_total(fullMat=fullMat,nsite=12),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="rare"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="rare"),na.rm=T))

FOri_Local_RARE_SE<-c(std.error(calcORI_total(fullMat=fullMat,nsite=12)),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="rare")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="rare")))

## Common

FOri_Local_COMM<-c(mean(calcORI_total(fullMat=fullMat,nsite=12),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="common"),na.rm=T),
                   mean(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="common"),na.rm=T))

FOri_Local_COMM_SE<-c(std.error(calcORI_total(fullMat=fullMat,nsite=12)),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="common")),
                      std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="common")))

# Scenario losing species randomly

FOri_Local_mean<-matrix(NA,nrow=10,ncol=1000) # ncol = n randomizations

for(i in 1:1000){
  FOri_Local_mean[,i]<-c(mean(calcORI_total(fullMat=fullMat,nsite=12),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="random"),na.rm=T),
                         mean(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="random"),na.rm=T))
}

FOri_Local_se<-matrix(NA,nrow=10,ncol=1000) # ncol = n randomizations

for(i in 1:1000){
  FOri_Local_se[,i]<-c(std.error(calcORI_total(fullMat=fullMat,nsite=12)),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.1,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.2,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.3,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.4,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.5,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.6,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.7,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.8,rule="random")),
                       std.error(calcORI(fullMat=fullMat,nsite=12,prob=0.9,rule="random")))
}

FOri_Local_NULL_mean<-rowMeans(FOri_Local_mean)
FOri_Local_NULL_SE<-rowMeans(FOri_Local_se)

# Plot local species loss scenarios - FOri
erosion<-seq(0,90,10)
plot(erosion,FOri_Local_NULL_mean,type="l",ylim=c(0.05,0.245),
     xlim=c(0,100),col="gray",pch=16,cex=1.5,xlab="Local species loss (%)",ylab="FOri",
     font.lab=2,cex.lab=1.4,cex.axis=1)
points(erosion[-1],FOri_Local_NULL_mean[-1],type="p",col="gray",pch=21,bg="gray",cex=1)
arrows(x0 = erosion, x1=erosion, y0=FOri_Local_NULL_mean-FOri_Local_NULL_SE, y1=FOri_Local_NULL_mean
       +FOri_Local_NULL_SE,code=0,angle=90,length = 0.05, col="gray")
points(erosion,FOri_Local_RARE,type="l",ylim=c(0,1),col="black")
points(erosion[-1],FOri_Local_RARE[-1],type="p",col="black",pch=21,bg="black",cex=1)
points(erosion,FOri_Local_COMM,type="l",ylim=c(0,1),col="black",lty=2)
points(erosion[-1],FOri_Local_COMM[-1],type="p",col="black",pch=21,bg="white",cex=1)
arrows(x0 = erosion, x1=erosion, y0=FOri_Local_RARE-FOri_Local_RARE_SE, y1=FOri_Local_RARE
       +FOri_Local_RARE_SE,code=0,angle=90,length = 0.05)
arrows(x0 = erosion, x1=erosion, y0=FOri_Local_COMM-FOri_Local_COMM_SE, y1=FOri_Local_COMM
       +FOri_Local_COMM_SE,code=0,angle=90,length = 0.05)

legend(0,0.20,legend=c("common","random","rare"),col=c("black","gray","black"),
       cex=1,lty=c(3,1,1),lwd=1.3,xjust=-0.01,yjust=0.5,horiz=F,text.font=0.5,x.intersp=0,bty="n")
