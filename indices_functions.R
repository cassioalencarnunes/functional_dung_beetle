## Functions designed by Leitão et al. (2016) to calculate Functional Richness (FRic)
## Functional Specialization (FSpe) and Functional Originality (FOri) with different levels
## of species extinctions. The functions calculate functional indices in scenarios losing
## rarest species first, most common species first and losing species randomly.

# LEITÃO, R.P., J. ZUANON, S. VILLÉGER, S.E. WILLIAMS, C. BARALOTO, C. FORTUNEL, F.P. MENDONÇA,
# and D. MOUILLOT. 2016. Rare species contribute disproportionately to the functional structure
# of species assemblages. Proc. R. Soc. B Biol. Sci. 283: 20160084. 
# Available at: https://doi.org/10.1098/rspb.2016.0084


# 'fullmat' is a matrix with species in rows and sites in columns. In addition, we must include
# the Rarity Index for each species and also the PCoA axis that form the functional space in
# the columns of the matrix

# 'nsite' is the number of sites in the matrix

# 'prob' is the percentage of species that will be extinct (0.1, 0.2, 0.3 ... 0.9)

# In 'rule' we can set the type of extinction scenario we want to simulate
#   "rare", "common" and "random"


#1)#Function to compute Functional Richness (FRic)
###############RUN FUNCTION calcFRic##########################
calcFRic<-function(fullMat,nsite,prob,rule){
  colCom<-1:nsite
  colPCOA<-grep("Axis",names(fullMat))
  siteNames<-names(fullMat)[colCom]
  FRic_sites<-sapply(siteNames,function(s,rule){
    cat(s,"\n")
    site<-fullMat[,s]
    pcoaMat<-fullMat[(c(FALSE,TRUE)[site+1]),colPCOA]
    if(rule=="common"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted,decreasing=TRUE)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      matFinal<-pcoaMat[names(riSiteSorted),]
    }
    if(rule=="rare"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      matFinal<-pcoaMat[names(riSiteSorted),]
    }
    if(rule=="random"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      toRem<-sample(1:length(riSite),nRem)
      matFinal<-pcoaMat[-toRem,]
    }
    hullSub<-convhulln(matFinal,"FA")$vol
    hullCom<-convhulln(pcoaMat,"FA")$vol
    prop<-hullSub/hullCom
  },rule=rule)
  return(FRic_sites)
}
##########END OF FUNCTION calcFRic##########################



#2)#Function to compute Functional Specialization (FSpe)
################RUN FUNCTION calcSPE##########################
calcSPE<-function(fullMat,nsite,prob,rule){
  colCom<-1:nsite
  siteNames<-names(fullMat)[colCom]
  SPE_sites<-sapply(siteNames,function(s,rule){
    cat(s,"\n")
    site<-fullMat[,s]
    SPE<-fullMat[(c(FALSE,TRUE)[site+1]),"speS"]
    names(SPE)<-
      rownames(fullMat[(c(FALSE,TRUE)[site+1]),])
    if(rule=="common"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted,decreasing=TRUE)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      vecFinal<-SPE[names(riSiteSorted)]
    }
    if(rule=="rare"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      vecFinal<-SPE[names(riSiteSorted)]
    }
    if(rule=="random"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      toRem<-sample(1:length(riSite),nRem)
      riSite<-riSite[-toRem]
      vecFinal<-SPE[names(riSite)]
    }
    mean(vecFinal)
  },rule=rule)
  return(SPE_sites)
}
#Computing total FSpe#
calcSPE_total <-function(fullMat,nsite){
  colCom<-1:nsite
  siteNames<-names(fullMat)[colCom]
  SPE_sites<-sapply(siteNames,function(s){
    cat(s,"\n")
    site<-fullMat[,s]
    SPE<-fullMat[(c(FALSE,TRUE)[site+1]),"speS"]
    names(SPE)<-rownames(fullMat[(c(FALSE,TRUE)[site+1]),])
    vecFinal<-SPE[names(SPE)]
    mean(vecFinal)
  })
  return(SPE_sites)
}
#############END OF FUNCTION calcSPE##########################

#3)#Function to compute Functional Originality (FOri)
##################RUN FUNCTION calcORI#######################
calcORI<-function(fullMat,nsite,prob,rule){
  colCom<-1:nsite
  siteNames<-names(fullMat)[colCom]
  ORI_sites<-sapply(siteNames,function(s,rule){
    cat(s,"\n")
    site<-fullMat[,s]
    ORI<-fullMat[(c(FALSE,TRUE)[site+1]),"oriS"]
    names(ORI)<-rownames(fullMat[(c(FALSE,TRUE)[site+1]),])
    if(rule=="common"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted,decreasing=TRUE)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      vecFinal<-ORI[names(riSiteSorted)]
    }
    if(rule=="rare"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      riSiteSorted<-riSite
      riSiteSorted<-sort(riSiteSorted)
      vecRem<-1:nRem
      riSiteSorted<-riSiteSorted[-vecRem]
      vecFinal<-ORI[names(riSiteSorted)]
    }
    if(rule=="random"){
      riSite<-fullMat$RI*site
      names(riSite)<-rownames(fullMat)
      riSite<-riSite[riSite!=0]
      nRem<-round(length(riSite)*prob)
      if(nRem==0) nRem<-1
      toRem<-sample(1:length(riSite),nRem)
      riSite<-riSite[-toRem]
      vecFinal<- ORI[names(riSite)]
    }
    mean(vecFinal)
  },rule=rule)
  return(ORI_sites)
}
#Computing total FOri#
calcORI_total<-function(fullMat,nsite){
  colCom<-1:nsite
  siteNames<-names(fullMat)[colCom]
  ORI_sites<-sapply(siteNames,function(s){
    cat(s,"\n")
    site<-fullMat[,s]
    ORI<-fullMat[(c(FALSE,TRUE)[site+1]),"oriS"]
    names(ORI)<-rownames(fullMat[(c(FALSE,TRUE)[site+1]),])
    vecFinal<-ORI[names(ORI)]
    mean(vecFinal)
  })
  return(ORI_sites)
}
############END OF FUNCTION calcORI##########################