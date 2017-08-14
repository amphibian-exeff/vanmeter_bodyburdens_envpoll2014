library(ggplot2)
library(plyr)

summary(allfrog)
colnames(allfrog)

# indices
#frogs
unique.frogs <- unique(allfrog$Species)  
for(frog in unique.frogs){print(frog)}
index.green <- which(allfrog$Species=="Green treefrog")
print(paste("# green =",length(index.green)))
index.barking <- which(allfrog$Species=="Barking treefrog")
print(paste("# barking =",length(index.barking)))
index.mole <- which(allfrog$Species=="Mole salamander")
print(paste("# mole salamanders =",length(index.mole)))
index.leopard <- which(allfrog$Species=="Leopard frog")
print(paste("# leopard =",length(index.leopard)))
index.fowlers <- which(allfrog$Species=="Fowlers toad")
print(paste("# fowlers =",length(index.fowlers)))
index.gray <- which(allfrog$Species=="Gray treefrog")
print(paste("# gray =",length(index.gray)))
index.cricket <- which(allfrog$Species=="Cricket frog")
print(paste("# cricket =",length(index.cricket)))
index.narrowmouth <- which(allfrog$Species=="Narrowmouth toad")
print(paste("# narrowmouth =",length(index.narrowmouth)))
print(paste("frog records =",length(allfrog$Species)))
count.frogs = length(index.green) + length(index.barking)+ length(index.mole)+ length(index.leopard) +
              length(index.fowlers)+ length(index.gray)+ length(index.cricket)+ length(index.narrowmouth)
print(paste("frog species records =",count.frogs))

#chemicals
unique.chemicals <- unique(allfrog$Chemical)
for(chemical in unique.chemicals){print(chemical)}
index.atrazine <- which(allfrog$Chemical=="Atrazine")
print(paste("# atrazine =",length(index.atrazine)))
index.fipronil <- which(allfrog$Chemical=="Fipronil")
print(paste("# fipronil =",length(index.fipronil)))
index.pendimethalin <- which(allfrog$Chemical=="Pendimethalin")
print(paste("# pendimethalin =",length(index.pendimethalin)))
index.triadimefon <- which(allfrog$Chemical=="Triadimefon")
print(paste("# triadimefon =",length(index.triadimefon)))
index.imidacloprid <- which(allfrog$Chemical=="Imidacloprid")
print(paste("# imidacloprid =",length(index.imidacloprid)))
Nchemicals = length(index.atrazine)+length(index.fipronil)+length(index.pendimethalin)+length(index.triadimefon)+length(index.imidacloprid)
print(paste("# chemicals =",Nchemicals))
#metabolites
index.sulfone <- which(allfrog$Chemical=="Fipronil-Sulfone")
print(paste("# sulfone =",length(index.sulfone)))
index.triadimenol <- which(allfrog$Chemical=="Triadimenol")
print(paste("# triadimenol =",length(index.triadimenol)))
index.deisopropyl <- which(allfrog$Chemical=="Deisopropyl Atrazine")
print(paste("# deisopropyl =",length(index.deisopropyl)))
index.desethyl <- which(allfrog$Chemical=="Desethyl Atrazine")
print(paste("# desethyl =",length(index.desethyl)))
Nmetabolites=length(index.sulfone)+length(index.triadimenol)+length(index.deisopropyl)+length(index.desethyl)
print(paste("# metabolites =",Nmetabolites))
#totals
index.totalatrazine <- which(allfrog$Chemical=="Total Atrazine")
print(paste("# total atrazine =",length(index.totalatrazine)))
index.totaltriadimefon <- which(allfrog$Chemical=="Total Triadimefon")
print(paste("# total triadimefon=",length(index.totaltriadimefon)))
index.totalfipronil <- which(allfrog$Chemical=="Total Fipronil")
print(paste("# total fipronil=",length(index.totalfipronil)))
Ntotals = length(index.totalatrazine)+length(index.totaltriadimefon)+length(index.totalfipronil)
print(paste("# totals =",Ntotals))

Ntotaltotal = Nchemicals + Nmetabolites+Ntotals
print(paste("# total chemical entries =",Ntotaltotal))
print(paste("frog species records =",count.frogs))

#instruments
unique.instruments <- unique(allfrog$Instrument)
for(instrument in unique.instruments){print(instrument)}
index.gcms <- which(allfrog$Instrument=="GCMS")
index.lcms <- which(allfrog$Instrument=="LCMS")
#applications
unique.applications <- unique(allfrog$Application)
for(application in unique.applications){print(application)}
index.soil <- which(allfrog$Application=="Soil")
index.overspray <- which(allfrog$Application=="Overspray")

#construct some factor fields as labels
attach(allfrog)
allfrog$ChemLabel <- paste("Log",allfrog$logKow,allfrog$Chemical,allfrog$Application,allfrog$Instrument)
allfrog$ChemLabel <- as.factor(allfrog$ChemLabel)
unique(paste(Chemical,Application,Instrument))

##############################
#basic histograms and test for normality
allsoil <- allfrog[index.soil,]
dim(allsoil)
#allsoil.lcms <- allsoil[which(allsoil$Instrument=="LCMS"),]
#allsoil.gcms <- allsoil[which(allsoil$Instrument=="GCMS"),]
#View(allsoil)
#View(allsoil.lcms)
#View(allsoil.gcms)
alloverspray <- allfrog[index.overspray,]
dim(alloverspray)
#View(alloverspray)
unique(alloverspray$Species)
index.allsoil.overspray <- which(allfrog$Species==unique(alloverspray$Species))
allsoil.overspray <- allsoil[index.allsoil.overspray,]
dim(allsoil.overspray)
#View(alloverspray)
#alloverspray.lcms <- alloverspray[which(alloverspray$Instrument=="LCMS"),]
#alloverspray.gcms <- alloverspray[which(alloverspray$Instrument=="GCMS"),]
#View(alloverspray)

## lump triademefons and fipronils and atrazines (tba)
## barkers and greens
#ignore frogs as a factor for distribution fitting
##LCMS
pdf(paste(frog_out,"hist_app_overspray.pdf",sep=""),width=11,height=8)
par(mfrow=c(2,2))
for(chemical in unique.chemicals){
    chem.soil <- allsoil.overspray$TissueConc[allsoil.overspray$Chemical==chemical]
    chem.overspray <- alloverspray$TissueConc[alloverspray$Chemical==chemical]
    this.instrument <- unique(allsoil.overspray$Instrument[allsoil.overspray$Chemical==chemical])
    #report out sample size
    print(paste(chemical,this.instrument, "soil samples = ", length(chem.soil)," overspray samples = ", length(chem.overspray)))
    if(length(chem.soil)>0 && length(chem.overspray)>0){
      histmin <- min(c(chem.soil,chem.overspray),na.rm=TRUE)
      histmax <- max(c(chem.soil,chem.overspray),na.rm=TRUE)
      t.p <- round(t.test(chem.soil,chem.overspray)$p.value,digits=5)
      hist(chem.soil,main=paste(this.instrument,chemical,"p=",t.p),xlab="Soil Application: Tissue Concentration",col="blue",xlim=c(histmin,histmax))
      hist(chem.overspray,main=paste(this.instrument,chemical,"p=",t.p),xlab="Overspray Application: Tissue Concentration",col="red",xlim=c(histmin,histmax))
    }
}
dev.off()

#lcms boxplots for barkers and greens that compare soil versus overspray
# for imidacloprid, total atrazine, total triadimefon, total fipronil, and pendimethalin
index.goodchems = c(index.imidacloprid,index.totalatrazine,index.totaltriadimefon,index.totalfipronil,index.pendimethalin)
spray.boxplot <- allfrog[index.goodchems,]
spray.boxplot <- spray.boxplot[spray.boxplot$Instrument=="LCMS",]
View(spray.boxplot[spray.boxplot$Species=="Barking treefrog"|spray.boxplot$Species=="Green treefrog",])
spray.boxplot <- spray.boxplot[spray.boxplot$Species=="Barking treefrog"|spray.boxplot$Species=="Green treefrog",]
dim(spray.boxplot)
View(spray.boxplot)
colnames(spray.boxplot)
spray.boxplot[which(spray.boxplot$Chemical=="Total Atrazine"),3] = "Atrazine"
spray.boxplot[which(spray.boxplot$Chemical=="Total Triadimefon"),3] = "Triadimefon"
spray.boxplot[which(spray.boxplot$Chemical=="Total Fipronil"),3] = "Fipronil"

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

dim(spray.boxplot)
View(spray.boxplot)
spray.boxplot[11:15,11] <- NA
pdf(paste(frog_out,"boxplot_soil_spray_bcf.pdf",sep=""),width=8.5,height=11)
    spray.barkingtreefrog <- na.omit(spray.boxplot[spray.boxplot=="Barking treefrog",])
    spray.factors <- reorder(spray.barkingtreefrog$Chemical, spray.barkingtreefrog$logKow)
    p1 <- qplot(spray.factors, BCF, fill=factor(Application), data=spray.barkingtreefrog,
                  geom="boxplot",xlab="",ylab="Barking treefrog BCF")+annotate("text", x=5, y=3.3, label="A")+
                  annotate("text", x=1, y=-0.25, label="***")+annotate("text", x=2, y=-0.25, label="*")+
                  annotate("text", x=3, y=-0.25, label="**")+annotate("text", x=4, y=-0.25, label="***")+
                  annotate("text", x=5, y=-0.25, label="***")+
                  theme_bw() +scale_fill_grey(start=0.5, end=1) + labs(fill="Application")
    spray.greentreefrog <- na.omit(spray.boxplot[spray.boxplot=="Green treefrog",])
    p2 <- qplot(reorder(Chemical,logKow), BCF, fill=factor(Application), data=spray.greentreefrog, 
                geom="boxplot",xlab="Pesticide",ylab="Green treefrog BCF")+annotate("text", x=4, y=1.2, label="B")+
                annotate("text", x=1, y=-0.25, label="***")+annotate("text", x=2, y=-0.25, label="**")+
                annotate("text", x=3, y=-0.25, label="**")+annotate("text", x=4, y=-0.25, label="***")+
                theme_bw()+scale_fill_grey(start=0.5, end=1) + labs(fill="Application")
    multiplot(p1, p2)
dev.off()

pdf(paste(frog_out,"boxplot_soil_spray_tissueconc.pdf",sep=""),width=8.5,height=11)
spray.barkingtreefrog <- na.omit(spray.boxplot[spray.boxplot=="Barking treefrog",])
spray.factors <- reorder(spray.barkingtreefrog$Chemical, spray.barkingtreefrog$logKow)
p1 <- qplot(spray.factors, TissueConc, fill=factor(Application), data=spray.barkingtreefrog,
            geom="boxplot",xlab="",ylab="Barking treefrog Tissue Concentration (ppm)")+annotate("text", x=5, y=17, label="A")+
  annotate("text", x=1, y=-1.25, label="***")+annotate("text", x=2, y=-1.25, label="*")+
  annotate("text", x=3, y=-1.25, label="**")+annotate("text", x=4, y=-1.25, label="***")+
  annotate("text", x=5, y=-1.25, label="***")+
  theme_bw() +scale_fill_grey(start=0.5, end=1) + labs(fill="Application")
spray.greentreefrog <- na.omit(spray.boxplot[spray.boxplot=="Green treefrog",])
p2 <- qplot(reorder(Chemical,logKow), TissueConc, fill=factor(Application), data=spray.greentreefrog, 
            geom="boxplot",xlab="Pesticide",ylab="Green treefrog Tissue Concentration (ppm)")+annotate("text", x=4, y=21, label="B")+
  annotate("text", x=1, y=-1.25, label="***")+annotate("text", x=2, y=-1.25, label="**")+
  annotate("text", x=3, y=-1.25, label="**")+annotate("text", x=4, y=-1.25, label="***")+
  theme_bw()+scale_fill_grey(start=0.5, end=1) + labs(fill="Application")
multiplot(p1, p2)
dev.off()



pdf(paste(frog_out,"barchart_soil_spray.pdf",sep=""),width=8.5,height=11)
  #create a data frame with averages and standard deviations
  bt <- spray.boxplot[spray.boxplot=="Barking treefrog",]
  bcf.avg<-ddply(bt, c("Chemical", "Application"), function(df)
    return(c(bcf.avg=mean(df$BCF), bcf.sd=sd(df$BCF),bcf.logKow=mean(df$logKow))))
  #create the barplot component
  dodge <- position_dodge(width=0.9)
  avg.plot<-qplot(reorder(Chemical,bcf.logKow), bcf.avg, fill=factor(Application), 
                  data=bcf.avg, xlab="",ylab="Barking treefrog BCF",geom="bar", position="dodge")
  #add error bars
  p1 <- avg.plot+geom_errorbar(aes(ymax=bcf.avg+bcf.sd, ymin=bcf.avg-bcf.sd),position="dodge")+
      annotate("text", x=5, y=3.3, label="A")+theme_bw()+ labs(fill="Application")

  gt <- spray.boxplot[spray.boxplot=="Green treefrog",]
  bcf.avg<-ddply(gt, c("Chemical", "Application"), function(df)
    return(c(bcf.avg=mean(df$BCF), bcf.sd=sd(df$BCF),bcf.logKow=mean(df$logKow))))
  bcf.avg[5,3]=NA
  bcf.avg[5,4]=NA
  #create the barplot component
  dodge <- position_dodge(width=0.9)
  avg.plot<-qplot(reorder(Chemical,bcf.logKow), bcf.avg, fill=factor(Application), 
                  data=bcf.avg, xlab="Pesticide",ylab="Green treefrog BCF",geom="bar", position="dodge")
  #add error bars
  p2 <- avg.plot+geom_errorbar(aes(ymax=bcf.avg+bcf.sd, ymin=bcf.avg-bcf.sd), position="dodge")+
      annotate("text", x=5, y=1.2, label="B")+theme_bw()+ labs(fill="Application")
  multiplot(p1, p2)
dev.off()

# ##GCMS
# pdf(paste(frog_out,"hist_app_gcms.pdf",sep=""),width=11,height=8)
# par(mfrow=c(2,2))
# for(chemical in unique.chemicals){
#     chem.soil <- allsoil.gcms$TissueConc[allsoil$Chemical==chemical]
#     chem.overspray <- alloverspray.gcms$TissueConc[alloverspray$Chemical==chemical]
#     #report out sample size
#     print(paste(chemical, "gcms ","soil samples = ", length(chem.soil)," overspray samples = ", length(chem.overspray)))    
#     if(length(chem.soil)>0 && length(chem.overspray)>0){
#       histmin <- min(c(chem.soil,chem.overspray),na.rm=TRUE)
#       histmax <- max(c(chem.soil,chem.overspray),na.rm=TRUE)
#       t.p <- round(t.test(chem.soil,chem.overspray)$p.value,digits=5)
#       hist(chem.soil,main=paste("GCMS:",chemical,"p=",t.p),xlab="Soil Application: Tissue Concentration",col="blue",xlim=c(histmin,histmax))
#       hist(chem.overspray,main=paste("GCMS:",chemical,"p=",t.p),xlab="Overspray Application: Tissue Concentration",col="red",xlim=c(histmin,histmax))
#     }
# }
# dev.off()

#################
#plot means
AppLabels <- c("Imidacloprod ()","Fipronil (1.43 mg)","Triadimefon (3.57 mg)","Pendimethalin (25.14 mg)","Atrazine (29.34 mg)")
pdf(paste(frog_out,"soil_means_w_CIs_AllChemicals.pdf",sep=""),width=11,height=8)
  par(mfrow=c(1,1))
  # soil concentrations
  plotmeans(log(allfrog$SoilConc)~allfrog$Chemical+allfrog$Application+allfrog$Instrument,xlab="Application Rate (mg/cm^2)", xaxt="n",
            ylab="ln(SoilConc)", main="Soil Concentrations for All Aquaria: Mean Plot with 95% CI",
            barwidth=2,col="dark green")
  axis(side=1,at=c(1,2,3,4),labels=AppLabels[1:4])
  plotmeans(log(allfrog$SoilConc)~allfrog$ChemLabel,xlab="Chemical", ylab="ln(SoilConc)", 
            main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")
dev.off()

plotmeans(log(permeability)~ChemLabel,xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")
plotmeans(permeability~ChemLabel,xlab="Chemical", ylab="Dermal Soil-Frog Transfer Coefficient", main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")

## plotmeans for 5 main chemicals- bcf
pdf(paste(frog_out,"soil_bcf_chemicals_lcms.pdf",sep=""),width=11,height=8)
#  index.somechemicals <- c(index.imidacloprid,index.totalatrazine,index.totaltriadimefon,index.totalfipronil,index.pendimethalin)
# 1,5,2,4,3
index.somechemicals <- c(index.imidacloprid,index.pendimethalin,index.totalatrazine,index.totalfipronil,index.totaltriadimefon)
  index.soil5chem <- intersect(index.somechemicals,index.soil)
  index.soil5chem.lcms <- intersect(index.soil5chem,index.lcms)
  #View(allfrog[index.soil5chemlcms,])
  KowLabels <- c("Imidacloprod \n Log Kow = 0.57","Atrazine \n Log Kow = 2.5","Triadimefon \n Log Kow = 3.11",
                 "Fipronil \n Log Kow = 4","Pendimethalin \n Log Kow = 5.18")
  plotmeans(allfrog[index.soil5chem.lcms,]$BCF~allfrog[index.soil5chem.lcms,]$ChemLabel,xlab="Chemical", xaxt="n",
            ylab="Soil-Frog BCF", main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")
  axis(side=1,at=c(1,2,3,4,5),labels=KowLabels[1:5])
dev.off()

## plotmeans for 5 main chemicals- ln(soil-frog dermal transfer coefficient)
pdf(paste(frog_out,"soil_dermaltransfer_chemicals_lcms.pdf",sep=""),width=11,height=8)
  index.somechemicals <- c(index.imidacloprid,index.totalatrazine,index.totaltriadimefon,index.totalfipronil,index.pendimethalin)
  index.soil5chem <- intersect(index.somechemicals,index.soil)
  index.soil5chem.lcms <- intersect(index.soil5chem,index.lcms)
  #View(allfrog[index.soil5chemlcms,])
  KowLabels <- c("Imidacloprod \n Log Kow = 0.57","Atrazine \n Log Kow = 2.5","Triadimefon \n Log Kow = 3.11",
                 "Fipronil \n Log Kow = 4","Pendimethalin \n Log Kow = 5.18")
  plotmeans(log(allfrog[index.soil5chem.lcms,]$AppFactor)~allfrog[index.soil5chem.lcms,]$ChemLabel,xlab="Chemical",xaxt="n", 
            ylab="ln(Soil-Frog Transfer Coefficient)", main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")
  axis(side=1,at=c(1,2,3,4,5),labels=KowLabels[1:5])
dev.off()

## plotmeans for 8 species- ln(soil-frog dermal transfer coefficient)
pdf(paste(frog_out,"soil_dermaltransfer_species_lcms.pdf",sep=""),width=14,height=8)
index.somechemicals <- c(index.imidacloprid,index.totalatrazine,index.totaltriadimefon,index.totalfipronil,index.pendimethalin)
index.soil5chem <- intersect(index.somechemicals,index.soil)
index.soil5chem.lcms <- intersect(index.soil5chem,index.lcms)
#View(allfrog[index.soil5chemlcms,])
KowLabels <- c("Imidacloprod \n Log Kow = 0.57","Atrazine \n Log Kow = 2.5","Triadimefon \n Log Kow = 3.11",
               "Fipronil \n Log Kow = 4","Pendimethalin \n Log Kow = 5.18")
plotmeans(log(allfrog[index.soil5chem.lcms,]$AppFactor)~allfrog[index.soil5chem.lcms,]$Species,xlab="Species",
          ylab="ln(Soil-Frog Transfer Coefficient)", main="All Species: Mean Plot with 95% CI",barwidth=2,col="dark green")
#axis(side=1,at=c(1,2,3,4,5),labels=KowLabels[1:5])
dev.off()

pdf(paste(frog_root,"means_w_CIs_AllSpecies.pdf",sep=""),width=11,height=8)
par(mfrow=c(1,1))
plotmeans(log(SoilConc)~Species,xlab="Species", ylab="ln(SoilConc)", main="All Chemicals: Mean Plot with 95% CI",barwidth=2,col="dark green")
plotmeans(log(permeability)~Species,xlab="Species", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="All Chemicals: Mean Plot with 95% CI",barwidth=2,col="dark green")
plotmeans(permeability~Species,xlab="Species", ylab="Dermal Soil-Frog Transfer Coefficient", main="All Chemicals: Mean Plot with 95% CI",barwidth=2,col="dark green")
dev.off()

pdf(paste(frog_root,"means_w_CIs_IndividualSpecies.pdf",sep=""),width=8,height=10)
par(mfrow=c(2,1))
plotmeans(log(permeability[index.barking])~ChemLabel[index.barking],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Barking treefrog: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.fowlers])~ChemLabel[index.fowlers],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Fowlers toad: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.gray])~ChemLabel[index.gray],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Gray treefrog: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.green])~ChemLabel[index.green],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Green treefrog: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.leopard])~ChemLabel[index.leopard],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Leopard frog: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.mole])~ChemLabel[index.mole],xlab="Chemical", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Mole salamander: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
dev.off()

pdf(paste(frog_root,"means_w_CIs_IndividualChemicals.pdf",sep=""),width=8,height=10)
par(mfrow=c(2,1))
plotmeans(log(permeability[index.atrazine])~Species[index.atrazine],xlab="Species", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Atrazine: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.fipronil])~Species[index.fipronil],xlab="Species", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Fipronil: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.pendimethalin])~Species[index.pendimethalin],xlab="Species", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Pendimethalin: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
plotmeans(log(permeability[index.triadimefon])~Species[index.triadimefon],xlab="Species", ylab="ln(Dermal Soil-Frog Transfer Coefficient)", main="Triadimefon: Mean Plot with 95% CI",ylim=c(-6.7,0),barwidth=2,col="dark green")
dev.off()

pdf(paste(frog_root,"checking_predictions.pdf",sep=""),width=11,height=8)
par(mfrow=c(1,1))
plot(logKow,permeability,col="dark red",xlim=c(2.3,5.7))
points(logKow+0.1,predict(lm.nospecies))
points(logKow+0.2,predict(ancova.species.1))
points(logKow+0.3,predict(ancova.species.2))
coplot(permeability~Species|molmass_gmol+Solat20C_mgL)
coplot(permeability~Species|Chemical+molmass_gmol)
coplot(permeability~Chemical|Species+molmass_gmol)
dev.off()


pdf(paste(frog_root,"allspecies_boxplots.pdf",sep=""),width=11,height=8)
par(mfrow=c(2,2))
boxplot(permeability~logKow,data=mlrfrog, xlab="Log Kow", ylab="Dermal Soil-Frog Transfer Coefficient",col="firebrick2")
boxplot(permeability~Solat20C_mgL,data=mlrfrog, xlab="Solubility", ylab="Dermal Soil-Frog Transfer Coefficient",col="firebrick2")
boxplot(permeability~molmass_gmol,data=mlrfrog, xlab="Molecular Mass", ylab="Dermal Soil-Frog Transfer Coefficient",col="firebrick2")
boxplot(permeability~Koc_gmL,data=mlrfrog, xlab="Koc", ylab="Dermal Soil-Frog Transfer Coefficient",col="firebrick2")

par(mfrow=c(2,2))
boxplot(Tissue~logKow,data=mlrfrog, xlab="Log Kow", ylab="Tissue Concentration",col="firebrick2")
boxplot(Tissue~Solat20C_mgL,data=mlrfrog, xlab="Solubility", ylab="Tissue Concentration",col="firebrick2")
boxplot(Tissue~molmass_gmol,data=mlrfrog, xlab="Molecular Mass", ylab="Tissue Concentration",col="firebrick2")
boxplot(Tissue~Koc_gmL,data=mlrfrog, xlab="Koc", ylab="Tissue Concentration",col="firebrick2")

par(mfrow=c(2,2))
boxplot(SoilConc~logKow,data=mlrfrog, xlab="Log Kow", ylab="Soil Concentration",col="firebrick2")
boxplot(SoilConc~Solat20C_mgL,data=mlrfrog, xlab="Solubility", ylab="Soil Concentration",col="firebrick2")
boxplot(SoilConc~molmass_gmol,data=mlrfrog, xlab="Molecular Mass", ylab="Soil Concentration",col="firebrick2")
boxplot(SoilConc~Koc_gmL,data=mlrfrog, xlab="Koc", ylab="Soil Concentration",col="firebrick2")

par(mfrow=c(2,2))
boxplot(AppFactor~logKow,data=mlrfrog, xlab="Log Kow", ylab="Uptake Proportion",col="firebrick2")
boxplot(AppFactor~Solat20C_mgL,data=mlrfrog, xlab="Solubility", ylab="Uptake Proportion",col="firebrick2")
boxplot(AppFactor~molmass_gmol,data=mlrfrog, xlab="Molecular Mass", ylab="Uptake Proportion",col="firebrick2")
boxplot(AppFactor~Koc_gmL,data=mlrfrog, xlab="Koc", ylab="Uptake Proportion",col="firebrick2")
dev.off()