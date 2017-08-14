library(reshape2)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(car)
library(plyr)

####load data############
frog_root<-"D:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\frog_reg13\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1 & raw_data$Application=="Soil"),] #select good data and Soil application
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon'))

raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog = subset(raw_data_t2, select=c(Species, app_rate_g_cm2, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa, Koc_gmL,
                                       HabFac, bodyweight2, SA_cm2_2, molmass_gmol, Solat20C_mgL, Density_gcm3, Chemical)) 
mlrfrog$Koc_gmL=as.numeric(as.character(mlrfrog$Koc_gmL)) #convert factors to numerics
mlrfrog$Koc_gmL=log10(mlrfrog$Koc_gmL)
names(mlrfrog)[9]='logkoc'
mlrfrog$AppFactor=log(mlrfrog$AppFactor)
names(mlrfrog)[3]='logAppFactor'
mlrfrog$SoilConc=log(mlrfrog$SoilConc)
names(mlrfrog)[5]='logSoilConc'
mlrfrog$TissueConc=(mlrfrog$TissueConc)
names(mlrfrog)[4]='logTissueConc'
mlrfrog$BCF=log(mlrfrog$BCF)
names(mlrfrog)[7]='logBCF'
mlrfrog$molmass_gmol=log(mlrfrog$molmass_gmol)
names(mlrfrog)[13]='logMol'
mlrfrog$Solat20C_mgL=log(mlrfrog$Solat20C_mgL)
names(mlrfrog)[14]='logSol'
mlrfrog$Density_gcm3=log(mlrfrog$Density_gcm3)
names(mlrfrog)[15]='logDen'

mlrfrog = mlrfrog[which(mlrfrog$Species != "Mole salamander"),] #remove Mole salamander
row.names(mlrfrog)<-NULL


###################Fig 1###############3
# c('#B8860B', '#FFFF00', '#FF4500', '#FB61D7', '#F08080', '#00FF7F', '#00FFFF')

mlrfrog=subset(mlrfrog, select=c(Species, Chemical, logTissueConc, HabFac))

Fig1_opts=theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line=element_line(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.title.x = element_blank(),
  axis.text.x=element_blank(),
  legend.position="none"
)

Fig1_opts_1=theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line=element_line(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.title.x = element_blank(),
  legend.position="none"
)

Fig1_opts_legend=theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line=element_line(),
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.title.x = element_blank(),
  legend.justification=c(0,1), 
  legend.position=c(0,1)
)

Fig1_opts_l=theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line=element_line(),
  plot.background = element_rect(fill = "transparent",colour = NA)
)

data_plot_Imidacloprid=mlrfrog[which(mlrfrog$Chemical=="Imidacloprid"), ]
p_Imidacloprid=ggplot(data_plot_Imidacloprid, aes(x= HabFac , y = logTissueConc,  fill = Species)) + geom_boxplot()+
  ylab("Tissue Concentrations (ppm)") + ggtitle("Imidacloprid") + Fig1_opts_1 
dd_Imidacloprid=ggplot_build(p_Imidacloprid)
dd_Imidacloprid$data[[1]]$xmin=c(0.60, 1.70, 2.00, 2.30, 3.30)
dd_Imidacloprid$data[[1]]$x=c(0.70, 1.800, 2.10, 2.40, 3.40)
dd_Imidacloprid$data[[1]]$xmax=c(0.80, 1.90, 2.20, 2.50,3.50)
dd_Imidacloprid$data[[1]]$fill=c('#B8860B', '#FF4500', '#FB61D7', '#F08080', '#00FFFF')
# p1=grid.draw(ggplot_gtable(dd_Imidacloprid)) 

data_plot_Pendimethalin=mlrfrog[which(mlrfrog$Chemical=="Pendimethalin"), ]
p_Pendimethalin=ggplot(data_plot_Pendimethalin, aes(x= HabFac , y = logTissueConc,  fill = Species)) + geom_boxplot()+
  ylab("Tissue Concentrations (ppm)") + ggtitle("Pendimethalin") + Fig1_opts# + coord_cartesian(ylim = c(0.0, 1.5)) 
dd_Pendimethalin=ggplot_build(p_Pendimethalin)
dd_Pendimethalin$data[[1]]$xmin=c(1.00, 1.70, 2.00, 2.30, 3.00)
dd_Pendimethalin$data[[1]]$x=c(1.10, 1.800, 2.10, 2.40, 3.10)
dd_Pendimethalin$data[[1]]$xmax=c(1.20, 1.90, 2.20, 2.50, 3.20)
dd_Pendimethalin$data[[1]]$fill=c('#FFFF00', '#FF4500', '#FB61D7', '#F08080', '#00FF7F')

# p2=grid.draw(ggplot_gtable(dd_Pendimethalin)) 


data_plot_TA=mlrfrog[which(mlrfrog$Chemical=="Total Atrazine"), ]
p_TA=ggplot(data_plot_TA, aes(x= HabFac , y = logTissueConc,  fill = Species)) + geom_boxplot()+
  ylab("Tissue Concentrations (ppm)") + ggtitle("Total Atrazine") + Fig1_opts 
dd_TA=ggplot_build(p_TA)
dd_TA$data[[1]]$xmin=c(1.00, 1.70, 2.00, 2.30, 3.00)
dd_TA$data[[1]]$x=c(1.10, 1.800, 2.10, 2.40, 3.10)
dd_TA$data[[1]]$xmax=c(1.20, 1.90, 2.20, 2.50, 3.20)
dd_TA$data[[1]]$fill=c('#FFFF00', '#FF4500', '#FB61D7', '#F08080', '#00FF7F')
# p3=grid.draw(ggplot_gtable(dd_TA)) 

data_plot_TF=mlrfrog[which(mlrfrog$Chemical=="Total Fipronil"), ]
p_TF=ggplot(data_plot_TF, aes(x= HabFac , y = logTissueConc,  fill = Species)) + geom_boxplot()+
  ylab("Tissue Concentrations (ppm)") + ggtitle("Total Fipronil") + Fig1_opts 
dd_TF=ggplot_build(p_TF)
dd_TF$data[[1]]$xmin=c(1.00, 1.70, 2.00, 2.30, 3.00)
dd_TF$data[[1]]$x=c(1.10, 1.800, 2.10, 2.40, 3.10)
dd_TF$data[[1]]$xmax=c(1.20, 1.90, 2.20, 2.50, 3.20)
dd_TF$data[[1]]$fill=c('#FFFF00', '#FF4500', '#FB61D7', '#F08080', '#00FF7F')
# p4=grid.draw(ggplot_gtable(dd_TF)) 

data_plot_TT=mlrfrog[which(mlrfrog$Chemical=="Total Triadimefon"), ]
# 
# limit_plot=c(levels(data_plot_TT$HabFac), "A", "B", "C", "D")
# break_plot = c(0.75, 1.25, 1.75, 2.00, 2.25, 2.75, 3.25)
# species_label=unique(mlrfrog$Species)
p_TT=ggplot(data_plot_TT, aes(x= HabFac , y = logTissueConc,  fill = Species)) + geom_boxplot()+
  ylab("Tissue Concentrations (ppm)") + ggtitle("Total Triadimefon") + Fig1_opts +
  theme(plot.margin=unit(c(0,0,1.0,0), "cm"))

dd_TT=ggplot_build(p_TT)
dd_TT$data[[1]]$xmin=c(0.60, 1.00, 1.70, 2.00, 2.30, 3.00, 3.30)
dd_TT$data[[1]]$x=c(0.70, 1.10, 1.800, 2.10, 2.40, 3.10, 3.40)
dd_TT$data[[1]]$xmax=c(0.80, 1.20, 1.90, 2.20, 2.50, 3.20, 3.50)
dd_TT$data[[1]]$fill[1:7]=c('#B8860B', '#FFFF00', '#FF4500', '#FB61D7', '#F08080', '#00FF7F', '#00FFFF')


# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# legend <- g_legend(p_TT)
# postscript("Fig1.eps", width = 10000, height = 20000, horizontal=TRUE)

tiff(file = "Fig1.tiff",  width = 8000, height = 14000, units = "px", res = 800, compression = "lzw") #
grid.arrange(ggplot_gtable(dd_Imidacloprid), ggplot_gtable(dd_Pendimethalin), ggplot_gtable(dd_TA),
             ggplot_gtable(dd_TF),ggplot_gtable(dd_TT), ncol=1)
grid.text(label="Cricket \n frog", x = 0.135, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Leopard \n frog", x = 0.255, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Barking \n treefrog", x = 0.46, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Gray \n treefrog", x = 0.555, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Green \n treefrog", x = 0.65, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Fowlers \n toad", x = 0.85, y = 0.015, gp=gpar(fontsize=10, col="black"))
grid.text(label="Narrowmouth \n toad", x = 0.94, y = 0.015, gp=gpar(fontsize=10, col="black"))
dev.off()


#########Figure 1.2##############
###Create a new data_frame for plotting####
mlrfrog=subset(mlrfrog, select=c(Species, Chemical, logTissueConc, HabFac))

Fig1_opts=theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line=element_line(),
  plot.background = element_rect(fill = "transparent",colour = NA)
)


data_plot_Arboreal=mlrfrog[which(mlrfrog$HabFac=="Arboreal"), ]
p1=ggplot(data_plot_Arboreal, aes(x=Species , y = logTissueConc,  fill = Chemical)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Arboreal") + ylab("Tissue Concentrations") + Fig1_opts

data_plot_Aquatic=mlrfrog[which(mlrfrog$HabFac=="Aquatic"), ]
data_plot_Aquatic_fake_1=cbind(expand.grid(Chemical=c('Pendimethalin','Total Atrazine','Total Fipronil'), Species="Cricket frog"), logTissueConc=5, HabFac="Aquatic")
data_plot_Aquatic_fake_2=cbind(expand.grid(Chemical=c('Imidacloprid'), Species="Leopard frog"), logTissueConc=5, HabFac="Aquatic")
data_plot_Aquatic=rbind(data_plot_Aquatic, data_plot_Aquatic_fake_1, data_plot_Aquatic_fake_2)
p2=ggplot(data_plot_Aquatic, aes(x=Species , y = logTissueConc,  fill = Chemical)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Aquatic") + ylab("Tissue Concentrations") + Fig1_opts

data_plot_Terrestrial=mlrfrog[which(mlrfrog$HabFac=="Terrestrial"), ]
data_plot_Terrestrial_fake_1=cbind(expand.grid(Chemical=c('Imidacloprid'), Species="Fowlers toad"), logTissueConc=5, HabFac="Aquatic")
data_plot_Terrestrial_fake_2=cbind(expand.grid(Chemical=c('Pendimethalin','Total Atrazine','Total Fipronil'), Species="Narrowmouth toad"), logTissueConc=5, HabFac="Aquatic")
data_plot_Terrestrial=rbind(data_plot_Terrestrial, data_plot_Terrestrial_fake_1, data_plot_Terrestrial_fake_2)
p3=ggplot(data_plot_Terrestrial, aes(x=Species , y = logTissueConc,  fill = Chemical)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Terrestrial") + ylab("Tissue Concentrations") + Fig1_opts


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(p1)

tiff(file = "Fig1.2.tiff",  width = 10000, height = 8000, units = "px", res = 800, compression = "lzw") #
grid.arrange(p1+ theme(legend.position="none"), p2+ theme(legend.position="none"), p3+ theme(legend.position="none"), legend, nrow=2)
dev.off()


#########Figure 1.3##############
###g <- ggplot_build(p7)##
###unique(g$data[[1]]["fill"])#######

p4=ggplot(data_plot_Arboreal, aes(x=Chemical, y = logTissueConc,  fill = Species)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Arboreal") + ylab("Tissue Concentrations") + Fig1_opts +  scale_fill_manual(values = c("#F8766D", "#C49A00", "#00C094"))

p5=ggplot(data_plot_Aquatic, aes(x=Chemical, y = logTissueConc,  fill = Species)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Aquatic") + ylab("Tissue Concentrations") + Fig1_opts +  scale_fill_manual(values = c("#00B6EB", "#FB61D7"))

p6=ggplot(data_plot_Terrestrial, aes(x=Chemical, y = logTissueConc,  fill = Species)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-4.5, 3.0)) + xlab("Terrestrial") + ylab("Tissue Concentrations") + Fig1_opts +  scale_fill_manual(values = c("#53B400", "#A58AFF"))

p7=ggplot(mlrfrog, aes(x=Chemical, y = logTissueConc,  fill = Species)) + geom_boxplot()

legend2 <- g_legend(p7)

tiff(file = "Fig1.3.tiff",  width = 10000, height = 8000, units = "px", res = 800, compression = "lzw") #
grid.arrange(p4+ theme(legend.position="none"), p5+ theme(legend.position="none"), p6+ theme(legend.position="none"), legend2, nrow=2)
dev.off()

