library(reshape2)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(car)

####load data############
#frog_root<-"D:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\frog_reg13\\" #PC path
frog_root<-"C:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\frog_reg13\\" #PC path
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
mlrfrog$TissueConc=log(mlrfrog$TissueConc)
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


#########Figure 2.1#################
###Regression BCF################
#######################################

data_plot=mlrfrog
lm_BCF=lm(mlrfrog$logBCF~mlrfrog$HabFac+mlrfrog$logkoc+mlrfrog$logSol)
lm_APP=lm(mlrfrog$logAppFactor~mlrfrog$HabFac+mlrfrog$logkoc+mlrfrog$logSol)

jpeg(file = "Fig2_color.jpg",  width = 10000, height = 8000, units = "px", res = 800) #

par(mfrow=c(3,2))

########################Arboreal#############################
data_Arboreal_index=which(mlrfrog$HabFac=="Arboreal")
data_Arboreal=mlrfrog[data_Arboreal_index, ]
data_Arboreal$fitted.values=lm_BCF$fitted.values[data_Arboreal_index]
data_Arboreal$residuals=lm_BCF$residuals[data_Arboreal_index]
pos_Arboreal_index=which(data_Arboreal$residuals>0)

lm_BCF_Arboreal=lm(data_Arboreal$logBCF~data_Arboreal$logkoc+data_Arboreal$logSol)
lm_BCF_Arboreal$coefficients[1]=2.4848
lm_BCF_Arboreal$coefficients[2]=-1.3730
lm_BCF_Arboreal$coefficients[3]=-0.2837 

s_Arboreal=scatterplot3d(data_Arboreal$logkoc, data_Arboreal$logSol, data_Arboreal$logBCF, color="blue",
                         main='Arboreal', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(BCF)', xlim=c(2.0,6.0))
s_Arboreal$plane3d(lm_BCF_Arboreal)
observ2d=s_Arboreal$xyz.convert(data_Arboreal$logkoc, data_Arboreal$logSol, data_Arboreal$logBCF)
pred2d=s_Arboreal$xyz.convert(data_Arboreal$logkoc, data_Arboreal$logSol, fitted(lm_BCF_Arboreal))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Arboreal$points3d(data_Arboreal[pos_Arboreal_index,]$logkoc, data_Arboreal[pos_Arboreal_index,]$logSol, data_Arboreal[pos_Arboreal_index,]$logBCF, col="red")


########################Arboreal APP#############################
data_Arboreal$residuals_APP=lm_APP$residuals[data_Arboreal_index]
pos_Arboreal_index=which(data_Arboreal$residuals_APP>0)

lm_APP_Arboreal=lm(data_Arboreal$logAppFactor~data_Arboreal$logkoc+data_Arboreal$logSol)
lm_APP_Arboreal$coefficients[1]=17.8905
lm_APP_Arboreal$coefficients[2]=-1.7280
lm_APP_Arboreal$coefficients[3]=-0.4383

s_Arboreal_APP=scatterplot3d(data_Arboreal$logkoc, data_Arboreal$logSol, data_Arboreal$logAppFactor, color="blue", zlim=c(4,16),
                         main='Arboreal', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(SPF)', xlim=c(2.0,6.0))
s_Arboreal_APP$plane3d(lm_APP_Arboreal)
observ2d=s_Arboreal_APP$xyz.convert(data_Arboreal$logkoc, data_Arboreal$logSol, data_Arboreal$logAppFactor)
pred2d=s_Arboreal_APP$xyz.convert(data_Arboreal$logkoc, data_Arboreal$logSol, fitted(lm_APP_Arboreal))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Arboreal_APP$points3d(data_Arboreal[pos_Arboreal_index,]$logkoc, data_Arboreal[pos_Arboreal_index,]$logSol, data_Arboreal[pos_Arboreal_index,]$logAppFactor, col="red")



########################Aquatic#############################
data_Aquatic_index=which(mlrfrog$HabFac=="Aquatic")
data_Aquatic=mlrfrog[data_Aquatic_index, ]
data_Aquatic$fitted.values=lm_BCF$fitted.values[data_Aquatic_index]
data_Aquatic$residuals=lm_BCF$residuals[data_Aquatic_index]
pos_Aquatic_index=which(data_Aquatic$residuals>0)

lm_BCF_Aquatic=lm(data_Aquatic$logBCF~data_Aquatic$logkoc+data_Aquatic$logSol)
lm_BCF_Aquatic$coefficients[1]=2.6474
lm_BCF_Aquatic$coefficients[2]=-1.3730
lm_BCF_Aquatic$coefficients[3]=-0.2837 

s_Aquatic=scatterplot3d(data_Aquatic$logkoc, data_Aquatic$logSol, data_Aquatic$logBCF, color="blue", zlim=c(-7.0,1),
                        main='Aquatic', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(BCF)', xlim=c(2.0,6.0))
s_Aquatic$plane3d(lm_BCF_Aquatic)
observ2d=s_Aquatic$xyz.convert(data_Aquatic$logkoc, data_Aquatic$logSol, data_Aquatic$logBCF)
pred2d=s_Aquatic$xyz.convert(data_Aquatic$logkoc, data_Aquatic$logSol, fitted(lm_BCF_Aquatic))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Aquatic$points3d(data_Aquatic[pos_Aquatic_index,]$logkoc, data_Aquatic[pos_Aquatic_index,]$logSol, data_Aquatic[pos_Aquatic_index,]$logBCF, col="red")

########################Aquatic APP#############################
data_Aquatic$residuals_APP=lm_APP$residuals[data_Aquatic_index]
pos_Aquatic_index=which(data_Aquatic$residuals_APP>0)

lm_APP_Aquatic=lm(data_Aquatic$logAppFactor~data_Aquatic$logkoc+data_Aquatic$logSol)
lm_APP_Aquatic$coefficients[1]=18.2152
lm_APP_Aquatic$coefficients[2]=-1.7280 
lm_APP_Aquatic$coefficients[3]=-0.4383

s_Aquatic_APP=scatterplot3d(data_Aquatic$logkoc, data_Aquatic$logSol, data_Aquatic$logAppFactor, color="blue", zlim=c(6.0,16),
                        main='Aquatic', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(SPF)', xlim=c(2.0,6.0))
s_Aquatic_APP$plane3d(lm_APP_Aquatic)
observ2d=s_Aquatic_APP$xyz.convert(data_Aquatic$logkoc, data_Aquatic$logSol, data_Aquatic$logAppFactor)
pred2d=s_Aquatic_APP$xyz.convert(data_Aquatic$logkoc, data_Aquatic$logSol, fitted(lm_APP_Aquatic))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Aquatic_APP$points3d(data_Aquatic[pos_Aquatic_index,]$logkoc, data_Aquatic[pos_Aquatic_index,]$logSol, data_Aquatic[pos_Aquatic_index,]$logAppFactor, col="red")



########################Terrestrial#############################
data_Terrestrial_index=which(mlrfrog$HabFac=="Terrestrial")
data_Terrestrial=mlrfrog[data_Terrestrial_index, ]
data_Terrestrial$fitted.values=lm_BCF$fitted.values[data_Terrestrial_index]
data_Terrestrial$residuals=lm_BCF$residuals[data_Terrestrial_index]
pos_Terrestrial_index=which(data_Terrestrial$residuals>0)

lm_BCF_Terrestrial=lm(data_Terrestrial$logBCF~data_Terrestrial$logkoc+data_Terrestrial$logSol)
lm_BCF_Terrestrial$coefficients[1]=3.3128
lm_BCF_Terrestrial$coefficients[2]=-1.3730
lm_BCF_Terrestrial$coefficients[3]=-0.2837 

s_Terrestrial=scatterplot3d(data_Terrestrial$logkoc, data_Terrestrial$logSol, data_Terrestrial$logBCF, color="blue", zlim=c(-7.0,1),
                            main='Terrestrial', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(BCF)', xlim=c(2.0,6.0))
s_Terrestrial$plane3d(lm_BCF_Terrestrial)
observ2d=s_Terrestrial$xyz.convert(data_Terrestrial$logkoc, data_Terrestrial$logSol, data_Terrestrial$logBCF)
pred2d=s_Terrestrial$xyz.convert(data_Terrestrial$logkoc, data_Terrestrial$logSol, fitted(lm_BCF_Terrestrial))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Terrestrial$points3d(data_Terrestrial[pos_Terrestrial_index,]$logkoc, data_Terrestrial[pos_Terrestrial_index,]$logSol, data_Terrestrial[pos_Terrestrial_index,]$logBCF, col="red")


########################Terrestrial#############################
data_Terrestrial$residuals_APP=lm_APP$residuals[data_Terrestrial_index]
pos_Terrestrial_index=which(data_Terrestrial$residuals_APP>0)

lm_APP_Terrestrial=lm(data_Terrestrial$logAppFactor~data_Terrestrial$logkoc+data_Terrestrial$logSol)
lm_APP_Terrestrial$coefficients[1]=18.5488
lm_APP_Terrestrial$coefficients[2]=-1.7280 
lm_APP_Terrestrial$coefficients[3]=-0.4383

s_Terrestrial_APP=scatterplot3d(data_Terrestrial$logkoc, data_Terrestrial$logSol, data_Terrestrial$logAppFactor, color="blue", zlim=c(7.0,16),
                            main='Terrestrial', xlab='log10(Koc)', ylab='log10(Sol)', zlab='log10(SPF)', xlim=c(2.0,6.0))
s_Terrestrial_APP$plane3d(lm_APP_Terrestrial)
observ2d=s_Terrestrial_APP$xyz.convert(data_Terrestrial$logkoc, data_Terrestrial$logSol, data_Terrestrial$logAppFactor)
pred2d=s_Terrestrial_APP$xyz.convert(data_Terrestrial$logkoc, data_Terrestrial$logSol, fitted(lm_APP_Terrestrial))
segments(observ2d$x, observ2d$y, pred2d$x, pred2d$y, lty=4)
s_Terrestrial_APP$points3d(data_Terrestrial[pos_Terrestrial_index,]$logkoc, data_Terrestrial[pos_Terrestrial_index,]$logSol, data_Terrestrial[pos_Terrestrial_index,]$logAppFactor, col="red")

dev.off()


par(mfrow=c(2,2))

qqPlot(lm_BCF_Arboreal, main="QQ Plot Arboreal", ylab="Studentized Residuals")
qqPlot(lm_BCF_Aquatic, main="QQ Plot Aquatic", ylab="Studentized Residuals")
qqPlot(lm_BCF_Terrestrial, main="QQ Plot Terrestrial", ylab="Studentized Residuals")





