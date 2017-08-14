library(mvoutlier)

frog_root<-"D:\\Dropbox\\Robin_data_test\\frog_reg6\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1),] #select good data
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon'))

raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog = subset(raw_data_t2, select=c(Species, app_rate_g_cm2, Chemical, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa, Koc_Lab,
                                       HabFac, Hlaw, bodyweight2, SA_cm2_2, molmass_gmol, Solat20C_mgL, Density_gcm3)) 
mlrfrog$Koc_Lab=as.numeric(as.character(mlrfrog$Koc_Lab)) #convert factors to numerics
mlrfrog$Koc_Lab=log(mlrfrog$Koc_Lab)
names(mlrfrog)[10]='logkoc'
mlrfrog$AppFactor=log(mlrfrog$AppFactor)
names(mlrfrog)[4]='logAppFactor'
mlrfrog$SoilConc=log(mlrfrog$SoilConc)
names(mlrfrog)[6]='logSoilConc'
mlrfrog$TissueConc=log(mlrfrog$TissueConc)
names(mlrfrog)[5]='logTissueConc'
mlrfrog$BCF=log(mlrfrog$BCF)
names(mlrfrog)[8]='logBCF'

data_Imi=mlrfrog[which(mlrfrog$Chemical=="Imidacloprid"),]
test_Imi= subset(data_Imi, select=c(logAppFactor, logBCF))
data_Imi$outliers=aq.plot(test_Imi, alpha=0.01)$outliers
covr <- covMcd(test_Imi, alpha = 0.5)
dist <- mahalanobis(test_Imi, center = covr$center, cov = covr$cov)
data_Imi$MD=dist
# s <- sort(dist, index = TRUE)
# q <- (0.5:length(dist))/length(dist)
# qchi <- qchisq(q, df = 2)
# data_Imi$qchi=qchi

data_Pen=mlrfrog[which(mlrfrog$Chemical=="Pendimethalin"),]
test_Pen= subset(data_Pen, select=c(logAppFactor, logBCF))
data_Pen$outliers=aq.plot(test_Pen, alpha=0.01)$outliers
covr <- covMcd(test_Pen, alpha = 0.5)
dist <- mahalanobis(test_Pen, center = covr$center, cov = covr$cov)
data_Pen$MD=dist

data_Atr=mlrfrog[which(mlrfrog$Chemical=="Total Atrazine"),]
test_Atr= subset(data_Atr, select=c(logAppFactor, logBCF))
data_Atr$outliers=aq.plot(test_Atr, alpha=0.01)$outliers
covr <- covMcd(test_Atr, alpha = 0.5)
dist <- mahalanobis(test_Atr, center = covr$center, cov = covr$cov)
data_Atr$MD=dist


data_Fip=mlrfrog[which(mlrfrog$Chemical=="Total Fipronil"),]
test_Fip= subset(data_Fip, select=c(logAppFactor, logBCF))
data_Fip$outliers=aq.plot(test_Fip, alpha=0.01)$outliers
covr <- covMcd(test_Fip, alpha = 0.5)
dist <- mahalanobis(test_Fip, center = covr$center, cov = covr$cov)
data_Fip$MD=dist

data_Tri=mlrfrog[which(mlrfrog$Chemical=="Total Triadimefon"),]
test_Tri= subset(data_Tri, select=c(logAppFactor, logBCF))
data_Tri$outliers=aq.plot(test_Tri, alpha=0.01)$outliers
covr <- covMcd(test_Tri, alpha = 0.5)
dist <- mahalanobis(test_Tri, center = covr$center, cov = covr$cov)
data_Tri$MD=dist

