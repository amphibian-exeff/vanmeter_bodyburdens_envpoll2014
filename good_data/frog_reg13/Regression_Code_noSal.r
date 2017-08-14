library(modeest)
library(QuantPsyc)

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
                                           HabFac, bodyweight2, SA_cm2_2, molmass_gmol, Solat20C_mgL, Density_gcm3)) 
mlrfrog$Koc_gmL=as.numeric(as.character(mlrfrog$Koc_gmL)) #convert factors to numerics
mlrfrog$Koc_gmL=log10(mlrfrog$Koc_gmL)
names(mlrfrog)[9]='logkoc'
mlrfrog$AppFactor=log10(mlrfrog$AppFactor)
names(mlrfrog)[3]='logAppFactor'
mlrfrog$SoilConc=log10(mlrfrog$SoilConc)
names(mlrfrog)[5]='logSoilConc'
mlrfrog$TissueConc=log(mlrfrog$TissueConc)
names(mlrfrog)[4]='logTissueConc'
mlrfrog$BCF=log10(mlrfrog$BCF)
names(mlrfrog)[7]='logBCF'
mlrfrog$molmass_gmol=log10(mlrfrog$molmass_gmol)
names(mlrfrog)[13]='logMol'
mlrfrog$Solat20C_mgL=log10(mlrfrog$Solat20C_mgL)
names(mlrfrog)[14]='logSol'
mlrfrog$Density_gcm3=log10(mlrfrog$Density_gcm3)
names(mlrfrog)[15]='logDen'

mlrfrog = mlrfrog[which(mlrfrog$Species != "Mole salamander"),] #remove Mole salamander
row.names(mlrfrog)<-NULL


#######################################
###Regression BCF################
#######################################
mlrfrog_BCF=subset(mlrfrog, select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 

nod=dim(mlrfrog_BCF)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_BCF)
regressors <- names(mlrfrog_BCF)
total_ita=dim(detail_select)[1]
R2_logBCF=array(dim=c(total_ita-1,4))
colnames(R2_logBCF)=c("R2", "Adj_R2", "AIC", "BIC")

###replace true and false in combination by variable names
detail_select_v=detail_select
d_names=names(detail_select)
for (i in d_names){
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
}
detail_select_v[is.na(detail_select_v)] = " "

for (i in 1:(total_ita-1)){
  vect=as.matrix(detail_select)[i,]
  lm_temp=lm(as.formula(paste("mlrfrog$logBCF", paste("mlrfrog_BCF$",regressors[vect], collapse=" + "), sep="~")))
  name_logBCF=paste("sum_logBCF", i, sep="")
  assign(name_logBCF, (lm_temp))
  R2_logBCF[i,1]=summary(lm_temp)$r.squared
  R2_logBCF[i,2]=summary(lm_temp)$adj.r.squared
  R2_logBCF[i,3]=AIC(lm_temp)
  R2_logBCF[i,4]=BIC(lm_temp)
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logBCF=c(apply(R2_logBCF,2,max)[1:2],apply(R2_logBCF,2,min)[3:4])
final_logBCF=c(which(R2_logBCF[,1]==col_best_logBCF[1]),which(R2_logBCF[,2]==col_best_logBCF[2]),
                     which(R2_logBCF[,3]==col_best_logBCF[3]),which(R2_logBCF[,4]==col_best_logBCF[4]))
#mfc_logBCF=mfv(final_logBCF)#find the most common combination ID
mfc_logBCF=final_logBCF[3]#find the most common combination ID
best_rg_logBCF=paste("sum_logBCF", mfc_logBCF, sep="")   #recommended regression
coef_logBCF=summary(eval(parse(text=best_rg_logBCF)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))



#######################################
###Regression Tissue conc################
#######################################
mlrfrog_tc=subset(mlrfrog, select=c(Species, HabFac, logkoc, app_rate_g_cm2, bodyweight2, SA_cm2_2, logMol, logSol, logDen))

nod=dim(mlrfrog_tc)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
detail_select=detail_select[-which(detail_select$Var5==TRUE & detail_select$Var6==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_tc)
regressors <- names(mlrfrog_tc)
total_ita=dim(detail_select)[1]
R2_logTissueConc=array(dim=c(total_ita-1,4))
colnames(R2_logTissueConc)=c("R2", "Adj_R2", "AIC", "BIC")

###replace true and false in combination by variable names
detail_select_v=detail_select
d_names=names(detail_select)
for (i in d_names){
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
}
detail_select_v[is.na(detail_select_v)] = " "

for (i in 1:(total_ita-1)){
  vect=as.matrix(detail_select)[i,]
  lm_temp=lm(as.formula(paste("mlrfrog$logTissueConc", paste("mlrfrog_tc$",regressors[vect], collapse=" + "), sep="~")))
  name_logTissueConc=paste("sum_logTissueConc", i, sep="")
  assign(name_logTissueConc, (lm_temp))
  R2_logTissueConc[i,1]=summary(lm_temp)$r.squared
  R2_logTissueConc[i,2]=summary(lm_temp)$adj.r.squared
  R2_logTissueConc[i,3]=AIC(lm_temp)
  R2_logTissueConc[i,4]=BIC(lm_temp)
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logTissueConc=c(apply(R2_logTissueConc,2,max)[1:2],apply(R2_logTissueConc,2,min)[3:4])
final_logTissueConc=c(which(R2_logTissueConc[,1]==col_best_logTissueConc[1]),which(R2_logTissueConc[,2]==col_best_logTissueConc[2]),
                     which(R2_logTissueConc[,3]==col_best_logTissueConc[3]),which(R2_logTissueConc[,4]==col_best_logTissueConc[4]))
#mfc_logTissueConc=mfv(final_logTissueConc)#find the most common combination ID
mfc_logTissueConc=final_logTissueConc[3]#find the most common combination ID
best_rg_logTissueConc=paste("sum_logTissueConc", mfc_logTissueConc, sep="")   #recommended regression
coef_logTissueConc=summary(eval(parse(text=best_rg_logTissueConc)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))




#######################################
###Regression appfactor################
#######################################
mlrfrog_af=subset(mlrfrog, select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 

nod=dim(mlrfrog_af)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_af)
regressors <- names(mlrfrog_af)
total_ita=dim(detail_select)[1]
R2_logAppFactor=array(dim=c(total_ita-1,4))
colnames(R2_logAppFactor)=c("R2", "Adj_R2", "AIC", "BIC")

###replace true and false in combination by variable names
detail_select_v=detail_select
d_names=names(detail_select)
for (i in d_names){
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
}
detail_select_v[is.na(detail_select_v)] = " "

for (i in 1:(total_ita-1)){
  vect=as.matrix(detail_select)[i,]
  lm_temp=lm(as.formula(paste("mlrfrog$logAppFactor", paste("mlrfrog_af$",regressors[vect], collapse=" + "), sep="~")))
  name_logAppFactor=paste("sum_logAppFactor", i, sep="")
  assign(name_logAppFactor, (lm_temp))
  R2_logAppFactor[i,1]=summary(lm_temp)$r.squared
  R2_logAppFactor[i,2]=summary(lm_temp)$adj.r.squared
  R2_logAppFactor[i,3]=AIC(lm_temp)
  R2_logAppFactor[i,4]=BIC(lm_temp)
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logAppFactor=c(apply(R2_logAppFactor,2,max)[1:2],apply(R2_logAppFactor,2,min)[3:4])
final_logAppFactor=c(which(R2_logAppFactor[,1]==col_best_logAppFactor[1]),which(R2_logAppFactor[,2]==col_best_logAppFactor[2]),
                     which(R2_logAppFactor[,3]==col_best_logAppFactor[3]),which(R2_logAppFactor[,4]==col_best_logAppFactor[4]))
#mfc_logAppFactor=mfv(final_logAppFactor)#find the most common combination ID
mfc_logAppFactor=final_logAppFactor[3]
best_rg_logAppFactor=paste("sum_logAppFactor", mfc_logAppFactor, sep="")   #recommended regression
coef_logAppFactor=summary(eval(parse(text=best_rg_logAppFactor)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))






































#######################################
###Regression soil conc################
#######################################
mlrfrog_sc_test=subset(mlrfrog,select=c(logSoilConc,app_rate_g_cm2,logKow,logkoc,Hlaw,VapPrs_mPa,logMol, logSol, logDen))
cor_matrix_sc_test=cor(subset(mlrfrog_sc_test))

#####Based on the correlation matrix, I can only keep a.logkoc, b.App rate
mlrfrog_sc=subset(mlrfrog_sc_test,select=-c(logSoilConc, logKow,Hlaw,VapPrs_mPa)) 
nod=dim(mlrfrog_sc)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_sc)
regressors <- names(mlrfrog_sc)
total_ita=dim(detail_select)[1]
R2_logSoilConc=array(dim=c(total_ita-1,4))
colnames(R2_logSoilConc)=c("R2", "Adj_R2", "AIC", "BIC")

###replace true and false in combination by variable names
detail_select_v=detail_select
d_names=names(detail_select)
for (i in d_names){
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
  eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
}
detail_select_v[is.na(detail_select_v)] = " "

for (i in 1:(total_ita-1)){
  vect=as.matrix(detail_select)[i,]
  lm_temp=lm(as.formula(paste("mlrfrog$logSoilConc", paste("mlrfrog_sc$",regressors[vect], collapse=" + "), sep="~")))
  name_logSoilConc=paste("sum_logSoilConc", i, sep="")
  assign(name_logSoilConc, (lm_temp))
  R2_logSoilConc[i,1]=summary(lm_temp)$r.squared
  R2_logSoilConc[i,2]=summary(lm_temp)$adj.r.squared
  R2_logSoilConc[i,3]=AIC(lm_temp)
  R2_logSoilConc[i,4]=BIC(lm_temp)
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logSoilConc=c(apply(R2_logSoilConc,2,max)[1:2],apply(R2_logSoilConc,2,min)[3:4])
final_logSoilConc=c(which(R2_logSoilConc[,1]==col_best_logSoilConc[1]),which(R2_logSoilConc[,2]==col_best_logSoilConc[2]),
                     which(R2_logSoilConc[,3]==col_best_logSoilConc[3]),which(R2_logSoilConc[,4]==col_best_logSoilConc[4]))
#mfc_logSoilConc=mfv(final_logSoilConc)#find the most common combination ID
mfc_logSoilConc=final_logSoilConc[3]#find the most common combination ID
best_rg_logSoilConc=paste("sum_logSoilConc", mfc_logSoilConc, sep="")   #recommended regression
coef_logSoilConc=summary(eval(parse(text=best_rg_logSoilConc)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))










