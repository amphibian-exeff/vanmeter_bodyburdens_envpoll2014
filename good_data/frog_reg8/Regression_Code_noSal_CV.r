library(modeest)
library(QuantPsyc)

####load data############
frog_root<-"I:\\Dropbox\\Robin_data_test\\frog_reg8\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1 & raw_data$Application=="Soil"),] #select good data and Soil application
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon'))

raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog_raw = subset(raw_data_t2, select=c(Chemical, Species, app_rate_g_cm2, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa, Koc_Lab,
                                           HabFac, Hlaw, bodyweight2, SA_cm2_2, molmass_gmol, Solat20C_mgL, Density_gcm3)) 
mlrfrog_raw$Koc_Lab=as.numeric(as.character(mlrfrog_raw$Koc_Lab)) #convert factors to numerics
mlrfrog_raw$Koc_Lab=log(mlrfrog_raw$Koc_Lab)
names(mlrfrog_raw)[10]='logkoc'
mlrfrog_raw$AppFactor=log(mlrfrog_raw$AppFactor)
names(mlrfrog_raw)[4]='logAppFactor'
mlrfrog_raw$SoilConc=log(mlrfrog_raw$SoilConc)
names(mlrfrog_raw)[6]='logSoilConc'
mlrfrog_raw$TissueConc=log(mlrfrog_raw$TissueConc)
names(mlrfrog_raw)[5]='logTissueConc'
mlrfrog_raw$BCF=log(mlrfrog_raw$BCF)
names(mlrfrog_raw)[8]='logBCF'
mlrfrog_raw$molmass_gmol=log(mlrfrog_raw$molmass_gmol)
names(mlrfrog_raw)[15]='logMol'
mlrfrog_raw$Solat20C_mgL=log(mlrfrog_raw$Solat20C_mgL)
names(mlrfrog_raw)[16]='logSol'
mlrfrog_raw$Density_gcm3=log(mlrfrog_raw$Density_gcm3)
names(mlrfrog_raw)[17]='logDen'
mlrfrog_raw = mlrfrog_raw[which(mlrfrog_raw$Species != "Mole salamander"),] #remove Mole salamander
row.names(mlrfrog_raw)<-NULL

#########################################
######Subset the data#########################
#########################################
C_name_pool=c('Imidacloprid','Pendimethalin','Total Atrazine','Total Fipronil','Total Triadimefon')
C_name=C_name_pool[5]
mlrfrog_CV=mlrfrog_raw[which(mlrfrog_raw$Chemical == C_name),]
mlrfrog=mlrfrog_raw[which(mlrfrog_raw$Chemical != C_name),]

#######################################
###Regression BCF################
#######################################
#####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, 

mlrfrog_BCF_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
mlrfrog_BCF_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 

nod=dim(mlrfrog_BCF_ind)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_BCF_ind)
regressors <- names(mlrfrog_BCF_ind)
total_ita=dim(detail_select)[1]
R2_logBCF=array(dim=c(total_ita-1,9))
colnames(R2_logBCF)=c("R2", "Adj_R2", "AIC", "BIC", "MSE_t", "R2_t", "Adj_R2_t", "AIC_t", "BIC_t")

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
  lm_temp=lm(as.formula(paste("mlrfrog$logBCF", paste("mlrfrog_BCF_ind$",regressors[vect], collapse=" + "), sep="~")))
  name_logBCF=paste("sum_logBCF", i, sep="")
  assign(name_logBCF, (lm_temp))
  
  lm_temp_CV=lm(as.formula(paste("mlrfrog_CV$logBCF", paste("mlrfrog_BCF_ind_CV$",regressors[vect], collapse=" + "), sep="~")))
  name_logBCF_CV=paste("sum_logBCF_CV", i, sep="")
  assign(name_logBCF_CV, (lm_temp_CV))  
  
  R2_logBCF[i,1]=summary(lm_temp)$r.squared
  R2_logBCF[i,2]=summary(lm_temp)$adj.r.squared
  R2_logBCF[i,3]=AIC(lm_temp)
  R2_logBCF[i,4]=BIC(lm_temp)
  R2_logBCF[i,5]=sum(lm_temp_CV$residuals^2)/lm_temp_CV$df.residual
  R2_logBCF[i,6]=summary(lm_temp_CV)$r.squared
  R2_logBCF[i,7]=summary(lm_temp_CV)$adj.r.squared
  R2_logBCF[i,8]=AIC(lm_temp_CV)
  R2_logBCF[i,9]=BIC(lm_temp_CV)
  
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logBCF=c(apply(R2_logBCF,2,max)[1:2],apply(R2_logBCF,2,min)[3:5])
final_logBCF=c(which(R2_logBCF[,1]==col_best_logBCF[1]),which(R2_logBCF[,2]==col_best_logBCF[2]),
                     which(R2_logBCF[,3]==col_best_logBCF[3]),which(R2_logBCF[,4]==col_best_logBCF[4]))
final_logBCF_CV=which(R2_logBCF[,5]==col_best_logBCF[5])
final_select_CV=rbind(cbind(detail_select_v[final_logBCF,],R2_logBCF[final_logBCF,]), cbind(detail_select_v[final_logBCF_CV,],R2_logBCF[final_logBCF_CV,]))

final_logBCF_CV=which(R2_logBCF[,5]==col_best_logBCF[5])
#mfc_logBCF=mfv(final_logBCF)#find the most common combination ID
mfc_logBCF=final_logBCF[3]#find the most common combination ID
best_rg_logBCF=paste("sum_logBCF", mfc_logBCF, sep="")   #recommended regression
coef_logBCF=summary(eval(parse(text=best_rg_logBCF)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))



#######################################
###Regression Tissue conc################
#######################################
#####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, d.app_rate_g_cm2, e.bodyweight2, f.SA_cm2_2, g. molmass_gmol
mlrfrog_tc_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, app_rate_g_cm2, bodyweight2, SA_cm2_2, logMol, logSol, logDen))
mlrfrog_tc_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, app_rate_g_cm2, bodyweight2, SA_cm2_2, logMol, logSol, logDen))

nod=dim(mlrfrog_tc_ind)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
detail_select=detail_select[-which(detail_select$Var5==TRUE & detail_select$Var6==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_tc_ind)
regressors <- names(mlrfrog_tc_ind)
total_ita=dim(detail_select)[1]
R2_logTissueConc=array(dim=c(total_ita-1,9))
colnames(R2_logTissueConc)=c("R2", "Adj_R2", "AIC", "BIC", "MSE_t", "R2_t", "Adj_R2_t", "AIC_t", "BIC_t")

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
  lm_temp=lm(as.formula(paste("mlrfrog$logTissueConc", paste("mlrfrog_tc_ind$",regressors[vect], collapse=" + "), sep="~")))
  name_logTissueConc=paste("sum_logTissueConc", i, sep="")
  assign(name_logTissueConc, (lm_temp))
  
  
  lm_temp_CV=lm(as.formula(paste("mlrfrog_CV$logTissueConc", paste("mlrfrog_tc_ind_CV$",regressors[vect], collapse=" + "), sep="~")))
  name_logTissueConc_CV=paste("sum_logTissueConc_CV", i, sep="")
  assign(name_logTissueConc_CV, (lm_temp_CV))    
  
  R2_logTissueConc[i,1]=summary(lm_temp)$r.squared
  R2_logTissueConc[i,2]=summary(lm_temp)$adj.r.squared
  R2_logTissueConc[i,3]=AIC(lm_temp)
  R2_logTissueConc[i,4]=BIC(lm_temp)
  R2_logTissueConc[i,5]=sum(lm_temp_CV$residuals^2)/lm_temp_CV$df.residual
  R2_logTissueConc[i,6]=summary(lm_temp_CV)$r.squared
  R2_logTissueConc[i,7]=summary(lm_temp_CV)$adj.r.squared
  R2_logTissueConc[i,8]=AIC(lm_temp_CV)
  R2_logTissueConc[i,9]=BIC(lm_temp_CV)  
  
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logTissueConc=c(apply(R2_logTissueConc,2,max)[1:2],apply(R2_logTissueConc,2,min)[3:5])
final_logTissueConc=c(which(R2_logTissueConc[,1]==col_best_logTissueConc[1]),which(R2_logTissueConc[,2]==col_best_logTissueConc[2]),
                     which(R2_logTissueConc[,3]==col_best_logTissueConc[3]),which(R2_logTissueConc[,4]==col_best_logTissueConc[4]))
final_logTissueConc_CV=which(R2_logTissueConc[,5]==col_best_logTissueConc[5])
final_select_CV=rbind(cbind(detail_select_v[final_logTissueConc,],R2_logTissueConc[final_logTissueConc,]), cbind(detail_select_v[final_logTissueConc_CV,],R2_logTissueConc[final_logTissueConc_CV,]))


#mfc_logTissueConc=mfv(final_logTissueConc)#find the most common combination ID
mfc_logTissueConc=final_logTissueConc[3]#find the most common combination ID
best_rg_logTissueConc=paste("sum_logTissueConc", mfc_logTissueConc, sep="")   #recommended regression
coef_logTissueConc=summary(eval(parse(text=best_rg_logTissueConc)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))



#######################################
###Regression appfactor################
#######################################
#####Based on the correlation matrix, I can only keep a.logkoc, b.HabFac, c.Species
mlrfrog_af_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
mlrfrog_af_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, logMol, logSol, logDen))


nod=dim(mlrfrog_af_ind)[2]
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_af_ind)
regressors <- names(mlrfrog_af_ind)
total_ita=dim(detail_select)[1]
R2_logAppFactor=array(dim=c(total_ita-1,9))
colnames(R2_logAppFactor)=c("R2", "Adj_R2", "AIC", "BIC", "MSE_t", "R2_t", "Adj_R2_t", "AIC_t", "BIC_t")

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
lm_temp=lm(as.formula(paste("mlrfrog$logAppFactor", paste("mlrfrog_af_ind$",regressors[vect], collapse=" + "), sep="~")))
name_logAppFactor=paste("sum_logAppFactor", i, sep="")
assign(name_logAppFactor, (lm_temp))

lm_temp_CV=lm(as.formula(paste("mlrfrog_CV$logAppFactor", paste("mlrfrog_af_ind_CV$",regressors[vect], collapse=" + "), sep="~")))
name_logAppFactor_CV=paste("sum_logAppFactor_CV", i, sep="")
assign(name_logAppFactor_CV, (lm_temp_CV))    

R2_logAppFactor[i,1]=summary(lm_temp)$r.squared
R2_logAppFactor[i,2]=summary(lm_temp)$adj.r.squared
R2_logAppFactor[i,3]=AIC(lm_temp)
R2_logAppFactor[i,4]=BIC(lm_temp)
R2_logAppFactor[i,5]=sum(lm_temp_CV$residuals^2)/lm_temp_CV$df.residual
R2_logAppFactor[i,6]=summary(lm_temp_CV)$r.squared
R2_logAppFactor[i,7]=summary(lm_temp_CV)$adj.r.squared
R2_logAppFactor[i,8]=AIC(lm_temp_CV)
R2_logAppFactor[i,9]=BIC(lm_temp_CV)  
}


####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
col_best_logAppFactor=c(apply(R2_logAppFactor,2,max)[1:2],apply(R2_logAppFactor,2,min)[3:5])
final_logAppFactor=c(which(R2_logAppFactor[,1]==col_best_logAppFactor[1]),which(R2_logAppFactor[,2]==col_best_logAppFactor[2]),
              which(R2_logAppFactor[,3]==col_best_logAppFactor[3]),which(R2_logAppFactor[,4]==col_best_logAppFactor[4]))
final_logAppFactor_CV=which(R2_logAppFactor[,5]==col_best_logAppFactor[5])
final_select_CV=rbind(cbind(detail_select_v[final_logAppFactor,],R2_logAppFactor[final_logAppFactor,]), cbind(detail_select_v[final_logAppFactor_CV,],R2_logAppFactor[final_logAppFactor_CV,]))

#mfc_logAppFactor=mfv(final_logAppFactor)#find the most common combination ID
mfc_logAppFactor=final_logAppFactor[3]
best_rg_logAppFactor=paste("sum_logAppFactor", mfc_logAppFactor, sep="")   #recommended regression
coef_logAppFactor=summary(eval(parse(text=best_rg_logAppFactor)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))





# #######################################
# ###Regression soil conc################
# #######################################
# mlrfrog_sc_test=subset(mlrfrog,select=c(logSoilConc,app_rate_g_cm2,logKow,logkoc,Hlaw,VapPrs_mPa,logMol, logSol, logDen))
# cor_matrix_sc_test=cor(subset(mlrfrog_sc_test))
# 
# #####Based on the correlation matrix, I can only keep a.logkoc, b.App rate
# mlrfrog_sc=subset(mlrfrog_sc_test,select=-c(logSoilConc, logKow,Hlaw,VapPrs_mPa)) 
# nod=dim(mlrfrog_sc)[2]
# detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
# detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
# row.names(detail_select)=NULL
# names(detail_select)=names(mlrfrog_sc)
# regressors <- names(mlrfrog_sc)
# total_ita=dim(detail_select)[1]
# R2_logSoilConc=array(dim=c(total_ita-1,4))
# colnames(R2_logSoilConc)=c("R2", "Adj_R2", "AIC", "BIC")
# 
# ###replace true and false in combination by variable names
# detail_select_v=detail_select
# d_names=names(detail_select)
# for (i in d_names){
#   eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
#   eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
# }
# detail_select_v[is.na(detail_select_v)] = " "
# 
# for (i in 1:(total_ita-1)){
#   vect=as.matrix(detail_select)[i,]
#   lm_temp=lm(as.formula(paste("mlrfrog$logSoilConc", paste("mlrfrog_sc$",regressors[vect], collapse=" + "), sep="~")))
#   name_logSoilConc=paste("sum_logSoilConc", i, sep="")
#   assign(name_logSoilConc, (lm_temp))
#   R2_logSoilConc[i,1]=summary(lm_temp)$r.squared
#   R2_logSoilConc[i,2]=summary(lm_temp)$adj.r.squared
#   R2_logSoilConc[i,3]=AIC(lm_temp)
#   R2_logSoilConc[i,4]=BIC(lm_temp)
# }
# 
# 
# ####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
# col_best_logSoilConc=c(apply(R2_logSoilConc,2,max)[1:2],apply(R2_logSoilConc,2,min)[3:4])
# final_logSoilConc=c(which(R2_logSoilConc[,1]==col_best_logSoilConc[1]),which(R2_logSoilConc[,2]==col_best_logSoilConc[2]),
#                     which(R2_logSoilConc[,3]==col_best_logSoilConc[3]),which(R2_logSoilConc[,4]==col_best_logSoilConc[4]))
# #mfc_logSoilConc=mfv(final_logSoilConc)#find the most common combination ID
# mfc_logSoilConc=final_logSoilConc[3]#find the most common combination ID
# best_rg_logSoilConc=paste("sum_logSoilConc", mfc_logSoilConc, sep="")   #recommended regression
# coef_logSoilConc=summary(eval(parse(text=best_rg_logSoilConc)))$coefficients 
# #coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))
# 
# 





