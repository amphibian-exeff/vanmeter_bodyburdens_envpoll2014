library(modeest)
library(QuantPsyc)

####load data############
frog_root<-"I:\\Dropbox\\Robin_data_test\\frog_reg3\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1),] #select good data
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                        'Total Fipronil','Total Triadimefon'))

raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog = subset(raw_data_t2, select=c(Species, app_rate_g_cm2, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa,
                                           HabFac, Hlaw, bodyweight2, SA_cm2_2)) 

mlrfrog$AppFactor=log(mlrfrog$AppFactor)
names(mlrfrog)[3]='logAppFactor'
mlrfrog$SoilConc=log(mlrfrog$SoilConc)
names(mlrfrog)[5]='logSoilConc'
mlrfrog$TissueConc=log(mlrfrog$TissueConc)
names(mlrfrog)[4]='logTissueConc'
mlrfrog$BCF=log(mlrfrog$BCF)
names(mlrfrog)[7]='logBCF'

#######################################
###Regression BCF################
#######################################
mlrfrog_BCF_test=subset(mlrfrog,select=c(logBCF,Species,HabFac,logKow,Hlaw,VapPrs_mPa))
cor_matrix_BCF_test=cor(subset(mlrfrog_BCF_test,select=-c(HabFac,Species)))

#####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkow 
mlrfrog_BCF=subset(mlrfrog_BCF_test,select=-c(logBCF,Hlaw,VapPrs_mPa)) 

detail_select=expand.grid(rep(list(c(TRUE,FALSE)),3))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_BCF)
regressors <- names(mlrfrog_BCF)
total_ita=dim(detail_select)[1]
R2_logBCF=array(dim=c(total_ita-1,4))
colnames(R2_logBCF)=c("R2", "Adj_R2", "AIC", "BIC")

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
mfc_logBCF=final_logBCF[4]#find the most common combination ID
best_rg_logBCF=paste("sum_logBCF", mfc_logBCF, sep="")   #recommended regression
coef_logBCF=summary(eval(parse(text=best_rg_logBCF)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))




#######################################
###Regression Tissue conc################
#######################################
mlrfrog_tc_test=subset(mlrfrog,select=c(logTissueConc,Species,HabFac,logKow,Hlaw,VapPrs_mPa,app_rate_g_cm2,bodyweight2,SA_cm2_2))
cor_matrix_tc_test=cor(subset(mlrfrog_tc_test,select=-c(HabFac,Species)))

#####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkow, d.app_rate_g_cm2, e.bodyweight2, f.SA_cm2_2
mlrfrog_tc=subset(mlrfrog_tc_test,select=-c(logTissueConc,Hlaw,VapPrs_mPa)) 

detail_select=expand.grid(rep(list(c(TRUE,FALSE)),6))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_tc)
regressors <- names(mlrfrog_tc)
total_ita=dim(detail_select)[1]
R2_logTissueConc=array(dim=c(total_ita-1,4))
colnames(R2_logTissueConc)=c("R2", "Adj_R2", "AIC", "BIC")

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
mfc_logTissueConc=final_logTissueConc[4]#find the most common combination ID
best_rg_logTissueConc=paste("sum_logTissueConc", mfc_logTissueConc, sep="")   #recommended regression
coef_logTissueConc=summary(eval(parse(text=best_rg_logTissueConc)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))


#######################################
###Regression soil conc################
#######################################
mlrfrog_sc_test=subset(mlrfrog,select=c(logSoilConc,app_rate_g_cm2,logKow,Hlaw,VapPrs_mPa))
cor_matrix_sc_test=cor(subset(mlrfrog_sc_test))

#####Based on the correlation matrix, I can only keep a.logkow, b.App_rate
mlrfrog_sc=subset(mlrfrog_sc_test,select=-c(logSoilConc, Hlaw,VapPrs_mPa)) 
cor_matrix_sc=cor(subset(mlrfrog_sc))

detail_select=expand.grid(rep(list(c(TRUE,FALSE)),2))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_sc)
regressors <- names(mlrfrog_sc)
total_ita=dim(detail_select)[1]
R2_logSoilConc=array(dim=c(total_ita-1,4))
colnames(R2_logSoilConc)=c("R2", "Adj_R2", "AIC", "BIC")

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
mfc_logSoilConc=final_logSoilConc[4]#find the most common combination ID
best_rg_logSoilConc=paste("sum_logSoilConc", mfc_logSoilConc, sep="")   #recommended regression
coef_logSoilConc=summary(eval(parse(text=best_rg_logSoilConc)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))




#######################################
###Regression appfactor################
#######################################
mlrfrog_af_test=subset(mlrfrog,select=c(logAppFactor,Species,HabFac,logKow,Hlaw,VapPrs_mPa))
cor_matrix_af_test=cor(subset(mlrfrog_af_test,select=-c(HabFac,Species)))

#####Based on the correlation matrix, I can only keep a.logkow, b.HabFac, c.Species
mlrfrog_af=subset(mlrfrog_af_test,select=-c(logAppFactor,Hlaw,VapPrs_mPa)) 
detail_select=expand.grid(rep(list(c(TRUE,FALSE)),3))
detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
row.names(detail_select)=NULL
names(detail_select)=names(mlrfrog_af)
regressors <- names(mlrfrog_af)
total_ita=dim(detail_select)[1]
R2_logAppFactor=array(dim=c(total_ita-1,4))
colnames(R2_logAppFactor)=c("R2", "Adj_R2", "AIC", "BIC")

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
mfc_logAppFactor=final_logAppFactor[4]
best_rg_logAppFactor=paste("sum_logAppFactor", mfc_logAppFactor, sep="")   #recommended regression
coef_logAppFactor=summary(eval(parse(text=best_rg_logAppFactor)))$coefficients 
#coef_std_logAppFactor=lm.beta(eval(parse(text=best_rg_logAppFactor)))









