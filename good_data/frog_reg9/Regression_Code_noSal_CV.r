library(modeest)
library(QuantPsyc)

####load data############
frog_root<-"D:\\Dropbox\\Robin_data_test\\frog_reg9\\" #PC path
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


# ############create log-likelihood function################

#########################################
##  loglike_test(lm_temp$fit, mean(mlrfrog$logBCF), sd(mlrfrog$logBCF))=-110.8256
##  logLik(lm_temp)=-154
#########################################

loglike_test = function(data, mu, sigma) {
  loglike = 0
  for(obs in data){
    loglike = loglike +
      log(1/(sqrt(2*pi)*sigma) *
            exp(-1/2 * (obs - mu)^2/(sigma^2)))
  }
  return(loglike)
}
##############################################################

#########################################
######Subset the data#########################
#########################################
C_name_pool=c('Imidacloprid','Pendimethalin','Total Atrazine','Total Fipronil','Total Triadimefon')


#######################################
###Regression BCF################
#######################################

for (k in 1:5){
  
  C_name=C_name_pool[k]
  mlrfrog=mlrfrog_raw[which(mlrfrog_raw$Chemical != C_name),]
  mlrfrog_CV=mlrfrog_raw[which(mlrfrog_raw$Chemical == C_name),]

  #####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, 
  
  mlrfrog_BCF_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
  mlrfrog_BCF_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
  
  ####CHECK CORRELATIONS###########
  cor_matrix=cor(subset(mlrfrog_BCF_ind,select=-c(Species, HabFac)))
  name_final_logBCF=paste("cor_matrix", "_", k, sep="")
  assign(name_final_logBCF, (cor_matrix))  
  
  nod=dim(mlrfrog_BCF_ind)[2]
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
  detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
  row.names(detail_select)=NULL
  names(detail_select)=names(mlrfrog_BCF_ind)
  regressors <- names(mlrfrog_BCF_ind)
  total_ita=dim(detail_select)[1]
  R2_logBCF=array(dim=c(total_ita-1,8))
  colnames(R2_logBCF)=c("R2", "Adj_R2", "AIC", "BIC", "LL_training", "MSE_t", "LL_testing", "R2_t")
  
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
    lm_temp=lm(as.formula(paste("mlrfrog$logBCF", paste(regressors[vect], collapse=" + "), sep="~")), data=subset(mlrfrog_BCF_ind, select=regressors[vect]))
    
    name_logBCF=paste("sum_logBCF", i, "_", k, sep="")
    assign(name_logBCF, (lm_temp))  
    
    lm_temp_CV=predict(lm_temp, data.frame(mlrfrog_CV$logBCF, subset(mlrfrog_BCF_ind_CV, select=regressors[vect])), se.fit = TRUE, df=25)
    name_logBCF_CV=paste("sum_logBCF_CV", i, "_", k, sep="")
    assign(name_logBCF_CV, (lm_temp_CV))  
    
    R2_logBCF[i,1]=summary(lm_temp)$r.squared
    R2_logBCF[i,2]=summary(lm_temp)$adj.r.squared
    R2_logBCF[i,3]=AIC(lm_temp)
    R2_logBCF[i,4]=BIC(lm_temp)
#     R2_logBCF[i,5]=logLik(lm_temp)
    R2_logBCF[i,5]=sum((lm_temp$fit-mlrfrog$logBCF)^2)/(lm_temp$df.residual)
    R2_logBCF[i,6]=sum((lm_temp_CV$fit-mlrfrog_CV$logBCF)^2)/(dim(mlrfrog_CV)[1]-length(lm_temp$coefficients))
    
    R2_logBCF[i,7]=loglike_test((lm_temp_CV$fit-mlrfrog_CV$logBCF), 0, R2_logBCF[i,6]^0.5)
    
    R2_logBCF[i,8]=cor(cbind(lm_temp_CV$fit, mlrfrog_CV$logBCF))[1,2]^2
  }
  
  name_R2_logBCF=paste("R2_logBCF", "_", k, sep="")
  assign(name_R2_logBCF, (R2_logBCF))  
  
  
  ####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
  col_best_logBCF=c(apply(R2_logBCF,2,max)[1:2],apply(R2_logBCF,2,min)[3:7])
  final_logBCF=c(which(R2_logBCF[,1]==col_best_logBCF[1]),which(R2_logBCF[,2]==col_best_logBCF[2]),
                       which(R2_logBCF[,3]==col_best_logBCF[3]),which(R2_logBCF[,4]==col_best_logBCF[4]),
                       which(R2_logBCF[,5]==col_best_logBCF[5]),which(R2_logBCF[,6]==col_best_logBCF[6]),
                       which(R2_logBCF[,7]==col_best_logBCF[7]))
  name_final_logBCF=paste("final_logBCF", "_", k, sep="")
  assign(name_final_logBCF, (final_logBCF))  
  
}

LL_ovarall_logBCF_trainging=as.matrix(rowMeans(cbind(R2_logBCF_1[,5],R2_logBCF_2[,5],R2_logBCF_3[,5],R2_logBCF_4[,5],R2_logBCF_5[,5])))
find_LL_ovarall_logBCF_training=which(LL_ovarall_logBCF_trainging==max(LL_ovarall_logBCF_trainging))

LL_ovarall_logBCF_tesing=as.matrix(rowSums(cbind(R2_logBCF_1[,7],R2_logBCF_2[,7],R2_logBCF_3[,7],R2_logBCF_4[,7],R2_logBCF_5[,7])))
find_LL_ovarall_logBCF_tesing=which(LL_ovarall_logBCF_tesing==max(LL_ovarall_logBCF_tesing))
               
MSE_ovarall_logBCF=as.matrix(rowMeans(cbind(R2_logBCF_1[,6],R2_logBCF_2[,6],R2_logBCF_3[,6],R2_logBCF_4[,6],R2_logBCF_5[,6])))
find_MSE_ovarall_logBCF=which(MSE_ovarall_logBCF==min(MSE_ovarall_logBCF))

########plot##########
result_BCF=as.data.frame(cbind(c(1:47), LL_ovarall_logBCF_trainging, LL_ovarall_logBCF_tesing))

names(result_BCF)=c('ID', 'Train', 'Test')
dfm <- melt(result_BCF[,names(result_BCF)],id.vars = 1)
ggplot(dfm, aes(x = ID, y = value, fill=variable)) + geom_bar(stat = "identity", position = "dodge") + ylim(-1000, 0)+scale_x_continuous(breaks = seq(1,47,2))


#######################################
###Regression Tissue conc################
#######################################
mlrfrog_raw_TC=subset(mlrfrog_raw_TC, )

for (k in 1:5){
  
  C_name=C_name_pool[k]
  mlrfrog=mlrfrog_raw[which(mlrfrog_raw$Chemical != C_name),]
  mlrfrog_CV=mlrfrog_raw[which(mlrfrog_raw$Chemical == C_name),]

  #####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, d.app_rate_g_cm2, e.bodyweight2, f.SA_cm2_2, g. molmass_gmol
  mlrfrog_tc_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, app_rate_g_cm2, bodyweight2, SA_cm2_2, logMol, logSol, logDen))
  mlrfrog_tc_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, app_rate_g_cm2, bodyweight2, SA_cm2_2, logMol, logSol, logDen))

  ####cHECK CORRELATIONS###########
  cor_matrix=cor(subset(mlrfrog_tc_ind,select=-c(Species, HabFac)))
  name_final_logTissueConc=paste("cor_matrix", "_", k, sep="")
  assign(name_final_logTissueConc, (cor_matrix)) 


  nod=dim(mlrfrog_tc_ind)[2]
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
  detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
  detail_select=detail_select[-which(detail_select$Var5==TRUE & detail_select$Var6==TRUE),]
  row.names(detail_select)=NULL
  names(detail_select)=names(mlrfrog_tc_ind)
  regressors <- names(mlrfrog_tc_ind)
  total_ita=dim(detail_select)[1]
  R2_logTissueConc=array(dim=c(total_ita-1,7))
  colnames(R2_logTissueConc)=c("R2", "Adj_R2", "AIC", "BIC", "LL", "MSE_t", "R2_t")

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
    lm_temp=lm(as.formula(paste("mlrfrog$logTissueConc", paste(regressors[vect], collapse=" + "), sep="~")), data=subset(mlrfrog_tc_ind, select=regressors[vect]))
    name_logTissueConc=paste("sum_logTissueConc", i, "_", k, sep="")
    assign(name_logTissueConc, (lm_temp))
    
    lm_temp_CV=predict(lm_temp, data.frame(mlrfrog_CV$logTissueConc, subset(mlrfrog_tc_ind_CV, select=regressors[vect])), se.fit = TRUE)
    name_logTissueConc_CV=paste("sum_logTissueConc_CV", i, "_", k, sep="")
    assign(name_logTissueConc_CV, (lm_temp_CV))    
    
    R2_logTissueConc[i,1]=summary(lm_temp)$r.squared
    R2_logTissueConc[i,2]=summary(lm_temp)$adj.r.squared
    R2_logTissueConc[i,3]=AIC(lm_temp)
    R2_logTissueConc[i,4]=BIC(lm_temp)
    R2_logTissueConc[i,5]=loglike_test(lm_temp_CV$fit, mean(mlrfrog_CV$logTissueConc), sd(mlrfrog_CV$logTissueConc))
    R2_logTissueConc[i,6]=sum(lm_temp_CV$fit-mlrfrog_CV$logTissueConc)^2/(dim(mlrfrog_CV)[1]-length(lm_temp$coefficients))
    R2_logTissueConc[i,7]=corr(cbind(lm_temp_CV$fit, mlrfrog_CV$logTissueConc))^2
    
  }

  
  name_R2_logTissueConc=paste("R2_logTissueConc", "_", k, sep="")
  assign(name_R2_logTissueConc, (R2_logTissueConc))  
  
  ####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
  col_best_logTissueConc=c(apply(R2_logTissueConc,2,max)[1:2],apply(R2_logTissueConc,2,min)[3:6])
  final_logTissueConc=c(which(R2_logTissueConc[,1]==col_best_logTissueConc[1]),which(R2_logTissueConc[,2]==col_best_logTissueConc[2]),
                        which(R2_logTissueConc[,3]==col_best_logTissueConc[3]),which(R2_logTissueConc[,4]==col_best_logTissueConc[4]),
                        which(R2_logTissueConc[,5]==col_best_logTissueConc[5]),which(R2_logTissueConc[,6]==col_best_logTissueConc[6]))
  name_final_logTissueConc=paste("final_logTissueConc", "_", k, sep="")
  assign(name_final_logTissueConc, (final_logTissueConc))  

}

LL_ovarall_logTissueConc=as.matrix(rowSums(cbind(R2_logTissueConc_1[,5],R2_logTissueConc_2[,5],R2_logTissueConc_3[,5],R2_logTissueConc_4[,5],R2_logTissueConc_5[,5])))
find_LL_ovarall_logTissueConc=which(LL_ovarall_logTissueConc==max(LL_ovarall_logTissueConc))
               
MSE_ovarall_logTissueConc=as.matrix(rowMeans(cbind(R2_logTissueConc_1[,6],R2_logTissueConc_2[,6],R2_logTissueConc_3[,6],R2_logTissueConc_4[,6],R2_logTissueConc_5[,6])))
find_MSE_ovarall_logTissueConc=which(MSE_ovarall_logTissueConc==min(MSE_ovarall_logTissueConc))





#######################################
###Regression appfactor################
#######################################

for (k in 1:5){
  
  C_name=C_name_pool[k]
  mlrfrog=mlrfrog_raw[which(mlrfrog_raw$Chemical != C_name),]
  mlrfrog_CV=mlrfrog_raw[which(mlrfrog_raw$Chemical == C_name),]

  #####Based on the correlation matrix, I can only keep a.logkoc, b.HabFac, c.Species
  mlrfrog_af_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
  mlrfrog_af_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, logMol, logSol, logDen))


  nod=dim(mlrfrog_af_ind)[2]
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
  detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
  row.names(detail_select)=NULL
  names(detail_select)=names(mlrfrog_af_ind)
  regressors=names(mlrfrog_af_ind)
  total_ita=dim(detail_select)[1]
  R2_logAppFactor=array(dim=c(total_ita-1,7))
  colnames(R2_logAppFactor)=c("R2", "Adj_R2", "AIC", "BIC", "LL", "MSE_t", "R2_t")

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
  lm_temp=lm(as.formula(paste("mlrfrog$logAppFactor", paste(regressors[vect], collapse=" + "), sep="~")), data=subset(mlrfrog_af_ind, select=regressors[vect]))
  name_logAppFactor=paste("sum_logAppFactor", i, "_", k, sep="")
  assign(name_logAppFactor, (lm_temp))

  lm_temp_CV=predict(lm_temp, data.frame(mlrfrog_CV$logAppFactor, subset(mlrfrog_af_ind_CV, select=regressors[vect])), se.fit = TRUE)
  name_logAppFactor_CV=paste("sum_logAppFactor_CV", i, "_", k, sep="")
  assign(name_logAppFactor_CV, (lm_temp_CV))    

  R2_logAppFactor[i,1]=summary(lm_temp)$r.squared
  R2_logAppFactor[i,2]=summary(lm_temp)$adj.r.squared
  R2_logAppFactor[i,3]=AIC(lm_temp)
  R2_logAppFactor[i,4]=BIC(lm_temp)
  R2_logAppFactor[i,5]=loglike_test(lm_temp_CV$fit, mean(mlrfrog_CV$logAppFactor), sd(mlrfrog_CV$logAppFactor))
  R2_logAppFactor[i,6]=mean((lm_temp_CV$fit-mlrfrog_CV$logAppFactor)^2)
  R2_logAppFactor[i,7]=corr(cbind(lm_temp_CV$fit, mlrfrog_CV$logAppFactor))^2 
  }

  name_R2_logAppFactor=paste("R2_logAppFactor", "_", k, sep="")
  assign(name_R2_logAppFactor, (R2_logAppFactor))  


  ####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
  col_best_logAppFactor=c(apply(R2_logAppFactor,2,max)[1:2],apply(R2_logAppFactor,2,min)[3:5])
  final_logAppFactor=c(which(R2_logAppFactor[,1]==col_best_logAppFactor[1]),which(R2_logAppFactor[,2]==col_best_logAppFactor[2]),
                       which(R2_logAppFactor[,3]==col_best_logAppFactor[3]),which(R2_logAppFactor[,4]==col_best_logAppFactor[4]))
  name_final_logAppFactor=paste("final_logAppFactor", "_", k, sep="")
  assign(name_final_logAppFactor, (final_logAppFactor))  
}

LL_ovarall_logAppFactor=as.matrix(rowSums(cbind(R2_logAppFactor_1[,5],R2_logAppFactor_2[,5],R2_logAppFactor_3[,5],R2_logAppFactor_4[,5],R2_logAppFactor_5[,5])))
find_LL_ovarall_logAppFactor=which(LL_ovarall_logAppFactor==max(LL_ovarall_logAppFactor))
               
MSE_ovarall_logAppFactor=as.matrix(rowMeans(cbind(R2_logAppFactor_1[,6],R2_logAppFactor_2[,6],R2_logAppFactor_3[,6],R2_logAppFactor_4[,6],R2_logAppFactor_5[,6])))
find_MSE_ovarall_logAppFactor=which(MSE_ovarall_logAppFactor==min(MSE_ovarall_logAppFactor))



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





