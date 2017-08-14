library(modeest)
library(QuantPsyc)

####load data############
frog_root<-"I:\\Dropbox\\Robin_data_test\\frog_reg11\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1 & raw_data$Application=="Soil"),] #select good data and Soil application
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon'))

raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog_raw = subset(raw_data_t2, select=c(Chemical, Species, app_rate_g_cm2, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa, Koc_gmL,
                                           HabFac, bodyweight2, SA_cm2_2, molmass_gmol, Solat20C_mgL, Density_gcm3)) 
mlrfrog_raw$Koc_gmL=as.numeric(as.character(mlrfrog_raw$Koc_gmL)) #convert factors to numerics
mlrfrog_raw$Koc_gmL=log(mlrfrog_raw$Koc_gmL)
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
names(mlrfrog_raw)[14]='logMol'
mlrfrog_raw$Solat20C_mgL=log(mlrfrog_raw$Solat20C_mgL)
names(mlrfrog_raw)[15]='logSol'
mlrfrog_raw$Density_gcm3=log(mlrfrog_raw$Density_gcm3)
names(mlrfrog_raw)[16]='logDen'
mlrfrog_raw = mlrfrog_raw[which(mlrfrog_raw$Species != "Mole salamander"),] #remove Mole salamander
row.names(mlrfrog_raw)<-NULL


# ############create log-likelihood function################

#########################################
##  loglike_test(lm_temp$fit, mean(mlrfrog$logAPP), sd(mlrfrog$logAPP))=-110.8256
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
###Regression APP################
#######################################

for (k in 1:5){
  
  C_name=C_name_pool[k]
  mlrfrog=mlrfrog_raw[which(mlrfrog_raw$Chemical != C_name),]
  mlrfrog_CV=mlrfrog_raw[which(mlrfrog_raw$Chemical == C_name),]

  #####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, 
  
  mlrfrog_APP_ind=subset(mlrfrog,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
  mlrfrog_APP_ind_CV=subset(mlrfrog_CV,select=c(Species, HabFac, logkoc, logMol, logSol, logDen)) 
  
  ####CHECK CORRELATIONS###########
  cor_matrix=cor(subset(mlrfrog_APP_ind,select=-c(Species, HabFac)))
  name_final_logAPP=paste("cor_matrix", "_", k, sep="")
  assign(name_final_logAPP, (cor_matrix))  
  
  nod=dim(mlrfrog_APP_ind)[2]
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
  detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
  row.names(detail_select)=NULL
  names(detail_select)=names(mlrfrog_APP_ind)
  regressors <- names(mlrfrog_APP_ind)
  total_ita=dim(detail_select)[1]
  R2_logAPP=array(dim=c(total_ita-1,11))
  colnames(R2_logAPP)=c("R2", "Adj_R2", "AIC", "BIC", "MSE_T", "LL_T", "LL_R", "df_T", "df_R", "N_reg", "N_var")
  
  ###replace true and false in combination by variable names
  detail_select_v=detail_select
  d_names=names(detail_select)
  
  for (i in d_names){
    eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
    eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
  }
  detail_select_v[is.na(detail_select_v)] = " "

  
  
  NCV_res=c()
  CV_res=c()
  NCV_obs=c()
  CV_obs=c()
  
  for (i in 1:(total_ita-1)){
    vect=as.matrix(detail_select)[i,]
    lm_temp=lm(as.formula(paste("mlrfrog$logAppFactor", paste(regressors[vect], collapse=" + "), sep="~")), data=subset(mlrfrog_APP_ind, select=regressors[vect]))
    name_logAPP=paste("sum_logAPP", i, "_", k, sep="")
    assign(name_logAPP, (lm_temp))  
    
    lm_temp_CV=predict(lm_temp, data.frame(mlrfrog_CV$logAppFactor, subset(mlrfrog_APP_ind_CV, select=regressors[vect])), se.fit = TRUE)
    name_logAPP_CV=paste("sum_logAPP_CV", i, "_", k, sep="")
    assign(name_logAPP_CV, (lm_temp_CV))  
    
    ###creat a matrix to hold residules for final check##
    NCV_res=cbind(NCV_res, lm_temp$residuals)
    NCV_obs=cbind(NCV_obs, mlrfrog$logAppFactor)
    
    CV_res_temp=lm_temp_CV$fit-mlrfrog_CV$logAppFactor
    CV_res=cbind(CV_res, CV_res_temp)
    CV_obs=cbind(CV_obs, mlrfrog_CV$logAppFactor)
    
    R2_logAPP[i,1]=summary(lm_temp)$r.squared
    R2_logAPP[i,2]=summary(lm_temp)$adj.r.squared
    R2_logAPP[i,3]=AIC(lm_temp)
    R2_logAPP[i,4]=BIC(lm_temp)
    R2_logAPP[i,5]=sum(lm_temp_CV$fit-mlrfrog_CV$logAppFactor)^2/(dim(mlrfrog_CV)[1]-length(lm_temp$coefficients))
    R2_logAPP[i,6]=loglike_test((lm_temp_CV$fit-mlrfrog_CV$logAppFactor), 0, R2_logAPP[i,5]^0.5)
    R2_logAPP[i,7]=logLik(lm_temp)
    R2_logAPP[i,8]=dim(mlrfrog_CV)[1]-length(lm_temp$coefficients)
    R2_logAPP[i,9]=lm_temp$df.residual
    R2_logAPP[i,10]=length(lm_temp$coefficients)-1
    R2_logAPP[i,11]=sum(vect)
  }
  
#  (1 "R2", 2"Adj_R2", 3"AIC", 4"BIC", 5"MSE_T", 6"LL_T", 7"LL_R", 8"df_T", 9"df_R", 10"N_reg", 11"N_var")
  
  colnames(NCV_res)=seq(1:dim(NCV_res)[2])  
  name_NCV_res=paste("NCV_res", "_", k, sep="")
  assign(name_NCV_res, (NCV_res))    
  
  colnames(CV_res)=seq(1:dim(CV_res)[2])
  name_CV_res=paste("CV_res", "_", k, sep="")
  assign(name_CV_res, (CV_res))  

  colnames(NCV_obs)=seq(1:dim(NCV_obs)[2])  
  name_NCV_obs=paste("NCV_obs", "_", k, sep="")
  assign(name_NCV_obs, (NCV_obs))  

  colnames(CV_obs)=seq(1:dim(CV_obs)[2])  
  name_CV_obs=paste("CV_obs", "_", k, sep="")
  assign(name_CV_obs, (CV_obs)) 
  
  name_R2_logAPP=paste("R2_logAPP", "_", k, sep="")
  assign(name_R2_logAPP, (R2_logAPP))  

  
  ###The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
  col_best_logAPP=c(apply(R2_logAPP,2,max)[1:2],apply(R2_logAPP,2,min)[3:4])
  final_logAPP=c(which(R2_logAPP[,1]==col_best_logAPP[1]),which(R2_logAPP[,2]==col_best_logAPP[2]),
                       which(R2_logAPP[,3]==col_best_logAPP[3]),which(R2_logAPP[,4]==col_best_logAPP[4]))
  name_final_logAPP=paste("final_logAPP", "_", k, sep="")
  assign(name_final_logAPP, (final_logAPP))  
}

CV_res_all=rbind(CV_res_1,CV_res_2,CV_res_3,CV_res_4,CV_res_5)
NCV_res_all=rbind(NCV_res_1,NCV_res_2,NCV_res_3,NCV_res_4,NCV_res_5)

NCV_obs_all=rbind(NCV_obs_1,NCV_obs_2,NCV_obs_3,NCV_obs_4,NCV_obs_5)
CV_obs_all=rbind(CV_obs_1,CV_obs_2,CV_obs_3,CV_obs_4,CV_obs_5)

MC_t=1000

NCV_MC_RMSE_Final=data.frame(matrix(ncol = MC_t, nrow = 47))
NCV_MC_LL_Final=data.frame(matrix(ncol = MC_t, nrow = 47))
CV_LL_Final=data.frame(matrix(ncol = 1, nrow = 47))


for (kk in 1:MC_t){
  
  NCV_obs_MC_index=apply(NCV_obs_all, 2, function(t) sample(length(t), 131, replace=TRUE))  #draw 131 samples of each combinaions
  

  for (jj in 1:dim(NCV_obs_all)[2]){
    
    NCV_res_MC_temp=NCV_res_all[NCV_obs_MC_index[,jj],jj]
#     NCV_MC_RMSE_Final[jj,kk]=(sum(NCV_res_MC_temp^2)/131)^0.5
    NCV_MC_RMSE_Final[jj,kk]=(var(NCV_res_MC_temp)+mean(NCV_res_MC_temp)^2)^0.5
    
    NCV_obs_MC_temp=NCV_obs_all[NCV_obs_MC_index[,jj],jj]
    NCV_fit_MC_temp=NCV_obs_MC_temp+NCV_res_MC_temp
    NCV_MC_LL_Final[jj,kk]=loglike_test(NCV_fit_MC_temp, mean(NCV_obs_MC_temp), sd(NCV_obs_MC_temp))
    CV_LL_Final[jj,1]=loglike_test(CV_obs_all[,jj]+CV_res_all[,jj], mean(CV_obs_all[,jj]), sd(CV_obs_all[,jj]))
    
  }
  
}


CV_RMSE_Final_old=(colSums(CV_res_all^2)/131)^0.5
CV_RMSE_Final=(apply(CV_res_all,2,var)+apply(CV_res_all,2,mean)^2)^0.5

NCV_RMSE_Final_old=(colSums(NCV_res_all^2)/131)^0.5
NCV_RMSE_Final=(apply(NCV_res_all,2,var)+apply(NCV_res_all,2,mean)^2)^0.5

name_RMSE=c()
name_LL=c()

for (z in 1:MC_t){
  
  name_RMSE=cbind(name_RMSE, paste("RMSE_MC_",z,sep=""))
  name_LL=cbind(name_LL, paste("LL_MC_",z,sep=""))
}
      
      
Final_RMSE=data.frame(R2_logAPP[,10], R2_logAPP[,11], CV_RMSE_Final, NCV_RMSE_Final, NCV_MC_RMSE_Final)
names(Final_RMSE)=c('N_reg','N_var','RMSE_T','RMSE_R_A', name_RMSE)

#######Check######
rowMeans(NCV_MC_RMSE_Final)
NCV_RMSE_Final

Final_LL=data.frame(R2_logAPP[,10], R2_logAPP[,11], CV_LL_Final, NCV_MC_LL_Final)
names(Final_LL)=c('N_reg','N_var','LL_T', name_LL)


write.table(Final_RMSE, file = "Final_RMSE_APP.xls", sep="\t", row.names = F)
write.table(Final_LL, file = "Final_LL_APP.xls", sep="\t", row.names = F)


###############Compare###################################################
loglike_test((lm_temp$fit-mlrfrog$logAPP), 0, summary(lm_temp)$sigma)
loglike_test((lm_temp$fit), mean(mlrfrog$logAPP), sd(mlrfrog$logAPP))
logLik(lm_temp)
#########################################################################











########plot##########
library(reshape2)
library(ggplot2)
library(gridExtra)

###boxplot for training MC############
pdf(file = "Fig_app.pdf") #

Final_RMSE_MC=subset(Final_RMSE, select=-c(RMSE_T,RMSE_R_A))
dfm_MC <- melt(Final_RMSE_MC, id.vars = c("N_reg", "N_var"))

p1=ggplot(dfm_MC, aes(x=factor(N_reg), y = value)) + geom_boxplot() + ggtitle('RMSE_MC_Training_Reg')
p2=ggplot(dfm_MC, aes(x=factor(N_var), y = value)) + geom_boxplot() + ggtitle('RMSE_MC_Training_Var')


Final_LL_MC=subset(Final_LL, select=-c(LL_T))
dfm_LL_MC <- melt(Final_LL_MC, id.vars = c("N_reg", "N_var"))

p3=ggplot(dfm_LL_MC, aes(x=factor(N_reg), y = value)) + geom_boxplot() + ggtitle('LL_MC_Training_Reg')
p4=ggplot(dfm_LL_MC, aes(x=factor(N_var), y = value)) + geom_boxplot() + ggtitle('LL_MC_Training_Var')



###boxplot for testing############
Final_RMSE_T=subset(Final_RMSE, select=c(RMSE_T,RMSE_R_A, N_reg, N_var))
dfm_T <- melt(Final_RMSE_T, id.vars = c("N_reg", "N_var"))

p5=ggplot(dfm_T, aes(x=factor(N_reg), y = value)) + geom_boxplot() + ggtitle('RMSE_MC_Testing_Reg')
p6=ggplot(dfm_T, aes(x=factor(N_var), y = value)) + geom_boxplot() + ggtitle('RMSE_MC_Testing_Var')


Final_LL_T=subset(Final_LL, select=c(LL_T, N_reg, N_var))
dfm_LL_T <- melt(Final_LL_T, id.vars = c("N_reg", "N_var"))

p7=ggplot(dfm_LL_T, aes(x=factor(N_reg), y = value)) + geom_boxplot() + ggtitle('LL_MC_Testing_Reg')
p8=ggplot(dfm_LL_T, aes(x=factor(N_var), y = value)) + geom_boxplot() + ggtitle('LL_MC_Testing_Var')

############QQ#########
Final_RMSE_QQ=subset(Final_RMSE, select=c(RMSE_T,RMSE_R_A,(N_reg),(N_var)))

p9=ggplot(Final_RMSE_QQ, aes(RMSE_T, RMSE_R_A)) + geom_point(shape = 2) + facet_wrap(~ N_reg) + ggtitle('QQ_Reg')
p10=ggplot(Final_RMSE_QQ, aes(RMSE_T, RMSE_R_A)) + geom_point(shape = 2) + facet_wrap(~ N_var) + ggtitle('QQ_Var')

grid.arrange(p1,p2,p5,p6, ncol=2, nrow=2)
grid.arrange(p3,p4,p7,p8, ncol=2, nrow=2)
grid.arrange(p9)
grid.arrange(p10)

dev.off()





