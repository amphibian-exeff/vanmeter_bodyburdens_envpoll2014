library(modeest)
library(QuantPsyc)

####load data############
frog_root<-"D:\\Dropbox\\Robin_data_test\\frog_reg3\\" #PC path
setwd(frog_root)
raw_data<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") #Read data into a variable
raw_data_t1 = raw_data[which(raw_data$good==1),] #select good data
Chemical_selec_index=is.element(raw_data_t1$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon'))
raw_data_t2 = raw_data_t1[Chemical_selec_index,] #select only five chemicals
row.names(raw_data_t2)<-NULL
raw_data_t2$bodyweight2=raw_data_t2$bodyweight/raw_data_t2$X  #average the weight by number of frogs
raw_data_t2$SA_cm2_2=raw_data_t2$SA_cm2/raw_data_t2$X  #average the surface area by number of frogs

mlrfrog = subset(raw_data_t2, select=c(Species, app_rate_g_cm2, AppFactor, TissueConc, SoilConc, logKow, BCF, VapPrs_mPa, Koc_gmL,
                                           HabFac, Hlaw, bodyweight2, SA_cm2_2)) 
mlrfrog$Koc_gmL=log(mlrfrog$Koc_gmL)
names(mlrfrog)[9]='logkoc'
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
mlrfrog_BCF_test=subset(mlrfrog,select=c(logBCF,Species,HabFac,logKow,logkoc,Hlaw,VapPrs_mPa))
cor_matrix_BCF_test=cor(subset(mlrfrog_BCF_test,select=-c(HabFac,Species)))

#####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, 
mlrfrog_BCF=subset(mlrfrog_BCF_test,select=-c(logBCF,logKow,Hlaw,VapPrs_mPa)) 





func3<-function(response_var, input_var) {
  
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),dim(input_var)[2])) #Check the number of combinations
  regressors <- names(input_var)
  total_ita=dim(detail_select)[1]
  R2_regression=array(dim=c(total_ita-1,4))
  colnames(R2_regression)=c("R2", "Adj_R2", "AIC", "BIC")
  
  response_var_name=deparse(substitute(response_var))
  input_var_name=deparse(substitute(input_var))
  
  
  for (i in 1:(total_ita-1)){
  vect=as.matrix(detail_select)[i,]
  lm_temp=lm(as.formula(paste(response_var_name, paste(input_var_name, "$", regressors[vect], collapse=" + "), sep="~")))
  name_regression=paste("sum", i, sep="")
  assign(name_regression, (lm_temp))
  R2_regression[i,1]=summary(lm_temp)$r.squared
  R2_regression[i,2]=summary(lm_temp)$adj.r.squared
  R2_regression[i,3]=AIC(lm_temp)
  R2_regression[i,4]=BIC(lm_temp)
  }
  
  col_best_regression=c(apply(R2_regression,2,max)[1:2],apply(R2_regression,2,min)[3:4])
  final_regression=c(which(R2_regression[,1]==col_best_regression[1]),which(R2_regression[,2]==col_best_regression[2]),
                 which(R2_regression[,3]==col_best_regression[3]),which(R2_regression[,4]==col_best_regression[4]))
  mfc_regression=mfv(final_regression)#find the most common combination ID
  best_rg_regression=paste("sum_logBCF", mfc_logBCF, sep="")   #recommended regression
  coef_logBCF=summary(eval(parse(text=best_rg_logBCF)))$coefficients   
  
  
  
  
  return(R2_regression)
  
}
  
aaa=func3(mlrfrog$logBCF,mlrfrog_BCF)
  
  

  
  
  
  


