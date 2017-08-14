library(DAAG)

####load data############
frog_root<-"I:\\Dropbox\\Robin_data_test\\frog_reg9\\" #PC path
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

#############################################################
cor_matrix=cor(subset(mlrfrog_raw,select=c(logkoc, logMol, logSol, logDen)))



  #######################################
  ###Regression BCF################
  #######################################
  #####Based on the correlation matrix, I can only keep a.Species, b.HabFac, c.logkoc, 
  
  mlrfrog_raw=subset(mlrfrog_raw,select=c(logBCF, Species, HabFac, logkoc, logMol, logSol, logDen)) 
  mlrfrog_raw_list=c('Species', 'HabFac', 'logkoc', 'logMol', 'logSol', 'logDen')
  nod=dim(mlrfrog_raw)[2]-1
  detail_select=expand.grid(rep(list(c(TRUE,FALSE)),nod))
  detail_select=detail_select[-which(detail_select$Var1==TRUE & detail_select$Var2==TRUE),]
  row.names(detail_select)=NULL
  names(detail_select)=mlrfrog_raw_list
  regressors=mlrfrog_raw_list
  total_ita=dim(detail_select)[1]
  R2_logBCF=array(dim=c(total_ita-1,1))
  colnames(R2_logBCF)=c("SS")

  ###replace true and false in combination by variable names
  detail_select_v=detail_select
  d_names=names(detail_select)
  for (i in d_names){
    eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",TRUE,"]=","i", sep="")))
    eval(parse(text=paste("detail_select_v$",as.name(i),"[", "detail_select_v$",as.name(i),"==",FALSE,"]=",NA, sep="")))   
  }
  detail_select_v[is.na(detail_select_v)] = " "


#   ####cHECK CORRELATIONS###########
#   cor_matrix=cor(subset(mlrfrog_BCF_ind,select=-c(Species, HabFac)))
#   name_final_logBCF=paste("cor_matrix", "_", k, sep="")
#   assign(name_final_logBCF, (cor_matrix))  
#   
#######################################################################################################
cv.lm_tao = function (df = houseprices, form.lm = formula(sale.price ~ area), 
          m = 3, dots = FALSE, seed = 29, plotit = c("Observed", "Residual"), 
          main = "Small symbols show cross-validation predicted values", 
          legend.pos = "topleft", printit = TRUE) 
{
  gphtype <- ""
  if (is.logical(plotit)) {
    if (plotit) 
      gphtype <- "Observed"
  }
  else if (is.character(plotit)) {
    if (!(plotit[1] %in% c("Observed", "Residual", ""))) 
      stop(paste("Illegal argument plotit =", plotit[1]))
    gphtype <- plotit[1]
    if (plotit[1] %in% c("Observed", "Residual")) 
      plotit <- TRUE
  }
  else stop("Argument plotit must be logical or character")
  if (class(form.lm) == "formula") 
    form <- form.lm
  else if (class(form.lm) %in% c("call", "lm")) 
    form <- formula(form.lm)
  else stop("form.lm must be formula or call or lm object")
  formtxt <- deparse(form)
  mf <- model.frame(form, data = df)
  ynam <- attr(mf, "names")[attr(attr(mf, "terms"), "response")]
  df.lm <- lm(mf)
  tm <- terms(mf)
  xcolumns <- labels(tm)
  n <- nrow(df)
  df[, ynam] <- model.response(mf)
  df[, "Predicted"] <- predict(df.lm)
  df[, "cvpred"] <- numeric(n)
  yval <- mf[, ynam]
  if (gphtype == "Residual") 
    yval <- yval - df[, "Predicted"]
  if (!is.null(seed)) 
    set.seed(seed)
  n <- dim(df)[1]
  rand <- sample(n)%%m + 1
  foldnum <- sort(unique(rand))
  for (i in foldnum) {
    rows.in <- rand != i
    rows.out <- rand == i
    subs.lm <- lm(form, data = df[rows.in, ])
    df[rows.out, "cvpred"] <- predict(subs.lm, newdata = df[rows.out, 
                                                            ])
  }
  if (length(xcolumns) == 1) {
    stline <- TRUE
    xnam <- xcolumns
  }
  else {
    stline <- FALSE
    xnam <- "Predicted"
  }
  if (printit) {
    options(digits = 3)
    print(anova(df.lm))
    cat("\n")
  }
  if (plotit) {
    oldpar <- par(mar = par()$mar - c(1, 0, 2, 0))
    on.exit(par(oldpar))
    coltypes <- palette()[c(2, 3, 6, 1, 4:5, 7)]
    if (m > 7) 
      coltypes <- c(coltypes, rainbow(m - 7))
    ltypes <- 1:m
    ptypes <- 2:(m + 1)
    par(lwd = 2)
    if (stline) 
      xlab <- xnam
    else {
      xlab <- "Predicted (fit to all data)"
      cat("\n")
      warning(paste("\n\n As there is >1 explanatory variable, cross-validation\n", 
                    "predicted values for a fold are not a linear function\n", 
                    "of corresponding overall predicted values.  Lines that\n", 
                    "are shown for the different folds are approximate\n"))
    }
    ylab <- ynam
    if (gphtype == "Residual") 
      ylab <- paste(ynam, " (offset from predicted using all data)")
    plot(as.formula(paste("yval ~", xnam)), data = df, ylab = ylab, 
         type = "p", pch = ptypes[rand], col = coltypes[rand], 
         cex = 1.25, xlab = xlab)
    title(main = main, cex = 1.05)
    if (dots) {
      with(df, points(as.formula(paste("yval ~", xnam)), 
                      data = df, type = "p", pch = 16, col = coltypes[rand], 
                      cex = 1))
    }
  }
  if (printit | plotit) {
    sumss <- 0
    sumdf <- 0
    for (i in foldnum) {
      rows.in <- rand != i
      rows.out <- rand == i
      n.out <- sum(rows.out)
      resid <- df[rows.out, ynam] - df[rows.out, "cvpred"]
      ss <- sum(resid^2)
      sumss <- sumss + ss
      if (printit) {
        fold_data <- t(cbind(df[rows.out, c(xnam, "cvpred", 
                                            ynam)], resid))
        rownames(fold_data) = c(xnam, "cvpred", ynam, 
                                "CV residual")
        cat("\nfold", i, "\n")
        cat("Observations in test set:", n.out, "\n")
        print(fold_data, collab = rep("", n.out))
        cat("\nSum of squares =", round(ss, 2), "   Mean square =", 
            round(ss/n.out, 2), "   n =", n.out, "\n")
      }
      if (plotit) {
        xval <- df[rows.out, xnam]
        nminmax <- c(which.min(xval), which.max(xval))
        cvpred <- df[rows.out, "cvpred"]
        if (gphtype == "Residual") 
          cvpred <- cvpred - df[rows.out, "Predicted"]
        points(xval, cvpred, col = coltypes[i], pch = ptypes[i], 
               cex = 0.75, lwd = 1)
        n1 <- which.min(xval)
        n2 <- which.max(xval)
        fold.lm <- lm(cvpred ~ xval)
        fold.b <- coef(fold.lm)
        lines(xval[c(n1, n2)], fold.b[1] + fold.b[2] * 
                xval[c(n1, n2)], col = coltypes[i], lty = ltypes[i])
        topleft <- par()$usr[c(1, 4)]
        par(lwd = 1, col = 1)
        legend(x = legend.pos, legend = paste("Fold", 
                                              1:m), pch = ptypes, lty = ltypes, col = coltypes, 
               cex = 0.75)
      }
    }
  }
  sumdf <- sum(!is.na(df[, "Predicted"]))
  if (printit) {
    cat("\nOverall", "(Sum over all", n.out, "folds)", "\n")
    print(c(ms = sumss/sumdf))
  }
  attr(df, "ms") <- sumss/sumdf
  attr(df, "df") <- sumdf
  invisible(df)
  return (sumss/sumdf)
}
#####################################################################################################
  
  for (i in 1:(total_ita-1)){
    vect=as.matrix(detail_select)[i,]
    lm_temp=as.formula(paste("logBCF", paste(regressors[vect], collapse=" + "), sep="~"))

    cv_temp=cv.lm_tao(mlrfrog_raw, lm_temp, printit=FALSE)

    R2_logBCF[i,1]=cv_temp
    
  }
  
  
  ####The highest R2(1), adjusted R2(2), lowest AIC(3) and BIC(4)
  col_best_logBCF=c(apply(R2_logBCF,2,min)[1])
  final_logBCF=c(which(R2_logBCF[,1]==col_best_logBCF[1]))
  








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





