#what are units on AppFactor?
mlrfrog$SoilFrogTransferCoeff <- (mlrfrog$AppFactor/mlrfrog$SA_cm2)/1235

View(mlrfrog)
attach(mlrfrog)
summary(mlrfrog)
#summary(logmlrfrog)
names(mlrfrog)
dim(mlrfrog)
class(mlrfrog)
mlrfrog.matrix <- as.matrix(mlrfrog[,-c(1,12,13)])

View(cor(mlrfrog.matrix))
pdf(paste(frog_root,"corrplot_matrix.pdf",sep=""),width=8,height=8)
corrplot(cor(mlrfrog.matrix),type="upper",title="All")
dev.off()

#Fligner-Killeen test of homogeneity of variances
fligner.test(permeability ~ Species)
# very small p-value, reject homogeneity, so no one-way anova

lm.nospecies <- lm(permeability~logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog) #no frog factor
lm.nospecies.log <- lm(log(permeability)~logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog) #no frog factor
summary(lm.nospecies)
summary(lm.nospecies.log)
confint(lm.nospecies)
confint(lm.nospecies.log)

#soil
lm.soil <- lm(SoilConc~mg_aquarium+logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog)
summary(lm.soil)
lm.soil.log <- lm(log(SoilConc)~mg_aquarium+logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog)
summary(lm.soil.log)

ancova.species.1 <- lm(permeability~Species*logKow*Solat20C_mgL*molmass_gmol,data=mlrfrog) #different slopes and intercepts
ancova.species.1.log <- lm(log(permeability)~Species*logKow*Solat20C_mgL*molmass_gmol,data=mlrfrog) #different slopes and intercepts

summary(ancova.species.1)
summary(ancova.species.1.log)
confint(ancova.species.1)
confint(ancova.species.1.log)

ancova.species.2 <- lm(permeability~Species+logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog) #common slope and different intercepts
ancova.species.2.log <- lm(log(permeability)~Species+logKow+Solat20C_mgL+molmass_gmol,data=mlrfrog) #common slope and different intercepts
ancova.species.2.logall <- lm(log(permeability)~Species+logKow+log(Solat20C_mgL)+log(molmass_gmol),data=mlrfrog) #common slope and different intercepts

summary(ancova.species.2)
summary(ancova.species.2.log)
summary(ancova.species.2.logall)

coef(ancova.species.2)
fitted(ancova.species.2)
confint(ancova.species.2)

coef(ancova.species.2.log)
fitted(ancova.species.2.log)
confint(ancova.species.2.log)

coef(ancova.species.2.logall)
fitted(ancova.species.2.logall)
confint(ancova.species.2.logall)

#test residuals for normality
shapiro.test(ancova.species.2.log$resid)

#AICs
AIC(lm.nospecies)
AIC(ancova.species.1)
AIC(ancova.species.2)

AIC(lm.nospecies.log)
AIC(ancova.species.1.log)
AIC(ancova.species.2.log)
AIC(ancova.species.2.logall)

#compare anovas
anova(lm.nospecies,ancova.logKow.2,ancova.logKow.1)
plot(ancova.species.2)


step(ancova.species.2)

#regression plots
#1) residuals vs. fitted values (a good model will show no pattern); 
#2) the qqnorm plot we saw above (values should be on the dashed line); 
#3) scale-location graph, indicating heteroscedasticity ; and 
#4) standardized residuals vs. leverage and Cook's distance, which is handy for identifying outliers. 
pdf(paste(frog_root,"regression_plots.pdf",sep=""),width=11,height=8)
plot(lm.nospecies)
plot(lm.soil.log)
plot(ancova.species.1)
plot(ancova.species.2)
plot(lm.nospecies.log)
plot(ancova.species.1.log)
plot(ancova.species.2.log)
dev.off()



pdf(paste(frog_root,"permeability_v_explanatoryvariables.pdf",sep=""),width=8,height=11)
par(mfrow=c(3,1))
#######logKow
plot(logKow,log(permeability),pch=20)
#barking
points(logKow[index.barking],log(permeability)[index.barking],col="brown",pch=19)
pred.mat <- unique(cbind(logKow[index.barking],fitted(ancova.species.2.log)[index.barking]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="brown")
#fowlers
points(logKow[index.fowlers],log(permeability)[index.fowlers],col="red",pch=19)
pred.mat <- unique(cbind(logKow[index.fowlers],fitted(ancova.species.2.log)[index.fowlers]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="red")
#gray
points(logKow[index.gray],log(permeability)[index.gray],col="gray",pch=19)
pred.mat <- unique(cbind(logKow[index.gray],fitted(ancova.species.2.log)[index.gray]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark gray")
#green
points(logKow[index.green],log(permeability)[index.green],col="dark green",pch=19)
pred.mat <- unique(cbind(logKow[index.green],fitted(ancova.species.2.log)[index.green]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark green")
#leopard
points(logKow[index.leopard],log(permeability)[index.leopard],col="orange",pch=19)
pred.mat <- unique(cbind(logKow[index.leopard],fitted(ancova.species.2.log)[index.leopard]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="orange")
#mole
points(logKow[index.mole],log(permeability)[index.mole],col="purple",pch=19)
pred.mat <- unique(cbind(logKow[index.mole],fitted(ancova.species.2.log)[index.mole]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="purple")
#legend
legend(4, -3.7, c("Fowler's toad", "Leopard frog", "Gray treefrog", "Barking treefrog", "Green treefrog", "Mole salamander"), 
       col = c("red","orange","dark gray","brown","dark green","purple"),
       text.col = "green4", lty = c(1,1,1,1,1,1), pch = c(19,19,19,19,19,19),
       merge = TRUE, bg = 'gray90')
#title
title("ln(Dermal Soil-Frog Transfer Coefficient)~Species+logKow+Solat20C_mgL+molmass_gmol\n Adjusted R^2 = 0.63")

#######molmass_gmol
plot(log(molmass_gmol),log(permeability),pch=20)
#barking
points(log(molmass_gmol[index.barking]),log(permeability)[index.barking],col="brown",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.barking]),fitted(ancova.species.2.log)[index.barking]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="brown")
#fowlers
points(log(molmass_gmol[index.fowlers]),log(permeability)[index.fowlers],col="red",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.fowlers]),fitted(ancova.species.2.log)[index.fowlers]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="red")
#gray
points(log(molmass_gmol[index.gray]),log(permeability)[index.gray],col="gray",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.gray]),fitted(ancova.species.2.log)[index.gray]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark gray")
#green
points(log(molmass_gmol[index.green]),log(permeability)[index.green],col="dark green",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.green]),fitted(ancova.species.2.log)[index.green]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark green")
#leopard
points(log(molmass_gmol[index.leopard]),log(permeability)[index.leopard],col="orange",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.leopard]),fitted(ancova.species.2.log)[index.leopard]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="orange")
#mole
points(log(molmass_gmol[index.mole]),log(permeability)[index.mole],col="purple",pch=19)
pred.mat <- unique(cbind(log(molmass_gmol[index.mole]),fitted(ancova.species.2.log)[index.mole]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="purple")

#######Solat20C_mgL
plot(log10(Solat20C_mgL),log(permeability),pch=20)
#barking
points(log10(Solat20C_mgL[index.barking]),log(permeability)[index.barking],col="brown",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.barking]),fitted(ancova.species.2.log)[index.barking]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="brown")
#fowlers
points(log10(Solat20C_mgL[index.fowlers]),log(permeability)[index.fowlers],col="red",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.fowlers]),fitted(ancova.species.2.log)[index.fowlers]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="red")
#gray
points(log10(Solat20C_mgL[index.gray]),log(permeability)[index.gray],col="gray",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.gray]),fitted(ancova.species.2.log)[index.gray]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark gray")
#green
points(log10(Solat20C_mgL[index.green]),log(permeability)[index.green],col="dark green",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.green]),fitted(ancova.species.2.log)[index.green]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="dark green")
#leopard
points(log10(Solat20C_mgL[index.leopard]),log(permeability)[index.leopard],col="orange",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.leopard]),fitted(ancova.species.2.log)[index.leopard]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="orange")
#mole
points(log10(Solat20C_mgL[index.mole]),log(permeability)[index.mole],col="purple",pch=19)
pred.mat <- unique(cbind(log10(Solat20C_mgL[index.mole]),fitted(ancova.species.2.log)[index.mole]))
pred.sorted <- pred.mat[order(pred.mat[,1],pred.mat[,2]), ]
lines(pred.sorted,col="purple")

dev.off()

# ###############
# unique(mlrfrog$Species)
# unique(mlrfrog$Application.Type)
# mlr_mole <- mlrfrog[which(mlrfrog$Species == "Mole salamander" & mlrfrog$Application.Type == "Soil "),]
# mlr_leopard <- mlrfrog[which(mlrfrog$Species == "Leopard frog"),]
# mlr_fowlers <- mlrfrog[which(mlrfrog$Species == "Fowlers toad"),]
# 
# unique(mlrfrog$Chemical)
# mlr_atrazine<- mlrfrog[which(mlrfrog$Chemical == "Atrazine"),]
# mlr_pendimethalin<- mlrfrog[which(mlrfrog$Chemical == "Pendimethalin"),]
# mlr_fipronil<- mlrfrog[which(mlrfrog$Chemical == "Fipronil"),]
# mlr_triadimefon<- mlrfrog[which(mlrfrog$Chemical == "Triademefon"),]
# 
# pdf(paste(frog_root,"dermal_boxplots.pdf",sep=""),width=11,height=8)
#   par(mfrow=c(2,2))
#   boxplot(BCF~Chemical,data=mlrfrog, main="Amphibian Dermal Exposure", xlab="Pesticide", ylab="Frog-Soil Bioconcentration Factor",col="firebrick2",ylim=c(0,0.33))
#   boxplot(BCF~Chemical,data=mlr_mole, main="Mole salamander Dermal Exposure", xlab="Pesticide", ylab="Frog-Soil Bioconcentration Factor",col="forestgreen",ylim=c(0,0.33))
#   boxplot(BCF~Chemical,data=mlr_leopard, main="Leopard frog Dermal Exposure", xlab="Pesticide", ylab="Frog-Soil Bioconcentration Factor",col="forestgreen",ylim=c(0,0.33))
#   boxplot(BCF~Chemical,data=mlr_fowlers, main="Fowler's toad Dermal Exposure", xlab="Pesticide", ylab="Frog-Soil Bioconcentration Factor",col="forestgreen",ylim=c(0,0.33))
#   boxplot(BCF~Species,data=mlr_atrazine, main="Atrazine Dermal Exposure", xlab="Species", ylab="Frog-Soil Bioconcentration Factor",col="lightskyblue3",ylim=c(0,0.33))
#   boxplot(BCF~Species,data=mlr_pendimethalin, main="Pendimethalin Dermal Exposure", xlab="Species", ylab="Frog-Soil Bioconcentration Factor",col="lightskyblue3",ylim=c(0,0.33))
#   boxplot(BCF~Species,data=mlr_fipronil, main="Fipronil Dermal Exposure", xlab="Species", ylab="Frog-Soil Bioconcentration Factor",col="lightskyblue3",ylim=c(0,0.33))
#   boxplot(BCF~Species,data=mlr_triadimefon, main="Triadimefon Dermal Exposure", xlab="Species", ylab="Frog-Soil Bioconcentration Factor",col="lightskyblue3",ylim=c(0,0.33))
# dev.off()
# 
# #models
# ################################1 variable
# lm1 <- lm(BCF ~ logKow, data = mlrfrog) 
# summary(lm1)
# # Call:
# #   lm(formula = BCF ~ logKow, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.07383 -0.05032 -0.02561  0.02929  0.23159 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)
# # (Intercept) 0.033595   0.039309   0.855    0.397
# # logKow      0.010919   0.009979   1.094    0.279
# # 
# # Residual standard error: 0.06992 on 49 degrees of freedom
# # Multiple R-squared: 0.02385,  Adjusted R-squared: 0.003928 
# # F-statistic: 1.197 on 1 and 49 DF,  p-value: 0.2792 
# 
# lm2 <- lm(BCF ~ Solat20C_mgL, data = mlrfrog)
# summary(lm2)
# # Call:
# #   lm(formula = BCF ~ Solat20C_mgL, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.08412 -0.04129 -0.01397  0.01882  0.20955 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   1.005e-01  1.023e-02   9.826 3.60e-13 ***
# #   Solat20C_mgL -3.251e-04  7.479e-05  -4.348 6.94e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.06011 on 49 degrees of freedom
# # Multiple R-squared: 0.2784,  Adjusted R-squared: 0.2636 
# # F-statistic:  18.9 on 1 and 49 DF,  p-value: 6.937e-05
# 
# lm3 <- lm(BCF ~ molmass_gmol, data = mlrfrog)
# summary(lm3)
# # Call:
# #   lm(formula = BCF ~ molmass_gmol, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.09096 -0.04011 -0.01417  0.02727  0.17326 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -0.0810394  0.0334525  -2.423   0.0192 *  
# #   molmass_gmol  0.0004956  0.0001029   4.818 1.44e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05829 on 49 degrees of freedom
# # Multiple R-squared: 0.3214,  Adjusted R-squared: 0.3076 
# # F-statistic: 23.21 on 1 and 49 DF,  p-value: 1.442e-05 
# 
# lm4 <- lm(BCF ~ bodyweight, data = mlrfrog)
# summary(lm4)
# # Call:
# #   lm(formula = BCF ~ bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.07569 -0.04397 -0.02414  0.02885  0.23303 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)  
# # (Intercept) 0.0736915  0.0292156   2.522    0.015 *
# #   bodyweight  0.0007442  0.0131290   0.057    0.955  
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.07076 on 49 degrees of freedom
# # Multiple R-squared: 6.557e-05,  Adjusted R-squared: -0.02034 
# # F-statistic: 0.003213 on 1 and 49 DF,  p-value: 0.955
# 
# #################################2 variables
# lm5 <- lm(BCF ~ logKow + Solat20C_mgL, data = mlrfrog) 
# summary(lm5)
# # Call:
# #   lm(formula = BCF ~ logKow + Solat20C_mgL, data = mlrfrog)
# # 
# # Residuals:
# #   Min        1Q    Median        3Q       Max 
# # -0.075251 -0.035509 -0.009936  0.019446  0.207872 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   1.478e-01  4.282e-02   3.451  0.00117 ** 
# #   logKow       -1.135e-02  9.987e-03  -1.136  0.26155    
# # Solat20C_mgL -3.762e-04  8.704e-05  -4.322 7.76e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05994 on 48 degrees of freedom
# # Multiple R-squared: 0.2973,  Adjusted R-squared: 0.268 
# # F-statistic: 10.15 on 2 and 48 DF,  p-value: 0.0002104
# 
# lm6 <- lm(BCF ~ logKow + molmass_gmol, data = mlrfrog)
# summary(lm6)
# # Call:
# #   lm(formula = BCF ~ logKow + molmass_gmol, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.09095 -0.04017 -0.01414  0.02724  0.17327 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -8.115e-02  4.150e-02  -1.956   0.0563 .  
# # logKow        4.151e-05  8.735e-03   0.005   0.9962    
# # molmass_gmol  4.954e-04  1.080e-04   4.588 3.23e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.0589 on 48 degrees of freedom
# # Multiple R-squared: 0.3214,  Adjusted R-squared: 0.2931 
# # F-statistic: 11.37 on 2 and 48 DF,  p-value: 9.089e-05 
# 
# lm7 <- lm(BCF ~ logKow + bodyweight, data = mlrfrog)
# summary(lm7)
# # Call:
# #   lm(formula = BCF ~ logKow + bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.07376 -0.05029 -0.02541  0.02920  0.23137 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)
# # (Intercept) 0.0330289  0.0475797   0.694    0.491
# # logKow      0.0109118  0.0100880   1.082    0.285
# # bodyweight  0.0002833  0.0131132   0.022    0.983
# # 
# # Residual standard error: 0.07064 on 48 degrees of freedom
# # Multiple R-squared: 0.02386,  Adjusted R-squared: -0.01681 
# # F-statistic: 0.5866 on 2 and 48 DF,  p-value: 0.5602 
# 
# lm8 <- lm(BCF ~ Solat20C_mgL + molmass_gmol, data = mlrfrog)
# summary(lm8)
# # Call:
# #   lm(formula = BCF ~ Solat20C_mgL + molmass_gmol, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.09960 -0.02392 -0.00846  0.01999  0.16462 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -3.405e-02  3.175e-02  -1.072 0.288895    
# # Solat20C_mgL -2.574e-04  6.557e-05  -3.926 0.000275 ***
# #   molmass_gmol  4.101e-04  9.301e-05   4.409 5.83e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05124 on 48 degrees of freedom
# # Multiple R-squared: 0.4864,  Adjusted R-squared: 0.465 
# # F-statistic: 22.73 on 2 and 48 DF,  p-value: 1.136e-07
# 
# lm9 <- lm(BCF ~ Solat20C_mgL + bodyweight, data = mlrfrog)
# summary(lm9)
# # Call:
# #   lm(formula = BCF ~ Solat20C_mgL + bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.08352 -0.04063 -0.01386  0.01755  0.20665 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   9.314e-02  2.545e-02   3.660 0.000628 ***
# #   Solat20C_mgL -3.265e-04  7.561e-05  -4.319 7.83e-05 ***
# #   bodyweight    3.586e-03  1.128e-02   0.318 0.751843    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.06067 on 48 degrees of freedom
# # Multiple R-squared: 0.2799,  Adjusted R-squared: 0.2499 
# # F-statistic: 9.328 on 2 and 48 DF,  p-value: 0.0003781 
# 
# lm10 <- lm(BCF ~ molmass_gmol + bodyweight, data = mlrfrog)
# summary(lm10)
# # Call:
# #   lm(formula = BCF ~ molmass_gmol + bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.08750 -0.04180 -0.01253  0.03173  0.16777 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -0.0957405  0.0427233  -2.241   0.0297 *  
# #   molmass_gmol  0.0005015  0.0001041   4.816 1.51e-05 ***
# #   bodyweight    0.0061266  0.0109493   0.560   0.5784    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05871 on 48 degrees of freedom
# # Multiple R-squared: 0.3258,  Adjusted R-squared: 0.2977 
# # F-statistic:  11.6 on 2 and 48 DF,  p-value: 7.776e-05
# 
# ##################################3 variables
# lm11 <- lm(BCF ~ logKow + Solat20C_mgL + molmass_gmol, data = mlrfrog) 
# summary(lm11)
# # Call:
# #   lm(formula = BCF ~ logKow + Solat20C_mgL + molmass_gmol, data = mlrfrog)
# # 
# # Residuals:
# #   Min        1Q    Median        3Q       Max 
# # -0.106425 -0.026123 -0.003145  0.021648  0.157794 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   3.198e-02  4.229e-02   0.756    0.453    
# # logKow       -1.877e-02  8.335e-03  -2.252    0.029 *  
# #   Solat20C_mgL -3.357e-04  7.192e-05  -4.668 2.56e-05 ***
# #   molmass_gmol  4.470e-04  9.080e-05   4.924 1.09e-05 ***
# #   ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.0492 on 47 degrees of freedom
# # Multiple R-squared: 0.5364,  Adjusted R-squared: 0.5068 
# # F-statistic: 18.13 on 3 and 47 DF,  p-value: 5.91e-08 
# 
# lm12 <- lm(BCF ~ logKow + Solat20C_mgL + bodyweight, data = mlrfrog)
# summary(lm12)
# # Call:
# #   lm(formula = BCF ~ logKow + Solat20C_mgL + bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.07905 -0.03922 -0.01034  0.01650  0.20417 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   1.397e-01  4.768e-02   2.930  0.00522 ** 
# #   logKow       -1.164e-02  1.010e-02  -1.153  0.25490    
# # Solat20C_mgL -3.793e-04  8.815e-05  -4.302 8.49e-05 ***
# #   bodyweight    4.537e-03  1.127e-02   0.403  0.68905    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.06047 on 47 degrees of freedom
# # Multiple R-squared: 0.2997,  Adjusted R-squared: 0.255 
# # F-statistic: 6.704 on 3 and 47 DF,  p-value: 0.0007363 
# 
# lm13 <- lm(BCF ~ logKow + molmass_gmol + bodyweight, data = mlrfrog)
# summary(lm13)
# # Call:
# #   lm(formula = BCF ~ logKow + molmass_gmol + bodyweight, data = mlrfrog)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.08755 -0.04141 -0.01218  0.03159  0.16769 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -0.0950581  0.0487455  -1.950   0.0571 .  
# # logKow       -0.0002659  0.0088157  -0.030   0.9761    
# # molmass_gmol  0.0005024  0.0001095   4.588 3.33e-05 ***
# #   bodyweight    0.0061476  0.0110870   0.554   0.5819    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05933 on 47 degrees of freedom
# # Multiple R-squared: 0.3258,  Adjusted R-squared: 0.2828 
# # F-statistic: 7.572 on 3 and 47 DF,  p-value: 0.0003125 
# 
# lm14 <- lm(BCF ~ Solat20C_mgL + molmass_gmol + bodyweight, data = mlrfrog)
# summary(lm14)
# # Call:
# #   lm(formula = BCF ~ Solat20C_mgL + molmass_gmol + bodyweight, 
# #      data = mlrfrog)
# # 
# # Residuals:
# #   Min        1Q    Median        3Q       Max 
# # -0.095441 -0.030926 -0.006074  0.022232  0.157862 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  -5.165e-02  3.909e-02  -1.321 0.192780    
# # Solat20C_mgL -2.593e-04  6.588e-05  -3.935 0.000273 ***
# #   molmass_gmol  4.167e-04  9.378e-05   4.443 5.37e-05 ***
# #   bodyweight    7.473e-03  9.603e-03   0.778 0.440339    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.05145 on 47 degrees of freedom
# # Multiple R-squared: 0.4929,  Adjusted R-squared: 0.4605 
# # F-statistic: 15.23 on 3 and 47 DF,  p-value: 4.673e-07 
# 
# #####################################4 variables
# lm15 <- lm(BCF ~ logKow + Solat20C_mgL + molmass_gmol + bodyweight, data = mlrfrog) 
# summary(lm15)
# # Call:
# #   lm(formula = BCF ~ logKow + Solat20C_mgL + molmass_gmol + bodyweight, 
# #      data = mlrfrog)
# # 
# # Residuals:
# #   Min        1Q    Median        3Q       Max 
# # -0.101452 -0.026480 -0.007485  0.024500  0.148969 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)   1.250e-02  4.635e-02   0.270   0.7886    
# # logKow       -1.955e-02  8.365e-03  -2.337   0.0238 *  
# #   Solat20C_mgL -3.413e-04  7.209e-05  -4.735 2.13e-05 ***
# #   molmass_gmol  4.570e-04  9.127e-05   5.007 8.58e-06 ***
# #   bodyweight    9.445e-03  9.215e-03   1.025   0.3108    
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# # 
# # Residual standard error: 0.04917 on 46 degrees of freedom
# # Multiple R-squared: 0.5467,  Adjusted R-squared: 0.5073 
# # F-statistic: 13.87 on 4 and 46 DF,  p-value: 1.693e-07 
# 
# #log-likelihoods
# ll_lm1<-logLik(lm1)
# ll_lm2<-logLik(lm2)
# ll_lm3<-logLik(lm3)
# ll_lm4<-logLik(lm4)
# ll_lm5<-logLik(lm5)
# ll_lm6<-logLik(lm6)
# ll_lm7<-logLik(lm7)
# ll_lm8<-logLik(lm8)
# ll_lm9<-logLik(lm9)
# ll_lm10<-logLik(lm10)
# ll_lm11<-logLik(lm11)
# ll_lm12<-logLik(lm12)
# ll_lm13<-logLik(lm13)
# ll_lm14<-logLik(lm14)
# ll_lm15<-logLik(lm15)
# 
# #aic
# aic1<-AIC(ll_lm1)
# aic2<-AIC(ll_lm2)
# aic3<-AIC(ll_lm3)
# aic4<-AIC(ll_lm4)
# aic5<-AIC(ll_lm5)
# aic6<-AIC(ll_lm6)
# aic7<-AIC(ll_lm7)
# aic8<-AIC(ll_lm8)
# aic9<-AIC(ll_lm9)
# aic10<-AIC(ll_lm10)
# aic11<-AIC(ll_lm11)
# aic12<-AIC(ll_lm12)
# aic13<-AIC(ll_lm13)
# aic14<-AIC(ll_lm14)
# aic15<-AIC(ll_lm15)
# 
# aic1
# aic2
# aic3
# aic4
# aic5
# aic6
# aic7
# aic8
# aic9
# aic10
# aic11
# aic12
# aic13
# aic14
# aic15
# 
# #model probabilities
# best.aic <- min(aic1,aic2,aic3,aic4,aic5,aic6,aic7,aic8,aic9,aic10,aic11,aic12,aic13,aic14,aic15)
# mp1 <- exp( -0.5 * (aic1-best.aic))
# mp2 <- exp( -0.5 * (aic2-best.aic))
# mp3 <- exp( -0.5 * (aic3-best.aic))
# mp4 <- exp( -0.5 * (aic4-best.aic))
# mp5 <- exp( -0.5 * (aic5-best.aic))
# mp6 <- exp( -0.5 * (aic6-best.aic))
# mp7 <- exp( -0.5 * (aic7-best.aic))
# mp8 <- exp( -0.5 * (aic8-best.aic))
# mp9 <- exp( -0.5 * (aic9-best.aic))
# mp10 <- exp( -0.5 * (aic10-best.aic))
# mp11 <- exp( -0.5 * (aic11-best.aic))
# mp12 <- exp( -0.5 * (aic12-best.aic))
# mp13 <- exp( -0.5 * (aic13-best.aic))
# mp14 <- exp( -0.5 * (aic14-best.aic))
# mp15 <- exp( -0.5 * (aic15-best.aic))
# 
# mp1
# #[1] 4.195434e-08
# mp2
# #[1] 9.298554e-05
# mp3
# #[1] 0.000446225
# mp4
# #[1] 2.270836e-08
# mp5
# #[1] 6.729594e-05
# mp6
# #[1] 0.000164159
# mp7
# #[1] 1.543797e-08
# mp8
# #[1] 0.1993579
# mp9
# #[1] 3.609373e-05
# mp10
# #[1] 0.0001937578
# mp11
# #[1] 1
# mp12
# #[1] 2.702879e-05
# mp13
# #[1] 7.131468e-05
# mp14
# #[1] 0.1016554
# mp15
# #[1] 0.6542875
# 
# plot(lm2)
# 
# 
