library(reshape2)
library(ggplot2)
library(gridExtra)

R2_BCF_filter=as.data.frame(subset(R2_logBCF, select=c(R2,N_var,BIC)))
R2_BCF_melt <- melt(R2_BCF_filter, id.vars = c("R2","N_var","BIC"))
p1=ggplot(R2_BCF_melt, aes(x=factor(N_var), y = R2)) + geom_boxplot() + ggtitle('logBCF')+xlab("Number of variables")
p4=ggplot(R2_BCF_melt, aes(x=factor(N_var), y = BIC)) + geom_boxplot() + ggtitle('logBCF')+xlab("Number of variables")


R2_APP_filter=as.data.frame(subset(R2_logAPP, select=c(R2,N_var,BIC)))
R2_APP_melt <- melt(R2_APP_filter, id.vars = c("R2","N_var","BIC"))
p2=ggplot(R2_APP_melt, aes(x=factor(N_var), y = R2)) + geom_boxplot() + ggtitle('logAPP')+xlab("Number of variables")
p5=ggplot(R2_APP_melt, aes(x=factor(N_var), y = BIC)) + geom_boxplot() + ggtitle('logAPP')+xlab("Number of variables")


R2_TC_filter=as.data.frame(subset(R2_logTC, select=c(R2,N_var,BIC)))
R2_TC_melt <- melt(R2_TC_filter, id.vars = c("R2","N_var","BIC"))
p3=ggplot(R2_TC_melt, aes(x=factor(N_var), y = R2)) + geom_boxplot() + ggtitle('logTC')+xlab("Number of variables")
p6=ggplot(R2_TC_melt, aes(x=factor(N_var), y = BIC)) + geom_boxplot() + ggtitle('logTC')+xlab("Number of variables")


tiff(file = "FigR2.tiff",  width = 10000, height = 8000, units = "px", res = 800, compression = "lzw") #
grid.arrange(p1,p2,p3, ncol=2, nrow=2)
dev.off()

tiff(file = "FigBIC.tiff",  width = 10000, height = 8000, units = "px", res = 800, compression = "lzw") #
grid.arrange(p4,p5,p6, ncol=2, nrow=2)
dev.off()