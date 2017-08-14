library(package="ellipse")

#http://rais.ornl.gov/cgi-bin/tools/TOX_search
pchem_root <- path.expand("~/Dropbox/terr_models/pchem_properties/")
pchem_data <- read.csv(paste(pchem_root,"rais_pchem_properties.csv", sep=""))

dim(pchem_data)
View(pchem_data)
summary(pchem_data)
colnames(pchem_data)
pchem_corr <- pchem_data[,c(5,7,8,9,11,13,14,15,18,19,21,22)]
colnames(pchem_corr)
#[1] "Soil.to.Dry.Plant.Uptake.."                    "Density..g.cm3..."                            
#[3] "Diffusivity.in.Air..cm2.s..."                  "Diffusivity.in.Water..cm2.s..."               
#[5] "Unitless.Henry.s.Law.Constant.."               "Organic.Carbon.Partition.Coefficient..L.kg..."
#[7] "Skin.Permeability.Constant..cm.hr..."          "Log.of.Octanol.Water.Partition.Coefficient."  
#[9] "Molecular.Weight..g.mol..."                    "RAGS.Part.E.Dermal.Absorption.Factor."        
#[11] "Water.Solubility..mg.L..."                     "Vapor.Pressure..mm.Hg.."  
colnames(pchem_corr) <- c("SoPl","Dens","DiA","DiW","HLC","Koc","SPC","LKow","MW","DAF","Sol","VP")
corr.pchem_data <- cor(pchem_corr,use="na.or.complete")
ord <- order(corr.pchem_data[1,])
xc <- corr.pchem_data[ord, ord]
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
#View(xc)



panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * abs(r))
}


par(mfrow=c(2,1))
plotcorr(xc, col=colors[5*xc + 6])
pairs(xc, lower.panel=panel.smooth, upper.panel=panel.cor)
