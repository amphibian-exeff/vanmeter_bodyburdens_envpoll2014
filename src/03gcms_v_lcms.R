dim(allfrog)
colnames(allfrog)

unique.species <- unique(allfrog$Species)
length(unique.species)
unique.chemicals <- unique(allfrog$Chemical)
length(unique.chemicals)
unique.instruments <- unique(allfrog$Instrument)
length(unique.instruments)
unique.applications <- unique(allfrog$Application)
length(unique.applications)

col.burdens <- which(colnames(allfrog)=="TissueConc")

pdf(paste(frog_out,"hist_ttest_gcms_lcms.pdf",sep=""),width=11,height=8)
  par(mfrow=c(3,2))
  for(species in unique.species){
    for(chemical in unique.chemicals){
      for(application in unique.applications){
        these.gcms <- which(allfrog$Species == species & allfrog$Chemical == chemical & 
                allfrog$Application == "Soil" & allfrog$Application == application & allfrog$Instrument == "GCMS")
        these.lcms <- which(allfrog$Species == species & allfrog$Chemical == chemical & 
                              allfrog$Application == "Soil" & allfrog$Application == application & allfrog$Instrument == "LCMS")
        if(length(these.gcms)>0 & length(these.lcms)>0){
          burdens.gcms <- allfrog[these.gcms,col.burdens]
          burdens.lcms <- allfrog[these.lcms,col.burdens]
          this.t.test <- t.test(burdens.gcms,burdens.lcms)
          burdens.min <- min(c(burdens.gcms,burdens.lcms))
          burdens.max <- max(c(burdens.gcms,burdens.lcms))
          hist(burdens.gcms,xlim=c(burdens.min,burdens.max),col="green")
          hist(burdens.lcms,xlim=c(burdens.min,burdens.max),col="blue",main=paste("LCMS;",chemical,"; ", application, ";\n",species,"; p-value= ",this.t.test[3]))     
        }
      }
    }
  }
dev.off()
