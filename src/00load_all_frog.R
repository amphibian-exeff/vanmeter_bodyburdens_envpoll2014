#library(gplots)
#library(corrplot)

#http://plantecology.syr.edu/fridley/bio793/lm.html
#http://cw.routledge.com/textbooks/9780805861853/r/ch13.pdf

#surface area of aquarium is 1235 cm^2 = 25*49

.Platform
Sys.info()[4]

#data directory
#mac/unix
if(.Platform$OS.type=="unix"){
  if(Sys.info()[4]=="stp-air.local"){#Tom Mac
    frog_root<-path.expand("~/Dropbox/amphib_dermalexposure/DATA/good_data/") 
  }
}
#windows
if(.Platform$OS.type=="windows"){
  if(Sys.info()[4]=="DC2626UTPURUCKE"){#Tom's work pc
    frog_root<-"C:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\"
  }
  if(Sys.info()[4]=="DC2626URVANMETE"){#Robin's work pc
    frog_root<-"C:\\Dropbox\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\" 
  }
  if(Sys.info()[4]=="DC2626UMCYTERSK"){#Mike's work pc
    frog_root<-"C:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\" 
  }
  if(Sys.info()[4]=="TH-THINK"){#Tao's work pc
    frog_root<-"D:\\Dropbox\\amphib_dermalexposure\\DATA\\good_data\\" 
  }
}


#subdirectories
frog_src <- paste(frog_root,"src/",sep="")
frog_out <- paste(frog_root,"output/",sep="")

#read data
allfrog<-read.table(paste(frog_root,"good_data.csv",sep=""), header = TRUE, sep = ",") 
dim(allfrog)
colnames(allfrog)
unique(allfrog$good)
goodrows <- which(allfrog$good==1)
allfrog <- allfrog[goodrows,]
dim(allfrog)