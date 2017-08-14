library(nortest)
dim(allfrog)
colnames(allfrog)

unique.species <- unique(allfrog$Species)
unique.species = unique.species[which (unique.species!="Mole salamander")]
# unique.chemicals <- unique(allfrog$Chemical)

unique.chemicals=c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                       'Total Fipronil','Total Triadimefon')

# bodyweight across 8 species
for(i in unique.species)
{
  z <- subset(allfrog,Species==i&Application=="Soil")
  
  z_index=is.element(z$Chemical,c('Imidacloprid','Pendimethalin','Total Atrazine',
                                                          'Total Fipronil','Total Triadimefon'))
  z = z[z_index,] #select only five chemicals
#   z = z[which(z$Chemical==unique.chemicals),]
  x = z[,12]
  x2 <- x[!is.na(x)]
  x3 <- x[!is.na(x)]
  s = paste(i, "Bodyweight", sep = " ")
  print(s)
  print(ad.test(x2))
  print(length(x2))
  
  for (j in 1:length(x3))
  {
    if (x3[j] == 0)
    {
      x3[j] = 1
    }
  }
  logx = log(x3)
  print(ad.test(logx))
}

for(i in unique.chemicals)
{

  z <- subset(allfrog, Application=="Soil")
  z = z[which(z$Chemical==i),]

  #   z = z[which(z$Chemical==unique.chemicals),]
  x = z[,18]
  x2 <- x[!is.na(x)]
  x3 <- x[!is.na(x)]
  s = paste(i, "AppFactor", sep = " ")
  print(s)
  print(ad.test(x2))
  print(length(x2))
  
  for (j in 1:length(x3))
  {
    if (x3[j] == 0)
    {
      x3[j] = 1
    }
  }
  logx = log(x3)
  print(ad.test(logx))
}


# body burdens across 8 species and 5 pesticides
for(i in unique.species)
{
  z <- subset(allfrog,Species==i)
  for(j in unique.chemicals)
  {
    zz = subset(z,Chemical==j)
    x = zz[,8]
    x2 <- x[!is.na(x)]
    x3 <- x[!is.na(x)]
    if (length(x2)>7)
    {
      s = paste(i,j,"tissue concentration", sep=" ")
      print(s)
      print(ad.test(x2))
  
      for (j in 1:length(x2))
      {
        if (x2[j] == 0)
        {
          x2[j] = 1
        }
      }
      
      logx = log(x2)
      print(ad.test(logx))
    }
  }
}
    
# soil concentrations across chemicals
for(i in unique.chemicals)
{
  z <- subset(allfrog,Chemical==i)
  x = unique(z$SoilConc)
  logx = log(x)
  if (length(x)>7)
  {
    print(i)
    print(ad.test(x))
    print(ad.test(logx))
  }
}