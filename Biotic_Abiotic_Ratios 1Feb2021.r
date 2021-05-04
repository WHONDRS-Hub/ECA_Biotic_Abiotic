# use transformation database with each transformation assigned to be biotic or abiotic or unknown
# transformations will be parsed into biotic and abiotic groups
# then the relative contributions of those two groups will be calculated for each sample and ratio will, in turn, be calculated
# each sample should have 3 metrics: number of biotic transformations, number of abiotic transformations, and the ratio of those two
# those outcomes will then be combined with metadata, including spatial position and the basic climate variables previously extracted

#### start function to turn factors to character

fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}

# set the working directory
setwd("//PNL/Projects/ECA_Project/Biotic_Abiotic")

# read in the transformation database, take a look at the first few rows, and check the variable types
trans = fac.to.char.fun(read.csv("Biotic-abiotic-transfromation-classification.csv")); head(trans); str(trans)

unique(trans$Biotic.abiotic)

trans = trans [-which(is.na(trans$Biotic.abiotic)),]; head(trans)
trans$Biotic.abiotic = gsub(pattern = "Bitoic/abiotic",replacement = "Both",x = trans$Biotic.abiotic)

trans$Biotic.abiotic = gsub(pattern = "\\*", "", trans$Biotic.abiotic)
trans$Biotic.abiotic = gsub(pattern = "Biotic based on Kegg but we didn't use for the biotic/abiotic paper",replacement = "Biotic",x = trans$Biotic.abiotic)
unique(trans$Biotic.abiotic)

#read in transformation profiles
prof.tran = fac.to.char.fun(read.csv("S19S_Sed-Water_08.12_Trans_Profiles.csv")); prof.tran[1:5,1:5]; str(prof.tran)

#new matrix, sample name, abundance of transformation type, ratio
trans.comp = numeric()
for (samp.name in colnames(prof.tran)[grep(pattern = 'Sample',x = colnames(prof.tran))]) { # compute the abundances
  
  abiotic.abund = sum(prof.tran[which(prof.tran$Name %in% trans$Name[grep(pattern = 'Abiotic',x = trans$Biotic.abiotic)]  ),which(colnames(prof.tran) == samp.name)])
    
  biotic.abund = sum(prof.tran[which(prof.tran$Name %in% trans$Name[grep(pattern = 'Biotic',x = trans$Biotic.abiotic)]  ),which(colnames(prof.tran) == samp.name)])
    
  both.abund = sum(prof.tran[which(prof.tran$Name %in% trans$Name[grep(pattern = 'Both',x = trans$Biotic.abiotic)]  ),which(colnames(prof.tran) == samp.name)])

  trans.comp = rbind(trans.comp,c(samp.name,abiotic.abund,biotic.abund,both.abund))
  
}

# do some formatting
colnames(trans.comp) = c("Sample_ID","Abiotic.abund","Biotic.abund","Both.abund")
trans.comp = fac.to.char.fun(as.data.frame(trans.comp)); 
trans.comp$Abiotic.abund = as.numeric(trans.comp$Abiotic.abund)
trans.comp$Biotic.abund = as.numeric(trans.comp$Biotic.abund)
trans.comp$Both.abund = as.numeric(trans.comp$Both.abund)
trans.comp$Total.trans = trans.comp$Abiotic.abund + trans.comp$Biotic.abund + trans.comp$Both.abund
str(trans.comp)

#cal ratios

trans.comp$Abiotic.to.Biotic = trans.comp$Abiotic.abund/trans.comp$Biotic.abund
trans.comp$Abiotic.to.Both = trans.comp$Abiotic.abund/trans.comp$Both.abund
trans.comp$Biotic.to.Both = trans.comp$Biotic.abund/trans.comp$Both.abund

head(trans.comp)

# write out the file
write.csv(x = trans.comp,file = "trans_variables.csv", quote = F, row.names = F)


#analysis
pdf("abiotic_to_biotic_hist.pdf")
  hist(trans.comp$Abiotic.to.Biotic, breaks = 50)
dev.off()
#hist(log10(trans.comp$Abiotic.to.Biotic))
#hist(trans.comp$Abiotic.to.Both, breaks = 50)
#hist(trans.comp$Biotic.to.Both, breaks = 50)

### look at correlations with Latitude and Longitude

meta.data = fac.to.char.fun(read.csv("//PNL/Projects/ECA_Project/Biotic_Abiotic/metadata/WHONDRS_S19S_Metadata_v2.csv", header=F))
meta.data = meta.data[-1,]
colnames(meta.data) = meta.data[1,]
meta.data = meta.data[-1,]
colnames(meta.data) = gsub (" ", "_", colnames(meta.data))
str(meta.data)

#colnames(meta.data) = gsub (" ", "_", colnames(meta.data))

climate.all = read.csv("//PNL/Projects/ECA_Project/Biotic_Abiotic/metadata/WHONDRS_ClimVariables.csv")
metadata.all = merge(meta.data, climate.all, by = "Sample_ID")

relevant.columns = c("Sample_ID","Stream_Name","City","State_or_Province", "Country","US_Longitude_dec.deg","US_Latitude_dec.deg" ,"PET_mm_per_yr","AET_mm_per_yr","MAT_degrees_C", "MAP_mm_per_yr", "TAR_degrees_C","Distance_m")

trans.comp$Lat_dec.deg=-999
trans.comp$Long_dec.deg=-999
trans.comp$Sample.Set = -999
trans.comp$Sample.State = -999

#could create more for data we want to pull into the data frame
trans.comp$PET_mm_per_yr = -999
trans.comp$AET_mm_per_yr = -999
trans.comp$MAT_degrees_C = -999
trans.comp$MAP_mm_per_y = -999
trans.comp$TAR_degrees_C = -999

for (i in 1:nrow(metadata.all)) {
  trans.comp$Lat_dec.deg[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]=  as.numeric(metadata.all$`Latitude_of_upstream_site_(decimal_degrees)`[i])
  trans.comp$Long_dec.deg[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$`Longitude_of_upstream_site_(decimal_degrees)`[i])
  trans.comp$PET_mm_per_yr[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$PET_mm_per_year[i])
  trans.comp$AET_mm_per_yr[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$AET_mm_per_year[i])
  trans.comp$MAT_degrees_C[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$MAT_degrees_C[i])
  trans.comp$TAR_degrees_C[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$TAR_degree_C[i])
  trans.comp$MAP_mm_per_y[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$MAP_mm_per_year[i])
  
   if (metadata.all$`Sampling_location:_Country`[i] == 'USA' & metadata.all$`Sampling_location:_State/Province`[i] != 'Alaska' & metadata.all$`Sampling_location:_State/Province`[i] != 'Puerto Rico') {
    
    trans.comp$Sample.Set[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)] = "CONUS"
    trans.comp$Sample.State[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)] = metadata.all$`Sampling_location:_State/Province`[i]  
  } 
}
#filling columns in above for loop #leave grep as is
unique(trans.comp$Sample.State)
trans.comp[grep(pattern = "Florida", x = trans.comp$Sample.State),]

trans.comp.sed=trans.comp[grep(pattern = "Sed" , x= trans.comp$Sample_ID),]
trans.comp.sw=trans.comp[-grep(pattern = "Sed" , x= trans.comp$Sample_ID),]

sarah.plot.fun = function(data.in = trans.comp.sed,grep.pattern = "CONUS",y.var = "Abiotic.abund",x.var = "Long_dec.deg") {
  data.temp = data.in[grep(pattern = grep.pattern,x = data.in[,'Sample.Set']),]
  if (length(grep(pattern = 'Sed',x = data.temp$Sample_ID)) > 0) {
    main.temp = paste(grep.pattern,"Sediment",sep=" ")
    
  } else {
    
    main.temp = paste(grep.pattern,"Surface Water",sep=" ")
  
  }
  
  mod.to.plot = data.temp[,y.var] ~ data.temp[,x.var]
  mod.lm = summary(lm(mod.to.plot))
  plot(mod.to.plot,xlab=x.var,ylab=y.var,main = main.temp)
  p.val = round(mod.lm$coefficients[2,4],digits = 4)
  r.sq = round(mod.lm$r.squared,digits = 3)
  mtext(text = paste("p = ",p.val," ",sep=""),line = -2,adj = 1,side = 3)
  mtext(text = paste("R.sq = ",r.sq," ",sep=""),line = -3,adj = 1,side = 3)
  abline(mod.lm)
  
  return(mod.lm)
  
}


for (y.var.use in c("Abiotic.abund", "Biotic.abund", "Abiotic.to.Biotic","Total.trans")) {

  for (x.var.use in c("Lat_dec.deg", "Long_dec.deg","PET_mm_per_yr","AET_mm_per_yr","MAT_degrees_C","TAR_degrees_C","MAP_mm_per_y")) {

    pdf(paste(y.var.use,"_vs_",x.var.use,".pdf",sep=""))
    
    sarah.plot.fun(data.in = trans.comp.sed,grep.pattern = "CONUS", y.var = y.var.use, x.var = x.var.use )
    sarah.plot.fun(data.in = trans.comp.sw,grep.pattern = "CONUS", y.var = y.var.use, x.var = x.var.use )
    dev.off()##closes graphic device in R and saves files in directory
  }
  
}


##look at number of peaks
##manuscript figures

######
### Figure 1
######

# define density objects
biotic.sed.den = density(trans.comp.sed$Biotic.abund,from = 0)
biotic.sw.den = density(trans.comp.sw$Biotic.abund,from = 0)

abiotic.sed.den = density(trans.comp.sed$Abiotic.abund,from = 0)
abiotic.sw.den = density(trans.comp.sw$Abiotic.abund,from = 0)

abio.to.bio.sed.den = density(trans.comp.sed$Abiotic.to.Biotic)
abio.to.bio.sw.den = density(trans.comp.sw$Abiotic.to.Biotic)


# define functions
x.range.fun = function(sed.den,sw.den) {
  
  x.min = min(c(sed.den$x,sw.den$x))
  x.max = max(c(sed.den$x,sw.den$x))
  return(c(x.min,x.max))
  
}

y.range.fun = function(sed.den,sw.den) {
  
  y.min = 0
  y.max = max(c(sed.den$y,sw.den$y))
  return(c(y.min,y.max))
  
}

# do the plotting

pdf("Fig1_Density_Plots.pdf",height = 15)

  par(mfrow = c(3,1),pty="s") # 3 rows, 1 column from mfrow; pty = "s" makes the panels square
 
  plot(biotic.sed.den,xlim = x.range.fun(biotic.sed.den,biotic.sw.den),ylim=y.range.fun(biotic.sed.den,biotic.sw.den),main="",xlab="Number of Biotic Transformations",cex.lab=2,cex.axis=1.5,lwd=2,col="orange")
  points(biotic.sw.den,typ="l",lwd=2,col="blue")
  mtext(text = " A",line = -2,adj = 0,side = 3,cex = 1.5)
  legend(x = 42000,y = 6.5e-05,legend = c("Surface Water","Sediments"),col = c("blue","orange"),lty = 1,lwd = 2,bty = "n",cex=1.5)
  
  plot(abiotic.sed.den,xlim = x.range.fun(abiotic.sed.den,abiotic.sw.den),ylim=y.range.fun(abiotic.sed.den,abiotic.sw.den),main="",xlab="Number of Abiotic Transformations",cex.lab=2,cex.axis=1.5,lwd=2,col="orange")
  points(abiotic.sw.den,typ="l",lwd=2,col="blue")
  mtext(text = " B",line = -2,adj = 0,side = 3,cex = 1.5)

  plot(abio.to.bio.sed.den,xlim = x.range.fun(abio.to.bio.sed.den,abio.to.bio.sw.den),ylim=y.range.fun(abio.to.bio.sed.den,abio.to.bio.sw.den),main="",xlab="Ratio of Abiotic to Biotic Transformations",cex.lab=2,cex.axis=1.5,lwd=2,col="orange")
  points(abio.to.bio.sw.den,typ="l",lwd=2,col="blue")
  mtext(text = " C",line = -2,adj = 0,side = 3,cex = 1.5)

dev.off()

######
### Figure 2
######

# abiotic from water vs abiotic from sediment (each field site has one value for abiotic water and one value for abiotic sediment)
# repeat for biotic and the ratio. So 3 panels total.
# add regression lines and 1 to 1 lines
# add letters to panels
# if we have 1 to 1 and regression on each panel, then add legend differentiating those two

# step 1, collapse data into site level mean values for sw and sed, for each variable

trans.site = numeric()

head(trans.comp.sed)
unique.site = unique(sub(pattern = "_Sed.*",replacement = "" ,x=trans.comp.sed$Sample_ID ))

for (i in unique.site) {
  
  temp.sed.abiotic = mean(trans.comp.sed$Abiotic.abund[grep(pattern = i,x = trans.comp.sed$Sample_ID)])
  temp.sed.biotic = mean(trans.comp.sed$Biotic.abund[grep(pattern = i,x = trans.comp.sed$Sample_ID)])
  temp.sed.ratio = mean(trans.comp.sed$Abiotic.to.Biotic[grep(pattern = i,x = trans.comp.sed$Sample_ID)])
  
  temp.sw.abiotic = mean(trans.comp.sw$Abiotic.abund[grep(pattern = i,x = trans.comp.sw$Sample_ID)])
  temp.sw.biotic = mean(trans.comp.sw$Biotic.abund[grep(pattern = i,x = trans.comp.sw$Sample_ID)])
  temp.sw.ratio = mean(trans.comp.sw$Abiotic.to.Biotic[grep(pattern = i,x = trans.comp.sw$Sample_ID)])

  trans.site = rbind(trans.site,c(i,temp.sed.abiotic,temp.sed.biotic,temp.sed.ratio,temp.sw.abiotic,temp.sw.biotic,temp.sw.ratio))
  
}

colnames(trans.site) = c("Site_ID","sed.abiotic","sed.biotic","sed.ratio","sw.abiotic","sw.biotic","sw.ratio") 
trans.site = fac.to.char.fun(as.data.frame(trans.site))

for (i in 2:ncol(trans.site)) {
  
  trans.site[,i] = as.numeric(trans.site[,i])
  
}

head(trans.site)
str(trans.site)

# do plotting

## next thing is to use the function from above to do the plotting
## expand to biotic and the ratio for all 3 panels

mod.to.plot = trans.site$sw.abiotic ~ trans.site$sed.abiotic
mod.lm = summary(lm(mod.to.plot))
plot(mod.to.plot,xlab="",ylab="",main = "")
p.val = round(mod.lm$coefficients[2,4],digits = 4)
r.sq = round(mod.lm$r.squared,digits = 3)
mtext(text = paste("p = ",p.val," ",sep=""),line = -2,adj = 1,side = 3)
mtext(text = paste("R.sq = ",r.sq," ",sep=""),line = -3,adj = 1,side = 3)
abline(mod.lm,lwd=2)
abline(0,1,lty=2,lwd=2)




