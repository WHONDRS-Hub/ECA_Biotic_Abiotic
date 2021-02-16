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


for (i in 1:nrow(metadata.all)) {
  trans.comp$Lat_dec.deg[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]=  as.numeric(metadata.all$`Latitude_of_upstream_site_(decimal_degrees)`[i])
  trans.comp$Long_dec.deg[grep(pattern = metadata.all$Sample_ID[i],x = trans.comp$Sample_ID)]= as.numeric(metadata.all$`Longitude_of_upstream_site_(decimal_degrees)`[i])
  
}


#check why there are 3 samples per site ID
