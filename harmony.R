# 1. Preprocessing
# Read infiles
readInfiles = function(path){
  library(flowCore)
  fileNames = dir(path)
  filePath = sapply(fileNames, function(x){paste(path,x,sep='/')})
  files = lapply(filePath, function(x){read.FCS(x)}) 
  return(files)
}
files = readInfiles("~/master/infiles")

# Read marker names and remove empty channels
channelFilter = function(file){
  library(dplyr)
  dat = file@exprs[,-1] 
  colnames(dat) = markernames(file)
  nonEmptyChannels = file@parameters@data[c("name", "desc")] %>% filter(name != desc)
  dat = dat[,colnames(dat) %in% nonEmptyChannels$desc]
  return(dat)
}
dat = lapply(files, channelFilter)

# Remove outliers 
removeOutliers = function(dat){
  q1 = quantile(dat, 0.01)
  q99 = quantile(dat, 0.99)
  dat[dat<q1] = q1
  dat[dat>q99] = q99
  return(dat)
}
dat = lapply(dat, function(x){apply(x, 2, removeOutliers)})

# Scaling
dat = lapply(dat, scale)

# 2. Batch correction
library(harmony)
batch = as.factor(c(rep(1, nrow(dat[[1]])), rep(2, nrow(dat[[2]])))) 
dat.harmony = HarmonyMatrix(dat, batch)

# 3. Visualization
library(ggplot2)
