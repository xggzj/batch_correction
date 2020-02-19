# 1. Preprocessing
# Read infiles
# Read expression data (.fcs):
readInfiles = function(path){
  library(flowCore)
  fileNames = dir(path)
  filePath = sapply(fileNames, function(x){paste(path,x,sep='/')})
  files = lapply(filePath, function(x){read.FCS(x)}) 
  return(files)
}
files = readInfiles("~/master/infiles")

# Read labels (.csv):
readLabels = function(path){
  fileNames = dir(path)
  filePath = sapply(fileNames, function(x){paste(path,x,sep='/')})
  files = lapply(filePath, function(x){read.csv(x)}) 
  return(files)
}
labels = readLabels("~/master/labels")

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
# Bind all batches into a single expression matrix, Dat
datBind = function(dat){
  Dat = dat[[1]]
  n = length(dat)
  for (i in 1:(n-1)){
    Dat = rbind(Dat, dat[[i+1]])
  }
  return(Dat)
}
Dat = datBind(dat)

# Meta_data: a dataframe with variables to integrate including batch numbers and cell labels
Meta_data = data.frame(batch = factor(c(rep(1, nrow(labels[[1]])), rep(2, nrow(labels[[2]]))), 
                       label = factor(c(as.character(labels[[1]]$level1), as.character(labels[[2]]$level1)))))

Dat.harmony = HarmonyMatrix(Dat, Meta_data, "batch")

# 3. Visualization
library(ggplot2)
