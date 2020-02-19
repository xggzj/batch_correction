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
# Bind the expression data of all batches into a single matrix
datBind = function(dat){
  Dat = dat[[1]]
  n = length(dat)
  for (i in 1:(n-1)){
    Dat = rbind(Dat, dat[[i+1]])
  }
  return(Dat)
}
Dat = datBind(dat)

# Bind the cell type labels of all batches into a single dataframe
labelBind = function(labels){
  Label = labels[[1]]
  n = length(labels)
  for (i in 1:(n-1)){
    Label = rbind(Label, labels[[i+1]])
  }
  return(Label)
}
Label = labelBind(labels)

# Input the batch number:
Label$batch = c(rep(1,nrow(dat[[1]])), rep(2,nrow(dat[[2]])))

# MetaData: a dataframe with variables to integrate including batch numbers and cell type labels
MetaData = data.frame(batch = Label$batch, label = Label$level1)

#Meta_data = data.frame(batch = factor(c(rep(1, nrow(labels[[1]])), rep(2, nrow(labels[[2]]))), 
#                       label = factor(c(as.character(labels[[1]]$level1), as.character(labels[[2]]$level1)))))

library(harmony)
Dat.harmony = HarmonyMatrix(Dat, MetaData, "batch")

# Save the batch-corrected files separately
saveCorrectedFiles = function(Dat.harmony){
  n = length(files)
  r = c()
  for (i in 1:n){
    r[i] = nrow(files[[i]])
  }
  
  dat.harmony = list()
  r0 = 0
  for (i in 1:n){
    dat.harmony[[i]] = Dat.harmony[(r0+1):(r0+r[i]),]
    r0 = r0+r[i]
  }
  names(dat.harmony) = names(files)
  
  for (i in 1:n){
    filename = paste('harmony', names(dat.harmony[i]))
    write.FCS(flowFrame(dat.harmony[[i]]), "") names(dat.harmony[i])
  }
  return(dat.harmony)
}

dat.harmony = saveCorrectedFiles(Dat.harmony)



flowFrame()
write.FCS()

# 3. Visualization
library(ggplot2)
