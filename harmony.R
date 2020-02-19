# 1. Preprocessing
# Read infiles
readInfiles = function(path, filetype = c("fcs", "csv")){
  fileNames = dir(path)
  filePath = sapply(fileNames, function(x){paste(path,x,sep='/')})
  if (filetype == "fcs"){
    library(flowCore)
    files = lapply(filePath, function(x){read.FCS(x)})
  }
  if (filetype == "csv"){
    files = lapply(filePath, function(x){read.csv(x)})
  }
  return(files)
}
# Read expression data (.fcs):
files = readInfiles("~/master/infiles", "fcs")
# Read cell type data (.csv):
labels = readInfiles("~/master/labels", "csv")

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
bindlist = function(data){
  Data = data[[1]]
  n = length(data)
  for (i in 1:(n-1)){
    Data = rbind(Data, data[[i+1]])
  }
  return(Data)
}
# Bind the expression data of all batches into a single matrix
Dat = bindlist(dat)
# Bind the cell type labels of all batches into a single dataframe
Label = bindlist(labels)

# Input the batch number:
Label$batch = c(rep(1,nrow(dat[[1]])), rep(2,nrow(dat[[2]])))

# MetaData: a dataframe with variables to integrate including batch numbers and cell type labels
MetaData = data.frame(batch = Label$batch, label = Label$level1)

library(harmony)
Dat.harmony = HarmonyMatrix(Dat, MetaData, "batch")

# Save the batch-corrected files separately
saveCorrectedFiles = function(Dat.harmony){
  n = length(files)
  r = c()
  for (i in 1:n){r[i] = nrow(files[[i]])}
  
  dat.harmony = list()
  r0 = 0
  for (i in 1:n){
    dat.harmony[[i]] = Dat.harmony[(r0+1):(r0+r[i]),]
    r0 = r0+r[i]
  }
  names(dat.harmony) = names(files)
  
  for (i in 1:n){
    filename = gsub("fcs", "csv", names(dat.harmony[i]))
    write.csv(dat.harmony[[i]], filename, row.names = FALSE) 
  }

  return(dat.harmony)
}

dat.harmony = saveCorrectedFiles(Dat.harmony)

# 3. Visualization 
# UMAP
library(umap)
umap.harmony = umap(Dat.harmony)
umap.dat = umap(Dat)

# Plot
library(ggplot2)
df = data.frame(UMAP1 = umap.dat$layout[,1],
                UMAP2 = umap.dat$layout[,2],
                batch = Label$batch)
gp = ggplot(df, aes(UMAP1, UMAP2, color = ))
