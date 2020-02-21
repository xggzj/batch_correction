library(flowCore)
library(dplyr)

readInfiles = function(path, filetype = c("fcs", "csv")){
  fileNames = dir(path)
  filePath = sapply(fileNames, function(x){paste(path,x,sep='/')})
  if (filetype == "fcs"){    
    files = lapply(filePath, function(x){read.FCS(x)})
  }
  if (filetype == "csv"){
    files = lapply(filePath, function(x){read.csv(x)})
  }
  return(files)
}

channelFilter = function(file){
  dat = file@exprs[,-1] 
  colnames(dat) = markernames(file)
  nonEmptyChannels = file@parameters@data[c("name", "desc")] %>% filter(name != desc)
  dat = dat[,colnames(dat) %in% nonEmptyChannels$desc]
  return(dat)
}

removeOutliers = function(dat){
  q1 = quantile(dat, 0.01)
  q99 = quantile(dat, 0.99)
  dat[dat<q1] = q1
  dat[dat>q99] = q99
  return(dat)
}
