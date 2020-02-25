library(harmony)

bindlist = function(data){
  Data = data[[1]]
  n = length(data)
  for (i in 1:(n-1)){
    Data = rbind(Data, data[[i+1]])
  }
  return(Data)
}

saveCorrectedFiles = function(Dat.harmony, outpath){
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
    filename = names(dat.harmony[i])
    filepath = sapply(filename, function(x){paste(outpath, x, sep = '/')})
    write.FCS(flowFrame(dat.harmony[[i]]), filepath) 
  }
  
  return(dat.harmony)
}
