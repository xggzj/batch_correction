# Read infiles
library(flowCore)
file_1 = read.FCS("000027_d_C4.fcs") # Batch 1 raw data
file_2 = read.FCS("VG32_(tv.1)-PTc_d_C1.fcs") # Batch 2 raw data

# Preprocessing:
# Read and filter the marker names (short_name) from a flowFrame, remove vacant channels.
getsnames = function(file){
  batch = file@exprs
  colnames(batch) = file@parameters@data$desc
  batch = batch[, c('CD45',  'CD57', 'HLA-ABC', 'EQBeads', 'CD49d', 'CD19', 
                    'CD5', 'CD16', 'CD4', 'CD8a', 'CD11c', 'CD31', 
                    'CD25', 'CD64', 'CD123', 'gdTCR', 'CD3e', 'CD33',
                    'CD26', 'CD9', 'CD34', 'CD22', 'CD14', 'CD161', 
                    'CD29', 'HLA-DR', 'CD44', 'CD127', 'CD24', 'CD27', 
                    'CD38', 'CD45RA', 'CD20', 'CD7', 'IgD', 'CD56', 
                    'CD99', 'CD15', 'CD39', 'DNA-Ir191', 'DNA-Ir193', 'CD11b')]
  return(batch)
}

batch_1 = getsnames(file_1) 
batch_2 = getsnames(file_2)

# Remove outliers
removeOutliers = function(Dat){
  q1 = quantile(Dat, 0.001)
  q99 = quantile(Dat, 0.999)
  Dat[Dat<q1] = q1
  Dat[Dat>q99] = q99
  return(Dat)
}

batch_1 = apply(batch_1, 2, removeOutliers)
batch_2 = apply(batch_2, 2, removeOutliers)

# 1. Batch correction:
# For SAUCIE: 
# Save the pre-processed matrix of the two batches as .fcs files to run SAUCIE correction in the terminal
flowframe_1 = flowFrame(batch_1)
write.FCS(flowframe_1, "000027_d_C4_processed.fcs")
flowframe_2 = flowFrame(batch_2)
write.FCS(flowframe_1, "VG32_(tv.1)-PTc_d_C1_processed.fcs")

# Read infiles after SAUCIE-batch-corrected
file_3 = read.FCS("000027_d_C4_saucie.fcs") # Batch 1 SAUCIE-corrected data
file_4 = read.FCS("VG32_(tv.1)-PTc_d_C1_saucie.fcs") # Batch 2 SAUCIE-corrected data
# Concatenate two batch-corrected data into a single matrix 
dat.saucie = rbind(file_3@exprs, file_4@exprs)

# For ComBat, limma and Harmony:  
# Concatenate the raw data of two batches into a single matrix with cells in the rows and features in the columns
dat = rbind(batch_1, batch_2) # Raw data of two batches before batch-correction
batch = as.factor(c(rep(1, nrow(batch_1)), rep(2, nrow(batch_2)))) # Batch covariate

# The input for ComBat and limma should be a matrix with features in the rows and samples in the columns
# The input for Harmony should be a matrix with samples in the rows and features in the columns
# Scaling after correction 

# ComBat
library(sva)
dat.combat = ComBat(t(dat), batch)  
dat.combat = t(dat.combat)
dat.combat = scale(dat.combat)
# Limma
library(limma)
dat.limma = removeBatchEffect(t(dat), batch)
dat.limma = t(dat.limma)
dat.limma = scale(dat.limma)
# Harmony
library(harmony)
dat.harmony = HarmonyMatrix(dat, batch)
dat.harmony = scale(dat.harmony)

# 2. Visulization:
library(ggplot2)

# (1) UMAP
# Cells in rows, features in columns
library(umap)
Umap = function(Dat, Batch, method){
  dat.umap = umap(Dat)
  df = data.frame(x = dat.umap$layout[,1],
                  y = dat.umap$layout[,2],
                  batch = Batch)
  gp = ggplot(df, aes(x, y, colour = batch)) +
    geom_point(size=-0.1) +
    labs(x = "UMAP1", y = "UMAP2", title = method)
  return(gp)
}

# Before correction:
umap.dat = Umap(dat, batch, "Raw")

# After correction:
# ComBat
umap.combat = Umap(dat.combat, batch, "ComBat")
# Limma 
umap.limma = Umap(dat.limma, batch, "limma")
# Harmony
umap.harmony = Umap(dat.harmony, batch, "Harmony")
# SAUCIE
umap.saucie = Umap(dat.saucie, batch, "SAUCIE")

# (2) t-SNE
# Cells in rows, features in columns
library(Rtsne)
Tsne = function(Dat, Batch, method){
  dat.tsne = Rtsne(Dat)
  df = data.frame(x = dat.tsne$Y[,1],
                  y = dat.tsne$Y[,2],
                  batch = Batch)
  gp = ggplot(df, aes(x, y, colour = batch)) +
    geom_point(size=-0.1) +
    labs(x = "t-SNE1", y = "t-SNE2", title = method)
  return(gp)
}

# Before correction:
tsne.dat = Tsne(dat, batch, "Raw")

# After correction:
# ComBat
tsne.combat = Tsne(dat.combat, batch, "ComBat")
# Limma 
tsne.limma = Tsne(dat.limma, batch, "limma")
# Harmony
tsne.harmony = Tsne(dat.harmony, batch, "Harmony")
# SAUCIE
tsne.saucie = Tsne(dat.saucie, batch, "SAUCIE")

# 3. Metrics:
# PCA and select the top 20 PCs
# Features in rows, cells in columns
PC20 = function(Dat_correction){
  dat.prcomp = prcomp(Dat_correction, scale. = TRUE)
  dat.pca = as.data.frame(dat.prcomp$rotation)
  pc20 = dat.pca[,c(1:20)]
  return(pc20)
}

# ComBat
pc20.combat = PC20(t(dat.combat))
# Limma
pc20.limma = PC20(t(dat.limma))
# Harmony
pc20.harmony = PC20(t(dat.harmony))
# SAUCIE
pc20.saucie = PC20(t(dat.saucie))

# (1) kBET
library(kBET)
KBET = function(pc20, Batch){
  # Run kBET using a predefined list of k values (equal to 5%, 10%, 15%, 20%, and 25% of the sample size)
  # Get the median of all kBET rejection rates of each input of k value
  kBET.5 = kBET(pc20, Batch, k0 = 0.05*nrow(pc20), plot = FALSE)
  median.5 = median(kBET.5$stats$kBET.observed)
  kBET.10 = kBET(pc20, Batch, k0 = 0.1*nrow(pc20), plot = FALSE)
  median.10 = median(kBET.10$stats$kBET.observed)
  kBET.15 = kBET(pc20, Batch, k0 = 0.15*nrow(pc20), plot = FALSE)
  median.15 = median(kBET.15$stats$kBET.observed)
  kBET.20 = kBET(pc20, Batch, k0 = 0.2*nrow(pc20), plot = FALSE)
  median.20 = median(kBET.20$stats$kBET.observed)
  kBET.25 = kBET(pc20, Batch, k0 = 0.25*nrow(pc20), plot = FALSE)
  median.25 = median(kBET.25$stats$kBET.observed)
  # Calculate the median of the five different k values as the final kBET score
  median = median(c(median.5, median.10, median.15, median.20, median.25))
  result = list("5%" = median.5, "10%" = median.10, "15%" = median.15, "20%" = median.20, "25%" = median.25, "median" = median)
  return(result)
}

# ComBat
kBET.combat = KBET(pc20.combat, batch)
# Limma
kBET.limma = KBET(pc20.limma, batch)
# Harmony
kBET.harmony = KBET(pc20.harmony, batch)
# SAUCIE
kBET.saucie = KBET(pc20.saucie, batch)

# (2) LISI
library(lisi)
# iLISI: integration LISI
LISI = function(pc20, Batch){
  dat.lisi = compute_lisi(pc20, data.frame(batch = Batch), 'batch')
  # Use the maximum and minimum score to scale the lisi scores
  dat.lisi.scale = (dat.lisi$batch-min(dat.lisi$batch))/(max(dat.lisi$batch)-min(dat.lisi$batch)) 
  median = median(dat.lisi.scale)
  return(median)
}

# ComBat
lisi.combat = LISI(pc20.combat, batch)
# Limma
lisi.limma = LISI(pc20.limma, batch)
# Harmony
lisi.harmony = LISI(pc20.harmony, batch)
# SAUCIE
lisi.saucie = LISI(pc20.saucie, batch)

# (3) ASW
library(fpc)
ASW20 = function(pc20){
  # Determine the optimal clusters
  pamk.best = pamk(pc20)
  asw20 = c()
  # To get a stable result, repeat ASW for 20 times and calculate the median of 20 ASW scores
  for (i in 1:20){
    pamk.result = kmeans(pc20, pamk.best$nc)
    pamk.stats = cluster.stats(dist(pc20)^2, pamk.result$cluster)
    asw = pamk.stats$avg.silwidth
    asw20[i] = asw
  }
  # Normalize the ASW scores to 0-1
  asw20_scale = (asw20-min(asw20))/(max(asw20)-min(asw20)) 
  asw_median = median(asw20)
  return(asw_median)
}

# ComBat
ASW.combat = ASW20(pc20.combat)
# Limma
ASW.limma = ASW20(pc20.limma)
# Harmony
ASW.harmony = ASW20(pc20.harmony)
# SAUCIE
ASW.saucie = ASW20(pc20.saucie)

# (4) ARI
library(stats)
library(mclust)
ARI20 = function(pc20, Batch){
  ari20 = c()
  # To get a stable result, repeat ARI for 20 times
  for (i in 1:20){
    dat.kmeans = kmeans(pc20, 2)
    ari20[i] = adjustedRandIndex(Batch, dat.kmeans$cluster)
  }
  # Normalize the ASW scores to 0-1
  ari20_scale = (ari20-min(ari20))/(max(ari20)-min(ari20))
  ari_median = median(ari20_scale)
  return(ari_median)
}

# ComBat
ARI.combat = ARI20(pc20.combat, batch)
# Limma
ARI.limma = ARI20(pc20.limma, batch)
# Harmony
ARI.harmony = ARI20(pc20.harmony, batch)
# SAUCIE
ARI.saucie = ARI20(pc20.saucie, batch)
