library(umap)
source('preprocessing.R', chdir = TRUE)
source('harmony.R', chdir = TRUE)
source('plot.R', chdir = True)

# 1. Preprocessing
# Read infiles
# Read expression data (.fcs):
# Read cell type data (.csv):
# Read marker names and remove empty channels
# Remove outliers 
# Scaling

# 2. Batch correction
# Bind the expression data of all batches into a single matrix
# Bind the cell type labels of all batches into a single dataframe
# Input the batch number:
# MetaData: a dataframe with variables to integrate including batch numbers and cell type labels
# Save the batch-corrected files separately

# 3. Visualization 
# UMAP
# Plot


files = readInfiles("~/master/infiles", "fcs")
labels = readInfiles("~/master/labels", "csv")
dat = lapply(files, channelFilter)
dat = lapply(dat, function(x){apply(x, 2, removeOutliers)})
dat = lapply(dat, scale)

Dat = bindlist(dat)
Label = bindlist(labels)
Label$batch = c(rep(1,nrow(dat[[1]])), rep(2,nrow(dat[[2]])))
MetaData = data.frame(batch = Label$batch, label = Label$level1)
Dat.harmony = HarmonyMatrix(Dat, MetaData, "batch")
dat.harmony = saveCorrectedFiles(Dat.harmony)

umap.harmony = umap(Dat.harmony)
umap.dat = umap(Dat)

gp = plotForUmap(umap)