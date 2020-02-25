library(umap)
source('preprocessing.R', chdir = TRUE)
source('harmony.R', chdir = TRUE)
source('plot.R', chdir = TRUE)

# 1. Preprocessing
# Read infiles: expression data (.fcs), cell type data (.csv)
# Read marker names and remove empty channels
# Remove outliers 
# Scaling

# 2. Batch correction
# Bind the expression data of all batches into a single matrix
# Bind the cell type labels of all batches into a single dataframe
# Input the batch number
# MetaData: a dataframe with variables to integrate including batch numbers and cell type labels
# Harmony batch effect correction
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
Label$batch = as.factor(c(rep(1,nrow(dat[[1]])), rep(2,nrow(dat[[2]]))))
MetaData = data.frame(batch = Label$batch, level1 = Label$level1, level2 = Label$level2)
Dat.harmony = HarmonyMatrix(Dat, MetaData, "batch", do_pca = FALSE)
dat.harmony = saveCorrectedFiles(Dat.harmony, "~/master/output")

umap.harmony = umap(Dat.harmony)
umap.dat = umap(Dat)

gp1 = umapPlot(umap.dat, "Raw", "batch")
gp2 = umapPlot(umap.dat, "Raw", "level1")
gp3 = umapPlot(umap.dat, "Raw", "level2")
gp4 = umapPlot(umap.harmony, "Harmony", "batch")
gp5 = umapPlot(umap.harmony, "Harmony", "level1")
gp6 = umapPlot(umap.harmony, "Harmony", "level2")

