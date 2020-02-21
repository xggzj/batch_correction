library(ggplot2)



# Plot

plotForUmap = function(umap, Label) {
    df = data.frame(UMAP1 = umap.dat$layout[,1],
                UMAP2 = umap.dat$layout[,2],
                batch = Label$batch)
    gp = ggplot(df, aes(UMAP1, UMAP2, color = ))

    return(gp)
}


