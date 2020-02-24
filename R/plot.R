library(ggplot2)

umapPlot = function(umap, method, colored_by = c("labels", "batch")){
  if(colored_by == "labels"){
    df = data.frame(UMAP1 = umap$layout[,1],
                    UMAP2 = umap$layout[,2],
                    Labels = Label$level1)
    gp = ggplot(df, aes(UMAP1, UMAP2, color = Labels)) +
      geom_point(size = -0.1, alpha = 0.1) +
      ggtitle(method) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            legend.title = element_text(size = 15)) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 5)))
  }
  if(colored_by == "batch"){
    df = data.frame(UMAP1 = umap$layout[,1],
                    UMAP2 = umap$layout[,2],
                    Batch = Label$batch)
    gp = ggplot(df, aes(UMAP1, UMAP2, color = Batch)) +
      geom_point(size = -0.1, alpha = 0.1) +
      ggtitle(method) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            legend.title = element_text(size = 15)) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 5)))
  }
  return(gp)
}

