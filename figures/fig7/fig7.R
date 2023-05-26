library(tidyverse)
library(maps)
library(usmap)
library(viridis)
library(dendextend)

load("data/fig7.RData")

gplots::heatmap.2(ind.mat, distfun = distfun,
                  hclustfun = hclustfun,
                  Colv = NA,
                  Rowv = hcp,
                  dendrogram = "row",
                  trace = "none",
                  density.info="none",
                  key = F,
                  margins = c(5,10),
                  col = c("aliceblue", "cadetblue1"),
                  colRow = col_labels)


plot_usmap(data = us_map_merge_2, values = "Cluster", regions = "states") +
  labs(fill = "Clusters") +
  #scale_fill_discrete(na.translate = F) +
  scale_fill_manual(values = c("maroon", "chartreuse3", "mediumpurple2", "coral" , "darkgoldenrod1", "aquamarine3"
                               , "bisque", "grey"
  ), na.translate = F)
