
library(iSEE)
library(RColorBrewer)

quickClusterColor <- function(i=5L) {
    # col.cluster <- brewer.pal(12, "Set3")[c(1,12,3,4,9)]
    colors <- brewer.pal(8, "Dark2")[c(1,2,3,4,8)]
    names(colors) <- c(1,2,3,0,"NA")
    return(col.cluster)
}

ecm <- ExperimentColorMap(
    colData = list(
        quickCluster=quickClusterColor
    )
)
