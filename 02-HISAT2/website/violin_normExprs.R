library(scater)
library(scran)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

outDir <- "violin_normExprs"
if(!dir.exists(outDir)){
  dir.create(outDir)
}

sce.endo <- readRDS("rds/13_sce.endo_clusters.rds")

col.set3 <- brewer.pal(12, "Set3")
col.cluster <- col.set3[c(1,12,3,4,9)]; names(col.cluster) <- c(1:3,0,"NA")
col.status<-col.set3[c(3,5:6)];names(col.status)<-levels(sce.endo$Status)
col.infection<-col.set3[c(3,8,10)];names(col.infection)<-levels(sce.endo$Infection)
col.time<-col.set3[c(2,7:8)];names(col.time)<-levels(sce.endo$Time)

getGeneId <- function(geneName){
  subset(fData(sce.endo), gene_name == geneName, "gene_id", drop = TRUE)
}

normExprsByName.4h <- function(geneId){
  sc <- sce.endo[,with(pData(sce.endo), Infection != "Mock" & Time == "4h")]
  gdata <- data.frame(
    norm_exprs = norm_exprs(sc)[geneId,],
    pData(sc)[,c("Infection", "Status", "quickCluster.4h")],
    row.names = sampleNames(sc)
  )
  geneName <- with(fData(sc), gene_name[gene_id == geneId])
  ggplot(gdata, aes(Status, norm_exprs)) +
    geom_violin(aes(fill = Status)) +
    geom_jitter(
      aes_string(colour = "quickCluster.4h"),
      width = 0.25, size = 3) +
    facet_grid(~Infection) +
    ggtitle(sprintf("%s - %s", geneName, geneId)) +
    scale_y_continuous(limits=range(norm_exprs(sc))) +
    scale_fill_manual(values = col.status) +
    scale_colour_manual(values = col.cluster) +
    theme_minimal()
}

savePlot <- function(geneName){
  gg <- normExprsByName.4h(getGeneId(geneName))
  ggsave(
    file.path(outDir, sprintf("4h_stimulated_%s.pdf", geneName)), gg,
    height = 5, width = 7
  )
  gg
}

savePlot("EBI3")
savePlot("PLAT")

savePlot("MS4A4A")
savePlot("CTSL")
savePlot("CTSD")


