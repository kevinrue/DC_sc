
# Based on the manuscript,
# Highlight the higher expression of some genes in D2350 infected
# relative to all other groups

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
sce.norm

sce.6h <- sce.norm[,sce.norm$Time == "6h"]
dim(sce.6h)

IL1A.id <- with(fData(sce.6h), gene_id[match("IL1A", gene_name)])

df <- data.frame(
  ne = norm_exprs(sce.6h)[IL1A.id,],
  Treatment = sce.6h$Treatment
)
