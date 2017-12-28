o.ifm.2h.02 <- readRDS("02-HISAT2/website/rds/o.ifm.2h_v8.rds")
o.ifm.2h.10 <- readRDS("10_website_sce/rds/o.ifm_2h.rds")

View(o.ifm.2h.02)
View(o.ifm.2h.10)

compNamesTables <- data.frame(
  old = rownames(o.ifm.2h.02),
  new = rownames(o.ifm.2h.10)
)
View(compNamesTables)


sce.pass.02 <- readRDS("02-HISAT2/website/rds/sce.pass.rds")

head(sce.pass.02@phenoData@data)
colnames(sce.pass.02@phenoData@data)
subset(sce.pass.02@phenoData@data, Time == "2h")
rownames(subset(sce.pass.02@phenoData@data, Time == "2h"))

setdiff(paste("Cell", rownames(o.ifm.2h.02), sep = "_"), rownames(o.ifm.2h.10))
head(sort(paste("Cell", rownames(o.ifm.2h.02), sep = "_")))
head(sort(rownames(o.ifm.2h.10)))
# Conclusion: both models have the same cells... in different orders!
# The key problem here is that values
