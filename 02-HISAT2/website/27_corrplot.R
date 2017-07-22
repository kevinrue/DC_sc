
library(corrplot)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
sce.norm

sce.4h <- sce.norm[,sce.norm$Time == "4h"]
dim(sce.4h)

geneNames <- c("CASP8","TLR4","CCL22","LYZ","CST3","COMMD3","CLEC5A","CLEC4G")
geneIds <- with(fData(sce.4h), gene_id[match(geneNames, gene_name)])

sce.genes <- sce.4h[geneIds,]
dim(sce.genes)

featureNames(sce.genes) <- geneNames

sce.d23 <- sce.genes[,sce.genes$Treatment == "STM-D23580_Violet +"]
dim(sce.d23)

sce.lt2 <- sce.genes[,sce.genes$Treatment == "STM-LT2_Violet +"]
dim(sce.lt2)

cor.d23 <- cor(t(norm_exprs(sce.d23)))
cor.lt2 <- cor(t(norm_exprs(sce.lt2)))

col2 <- colorRampPalette(
  rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
    "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"
  ))
)

pdf("27_out/corplot_D23580.pdf", height = 10, width = 10)
corrplot(
  cor.d23, col = col2(200), method = "number", number.cex = 1.75,
  tl.pos = "lt", type = "lower", tl.cex = 2
)
corrplot(
  cor.d23, col = col2(200),
  type = "upper", tl.srt = 45, tl.pos = "n", add = TRUE
  )
dev.off()

pdf("27_out/corplot_LT2.pdf", height = 10, width = 10)
corrplot(
  cor.lt2, col = col2(200), method = "number", number.cex = 1.75,
  tl.pos = "lt", type = "lower", tl.cex = 2
)
corrplot(
  cor.lt2, col = col2(200),
  type = "upper", tl.srt = 45, tl.pos = "n", add = TRUE
)
dev.off()



