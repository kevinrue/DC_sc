
library(scater)
library(scde)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# load data ----

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

# phenotype colour code ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

# subset SCESet ----

sce.endo <- sce.norm[!isSpike(sce.norm),]
dim(sce.endo)

sce.4h <- sce.endo[,sce.endo$Time == '4h']
dim(sce.4h)

# scde pre-processing ----

o.ifm.4h <- readRDS("rds/o.ifm.4h_v8.rds")

filterCounts <- function(m, counts = 10, cells = 10){
  apply(m, 1, function(e){
    return(sum(e >= counts) >= cells)
  })
}

sg.4h <- droplevels(sce.4h$Group)
keep.4h <- filterCounts(counts(sce.4h))
cd.4h <- counts(sce.4h)[keep.4h,]; storage.mode(cd.4h) <- 'integer'
dim(cd.4h)

sce.ifm.4h <- sce.4h[,rownames(o.ifm.4h)]
cd.ifm.4h <- counts(sce.ifm.4h); storage.mode(cd.ifm.4h) <- "integer"

o.prior.4h <- scde.expression.prior(
  models = o.ifm.4h, counts = cd.ifm.4h, show.plot = TRUE
)

# Helper functions for scde ----

convert.z.score <- function(x, one.sided = NULL) {
  z <- x$Z
  if(is.null(one.sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one.sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);
  }
  x <- cbind(
    x,
    p.value = pval
  )
  return(x);
}

addGENENAME <- function(x){
  x <- cbind(
    GENENAME = with(fData(sce.endo), gene_name[match(rownames(x), gene_id)]),
    x
  )
  return(x)
}

orderResults <- function(x){
  x <- x[with(x, order(abs(Z), decreasing = TRUE)),]
  return(x)
}

normExprsById <- function(geneId){
  geneName <- subset(fData(sce.endo),gene_id==geneId,"gene_name",drop=TRUE)
  gdata <- data.frame(
    norm_exprs = norm_exprs(sce.endo)[geneId,],
    pData(sce.endo)[,c("Infection","Status","Time")],
    row.names = sampleNames(sce.endo)
  )
  ggplot(gdata, aes(gsub("-", "\n", Infection), norm_exprs)) +
    geom_violin(
      aes(fill = Infection, alpha = Status),
      draw_quantiles = c(0.25, 0.5, 0.75)
    ) +
    geom_jitter(width = 0.1) +
    facet_grid(Time ~ Status) +
    scale_fill_manual(values = col.infection) +
    scale_alpha_discrete(range = c(0.3, 0.6, 0.9)) +
    ggtitle(sprintf("%s - %s", geneId, geneName)) +
    theme_minimal() +
    labs(y = "Normalised expression", x = "Infection") +
    guides(alpha = "none")
}

# Identify APOE +/- cells ----

APOE.id <- with(fData(sce.ifm.4h), gene_id[match("APOE", gene_name)])

APOE.sce <- sce.ifm.4h[APOE.id,]

library(GOexpress)
ggplot(
  ggfy(
    APOE.sce,
    pheno = c("Infection", "Status", "Treatment"),
    assay = "norm_exprs"
  ),
  aes(norm_exprs, colour = Treatment)
) +
  geom_dotplot(
    method = "histodot", dotsize = 0.5,
    stackdir = "down", stackratio = 3/4, width = 0.5) +
  geom_density() +
  facet_wrap(~Infection + Status, ncol = 1)

# DE ----

sce.ifm.4h$APOE <- ifelse(
  norm_exprs(sce.ifm.4h)[APOE.id,] > 0, "Positive", "Negative"
)
with(pData(sce.ifm.4h), table(Treatment, APOE))

sg.test <-
  factor(pData(sce.ifm.4h)[,"APOE"], levels = c("Positive", "Negative"))
sg.test[sce.ifm.4h$Treatment != "STM-D23580_Violet +"] <- NA
table(sg.test, useNA = "ifany")

scde.d23.violet.apoe <- scde.expression.difference(
  o.ifm.4h, cd.ifm.4h, o.prior.4h, sg.test, n.cores = 4, verbose = 1
)

scde.res <- orderResults(addGENENAME(convert.z.score(scde.d23.violet.apoe)))
View(scde.res)

head(scde.res[,c("GENENAME", "mle", "Z", "p.value")], n = 10)

normExprsById("ENSG00000182566")

# Correlation plots ----

sce.d23.violet <- sce.ifm.4h[,sce.ifm.4h$Treatment == "STM-D23580_Violet +"]
dim(sce.d23.violet)

corAPOE <- function(geneId = "ENSG00000182566"){
  geneIds <- c("ENSG00000130203", geneId)
  geneExprs <- norm_exprs(sce.d23.violet)[geneIds,]
  geneExprs <- as.data.frame(t(geneExprs))
  genesNames <- with(fData(sce.d23.violet), gene_name[match(geneIds, gene_id)])
  colnames(geneExprs) <- genesNames
  ggplot(geneExprs, aes_string(x = genesNames[1], y = genesNames[2])) +
    geom_point() +
    scale_y_continuous(limits = range(norm_exprs(sce.ifm.4h))) +
    scale_x_continuous(limits = range(norm_exprs(sce.ifm.4h)))
}

corAPOE("ENSG00000110536")
corAPOE("ENSG00000169727")
corAPOE("ENSG00000102962")
corAPOE("ENSG00000064012")

scde.test.gene.expression.difference(
  "ENSG00000064012", # CASP8
  o.ifm.4h, cd.ifm.4h, o.prior.4h, sg.test, n.cores = 4, verbose = 1
)

write.csv(scde.res, "26_out/scde.D23580.violet_APOE.csv")

gene.res <- scde.res[
  match(
    c(
      "APOE","CASP8","TLR4","LIMA1","CCL22","APOC1","LYZ","AP1G2","AFTPH", # 9
      "HLA-DQA1","CST3","CLEC4G","COMMD3","CLEC5A","PTGER3","COPS7B" # 7
    ),
    scde.res$GENENAME),]

ggplot(scde.res, aes(Z, -log10(p.value))) +
  geom_point(aes(colour = p.value < 0.01)) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey") +
  geom_text_repel(
    aes(label = GENENAME), gene.res,
    size = 2, min.segment.length = unit(0, "mm"), alpha = 0.75,
    nudge_x = ifelse(gene.res$Z > 0, 1, -1), nudge_y = 0, segment.size = 0.2,
    fontface = "bold"
  ) +
  scale_x_continuous(limits = max(scde.res$Z) * c(-1,1)) +
  guides(colour = "none") +
  labs(y = expression(-log[10]*" ("*italic(P)*"-value)")) +
  theme_bw()

ggsave("26_out/volcano_APOE_PosNeg_Vogel.pdf", height = 5, width = 6.5)
