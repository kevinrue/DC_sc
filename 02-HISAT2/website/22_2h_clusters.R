
# setup ----

library(scater)
library(scran)
library(ggplot2)
library(scde)
library(ComplexHeatmap)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(goseq)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
o.ifm.2h <- readRDS("rds/o.ifm.2h_v8.rds")

scde.res <- readRDS("rds/13_cluster_scde.res_v2.rds")
goseq.res <- readRDS("rds/13_cluster_goseq.res_v2.rds")

scde.groupRes <- readRDS("rds/scde.res_v8.rds")

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

sce.endo <- sce.norm[!isSpike(sce.norm),]

# Subset by time point
sce.2h <- sce.endo[,sce.endo$Time == '2h']

# Filter detected genes
filterCounts <- function(m, counts = 10, cells = 10){
  apply(m, 1, function(e){
    return(sum(e >= counts) >= cells)
  })
}

sg.2h <- droplevels(sce.2h$Group)
keep.2h <- filterCounts(counts(sce.2h))
cd.2h <- counts(sce.2h)[keep.2h,]; storage.mode(cd.2h) <- 'integer'

# Reorder SCESet to match count matrix
sce.ifm.2h <- sce.2h[,rownames(o.ifm.2h)]
cd.ifm.2h <- counts(sce.ifm.2h); storage.mode(cd.ifm.2h) <- "integer"

# Colour map for heat map
exprsRange <- range(norm_exprs(sce.ifm.2h))
colorMap <- circlize::colorRamp2(
  c(min(exprsRange), median(exprsRange), max(exprsRange)),
  c("blue","white","red"))


# Prior of expression
o.prior.2h <- scde.expression.prior(
  models = o.ifm.2h, counts = cd.ifm.2h, show.plot = FALSE
)

# Functions
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
  ggplot(gdata, aes(Infection, norm_exprs)) +
    geom_violin(
      aes(fill = Infection, alpha = Status),
      draw_quantiles = c(0.25, 0.5, 0.75)
    ) +
    geom_jitter(height = 0, width = 0.1) +
    facet_grid(Time ~ Status) +
    scale_fill_manual(values = col.infection) +
    ggtitle(sprintf("%s - %s", geneId, geneName)) +
    theme_minimal() +
    labs(y = "Normalised expression") +
    guides(alpha = "none")
}
single.scde.2h <- function(gene, groupTarget, groupRef){
  sg.test <- factor(sce.ifm.2h$quickCluster, levels = c(groupTarget, groupRef))
  scde.test.gene.expression.difference(
    gene, o.ifm.2h, cd.ifm.2h, o.prior.2h, sg.test,
    n.cores = 4, verbose = 1
  )
}

sig.levels <- c(0.05, 0.01)
volcano.sig <- data.frame(
  P = sig.levels,
  level = as.character(sig.levels)
)
volcano.mle <- function(x, sub = NULL){
  varName <- deparse(substitute(x))
  x <- convert.z.score(x)
  gg <- ggplot(x, aes(mle, -log10(p.value))) +
    geom_point(aes(colour = (cZ != 0))) +
    geom_hline(aes(yintercept=-log10(p.value),linetype=level),volcano.sig) +
    ggtitle(varName, sub)
  print(gg)
  return(x)
}
showTopGO <- function(x, direction = "over", cutoff = 0.05, min = 10){
  goCols <- c(
    "category","over_represented_pvalue","under_represented_pvalue",
    "numDEInCat","numInCat","ontology","term")
  directionCol <- sprintf("%s_represented_pvalue", direction)
  tmp.go <- x[,goCols]
  tmp.sig <- (tmp.go[,directionCol] < cutoff); tmp.go <- tmp.go[tmp.sig,]
  tmp.min <- (tmp.go[,"numInCat"] >= min); tmp.go <- tmp.go[tmp.min,]
  tmp.go <- tmp.go[order(tmp.go[,directionCol]),]
  tmp.go <- head(tmp.go, 50)
  DT::datatable(
    tmp.go, list(pageLength = 10, searching = TRUE), filter = "top",
    rownames = FALSE
  )
}

# identify clusters ----

clusters <- quickCluster(sce.ifm.2h,min.size=min(table(sce.ifm.2h$Treatment)))
sce.ifm.2h$quickCluster <- clusters
qc.2h <- rep(NA, ncol(sce.endo)); names(qc.2h) <- sampleNames(sce.endo)
qc.2h[sampleNames(sce.ifm.2h)] <- as.numeric(levels(clusters)[clusters])
sce.endo$quickCluster.2h <- as.factor(qc.2h)

# show clusters ----

with(pData(sce.endo), table(Treatment, quickCluster.2h))

col.cluster <- brewer.pal(12, "Set3")[c(1,12,3,4,9)]
names(col.cluster) <- c(1:3,0,"NA")

tmpGG <- data.frame(
  reducedDimension(sce.endo),
  Cluster = sce.endo$quickCluster.2h,
  Treatment = gsub("_", " ", sce.endo$Treatment),
  Infection = sce.endo$Infection,
  Status = sce.endo$Status
)
ggplot(tmpGG, aes(X1, X2, shape = Treatment, colour = Cluster)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x="Dimension 1",y="Dimension 2",title="2h") +
  scale_color_manual(values = col.cluster, na.value = col.cluster["NA"])

ggplot(tmpGG, aes(X1, X2, colour = Cluster)) +
  geom_point(size = 3) +
  facet_grid(Infection ~ Status) +
  theme_bw() +
  labs(x="Dimension 1",y="Dimension 2",title="2h") +
  scale_color_manual(values = col.cluster, na.value = col.cluster["NA"])

ggplot(subset(pData(sce.ifm.2h), Time == "2h")) +
  facet_wrap(~ quickCluster) +
  geom_bar(aes(Treatment, fill = quickCluster)) +
  scale_fill_manual(values = col.cluster) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("22_out/facet_wrap_cluster_count.pdf", height = 8, width = 9)

# DE between challenged cluster 2 vs challenged cluster 1 ----

sg.test <- rep(NA, ncol(sce.ifm.2h))
sg.test[
  with(pData(sce.ifm.2h), which(Infection != "Mock" & quickCluster == 2))] <-
  "challenged_cluster2"
sg.test[
  with(pData(sce.ifm.2h), which(Infection != "Mock" & quickCluster == 1))] <-
  "challenged_cluster1"
sg.test <- factor(sg.test, paste0("challenged_cluster", 2:1))
table(sg.test)

names(sg.test) <- sampleNames(sce.ifm.2h)
contrastName <- "challenged_cluster2-cluster1"
scde.res[[contrastName]] <- scde.expression.difference(
  o.ifm.2h, cd.2h, o.prior.2h, sg.test, n.cores = 4, verbose = 1
)

assign(contrastName, orderResults(addGENENAME(convert.z.score(scde.res[[contrastName]]))))

View(get(contrastName))

write.csv(get(contrastName), sprintf("22_out/%s.csv", contrastName))

# CD1A
scde.test.gene.expression.difference("ENSG00000158477", o.ifm.2h, cd.2h, o.prior.2h, sg.test)
# CTSL
scde.test.gene.expression.difference("ENSG00000135047", o.ifm.2h, cd.2h, o.prior.2h, sg.test)

normExprsById("ENSG00000135047")

# randomForest ----

library(randomForest)

sce.rf <- sce.ifm.2h[
  , sce.ifm.2h$quickCluster %in% c(1, 2) & sce.ifm.2h$Infection != "Mock"]

rf <- randomForest(
  x = t(norm_exprs(sce.rf)),
  y = droplevels(sce.rf$quickCluster),
  importance = TRUE,
  do.trace = 100,
  ntree = 500
)
plot(rf)

View(addGENENAME(as.data.frame(round(importance(rf), 2))))


rf.topGenes <- head(rownames(importance(rf))[
  order(importance(rf)[,"MeanDecreaseGini"], decreasing = TRUE)
  ], 50)

mat.rf <- norm_exprs(sce.rf)[rf.topGenes,]
rownames(mat.rf) <- fData(sce.rf)[rf.topGenes, "gene_name"]

ha <- HeatmapAnnotation(
  df = data.frame(
    quickCluster = sce.rf$quickCluster,
    Infection = sce.rf$Infection,
    Status = sce.rf$Status
  ),
  col = list(
    quickCluster = col.cluster,
    Infection = col.infection,
    Status = col.status
  )
)

hm <- Heatmap(
  mat.rf,
  colorMap, "norm.\nexprs.",
  cluster_rows = TRUE, row_names_side = "left",
  cluster_columns = TRUE,
  show_column_names = FALSE, show_column_dend = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10), top_annotation = ha,
  column_title = "2h")

draw(hm)

# Monocle pseudotemporal trajectory ----

library(monocle)

fd <- featureData(sce.ifm.2h)
fd$gene_short_name <- fd$gene_name

# Make newCellDataSet
HSMM <- newCellDataSet(
  counts(sce.ifm.2h)[!isSpike(sce.ifm.2h),],
  phenoData = phenoData(sce.ifm.2h),
  featureData = fd[!isSpike(sce.ifm.2h),],
  lowerDetectionLimit = 0.1,
  expressionFamily=negbinomial.size()
)
HSMM$isMock <- factor(
  ifelse(HSMM$Infection == "Mock", "Mock", "Challenged"),
  c("Mock", "Challenged")
)

# Compulsory step
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Unsupervised ?!
diff_test_res <- differentialGeneTest(HSMM, fullModelFormulaStr="~isMock")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2)
HSMM <- orderCells(HSMM)

write.csv(subset(diff_test_res, qval < 0.01), "22_out/diff_test_res.csv")

# Find state with maximum uninfected
root_state <- function(cds, refPheno = "Infection", refState = "Mock"){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)[,"State"], pData(cds)[,refPheno])[,refState]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) }else { return (1) }
}
HSMM <- orderCells(HSMM, root_state = root_state(HSMM))

HSMM$State <- factor(HSMM$State, c(4,3,2,5,1))

plot_cell_trajectory(HSMM, color_by = "State")
ggsave("22_out/pseudotime_State.pdf", height = 4, width = 5)
plot_cell_trajectory(HSMM, color_by = "Status")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
ggsave("22_out/pseudotime_Pseudotime.pdf", height = 4, width = 5)
plot_cell_trajectory(HSMM, color_by = "quickCluster")
plot_cell_trajectory(HSMM, color_by = "Treatment")

plot_cell_trajectory(HSMM, color_by = "Status") + facet_wrap(~State, nrow=1)

HSMM$vState <- with(pData(HSMM), factor(paste("State", State), paste("State", levels(State))))

# HSMM$State <- factor(HSMM$State, c(4,3,2,5,1))

plot_cell_trajectory(HSMM, color_by = "Treatment") +
  facet_wrap(~vState, nrow=1)
ggsave("22_out/pseudotime_Treatment_facetOrderedStates.pdf", height = 4, width = 12)

plot_cell_trajectory(HSMM, color_by = "quickCluster") + facet_wrap(~State, nrow=1)
ggsave("22_out/pseudotime_quickCluster_facetOrderedStates.pdf", height = 4, width = 12)

# BEAM branching point 1 ----

# T <- minSpanningTree(HSMM)

bds <- buildBranchCellDataSet(HSMM, branch_point = 1)

BEAM_1 <- BEAM(HSMM, branch_point = 1, cores = 4)
# BEAM_1 <- BEAM(
#   HSMM, cores = 4,
#   branch_states = c(2, 3), branch_labels = c("Activated", "Mock-like")
# )
BEAM_1 <- BEAM_1[order(BEAM_1$qval),]
BEAM_1 <- BEAM_1[,c("gene_short_name", "pval", "qval")]

View(BEAM_1)

pdf("22_out/heatmap_point1.pdf", height = 12, width = 6)
plot_genes_branched_heatmap(
  HSMM[row.names(subset(BEAM_1, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4,
  cores = 1,
  use_gene_short_name = T,
  show_rownames = T)
dev.off()

point1_genes <- row.names(
  subset(fData(HSMM), gene_short_name %in% c("LIPA", "APOE", "APOC1", "LAMP3")))

gg <- plot_genes_branched_pseudotime(
  HSMM[point1_genes,],
  branch_point=1,
  color_by = "Status",
  ncol=1, cell_size = 3
)
gg

# BEAM branching point 2 ----

BEAM_2 <- BEAM(HSMM, branch_point = 1, cores = 4)
# BEAM_2 <- BEAM(
#   HSMM, cores = 4,
#   branch_states = c(2, 3), branch_labels = c("Activated", "Mock-like")
# )
BEAM_2 <- BEAM_2[order(BEAM_2$qval),]
BEAM_2 <- BEAM_2[,c("gene_short_name", "pval", "qval")]
View(BEAM_2)

pdf("22_out/heatmap_point2.pdf", height = 12, width = 6)
plot_genes_branched_heatmap(
  HSMM[row.names(subset(BEAM_2, qval < 1e-4)),],
  branch_point = 2,
  num_clusters = 4,
  cores = 4,
  use_gene_short_name = T,
  show_rownames = T)
dev.off()

point1_genes <- row.names(
  subset(fData(HSMM), gene_short_name %in% c("APOE", "APOE", "APOC1", "LAMP3")))

gg <- plot_genes_branched_pseudotime(
  HSMM[point1_genes,],
  branch_point=1,
  color_by = "Status",
  ncol=1, cell_size = 3
)
gg

# Profile of states ----

sce.ifm.2h$State <- pData(HSMM)[sampleNames(sce.ifm.2h),"State"]
sce.ifm.2h$Pseudotime <- pData(HSMM)[sampleNames(sce.ifm.2h),"Pseudotime"]

ggplot(
  cbind(
    n.e. = norm_exprs(sce.ifm.2h)[
      subset(fData(sce.ifm.2h), gene_name == "CD52", gene_id, drop = TRUE)
    ,],
    pData(sce.ifm.2h)[,c("Time","Infection","Status","State")]
  ),
  aes(x = State, y = n.e.)
) + geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  geom_jitter(aes(colour = Status), height = 0, width = 0.3) +
  scale_y_continuous(limits = range(norm_exprs(sce.ifm.2h))) +
  labs(y = "Normalised expression", x = NULL) +
  theme_minimal()


assayData(bds)[["log2"]] <- log2(exprs(bds) + 1)
lipa.id <- subset(fData(bds), gene_name == "IL1B", gene_id, drop = TRUE)
lipa.df <- data.frame(
  log2 = assayData(bds)[["log2"]][lipa.id,],
  State = pData(bds)[,"vState"],
  Branch = pData(bds)[,"Branch"]
)
# ggplot(lipa.df) + geom_jitter(aes(as.factor(State), log2, colour = State))
ggplot(lipa.df) +
  geom_violin(aes(as.factor(State), log2, colour = State)) +
  geom_jitter(aes(as.factor(State), log2, colour = State), alpha = 0.5)

varLabels(bds)
bds$Branch
dim(bds)

with(pData(bds), table(Branch, Infection))

# Volcano plots ----

mleRange <-
  max(abs(do.call("c", lapply(scde.res, function(x){return(x$mle)}))))*c(-1,1)
