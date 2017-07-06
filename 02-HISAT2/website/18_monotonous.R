library(plyr)
library(dplyr)
library(broom)
library(monocle)

# Normalised
sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

# Subset to endogenous genes
sce.endo <- sce.norm[!isSpike(sce.norm),]; dim(sce.endo)
# Subset to stimulated DCs
sce.stim <- sce.endo[,sce.endo$Infection != "Mock"]; dim(sce.stim)

p.time <- as.numeric(gsub("h", "", sce.stim$Time))

# g1 <- data.frame(
#   Gene = norm_exprs(sce.stim)[1,],
#   Time = p.time
# )
#
# ggplot(g1) +
#   geom_jitter(aes(Time, Gene, colour = Time), height = 0, width = 0.2)
#
# glm.out <- glm(formula = Gene ~ Time, data = g1)
# lm.out <- lm(formula = Gene ~ Time, data = g1)
#
# summary(glm.out)
# tidy(glm.out)
#
# summary(lm.out)
# tidy(lm.out)
#
# hist(glm.out$residuals)

glm.out <- do.call(
  "rbind",
  apply(
    norm_exprs(sce.stim),
    1,
    function(x){
      myData <- data.frame(
        Gene = x,
        Time = p.time
      )
      glm.out <- glm(formula = Gene ~ Time, data = myData)
      tidy.out <- tidy(glm.out)
      subset(tidy.out, term == "Time", select = -term)
    }
  )
)


# Count of genes positively correlated with time in stimulated cells
dim(
  subset(glm.out, estimate > 0 & p.value < 0.01)
)

# Count of genes by significance and direction of time effect
ggplot(glm.out) +
  geom_histogram(aes(p.value, fill = estimate > 0), binwidth = 0.01)

glm.out$BH <- p.adjust(glm.out$p.value, "BH")

dim(
  subset(glm.out, estimate > 0 & BH < 0.01)
)

glm.out$gene_name <- with(fData(sce.stim), gene_name[match(rownames(glm.out), gene_id)])

write.csv(
  glm.out[order(-glm.out$estimate),],
  "18_out/glm.out.csv"
)

View(
  subset(glm.out, estimate > 0 & BH < 0.05)
)

head(
  subset(glm.out[order(glm.out$p.value),], estimate > 0 & BH < 0.05)
)

showGene <- function(geneid){
  ggData <- data.frame(
    NormExprs = norm_exprs(sce.endo)[geneid,],
    Time = sce.endo$Time,
    Challenged = factor(ifelse(sce.endo$Infection == "Mock", "Mock", "Challenged"), c("Mock", "Challenged")),
    Infection = sce.endo$Infection,
    Status = sce.endo$Status
  )
  ggplot(ggData) +
    facet_grid(~ Challenged) +
    geom_jitter(aes(Time, NormExprs, colour = Status), height = 0, width = 0.3) +
    labs(
      y = "Normalised expression",
      title = sprintf("%s - %s", geneid, with(fData(sce.endo), gene_name[match(geneid, gene_id)]))
    )
}

showGene("ENSG00000172183")

dim(
  subset(glm.out, estimate > 2 & BH < 0.01) # more stringent estimate
)

dim(
  subset(glm.out, estimate < -2 & BH < 0.01) # more stringent estimate
)

head(dplyr::arrange(
  subset(glm.out, estimate < 0 & BH < 0.01),
  estimate
))

fd <- featureData(sce.endo)
fd$gene_short_name <- fd$gene_name

# Make newCellDataSet
HSMM <- newCellDataSet(
  counts(sce.endo)[!isSpike(sce.endo),],
  phenoData = phenoData(sce.endo),
  featureData = fd[!isSpike(sce.endo),],
  lowerDetectionLimit = 0.1,
  expressionFamily=negbinomial.size()
)
HSMM$isMock <- factor(ifelse(HSMM$Infection == "Mock", "Mock", "Challenged"), c("Mock", "Challenged"))

# Compulsory step
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Select genes that increase monotonously in challenged DCs
ordering_genes <- rownames(subset(glm.out, estimate > 2 & BH < 0.01))

# Using genes to order the cells as above yields the following trajectory:
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2)
HSMM <- orderCells(HSMM)

pdf("18_out/plot_cell_trajectory_Time.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by = "Time")
dev.off()

pdf("18_out/plot_cell_trajectory_Status.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by = "Status")
dev.off()

pdf("18_out/plot_cell_trajectory_Infection.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by = "Infection")
dev.off()

pdf("18_out/plot_cell_trajectory_State.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

pdf("18_out/plot_cell_trajectory_Pseudotime.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by="Pseudotime")
dev.off()

pdf("18_out/plot_cell_trajectory_isMock.pdf", width = 6, height = 4)
plot_cell_trajectory(HSMM, color_by="isMock")
dev.off()

# Count of cells by experimental group and state
with(pData(HSMM), table(State, Time))
with(pData(HSMM), table(State, Infection))

# Finding genes that change in  as a function of pseudotime

HSMM.stim <- HSMM[,HSMM$Infection != "Mock"]

diff_test_res <- differentialGeneTest(HSMM.stim, fullModelFormulaStr="~sm.ns(Pseudotime)")

write.csv(
  diff_test_res[order(diff_test_res$qval),],
  "18_out/diff_test_res_Pseudotime.csv")

head(dplyr::arrange(
  diff_test_res[,c("gene_short_name", "pval", "qval")],
  pval
))

View(diff_test_res)

table(diff_test_res$qval < 0.01)

# Bug fix: avoid all-numeric cell names
sampleNames(HSMM) <- gsub("^", "C.", sampleNames(HSMM))

pdf("18_out/plot_genes_in_pseudotime_Time.pdf", width = 6, height = 10)
plot_genes_in_pseudotime(
  HSMM[head(rownames(diff_test_res)[order(diff_test_res$qval)], 10),],
  color_by = "Time") + theme_bw()
dev.off()

pdf("18_out/plot_genes_in_pseudotime_Infection.pdf", width = 6, height = 10)
plot_genes_in_pseudotime(
  HSMM[head(rownames(diff_test_res)[order(diff_test_res$qval)], 10),],
  color_by = "Infection") + theme_bw()
dev.off()

HSMM.top <- HSMM[head(rownames(diff_test_res)[order(diff_test_res$qval)]),]
plotExprs <- reshape2::melt(exprs(HSMM.stim))

pdf("18_out/plot_pseudotime_heatmap_30.pdf", width = 7, height = 4)
plot_pseudotime_heatmap(
  HSMM[head(rownames(diff_test_res)[order(diff_test_res$qval)], 30),],
  num_clusters = 3,
  cores = 1,
  show_rownames = T)
dev.off()
