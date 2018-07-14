library(monocle)
library(reshape2)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

fd <- rowData(sce.norm)
fd$gene_short_name <- fd$gene_name

HSMM <- newCellDataSet(
  counts(sce.norm)[!isSpike(sce.norm),],
  phenoData = AnnotatedDataFrame(as.data.frame(colData(sce.norm))),
  featureData = AnnotatedDataFrame(as.data.frame(fd[!isSpike(sce.norm),])),
  lowerDetectionLimit = 0.1,
  expressionFamily=negbinomial.size()
)

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 10)

# Same definition as for SCDE: >= 10 counts in >= 10 cells
expressed_genes <- row.names(subset(
  fData(HSMM), num_cells_expressed >= 10)
)

print(head(fData(HSMM)))
print(head(pData(HSMM)))

valid_cells <- row.names(pData(HSMM))

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

upper_bound <-
  10^(
    mean(log10(pData(HSMM)$Total_mRNAs)) +
      2*sd(log10(pData(HSMM)$Total_mRNAs))
  )
lower_bound <-
  10^(
    mean(log10(pData(HSMM)$Total_mRNAs)) -
      2*sd(log10(pData(HSMM)$Total_mRNAs))
  )

qplot(Total_mRNAs, data=pData(HSMM), color=Time, geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_minimal()

qplot(Total_mRNAs, data=pData(HSMM), color=Status, geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_minimal()

ggplot(pData(HSMM)) +
  facet_grid(Treatment ~ Time) +
  geom_histogram(aes(Total_mRNAs), binwidth = 1E5) +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound) +
  scale_x_continuous(labels = scales::scientific) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]) + 1)

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

qplot(value, geom="density", data=melted_dens_df) +
  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(counts)") +
  ylab("Density")

# subset the CellDataSet object to create one that includes only stimulated
# cells

HSMM_stim <- HSMM[,pData(HSMM)$Infection != "Mock"]
HSMM_stim <- estimateDispersions(HSMM_stim)

# Unsupervised ordering

# One effective way to isolate a set of ordering genes is to simply compare the
# cells collected at the beginning of the process to those at the end and find
# the differentially expressed genes, as described above.

diff_test_res <- differentialGeneTest(
  HSMM_stim[expressed_genes,],
  fullModelFormulaStr="~Time"
)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM_stim <- setOrderingFilter(HSMM_stim, ordering_genes)
plot_ordering_genes(HSMM_stim)

# For a number of reasons, Monocle works better if we can reduce the
# dimensionality of that space before we try to put the cells in order.

HSMM_stim <- reduceDimension(HSMM_stim, max_components=2)

HSMM_stim <- orderCells(HSMM_stim)

plot_cell_trajectory(HSMM_stim, color_by = "Time")

plot_cell_trajectory(HSMM_stim, color_by = "Status")

plot_cell_trajectory(HSMM_stim, color_by = "Infection")

plot_cell_trajectory(HSMM_stim, color_by = "Treatment")
