
names(assayData(sce.norm))

plot(sort(log10(rowSums(counts(sce.norm)) + 1 )))

sce.2h <- sce.norm[!isSpike(sce.norm), sce.norm$Time == "2h"]
dim(sce.2h)
summary(sce.2h$Time)

var.norm <- apply(norm_exprs(sce.2h), 1, var)
summary(var.norm)

# Remove feature @ 0 variance
exprs_to_plot <- scale(t(norm_exprs(sce.2h)))

var.scaled <- apply(exprs_to_plot, 2, var)

which.min(var.norm)
var.scaled[which.min(var.norm)]
plot(exprs_to_plot[,which.min(var.norm)])
plot(norm_exprs(sce.2h)[which.min(var.norm),])
plot(counts(sce.2h)[which.min(var.norm),])

plot(apply(exprs_to_plot, 2, var))
plot(matrixStats::colVars(exprs_to_plot))

keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
keep_feature[is.na(keep_feature)] <- FALSE
exprs_to_plot <- exprs_to_plot[, keep_feature]

pca.2h <- prcomp()
dim(pca.2h$x)
dim(pca.2h$rotation)

plot(pca.2h)

x.2h <- pca.2h$x
all(rownames(x.2h) == colnames(sce.2h))
x.2h <- cbind(x.2h, pData(sce.2h))

names(x.2h)

ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = Infection))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = Status))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = Lane))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = Plate))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = log10_total_counts))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = total_features))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = pct_exprs_feature_controls_MT))
ggplot(x.2h) +
  geom_point(aes(PC1, PC2, colour = pct_dropout))

pd <- pData(sce.2h)

# Identify of highest loading on PC2 ---

# str(pca.2h)
# dimnames(pca.2h$rotation)
# order(pca.2h$rotation[,"PC2"])
# order(abs(pca.2h$rotation[,"PC2"]))
# rownames(pca.2h$rotation)[order(abs(pca.2h$rotation[,"PC2"]))]

pc2.load <- pca.2h$rotation[,"PC2"]
pc2.load.order <- order(abs(pc2.load), decreasing = TRUE)

plot(sort(abs(pca.2h$rotation[,"PC2"]), decreasing = TRUE))

pc2.geneIds <- rownames(pca.2h$rotation)[pc2.load.order]
head(pc2.geneIds)

library(EnsDb.Hsapiens.v79)
# columns(EnsDb.Hsapiens.v79)
pc2.geneNames <- mapIds(EnsDb.Hsapiens.v79, pc2.geneIds, "GENENAME", "GENEID")

pc2.df <- data.frame(
  id = pc2.geneIds,
  name = pc2.geneNames,
  pca.2h$rotation[pc2.load.order,]
)


# Expression of most important genes on PC2 ----

names(assayData(sce.2h))
pc2.norm <- norm_exprs(sce.2h)[geneIds,]
pc2.count <- counts(sce.2h)[geneIds,]

plot(sort(log10(rowSums(pc2.count) + 1)))
pc2.noZero <- pc2.count[rowSums(pc2.count) > 0,]
dim(pc2.noZero)

pc2.count.noZero <- pc2.count[rownames(pc2.noZero),]
