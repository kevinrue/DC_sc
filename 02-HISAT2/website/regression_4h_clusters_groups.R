library(broom)
library(reshape2)
library(RColorBrewer)

# Load previous data ----

scde.clusters <- readRDS("rds/13_cluster_scde.res_v2.rds")
names(scde.clusters)

scde.groups <- readRDS("rds/scde.res_v8.rds")
names(scde.groups)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
sce.norm

# Colour code ----

colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- c("STM-D23580", "STM-LT2")

# Extract DE tables ----

clusters.de <- scde.clusters[["4h_cluster2-cluster1"]]
head(clusters.de, n=2)

groups.de.d23 <- scde.groups[["4h_STM-D23580_Violet +-4h_STM-D23580_Exposed"]]
head(groups.de.d23, n=2)

groups.de.lt2 <- scde.groups[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]]
head(groups.de.lt2, n=2)

# Convert Z-score to P-value ----

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

clusters.de <- convert.z.score(clusters.de)

groups.de.d23 <- convert.z.score(groups.de.d23)
groups.de.lt2 <- convert.z.score(groups.de.lt2)

# Annotate with gene name (optional) ----

# clusters.de$GENENAME <-
#   with(fData(sce.norm), gene_name[match(rownames(clusters.de), gene_id)])
#
# groups.de.d23$GENENAME <-
#   with(fData(sce.norm), gene_name[match(rownames(groups.de.d23), gene_id)])
#
# groups.de.lt2$GENENAME <-
#   with(fData(sce.norm), gene_name[match(rownames(groups.de.lt2), gene_id)])

# Prepare MLE data frame for ggplot ----

gg.mle <- data.frame(
  gene_id = c(
    rownames(groups.de.d23),
    rownames(groups.de.lt2)
  ),
  Group = rep(
    c("STM-D23580", "STM-LT2"),
    c(nrow(groups.de.d23), nrow(groups.de.lt2))
  ),
  MLE.group = c(
    groups.de.d23$mle,
    groups.de.lt2$mle
  ),
  MLE.cluster = clusters.de[
    c(
      rownames(groups.de.d23),
      rownames(groups.de.lt2)
    ),
    "mle"
  ],
  P.group = c(
    groups.de.d23$p.value,
    groups.de.lt2$p.value
  ),
  P.cluster = clusters.de[
    c(
      rownames(groups.de.d23),
      rownames(groups.de.lt2)
    ),
    "p.value"
  ],
  gene_name = with(
    fData(sce.norm),
    gene_name[
      match(c(
        rownames(groups.de.d23),
        rownames(groups.de.lt2)
      ),
      gene_id)
    ]
  )
)

ggplot(gg.mle) +
  geom_point(aes(MLE.cluster, -MLE.group, colour = Group), alpha = 0.2) +
  theme_minimal()

ggplot(gg.mle) +
  stat_density2d(
    aes(MLE.cluster, -MLE.group, fill = ..density..^0.25),
    geom = "tile", contour = FALSE, n = 200) +
  facet_grid(~ Group) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Prepare regression data for ggplot ----

d23.reg <- with(
  data.frame(
    D23580 = groups.de.d23$mle,
    Cluster1v2 = -clusters.de[rownames(groups.de.d23), 'mle']
  ),
  tidy(lm(D23580 ~ Cluster1v2))
)

lt2.reg <- with(
  data.frame(
    LT2 = groups.de.lt2$mle,
    Cluster1v2 = -clusters.de[rownames(groups.de.lt2), 'mle']
  ),
  tidy(lm(LT2 ~ Cluster1v2))
)

all.reg <- rbind(
  cbind(d23.reg, Group = "STM-D23580"),
  cbind(lt2.reg, Group = "STM-LT2")
)
all.reg$term <- rep(c("intercept", "slope"), 2)

gg.regression <- dcast(all.reg, Group ~ term, value.var = "estimate")
gg.regression$x <- 5
gg.regression$y <- -5:-6

ggplot(gg.mle) +
  # stat_density2d(
  #   aes(MLE.cluster, -MLE.group, fill = ..density..^0.25),
  #   geom = "tile", contour = FALSE, n = 200) +
  geom_point(aes(MLE.cluster, -MLE.group, colour = Group), alpha = 0.2) +
  geom_abline(
    aes(slope = slope, intercept = intercept, colour = Group),
    gg.regression
  ) +
  # facet_grid(~ Group) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Subset to genes DE between clusters ----

gg.sig.cluster <- subset(gg.mle, P.cluster < 0.01)
dim(gg.sig.cluster)

ggplot(gg.sig.cluster) +
  geom_point(aes(MLE.cluster, -MLE.group, colour = Group), alpha = 0.2) +
  theme_minimal()

# Regression on cluster DE genes only ----

genes.sig.cluster <- unique(as.character(gg.sig.cluster$gene_id))
length(genes.sig.cluster)

d23.reg <- with(
  data.frame(
    D23580 = groups.de.d23[genes.sig.cluster, "mle"],
    Cluster1v2 = -clusters.de[genes.sig.cluster, "mle"]
  ),
  tidy(lm(D23580 ~ Cluster1v2))
)

lt2.reg <- with(
  data.frame(
    LT2 = groups.de.lt2[genes.sig.cluster, "mle"],
    Cluster1v2 = -clusters.de[genes.sig.cluster, "mle"]
  ),
  tidy(lm(LT2 ~ Cluster1v2))
)

all.reg <- rbind(
  cbind(d23.reg, Group = "STM-D23580"),
  cbind(lt2.reg, Group = "STM-LT2")
)
all.reg$term <- rep(c("intercept", "slope"), 2)

gg.regression <- dcast(all.reg, Group ~ term, value.var = "estimate")
gg.regression$x <- 5
gg.regression$y <- -5:-6

# Correlation (Pearson) on genes DE between clusters ----

# cor.test(groups.de.d23[genes.sig.cluster, "mle"], -clusters.de[genes.sig.cluster, "mle"])
# cor.test(groups.de.lt2[genes.sig.cluster, "mle"], -clusters.de[genes.sig.cluster, "mle"])
#
# cor.test(groups.de.d23[genes.sig.cluster, "mle"], -clusters.de[genes.sig.cluster, "mle"])$estimate^2
# cor.test(groups.de.lt2[genes.sig.cluster, "mle"], -clusters.de[genes.sig.cluster, "mle"])$estimate^2

cor.d23 <- tidy(cor.test(-groups.de.d23[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"]))
cor.lt2 <- tidy(cor.test(-groups.de.lt2[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"]))

all.cor <- rbind(
  cbind(cor.d23, Group = "STM-D23580"),
  cbind(cor.lt2, Group = "STM-LT2")
)
all.cor$x <- -5
all.cor$y <- 7:6
all.cor$r2 <- sprintf("italic(R)^2~ %.2f", all.cor$estimate^2)

range.mle <- range(
  c(gg.sig.cluster$MLE.group, gg.sig.cluster$MLE.cluster)
)

gg <- ggplot(gg.sig.cluster) +
  geom_point(aes(-MLE.cluster, MLE.group, colour = Group), alpha = 0.25) +
  geom_abline(
    aes(slope = slope, intercept = intercept, colour = Group),
    gg.regression
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(x, y, label = r2, colour = Group), all.cor, parse = TRUE) +
  geom_text(aes(
    x, y,
    label = sprintf("slope: %.2f", slope), colour = Group),
    gg.regression
  ) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  theme_minimal() +
  scale_x_continuous(limits = range.mle) +
  scale_y_continuous(limits = range.mle) +
  scale_colour_manual(values = colours.infection) +
  labs(y = "MLE\n(Violet + / Exposed)", x = "MLE\n(Cluster 1 / Cluster 2)")

ggsave("figures/regression_4h_clusters_groups.pdf", height = 4, width = 5.5)

gg

# Subset D23 and LT2 to their DE genes ----

# gg.sig.groups <- subset(gg.mle, P.cluster < 0.01 & P.group < 0.01)
# dim(gg.sig.groups)
#
# head(gg.sig.groups)
# genes.sig.23 <- as.character(subset(
#   gg.sig.groups,
#   Group == "STM-D23580" & P.group < 0.01,
#   "gene_id", drop = TRUE))
#
# d23.reg <- tidy(lm(
#   -groups.de.d23[genes.sig.cluster, "mle"] ~ clusters.de[genes.sig.cluster, "mle"]
# ))
# lt2.reg <- tidy(lm(
#   -groups.de.lt2[genes.sig.cluster, "mle"] ~ clusters.de[genes.sig.cluster, "mle"]
# ))
#
# cor.test(-groups.de.d23[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"])
# cor.test(-groups.de.lt2[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"])
#
# cor.test(-groups.de.d23[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"])$estimate^2
# cor.test(-groups.de.lt2[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"])$estimate^2
#
# tidy(cor.test(-groups.de.d23[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"]))
# tidy(cor.test(-groups.de.lt2[genes.sig.cluster, "mle"], clusters.de[genes.sig.cluster, "mle"]))
#
# all.reg <- rbind(
#   cbind(d23.reg, Group = "STM-D23580"),
#   cbind(lt2.reg, Group = "STM-LT2")
# )
# all.reg$term <- rep(c("intercept", "slope"), 2)
#
# gg.regression <- dcast(all.reg, Group ~ term, value.var = "estimate")
#
# ggplot(gg.sig.cluster) +
#   # stat_density2d(
#   #   aes(MLE.cluster, -MLE.group, fill = ..density..^0.25),
#   #   geom = "tile", contour = FALSE, n = 200) +
#   geom_point(aes(MLE.cluster, -MLE.group, colour = Group), alpha = 0.2) +
#   geom_abline(
#     aes(slope = slope, intercept = intercept, colour = Group),
#     gg.regression
#   ) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   # facet_grid(~ Group) +
#   scale_fill_continuous(low = "white", high = "dodgerblue4") +
#   theme_minimal() +
#   theme(panel.grid = element_blank())
#
