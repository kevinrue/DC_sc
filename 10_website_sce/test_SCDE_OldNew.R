sce.endo <- readRDS("rds/sce.norm.rds")
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
    GENENAME = with(
      fData(sce.endo),
      gene_name[match(rownames(x), gene_id)]
    ),
    x
  )
  return(x)
}

old_scde_cluster <- readRDS("rds/13_cluster_scde.res_v2.rds")
names(old_scde_cluster)
# "4h_cluster2-cluster1" "6h_cluster2-cluster1"

old_scde_groups <- readRDS("rds/scde.res_v8.rds")
names(old_scde_groups)
# 4 contrasts

old4h.cluster <- addGENENAME(convert.z.score(old_scde_cluster[["4h_cluster2-cluster1"]]))
old4h.d23 <- addGENENAME(convert.z.score(old_scde_groups[["4h_STM-D23580_Violet +-4h_STM-D23580_Exposed"]]))
old4h.lt2 <- addGENENAME(convert.z.score(old_scde_groups[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]]))

venn::venn(
  x = list(
    d23 = rownames(old4h.d23)[old4h.d23$p.value < 0.01],
    cluster = rownames(old4h.cluster)[old4h.cluster$p.value < 0.01],
    lt2 = rownames(old4h.lt2)[old4h.lt2$p.value < 0.01]
  ), cexil = 2, cexsn = 2
)
