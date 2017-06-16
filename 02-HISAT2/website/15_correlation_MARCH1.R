
library(scater)
library(broom)
sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

sce.endo <- sce.norm[!isSpike(sce.norm),]; dim(sce.endo)

sce.6h.d23.inf <- sce.endo[,with(
  pData(sce.endo), Time == "6h" & Infection == "D23580" & Status == "infected")
]

geneId <- subset(fData(sce.6h.d23.inf), gene_name == "MARCH1", "gene_id", drop = TRUE)

cor.test(norm_exprs(sce.6h.d23.inf)[geneId,], norm_exprs(sce.6h.d23.inf)[1,], method = "pearson")


MARCH1.cor <- do.call("rbind", apply(
  norm_exprs(sce.6h.d23.inf)[which(featureNames(sce.6h.d23.inf) != geneId),],
  1,
  function(x){
    tidy(cor.test(
      norm_exprs(sce.6h.d23.inf)[geneId,], # MARCH1
      x,
      method = "pearson"))
  }
))

MARCH1.cor <- cbind(
  gene_name = fData(sce.6h.d23.inf)[rownames(MARCH1.cor), "gene_name"],
  BH = p.adjust(MARCH1.cor$p.value, "BH"),
  MARCH1.cor
)

write.csv(MARCH1.cor, "15_out/MARCH1_correlation.csv", quote = FALSE)

geneId2 <- "ENSG00000138246" # DNAJC13
geneId2 <- "ENSG00000130208" # APOC1
geneId2 <- "ENSG00000124222" # STX16
geneId2 <- "ENSG00000162924" # REL

ggdata <- data.frame(
  t(norm_exprs(sce.endo)[c(geneId, geneId2),]),
  pData(sce.endo)[,c("Time","Treatment", "Infection", "Status")]
)

geneName1 <- fData(sce.endo)[geneId, "gene_name"]
geneName2 <- fData(sce.endo)[geneId2, "gene_name"]
colnames(ggdata)[1:2] <- c(geneName1, geneName2)

ggplot(ggdata, aes_string(geneName1, geneName2)) +
  geom_point() +
  geom_smooth(aes(fill = Status), method = "lm") +
  theme_minimal() +
  facet_grid(Status + Infection ~ Time)
