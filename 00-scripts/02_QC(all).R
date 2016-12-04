source("settings.R")

library(scater)

sc_DC_gene <- readRDS(file.path(folder.rds, "sc_DC_gene.rds"))

# Filter detected features ----

library(S4Vectors)

# Keep only genes:
# - detected in >= 6 of 9 bulks
# - variance > 0
fr <- FilterRules(
    exprs = list(
        detected = function(x){
            rowSums(exprs(x[,x$Status == "BULK"]) > 0) >= 6
        },
        variance = function(x){
            apply(exprs(x), 1, function(x) {var(x) > 0})
        }
    )
)

fm <- evalSeparately(fr, sc_DC_gene)
summary(fm)

# Remove undetected features and zero variance
sc_DC_gene <- subsetByFilter(sc_DC_gene, fr)
dim(sc_DC_gene)

# Identify ERCC features ----

ercc.idx <- which(grepl("ERCC-[[:digit:]]", featureNames(sc_DC_gene)))
length(ercc.idx)

# Identify control cells ----

bulk.idx <- which(sc_DC_gene$Status == "BULK")
length(bulk.idx)
blank.idx <- which(sc_DC_gene$Status == "Blank")
length(blank.idx)

# Calculate QC metrics ----

sc_DC_gene <- calculateQCMetrics(
    sc_DC_gene,
    list(ERCC = ercc.idx),
    list(
        Bulk = bulk.idx,
        Blank = blank.idx
    )
)

saveRDS(sc_DC_gene, file.path(folder.rds, "sc_DC_gene_QC(384).rds"))

# scater_gui(sc_DC_gene)

# Plot SCESet ----

varLabels(sc_DC_gene)

p.SCESet <- plot(sc_DC_gene, block1 = "Status", block2 = "Infection",
                 colour_by = "Time", nfeatures = 300, exprs_values = "counts")
ggsave(
    file.path(folder.qc, "SCESet_Status-Infection-Time.pdf"),
    p.SCESet,
    width = 6, height = 4)
rm(p.SCESet)

# Plot of expression values ----

p.expr <- plotExpression(
    sc_DC_gene, rownames(sc_DC_gene)[
        sample(
            grep("CCL|CXCL|(IL^)", rownames(sc_DC_gene)),
            6
        )
        ],
    x = "Status", exprs_values = "exprs", colour = "Infection"
) + theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
)
ggsave(
    file.path(folder.qc, "Expression_Status-Time(6).pdf"),
    p.expr,
    width = 6, height = 4)

# PCs vs Well ----

phenoData(sc_DC_gene)
varLabels(sc_DC_gene)[1:10]

head(pData(sc_DC_gene))[,1:10]

p.pca <- plotPCA(
    sc_DC_gene,
    ncomponents = 3,
    colour_by = "Status",
    size_by = "Plate",
    shape_by = "Time"
) +
    theme(legend.position = "bottom")
ggsave(
    file.path(folder.qc, "PCA_Status-Plate-Time.pdf"),
    p.pca,
    width = 11,
    height = 12)

head(p.pca$data)
varLabels(sc_DC_gene)[1:10]

plate.exprData <- data.frame(
    Sample = colnames(sc_DC_gene),
    Row = factor(
        gsub("([[:upper:]*])[[:digit:]]*", "\\1", sc_DC_gene$Well),
        LETTERS[8:1]
    ),
    Column = factor(
        gsub("[[:upper:]*]*([[:digit:]*])", "\\1", sc_DC_gene$Well),
        as.character(1:12)
    ),
    Plate = sc_DC_gene$Plate,
    p.pca$data[colnames(sc_DC_gene),paste0("PC", 1:3)]
)
dim(plate.exprData)
table(duplicated(plate.exprData[,c("Row", "Column", "Plate")]))
head(plate.exprData)

p.platePC1 <- ggplot(plate.exprData, aes(x=Column, y=Row)) +
    geom_tile(aes(fill=PC1)) +
    facet_wrap(~Plate)
ggsave(
    file.path(folder.qc, "Plates_PC1.pdf"),
    p.platePC1,
    width = 12,
    height = 8)
p.platePC2 <- ggplot(plate.exprData, aes(x=Column, y=Row)) +
    geom_tile(aes(fill=PC2)) +
    facet_wrap(~Plate)
ggsave(
    file.path(folder.qc, "Plates_PC2.pdf"),
    p.platePC2,
    width = 12,
    height = 8)
p.platePC3 <- ggplot(plate.exprData, aes(x=Column, y=Row)) +
    geom_tile(aes(fill=PC3)) +
    facet_wrap(~Plate)
ggsave(
    file.path(folder.qc, "Plates_PC3.pdf"),
    p.platePC3,
    width = 12,
    height = 8)

# Detected features ----

# Cell QC metrics
varLabels(sc_DC_gene)
# Feature QC metrics
names(fData(sc_DC_gene))

table(sceset.detected.1$Status)

# because 9 Bulks (although we don't check which 9)
keep_feature <- rowSums(counts(sceset.detected.1) > 0) >= 9
sceset.detected.9 <- sceset.detected.1[keep_feature,]

saveRDS(sceset.detected.9, file.path(folder.rds, "sceset.detected.9.rds"))
# sceset.detected.9 <- readRDS(file.path(folder.rds, "sceset.detected.9.rds"))

rm(sceset.detected.1)

detected.9.featureCount <- nrow(sceset.detected.9)

detected.9.featureCount # 70525 transcripts left
raw.featureCount - detected.9.featureCount # 107703 transcripts removed
(raw.featureCount - detected.9.featureCount) / raw.featureCount # 60% removed

# Plot QC (highest-expression) ----

sceset.detected.9
get_exprs(sceset.detected.9, "counts")

qc.highest <- plotQC(
    sceset.detected.9, type = "highest-expression", exprs_values = "counts")
ggsave(
    file.path(folder.qc, "highest-expression.pdf"),
    qc.highest,
    width = 4,
    height = 8)

library(dplyr)
library(EnsDb.Hsapiens.v79) # Ensembl 79: Mar 2015 (GRCh38.p2)

write.table(
    data.frame(Name = mapIds(
        EnsDb.Hsapiens.v79,
        gsub( # Strip transcript version number
            "\\.[[:digit:]]*$",
            "",
            dplyr::arrange(
                fData(sceset.detected.9),
                desc(pct_total_counts)
            )[1:50,"feature_id"]),
        "GENENAME",
        "TXID"
    )),
    file.path(folder.qc, "highest-expression.txt"),
    sep = "\t",
    quote = FALSE, col.names = FALSE
)

# Plot QC (highest-expression) by Status ----

varLabels(sceset.detected.9)

p.DC <- plotQC(
    sceset.detected.9[, !sceset.detected.9$is_cell_control],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_DC.pdf"),
    p.DC,
    width = 4, height = 8
)

p.Bulk <- plotQC(
    sceset.detected.9[, sceset.detected.9$is_cell_control_Bulk],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_Bulks.pdf"),
    p.Bulk,
    width = 4, height = 8
)

p.Blank <- plotQC(
    sceset.detected.9[, sceset.detected.9$is_cell_control_Blank],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_Blanks.pdf"),
    p.Blank,
    width = 4, height = 8
)

# multiplot(p.DC, p.Bulk, p.Blank, cols = 3) # Lengthy, packed

p.uninfected <- plotQC(
    sceset.detected.9[, sceset.detected.9$Status == "uninfected"],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_Uninfected.pdf"),
    p.uninfected,
    width = 4, height = 8
)

p.exposed <- plotQC(
    sceset.detected.9[, sceset.detected.9$Status == "exposed"],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_exposed.pdf"),
    p.exposed,
    width = 4, height = 8
)

p.infected <- plotQC(
    sceset.detected.9[, sceset.detected.9$Status == "infected"],
    type = "highest-expression")
ggsave(
    file.path(folder.qc, "highest-expression_infected.pdf"),
    p.infected,
    width = 4, height = 8
)

p.dropout <- ggplot(pData(sceset.detected.9), aes(x=Status, y=pct_dropout)) +
    geom_boxplot(aes(fill=Status)) +
    geom_jitter(alpha = 0.3, width = 0.4, height = 0)
ggsave(
    file.path(folder.qc, "pct_dropout.pdf"),
    p.dropout,
    width = 6, height = 6)

# Plot QC (freq-vs-mean) ----

p.freqVsMean <- plotQC(sceset.detected.9, type = "exprs-freq-vs-mean")

ggsave(
    file.path(folder.qc, "Freq-vs-Mean.pdf"),
    p.freqVsMean,
    width = 6, height = 6)
