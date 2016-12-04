
stopifnot(
    requireNamespace("dplyr"),
    requireNamespace("EnsDb.Hsapiens.v79")
)

source("settings.R")

library(scater)

# Create necessary folders ----

sapply(
    c(folder.rds, folder.qc),
    function(x){
    if (!dir.exists(x)){
        dir.create(x)
    } else {
        return(FALSE)
    }
})

# Import data ----

kallisto.dirs <- list.files(folder.quant)

sc_DC <- readKallistoResults(
    samples = gsub("WTCHG_", "", kallisto.dirs),
    directories = file.path(folder.quant, kallisto.dirs)
)

# saveRDS(sc_DC, file.path(folder.rds, "sc.raw.rds"))
# Note: 32.2 MB in memory / 297.2 MB on disk
# sc.raw <- readRDS(file.path(folder.rds, "sc.raw.rds"))

samplesData <- read.csv(
    file.path(folder.expData, "samples.csv"),
    row.names = 1)

head(pData(sc_DC))
pData(sc_DC) <- cbind(
    samplesData,
    pData(sc_DC)
)

rm(samplesData)

# Collapse transcript to genes ----

fData(sc_DC)$feature_id <- gsub(
    "\\.[[:digit:]]$", "", fData(sc_DC)$feature_id
)

sc_DC <- getBMFeatureAnnos(
    sc_DC,
    attributes = c(
        "ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol",
        "chromosome_name", "transcript_biotype", "transcript_start",
        "transcript_end", "transcript_count"),
    feature_symbol = "hgnc_symbol",
    dataset = "hsapiens_gene_ensembl")

## Summarise expression at the gene level
sc_DC_gene <- summariseExprsAcrossFeatures(
    sc_DC, exprs_values="tpm", summarise_by="feature_id")

rm(sc_DC)

# Add gene symbols to featureNames to make them more intuitive ----

# Replace by
colnames(fData(sc_DC_gene)) <-
    gsub("exprs_collapsed_to", "ensembl_gene_id", colnames(fData(sc_DC_gene)))

library(biomaRt)
mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org"
)
bm.map <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart = mart)

fData(sc_DC_gene)$hgnc_symbol <- bm.map[
    match(fData(sc_DC_gene)$ensembl_gene_id, bm.map$ensembl_gene_id),
    "hgnc_symbol"
]
fData(sc_DC_gene)$hgnc_symbol[fData(sc_DC_gene)$hgnc_symbol == ""] <- NA
sum(is.na(fData(sc_DC_gene)$hgnc_symbol))
head(fData(sc_DC_gene))

rm(bm.map, mart)

## Add gene symbols to featureNames to make them more intuitive

table(grepl("^ENSG", featureNames(sc_DC_gene)))

featureNames(sc_DC_gene) <- ifelse(
    is.na(fData(sc_DC_gene)$hgnc_symbol),
    as.character(fData(sc_DC_gene)$ensembl_gene_id),
    paste(
        fData(sc_DC_gene)$hgnc_symbol,
        fData(sc_DC_gene)$ensembl_gene_id,
        sep = "_"
    )
)

saveRDS(sc_DC_gene, file.path(folder.rds, "sc_DC_gene.rds"))
