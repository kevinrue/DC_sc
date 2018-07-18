
require(readxl)
require(org.Hs.eg.db)

googleDrive <- "~/Google Drive/DC_Salmonella_2018/4th Nat Comm revision/tables"
# list.files(googleDrive, "xlsx$")
googleOut <- "~/Google Drive/DC_Salmonella_2018/4th Nat Comm revision/GO_tables_with_genes/BH_05"
dir.create(googleOut)

# SCDE (24 tables) ----

list.files(googleDrive, "xlsx$")
scdeMapFile <- file.path(googleDrive, "1.Table_Kevin_SCDE_GO.xlsx")
scdeMapFile <- file.path(googleDrive, "1.Table_Kevin_SCDE_GO_v2.xlsx")
file.exists(scdeMapFile)

excel_sheets(scdeMapFile)

scdeMapTable <- read_xlsx(scdeMapFile, "Sheet1", col_names = c("scde", "go"))

for (index in seq_len(nrow(scdeMapTable))) {
    message("Processing index: ", index)
    scdeFile <- file.path(googleDrive, paste0(scdeMapTable$scde[index], ".csv"))
    file.exists(scdeFile)
    goFile <- file.path(googleDrive, paste0(scdeMapTable$go[index], ".csv"))
    file.exists(goFile)

    scdeTable <- read.csv(scdeFile, as.is = TRUE, row.names = 1)
    deGeneIds <- rownames(subset(scdeTable, p.value < 0.01))
    goTable <- read.csv(goFile, as.is = TRUE, row.names = 1)
    rownames(goTable) <- goTable$category
    #BH	bonferroni	Genes
    goTable$BH <- p.adjust(goTable$over_represented_pvalue, method = "BH")
    goTable$bonferroni <- p.adjust(goTable$over_represented_pvalue, method = "bonferroni")
    GOtoAnnotate <- subset(goTable, BH < 0.1, "category", drop=TRUE)
    message("Annotating BH < 0.05: ", length(GOtoAnnotate))
    for (tmp_go in GOtoAnnotate) {
            message("Processing # ", match(tmp_go, goTable$category))
            message("numDEInCat: ", goTable[tmp_go, "numDEInCat"], appendLF = FALSE)
            GO1geneInfo <- suppressMessages(
                select(org.Hs.eg.db, tmp_go, c("SYMBOL", "ENSEMBL"), "GOALL")
            )
            GO1geneInfoDE <- subset(GO1geneInfo, ENSEMBL %in% deGeneIds)
            tmpSymbolsDE <- unique(GO1geneInfoDE$SYMBOL)
            message(" | Found in org.Hs.eg.db: ", length(tmpSymbolsDE))
            goTable[tmp_go, "Genes"] <- paste(tmpSymbolsDE, collapse = ",")
    }
    write.csv(goTable, file.path(googleOut, basename(goFile)), row.names = FALSE)
}

# SCDE of clusters (24 tables) ----

scdeMapFile <- file.path(googleDrive, "2.Table_Kevin_Cluster_GO.xlsx")
file.exists(scdeMapFile)

excel_sheets(scdeMapFile)

scdeMapTable <- read_xlsx(scdeMapFile, "Sheet1", col_names = c("scde", "go"))

for (index in seq_len(nrow(scdeMapTable))) {
    message("Processing index: ", index)
    scdeFile <- file.path(googleDrive, paste0(scdeMapTable$scde[index], ".csv"))
    file.exists(scdeFile)
    goFile <- file.path(googleDrive, paste0(scdeMapTable$go[index], ".csv"))
    file.exists(goFile)

    scdeTable <- read.csv(scdeFile, as.is = TRUE, row.names = 1)
    deGeneIds <- rownames(subset(scdeTable, p.value < 0.01))
    goTable <- read.csv(goFile, as.is = TRUE, row.names = 1)
    rownames(goTable) <- goTable$category
    #BH	bonferroni	Genes
    goTable$BH <- p.adjust(goTable$over_represented_pvalue, method = "BH")
    goTable$bonferroni <- p.adjust(goTable$over_represented_pvalue, method = "bonferroni")
    GOtoAnnotate <- subset(goTable, BH < 0.1, "category", drop=TRUE)
    message("Annotating BH < 0.05: ", length(GOtoAnnotate))
    for (tmp_go in GOtoAnnotate) {
            message("Processing # ", match(tmp_go, goTable$category))
            message("numDEInCat: ", goTable[tmp_go, "numDEInCat"], appendLF = FALSE)
            GO1geneInfo <- suppressMessages(
                select(org.Hs.eg.db, tmp_go, c("SYMBOL", "ENSEMBL"), "GOALL")
            )
            GO1geneInfoDE <- subset(GO1geneInfo, ENSEMBL %in% deGeneIds)
            tmpSymbolsDE <- unique(GO1geneInfoDE$SYMBOL)
            message(" | Found in org.Hs.eg.db: ", length(tmpSymbolsDE))
            goTable[tmp_go, "Genes"] <- paste(tmpSymbolsDE, collapse = ",")
    }
    write.csv(goTable, file.path(googleOut, basename(goFile)), row.names = FALSE)
}
