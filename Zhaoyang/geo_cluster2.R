
intable <- read.csv("GSE111543_info.csv", as.is = TRUE)

outtable <- intable[, c("geo_full_name", "Time", "Infection", "Status", "Treatment","DC1_2h", "DC2_2h", "Cluster_2h", "DC1_4h", "DC2_4h",  "Cluster_4h", "DC1_6h", "DC2_6h", "Cluster_6h")]

write.csv(x = outtable, file = "GSE111543_info_table.csv", row.names = FALSE)
