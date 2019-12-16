source("wtchg_geo.R")

sce.endo <- readRDS("../12_sce_infected/rds/20180707_sce.endo_clusters.rds")

geo_data <- data.frame(
    geo_full_name = wtchg_geo$geo,
    cell_index = gsub(".*_([[:digit:]]{6})", "\\1", wtchg_geo$geo),
    stringsAsFactors = FALSE
)

# Add metadata ----

sce.endo_data <- data.frame(
    cell_index = gsub("Cell_", "", colnames(sce.endo)),
    colData(sce.endo)[, c("Time", "Infection", "Status", "Treatment")],
    stringsAsFactors = FALSE
)

sce.endo_data <- merge(geo_data, sce.endo_data, by="cell_index")

# Add Diffusion map ----

fig1_data <- read.csv("../12_sce_infected/15b_out/Fig1a_data.csv", as.is = TRUE)

sce.endo_data <- merge(sce.endo_data, fig1_data[, c("DC1", "DC2", "cell_index")], by="cell_index")

# Add Diffusion map ----

fig2_data <- read.csv("../12_sce_infected/16b_chinese_email_out/Fig2a_data.csv", as.is = TRUE)

sce.endo_data <- merge(sce.endo_data, fig2_data, by="cell_index")

sce.endo_data <- sce.endo_data[, c("geo_full_name", "Time", "Infection", "Status",  "Treatment", "DC1", "DC2", "DC1_2h", "DC2_2h",  "Cluster_2h", "DC1_4h", "DC2_4h", "Cluster_4h", "DC1_6h", "DC2_6h",  "Cluster_6h")]
write.csv(x = sce.endo_data, file = "GSE111543_info.csv", row.names = FALSE)

# Sanity checks ----

ggplot(sce.endo_data) +
    geom_point(aes(DC1_2h, DC2_2h, color=as.factor(Cluster_2h)))

ggplot(sce.endo_data) +
    geom_point(aes(DC1_4h, DC2_4h, color=as.factor(Cluster_4h)))

ggplot(sce.endo_data) +
    geom_point(aes(DC1_6h, DC2_6h, color=as.factor(Cluster_6h)))

ggplot(sce.endo_data) +
    geom_point(aes(DC1, DC2, color=as.factor(Treatment))) +
    scale_color_brewer(palette = "Paired")
