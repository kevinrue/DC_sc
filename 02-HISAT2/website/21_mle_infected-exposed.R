library(ggplot2)
library(dplyr)
library(RColorBrewer)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
scde.res <- readRDS("rds/scde.res_v8.rds")

# Colours ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

names(scde.res)

# feature.all <- featureNames(sce.norm)
#
# comp.df <- data.frame(
#   D23580 = c(
#     scde.res[["2h_STM-D23580_Violet +-2h_STM-D23580_Exposed"]][feature.2h, "mle"],
#
#   )
# )

# MLE ----

feature.2h <- rownames(scde.res[["2h_STM-LT2_Violet +-2h_STM-LT2_Exposed"]])

comp.2h <- data.frame(
  D23580 = scde.res[["2h_STM-D23580_Violet +-2h_STM-D23580_Exposed"]][feature.2h, "mle"],
  LT2 = scde.res[["2h_STM-LT2_Violet +-2h_STM-LT2_Exposed"]][feature.2h, "mle"],
  row.names = feature.2h
)

ggplot(comp.2h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "2h")

feature.4h <- rownames(scde.res[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]])

comp.4h <- data.frame(
  D23580 = scde.res[["4h_STM-D23580_Violet +-4h_STM-D23580_Exposed"]][feature.4h, "mle"],
  LT2 = scde.res[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]][feature.4h, "mle"],
  row.names = feature.4h
)

ggplot(comp.4h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "4h")

feature.6h <- rownames(scde.res[["6h_STM-LT2_Violet +-6h_STM-LT2_Exposed"]])

comp.6h <- data.frame(
  D23580 = scde.res[["6h_STM-D23580_Violet +-6h_STM-D23580_Exposed"]][feature.6h, "mle"],
  LT2 = scde.res[["6h_STM-LT2_Violet +-6h_STM-LT2_Exposed"]][feature.6h, "mle"],
  row.names = feature.6h
)

ggplot(comp.6h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "6h")

comp.df <- rbind(
  cbind(comp.2h, Time = "2h"),
  cbind(comp.4h, Time = "4h"),
  cbind(comp.6h, Time = "6h")
)

comp.df$product <- with(comp.df, D23580 * LT2)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_point(aes(colour = product), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(MLE[STM-D23580]%*%MLE[STM-LT2])
  )

ggsave("21_out/mle_violet-exposed_Time.pdf", height = 4, width = 12)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_density2d() +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(Z[STM-D23580]%*%Z[STM-LT2])
  )

ggsave("21_out/mle_violet-exposed_Time_geom_density2d.pdf", height = 4, width = 12)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(Z[STM-D23580]%*%Z[STM-LT2])
  )

ggsave("21_out/mle_violet-exposed_Time_hex.pdf", height = 4, width = 12)


# ggplot(comp.df) +
#   facet_grid(~ Time) +
#   geom_hex(aes(LT2, D23580)) +
#   theme_minimal() + scale_fill_continuous()

# ggplot(comp.2h) +
#   geom_density_2d(aes(LT2, D23580)) +
#   theme_minimal()

# Z-score ----

feature.2h <- rownames(scde.res[["2h_STM-LT2_Violet +-2h_STM-LT2_Exposed"]])

comp.2h <- data.frame(
  D23580 = scde.res[["2h_STM-D23580_Violet +-2h_STM-D23580_Exposed"]][feature.2h, "Z"],
  LT2 = scde.res[["2h_STM-LT2_Violet +-2h_STM-LT2_Exposed"]][feature.2h, "Z"],
  row.names = feature.2h
)

ggplot(comp.2h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "2h")

feature.4h <- rownames(scde.res[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]])

comp.4h <- data.frame(
  D23580 = scde.res[["4h_STM-D23580_Violet +-4h_STM-D23580_Exposed"]][feature.4h, "Z"],
  LT2 = scde.res[["4h_STM-LT2_Violet +-4h_STM-LT2_Exposed"]][feature.4h, "Z"],
  row.names = feature.4h
)

ggplot(comp.4h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "4h")

feature.6h <- rownames(scde.res[["6h_STM-LT2_Violet +-6h_STM-LT2_Exposed"]])

comp.6h <- data.frame(
  D23580 = scde.res[["6h_STM-D23580_Violet +-6h_STM-D23580_Exposed"]][feature.6h, "Z"],
  LT2 = scde.res[["6h_STM-LT2_Violet +-6h_STM-LT2_Exposed"]][feature.6h, "Z"],
  row.names = feature.6h
)

ggplot(comp.6h) +
  geom_point(aes(LT2, D23580)) +
  theme_minimal() +
  labs(title = "6h")

comp.df <- rbind(
  cbind(comp.2h, Time = "2h"),
  cbind(comp.4h, Time = "4h"),
  cbind(comp.6h, Time = "6h")
)

comp.df$product <- with(comp.df, D23580 * LT2)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_point(aes(colour = product), size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(Z[STM-D23580]%*%Z[STM-LT2])
  )

ggsave("21_out/Z_violet-exposed_Time.pdf", height = 4, width = 12)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_density2d() +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(Z[STM-D23580]%*%Z[STM-LT2])
  )

ggsave("21_out/Z_violet-exposed_Time_geom_density2d.pdf", height = 4, width = 12)

ggplot(comp.df, aes(LT2, D23580)) +
  facet_grid(~ Time) +
  coord_fixed() +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "STM-LT2",
    y = "STM-D23580",
    colour = expression(Z[STM-D23580]%*%Z[STM-LT2])
  )

ggsave("21_out/Z_violet-exposed_Time_hex.pdf", height = 4, width = 12)

# Count of DE genes ----
