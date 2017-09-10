stopifnot(
  require(readxl),
  require(reshape2),
  require(ggplot2),
  require(broom)
)

# Input files ----

excelFile <- "IL1b.xlsx"

# Load data ----

# Load 24h data
data_24h <- read_excel(
  path = "IL1b.xlsx", range = "R_24!B1:F16")
if (interactive()){
  View(data_24h)
}

# Load 48h data
data_48h <- read_excel(
  path = "~/git/DC_sc/09_IL1b_ELISA/IL1b.xlsx", range = "R_48!B1:F16")
if (interactive()){
  View(data_48h)
}

# Melt data ----

long_24h <- melt(data_24h, variable.name = "Donor", value.name = "IL1b_pg.mL")
if (interactive()){
  View(long_24h)
}

long_48h <- melt(data_48h, variable.name = "Donor", value.name = "IL1b_pg.mL")
if (interactive()){
  View(long_48h)
}

# Reshape data for t-test ----

t_24h <- dcast(
  subset(long_24h, Stimulus %in% paste("STM", c("LT2","D23580"), sep = "-")), formula = Concentration + Donor ~ Stimulus, value.var = "IL1b_pg.mL"
)
if (interactive()){
  View(t_24h)
}

t_48h <- dcast(
  subset(long_48h, Stimulus %in% paste("STM", c("LT2","D23580"), sep = "-")), formula = Concentration + Donor ~ Stimulus, value.var = "IL1b_pg.mL"
)
if (interactive()){
  View(t_48h)
}

# Test 24h and plot ----

with(t_24h, t.test(x = `STM-D23580`, y = `STM-LT2`, paired = TRUE))
with(t_24h, t.test(x = log(`STM-D23580` / `STM-LT2`)))

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_point(aes(Concentration, `STM-D23580`), colour = "blue") +
  geom_point(aes(Concentration, `STM-LT2`), colour = "green") +
  theme_bw()

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 1, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, `STM-D23580`/`STM-LT2`), colour = "red") +
  scale_y_continuous(limits=with(
    t_24h,
    c(0,ceiling(max(`STM-D23580`/`STM-LT2`)))
  )) +
  theme_bw()

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, log2(`STM-D23580`/`STM-LT2`)), colour = "red") +
  scale_y_continuous(limits=with(
    t_24h,
    ceiling(max(`STM-D23580`/`STM-LT2`)) * c(-1,1)
  )) +
  theme_bw()

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, log2(`STM-D23580`/`STM-LT2`)), colour = "red") +
  scale_y_continuous(limits=with(
    t_24h,
    max(ceiling(abs(log2(`STM-D23580`/`STM-LT2`))))*c(-1,1)
  )) +
  theme_bw()

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_density(aes(`STM-D23580`/`STM-LT2`)) +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_24h,
    max(ceiling(abs(`STM-D23580`/`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

ggplot(t_24h) +
  geom_density(aes(`STM-D23580`/`STM-LT2`, colour=Donor), size=1.5) +
  geom_density(aes(`STM-D23580`/`STM-LT2`, colour="Combined"), size=2, alpha=1/5, linetype="dashed") +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_24h,
    max(ceiling(abs(`STM-D23580`/`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

ggplot(t_24h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 1, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, `STM-D23580`-`STM-LT2`), colour = "red") +
  scale_y_continuous(limits=with(
    t_24h,
    max(ceiling(abs(`STM-D23580`-`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()
ggsave("24h_difference.pdf", height = 4, width = 10)

ggplot(t_24h) +
  geom_density(aes(`STM-D23580`-`STM-LT2`, colour=Donor), size=1.5) +
  geom_density(aes(`STM-D23580`-`STM-LT2`, colour="Combined"), size=2, alpha=1/5, linetype="dashed") +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_24h,
    max(ceiling(abs(`STM-D23580`-`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

# Test 48h and plot ----

with(t_48h, t.test(x = `STM-D23580`, y = `STM-LT2`, paired = TRUE))
with(t_48h, t.test(x = log(`STM-D23580` / `STM-LT2`)))

ggplot(t_48h) +
  facet_grid(~ Donor) +
  geom_point(aes(Concentration, `STM-D23580`), colour = "blue") +
  geom_point(aes(Concentration, `STM-LT2`), colour = "green") +
  theme_bw()

ggplot(t_48h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 1, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, `STM-D23580`/`STM-LT2`), colour = "red") +
  scale_y_continuous(limits=with(
    t_48h,
    c(0,ceiling(max(`STM-D23580`/`STM-LT2`)))
  )) +
  theme_bw()

ggplot(t_48h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, log2(`STM-D23580`/`STM-LT2`)), colour = "red") +
  scale_y_continuous(limits=with(
    t_48h,
    max(ceiling(abs(log2(`STM-D23580`/`STM-LT2`))))*c(-1,1)
  )) +
  theme_bw()

ggplot(t_48h) +
  facet_grid(~ Donor) +
  geom_density(aes(`STM-D23580`/`STM-LT2`)) +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_48h,
    max(ceiling(abs(`STM-D23580`/`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

ggplot(t_48h) +
  geom_density(aes(`STM-D23580`/`STM-LT2`, colour=Donor), size=1.5) +
  geom_density(aes(`STM-D23580`/`STM-LT2`, colour="Combined"), size=2, alpha=1/5, linetype="dashed") +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_48h,
    max(ceiling(abs(`STM-D23580`/`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

ggplot(t_48h) +
  facet_grid(~ Donor) +
  geom_hline(yintercept = 1, colour="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(Concentration, `STM-D23580`-`STM-LT2`), colour = "red") +
  scale_y_continuous(limits=with(
    t_48h,
    max(ceiling(abs(`STM-D23580`-`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()
ggsave("48h_difference.pdf", height = 4, width = 10)

ggplot(t_48h) +
  geom_density(aes(`STM-D23580`-`STM-LT2`, colour=Donor), size=1.5) +
  geom_density(aes(`STM-D23580`-`STM-LT2`, colour="Combined"), size=2, alpha=1/5, linetype="dashed") +
  geom_vline(xintercept = 0, colour="black", linetype="dashed", alpha=0.5) +
  scale_x_continuous(limits=with(
    t_48h,
    max(ceiling(abs(`STM-D23580`-`STM-LT2`)))*c(-1,1)*1.5
  )) +
  theme_bw()

# Summary of t-tests ----

tidy(
  with(t_24h, t.test(x = `STM-D23580`, y = `STM-LT2`, paired = TRUE))
)

tidy(
  with(t_48h, t.test(x = `STM-D23580`, y = `STM-LT2`, paired = TRUE))
)

# tidy(
#   with(t_24h, t.test(x = `STM-D23580` / `STM-LT2`))
# )
#
# tidy(
#   with(t_48h, t.test(x = `STM-D23580` / `STM-LT2`))
# )

# Summary plot ----
head(t_24h)

with(t_24h, `STM-D23580`-`STM-LT2`)

ggplot(t_24h) +
  geom_histogram(aes(`STM-D23580`-`STM-LT2`))

ggplot(t_24h) +
  geom_bar(aes(`STM-D23580`-`STM-LT2`))

ggplot(t_24h) +
  geom_density(aes(`STM-D23580`-`STM-LT2`))
