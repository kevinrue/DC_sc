
rds <- readRDS("rds/o.ifm_2h.rds")
head(rds)

oldRds <- readRDS("~/git/DC_sc/02-HISAT2/website/rds/o.ifm.2h_v8.rds")
head(oldRds)

# Check that names match (not necessarily in order)
all(sort(rownames(rds)) == sort(paste("Cell", rownames(oldRds), sep = "_")))
# Check that names do _not_ match (in order)
all(rownames(rds) == paste("Cell", rownames(oldRds), sep = "_"))
# Check that values match (in order)
all(rds == oldRds)

# What I think happened here is that I messed up as follows:
# I took the oldRds file
# I renamed the rows taking the rownames of the SCE object,
# which don't match the order of the SCDE error model,
# as the latter reorders cells by grouping level

# Confirmed by running scde.error.models again on the newly named data

# Fix:
# 1) either run scde.error.models on the three time points
# 2) do a proper renaming of the old error models

for (tp in paste0(seq(2,6,2), "h")){
  oldRds <- readRDS(sprintf("~/git/DC_sc/02-HISAT2/website/rds/o.ifm.%s_v8.rds", tp))
  rownames(oldRds) <- paste("Cell", rownames(oldRds), sep = "_")
  saveRDS(oldRds, sprintf("~/git/DC_sc/10_website_sce/rds/o.ifm.%s.rds", tp))
}
