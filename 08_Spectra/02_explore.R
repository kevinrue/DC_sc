require(ggplot2)

library(RColorBrewer)
col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(2,4,7)]
# col.status <- brewer.pal(12, "Paired")[c(7,1,7)]
names(col.time) <- paste0(seq(2,6,2), "H")
names(col.infection) <- c("Ty2","STM-D23580","STM-LT2")
# names(col.status) <- levels(sce.norm$Status)

dirIn <- "~/Downloads/Spectra 2/"

filesIn <- list.files(dirIn)
length(filesIn)

# data1 <- read.table(file.path(dirIn, filesIn[1]))
# head(data1)
# plot(data1)

# filesIn
# length(grep("d23", filesIn))
# length(grep("lt2", filesIn))
# length(grep("ty2", filesIn))


read_spectrum <- function(file){
  rawData <- read.table(file)
  signal <- rawData
  return(signal)
}

read_intensity <- function(file){
  return(read_spectrum(file)[,2])
}

read_wavelength <- function(file){
  return(read_spectrum(file)[,1])
}

read_spectra <- function(files, dir = "."){
  # assumption: filename = (group)_inoculum_t=0_(replicate).txt
  # sample <- group_replicate
  rawList <- lapply(file.path(dir, files), "read_intensity")
  rawDataFrame <- data.frame(rawList)
  samples <- gsub("([[:alnum:]]+)_.*_([[:digit:]]+)\\.txt", "\\1_\\2", files)
  samples <- gsub("([[:alnum:]]+)_.*=0\\.txt", "\\1_0", samples)
  colnames(rawDataFrame) <- samples
  return(rawDataFrame)
}



# file.path(dirIn, filesIn)
# head(file.path(dirIn, filesIn))
# read_spectrum(head(file.path(dirIn, filesIn), 1))
# read_spectra(filesIn, "Spectra")
dataIn <- read_spectra(filesIn, dirIn)
dim(dataIn)

waveIn <- read_wavelength(file.path(dirIn, filesIn[1]))

dataWithWave <- cbind(
  Wavelength = waveIn,
  dataIn
)

data_300_1900 <- subset(dataWithWave, Wavelength > 300 & Wavelength < 1900)
dim(data_300_1900)

dataLong <- reshape2::melt(data_300_1900, id.vars="Wavelength")
dataLong <- dplyr::mutate(dataLong, group = gsub("_.*", "", variable))
dataLong <- dplyr::mutate(
  dataLong,
  Group = factor(
    group,
    levels = c("d23", "lt2", "ty2"),
    labels = c("STM-D23580", "STM-LT2", "Ty2")
  )
)
ggplot(dataLong) +
  geom_line(aes(Wavelength, value, colour = variable)) +
  guides(colour = "none") +
  theme_bw()

ggplot(dataLong) +
  geom_smooth(aes(Wavelength, value, colour = Group)) +
  scale_colour_manual(values = col.infection) +
  theme_bw()
ggsave("08_Spectra/02_smooth.pdf", width = 6, height = 4)

ggplot(dataLong) +
  geom_point(aes(Wavelength, value, colour = group), size=0.2) +
  theme_bw()

pca <- prcomp(t(dataIn[waveIn > 300 & waveIn < 1900,]), scale. = FALSE)
plot(pca)

pcaData <- data.frame(
  pca$x,
  group = gsub("_.*", "", rownames(pca$x))
)

pcaData <- dplyr::mutate(
  pcaData,
  Group = factor(
    group,
    levels = c("d23", "lt2", "ty2"),
    labels = c("STM-D23580", "STM-LT2", "Ty2")
  )
)

ggplot(pcaData) +
  geom_point(aes(PC1, PC2, colour = Group, shape = Group), size = 2) +
  scale_colour_manual(values = col.infection) +
  theme_bw()

ggsave("08_Spectra/02_PCA_unscaled.pdf", height = 4, width = 6)
