require(ggplot2)

dirIn <- "~/Downloads/Spectra"

filesIn <- list.files(dirIn)

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
read_spectrum(head(file.path(dirIn, filesIn), 1))
# read_spectra(filesIn, "Spectra")
dataIn <- read_spectra(filesIn, "Spectra")

dataWithWave <- cbind(
  Wavelength = read_wavelength(file.path(dirIn, filesIn[1])),
  dataIn
)

dataLong <- reshape2::melt(dataWithWave, id.vars="Wavelength")
dataLong <- dplyr::mutate(dataLong, group = gsub("_.*", "", variable))
ggplot(dataLong) +
  geom_line(aes(Wavelength, value, colour = variable)) +
  guides(colour = "none")

ggplot(dataLong) +
  geom_smooth(aes(Wavelength, value, colour = group)) +
  theme_bw()
ggsave("01_smooth.pdf", width = 6, height = 4)

ggplot(dataLong) +
  geom_point(aes(Wavelength, value, colour = group), size=0.2) +
  theme_bw()

pca <- prcomp(t(dataIn), scale. = FALSE)
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
  theme_bw()

ggsave("01_PCA_unscaled.pdf", height = 4, width = 6)
