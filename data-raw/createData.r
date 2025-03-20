## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
data <- readMat('data-raw/pt01epochdata.mat')
pt01EpochRaw <- data$a

## add channel names to the rows
goodChannels <- c(1:4, 7:36, 42:43, 46:69, 72:95)
sozChannels <- c(33:34, 62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01EpochRaw) <- channelNames$name[goodChannels]
sozIndex <- which(goodChannels %in% sozChannels == TRUE)
sozNames <- channelNames$name[sozChannels]

## Add time stamps to the columns
times <- -10 + (seq_len(ncol(pt01EpochRaw)) - 1) * 1e-3
colnames(pt01EpochRaw) <- sprintf("%+.3f", times)

pt01EcoG <- pt01EpochRaw[ , 9001L:12000L]
attr(pt01EcoG, "sozIndex") <- sozIndex
attr(pt01EcoG, "sozNames") <- sozNames
usethis::use_data(pt01EcoG, overwrite = TRUE)

## load fragility matrix
pt01Frag <- calcAdjFrag(epoch = pt01EcoG, window = 250, step = 125)
usethis::use_data(pt01Frag, overwrite = TRUE)


