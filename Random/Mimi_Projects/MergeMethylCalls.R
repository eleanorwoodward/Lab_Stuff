## Noah Greenwald

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

histones <- read.table("C:/Users/Noah/OneDrive/Work/Coding/short.txt", stringsAsFactors = F)
colnames(histones) <- c("PeakID", "chr", "start", "end", "strand", "normalized.tag.count", "region size", "findpeaks.score", "clonal.fold.change")

rbind(existingDF[1:r,],newrow,existingDF[-(1:r),])