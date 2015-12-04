## Quickie For Zuhana

folder <- ("C:/Users/Noah/OneDrive/Work/temp")
files <- list.files(folder)

indels.1783 <- read.delim(paste(folder, files[1], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
indels.46 <- read.delim(paste(folder, files[3], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
snps.46 <- read.delim(paste(folder, files[4], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
snps.1783 <- read.delim(paste(folder, files[2], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)

exac.afs <- c(.000005, .00005, .0005, .005, .01, .05)

maf <- indels.46
count <- numeric()
for (i in 1:length(exac.afs)){
  filtered.version <- run.exac(maf, exac.afs[i])
  count <- c(count, sum(filtered.version$germline))
}

barplot(count, names.arg = exac.afs, main = "Indels in Sample 46 removed")
