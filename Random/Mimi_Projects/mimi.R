## Mimi's oncopanel paper
gct.files <- read.delim("C:/Users/Noah/Onedrive/Work/temp/broad_values_by_arm.gct", stringsAsFactors = F)
gct.transformed <- read.delim("C:/Users/Noah/Onedrive/Work/temp/cghtransformed.txt", stringsAsFactors = F)


## to median center data
for (i in 3:nrow(gct.files)){
    #val <- median(as.numeric(gct.files[10, 3:148]))
    gct.files[i, 3:148] <- as.numeric(gct.files[i, 3:148]) + 10
}

write.table(gct.files, "C:/Users/Noah/OneDrive/Work/temp/transformed10.gct", sep = "\t", quote = F, row.names = F, col.names = F)
