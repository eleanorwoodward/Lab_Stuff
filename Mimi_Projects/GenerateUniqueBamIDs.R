## returns a list of unique files that match a given pattern
files <- read.delim("/xchip/beroukhimlab02/mimi/EGAD00001000807/files.txt", stringsAsFactors = F, header = F)

for (i in 1:nrow(files)){
    files[i,1] <- strsplit(files[i,1], "EGAR[0-9]*_")[[1]][2]
    files[i,1] <- strsplit(files[i, 1], "-")[[1]][1]
    write.table(unique(files[, 1]), "/xchip/beroukhimlab02/mimi/EGAD00001000807/unique.txt", row.names = F, quote = F, col.names = F)
    
}
