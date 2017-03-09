## Converts from gene name to chromosomal position
input.file <- read.delim(## your file here)
input.file$chr <- NA
input.file$pos <- NA

hg.reference <- read.delim("path/to/all_genes_hg38.txt", stringsAsFactors = F, comment.char = "", header = T)
reference <- hg.reference[, c(3, 5,6,13)]

for (i in 1:nrow(input.file)){
    current.name <- input.file$your_column_name_here[i]
    hits <- reference[reference$name2 == current.name, ]
    if (nrow(hits) == 0){
        warning(paste("Did not find gene ", current.name, " in database", sep = ""))
    }
    input.file$chr[i] <- hits$chrom[1]
    input.file$pos[i] <- hits$txStart[1]
}