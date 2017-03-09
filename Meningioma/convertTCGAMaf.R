x <- read.delim("/xchip/cga_home/galengao/finishedProducts/ABSOLUTE/MAFS/ACC.abs_mafs_2colfixed.txt", stringsAsFactors = F)
newfile <- x[, c("Hugo_Symbol", "Chromosome", "Start_position", "Variant_Classification", "ccf_hat")]
newfile <- newfile[newfile$Variant_Classification %in% c("Missense_Mutation", "Silent", "Splice_Site", "Frame_Shift_Del", "In_Frame_Del", "Nonsense_Mutation", "Frame_Shift_Ins", 
                                                         "Nonstop_Mutation", "Translation_Start_Site")]
saveRDS(newfile, "/xchip/beroukhimlab/noah/Meningioma/absolute_small.rds")