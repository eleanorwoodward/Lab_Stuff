## Noah Greenwald
## 9/17/15

## Helpful R functions for manipulating MAF output files for data analysis. 
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering/ExacFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//ESPFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//PoNFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//CallCompare.R")

source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("devtools")
# install.packages("R.utils")
# install.packages("NMF")
#install.packages("kohonen")
# install.packages("nycflights13")
#install.packages("randomcoloR")
#install.packages("NMF")
#install.packages("gridExtra")
#install.packages("cowplot")

## Nice packages
#library(grDevices)
library(RColorBrewer)
library(dplyr)
#library(rtracklayer)
library(ggplot2)
#library(curl)
library(R.utils)
#library(Biobase)
#library(kohonen)
#library(nycflights13)
library(reshape2)
library(randomcoloR)
library(NMF)
library(gridExtra)

## ggplot settings
rameen_theme <- theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), axis.line.x = element_line(size = 1.25, color = "black"), axis.line.y = element_line(size = 1.25, color = "black"), 
                      axis.text.x = element_text(angle = 90, hjust = 1))

FilterMaf <- function(maf, list, column.name, keep = TRUE){
  # Takes a maf, and returns a maf with only those rows who have a value in given column that matches an entry in the list. 
  
  # Args:
  #     maf: a maf
  #     list: list of column values to keep
  #     column.name: column in the maf which will be checked
  #     keep: optional, logical value. When false, discards selected rows
  
  # Returns: A maf matching input criteria
  
  ## Error Checking
  if (is.null(maf[[column.name]])){
    stop("Invalid column name")
  }
  
  attendance <- list %in% maf[[column.name]]
  if (sum(attendance) < length(attendance)){
    warning("List element ", list[which(!attendance)], " not found")
  }
  
  if (keep == TRUE){
    idx <- which(maf[[column.name]] %in% list)
    new.maf <- maf[idx, ]
    return (new.maf)
  }else{
    idx <- which(!(maf[[column.name]] %in% list))
    new.maf <- maf[idx, ]
    return (new.maf)
  }
}

PlotMaf <- function(maf, column.name, percent = 100, lower.margin = 1, title = ""){
    ## Takes a Maf file and produces a bar plot 
    
    # Args: 
    #     Maf: a maf file in array format. 
    #     column.name: a string with the name of the column that the barplot will be produced from. 
    #     lower.margin: optional, a multiplier that adjusts the marign at the bottom of the graph. 
    #     title: optional, gives the graph a title
    #     percent: optional, allows for plotting of only top x% of hits ranked by frequency
    
    ## Error checking
    if (is.null(maf[[column.name]])){
        stop("Invalid column name")
    }
    
    if (percent <= 0 | percent > 100){
        stop("Learn how to calculate a percent")
    }
    data <- maf[, column.name]
    tbl <- table(data)
    tbl <- sort(tbl, decreasing = T)
    tbl <- tbl[1:floor(length(tbl)*percent/100)]
    par(mar= c(5*lower.margin,4,4,2) + .1)
    size = (length(unique(data)))
    barplot(tbl, las = 2, main = title)
    
    #   if (size > 60){
    #       x = ceiling(size / 60)
    #       for (i in 1:x){
    #         end
    #         barplot(tbl, las = 2, main = c(title," ", i, " of ", x))
    #       }
    #   }
    #   else{
    #   barplot(tbl, las = 2, main = title)
    #   }
}

AverageMaf <- function(maf, list, list.column, int.column){
  # Takes a maf and produces a graph with the average value of a number of categories within the maf
  
  # Args:
  #     maf: a maf file in data.table format
  #     list: a list of the groupings the function will use to calcuate averages
  #     list.column: column of the maf which the list is drawn from
  #     int.column: column whose values will be summed for each entry
  
  # Returns: a list of averages corresponding to each item in input list, values of which are found in int.column
  
  ## Error Checking
  if (is.null(maf[[list.column]]) | is.null(maf[[int.column]])){
    stop("Invalid column name")
  }
  
  attendance <- list %in% maf[[list.column]]
  if (sum(attendance) < length(attendance)){
    warning("List element ", list[which(!attendance)], " not found")
  }
  
  ## Creates a vector to store average for each entry
  values <- c()
  
  ## loop through list, computing average for each entry and storing in vector
  for(i in 1:length(list)){
    temp.frame <- maf[which(maf[[list.column]] == list[i]), ]
    x <-  mean(temp.frame[[int.column]])
    values[i] <- x
  }
  return(values)
}


ReccurentMaf <- function(maf, column.name, threshold = 1){
  # Takes a maf file, and returns a filtered maf with only those entries in designated column with
  # a frequency greater than the threshold value
  
  # Args:
  #     maf: a maf file
  #     column.name: the designated column in maf file
  #     threshold: cutoff value for # of times column entry must be seen to be included
  
  # Returns: a filtered maf
  
  ## Error Checking
  if (is.null(maf[[column.name]])){
    stop("Invalid column name")
  }
  
  tbl <- table(maf[[column.name]])
  approved.list <- c()
  for (i in 1:length(tbl)){
    if (tbl[[i]] > threshold){
      approved.list <- c(approved.list, names(tbl[i]))
    }
  }
  FilterMaf(maf, approved.list, column.name)
}

RatioMaf <- function (maf, column.name, list1, list2){
  ## Takes a maf file, and returns the number of occurences of variables found in list 1 and 2
  
  ## Args:
  #       maf: a maf file
  #       column.name: the column in maf file that list elements are drawn from
  #       list1: list of possible values for given column
  #       list2: second list of possible values for given column
  
  # Returns: summary of occurences of selected elements in each list. 
  
  
  temp1 <- FilterMaf(maf, list1, column.name)
  temp2 <- FilterMaf(maf, list2, column.name)
  val1 <- nrow(temp1)
  val2 <- nrow(temp2)
  valtot <- nrow(maf)
  paste("There are a total of ", valtot, "mutations in your maf. Of these, ", val1, 
        " were in the first subset, and ", val2, " were in the second subset.")
}

PerSampleMaf <- function(maf, column.name, identifier.column = "Tumor_Sample_Barcode", list = maf[[column.name]]){
  ## Takes a maf, and returns a filtered maf that limits occurence to one per sample
  
  ## Args:
  #       maf: a maf file
  #       column.name: name of the column to be checked for uniqueness
  #       identifier.column: name of the column with sample IDs            
  #       list: values that will be kept in filtered maf. Default is all values in column
  
  # Returns: a maf with only those rows which contain an element from list in selected column, 
  # limited to once per sample. 
  filtered <- FilterMaf(maf, list, column.name, TRUE)
  cols <- c(identifier.column, column.name)
  dups <- filtered[, cols]
  idx <- duplicated(dups)
  filtered <- filtered[!idx, ]
}            

CompareMutsMaf <- function(maf, cutoff = 0, title = ""){
  ## Takes a maf and a minimum value for # of mutations, and returns a comparison chart of total 
  ## number of mutations in given genes, as well as # of samples with given mutation
  
  ## Args:
  #       maf: a maf file
  #       cutoff: minimum number of times a mutation is seen in dataset in order to be kept
  
  # Returns: a side-by-side barplot comparing frequencies of total vs per sample mutations per gene
  
  filtered <- ReccurentMaf(maf, "Hugo_Symbol", cutoff)
  per.gene <- table(filtered[["Hugo_Symbol"]])
  per.sample.maf <- PerSampleMaf(filtered, "Hugo_Symbol")
  per.sample <- table(per.sample.maf[["Hugo_Symbol"]])
  data <- rbind(per.gene, per.sample)
  par(mar=c(10, 4, 4, 2) + .1)
  barplot(data, beside = TRUE, legend.text = c("Total mutations per gene", "Samples with gene mutated"), 
          las = 2, main = title )
}

## Compares different cutoffs for filtering


FilterCutoffMaf <- function(maf1, maf2, cutoff, cut1 = "Filter 1", cut2 = "Filter2",  title = ""){
  maf1 <- PerSampleMaf(maf1, "Hugo_Symbol")
  maf1 <- ReccurentMaf(maf1, "Hugo_Symbol", cutoff)
  maf2 <- PerSampleMaf(maf2, "Hugo_Symbol")
  maf2 <- ReccurentMaf(maf2, "Hugo_Symbol", cutoff)
  data1 <- table(maf1[["Hugo_Symbol"]])
  data2 <- table(maf2[["Hugo_Symbol"]])
  data <- EqualizeTable(data1, data2)
  barplot(data, beside = TRUE, legend.text = c(cut1, cut2),
          las = 2, main = title)
}


## Takes two tables of values, makes sure all values from each table appear in the other, and orders the table
## so that all entries are in same position. 

EqualizeTable <- function(table1, table2){
  
  ## Takes all entries not found in table 1, sets them to zero value
  table1.names <- names(table1)
  table2.names <- names(table2)
  free <- !(table2.names %in% table1.names)
  vals <- c(rep(0, sum(free)))
  names(vals) <- table2.names[free]
  table1 <- c(table1, vals)
  table1 <- table1[order(names(table1))]
  
  ## Does the same thing, vice aversa
  free <- !(table1.names %in% table2.names)
  vals <- c(rep(0, sum(free)))
  names(vals) <- table1.names[free]
  table2 <- c(table2, vals)
  table2 <- table2[order(names(table2))]

  table3 <- rbind(table1, table2)
  table3
}


CleanYourRoom <- function(input.data, naughty.list = c()){
## Takes a data frame or matrix, and removes rows with problemtic entries

## Args:
##      input.data: the data to be input
##      naughty.list: additional parameters to exclude

    for(i in 1:ncol(input.data)){
        bad.boys <- input.data[, i] %in% c("", "-", "N/A", "na", NA, naughty.list)
        input.data <- input.data[!bad.boys, ]
    }
    return (input.data)
}

TestStatistic <- function(vector, mu){
## takes a vector of values, determines if average is significantly different from input
    
## Args:
##      vector: vector of numeric values
##      mu: value to compare vector values to    
    deviation <- sd(vector)
    average <- mean(vector)
    denom <- deviation / sqrt(length(vector))
    return((mu - average) / denom)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    ## checks if a number is an integer
    ## taken from Integer help page documentation examples
    abs(x - round(x)) < tol
}

PlotComut <- function(mut.maf1 = NULL, mut.maf2 = NULL, focal.cnv.maf = NULL, broad.cnv.maf = NULL, samples, input.samples, input.genes = "Hugo_Symbol", input.mutations = "Variant_classification", 
                      gene.cutoff = 1, file.path = NULL, unique.muts = NULL, col.vector = NULL, phenotypes = NULL, manual.order = NULL, fixed.order = NULL, return.matrix = FALSE, 
                      return.plot = FALSE, y.axis.font.size = 8, legend.font.size = 8, dimensions = c(7, 7), title = "Comut Plot"){
## Takes in a maf file, and generates a comut plot of common mutations in each of the samples
## adapted from code at https://benchtobioinformatics.wordpress.com/2015/05/25/how-to-make-a-co-mutation-plot/
    

## Args:
    ## mut.maf1: a maf file with mutation calls
    ## mut.maf2: an optional maf file with mutation calls, all of which will sorted below the first and labeled tier 4
    ## focal.cnv.maf: an optional maf (long format) of CNV calls with sample identifier as first column, gene as second, and type of CNV as third
    ## broad.cnv.maf: an optional maf
    ## samples: a vector with all individuals to be plotted (don't take from MAF, could be absent due to no mutations)
    ## input.samples: a column wthin the maf which specifies samples
    ## input.genes: a column within the maf which specifies genes
    ## input.mutations: a column within the maf describing the type of mutation
    ## gene.cutoff: a threshold for deciding how many genes to include in comut. if beteen 0 and 1, treated as percentage of total. otherwise, number of genes
    ## file.path: optional path specifying where the comut will be saved
    ## unique.muts: an optional vector with all the unique mutation classes present in dataset. Useful for consistent colors across multiple comuts
    ## col.vector
    ## phenotypes: an optional dataframe with samples as the first column, and subsequent columns taken as additional information to plot in given order
    ## manual.order: a vector of variables names which will overide the default ordering of rows by frequency for alternate sorting
    ## fixed.order: a vector of sample identifiers which determines ordering of all samples
    ## return.matrix: if TRUE, returns the wide dataframe used to sort the comut to enable statistical testing
    ## title: a title for the plot
    
    
    ## Error Checking
    if(!missing(mut.maf1)){
        if (is.null(mut.maf1[[input.samples]]) |is.null(mut.maf1[[input.genes]]) | is.null(mut.maf1[[input.mutations]])){
            stop("Invalid column name for mut.maf")
        }
        
        if (nrow(mut.maf1) == 0){
            stop("Invalid mut.maf, doesn't contain any entries")
        }
    }
    
    if(!missing(mut.maf2)){
        if (is.null(mut.maf2[[input.samples]]) |is.null(mut.maf2[[input.genes]]) | is.null(mut.maf2[[input.mutations]])){
            stop("Invalid column name for mut.maf2")
        }   
        
        if(nrow(mut.maf2) == 0){
            stop("Invalid mut.maf2, doesn't contain any entries")
        }
    }
    
    
    if (gene.cutoff < 0){
        stop("Invalid value for gene cutoff: must be an integer greater than 1 or a real number between 0 and 1")
    }
    
    if (gene.cutoff > 1){
        if(!is.wholenumber(gene.cutoff)){
            stop("Invalid value for gene cutoff: must be an integer greater than 1 or a real number between 0 and 1")
        }
    }
    
    if (!is.vector(samples)){
        warning("samples not passed in as vector, converting first column")
        samples <- samples[, 1]
    }
    
    if (!missing(fixed.order)){
        if (sum(!(fixed.order %in% samples) != 0)){
            stop("Some samples in fixed ordering variable aren't present in samples list")
        }
    }
    
    ## initailize long dataframe with all possible gene-sample combinations
    if(!missing(mut.maf1)){
        all.samples <- unique(samples)
        all.genes <- unique(mut.maf1[, input.genes])
        df <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
        colnames(df) <- c("samples", "genes")
        
        ## optimized code: use merge + melt to generate df with all called mutations
        print("Initializing mut.maf1")
        mut.maf1.long <- melt(mut.maf1[, c(input.samples, input.genes, input.mutations)], id = c(input.samples, input.genes))
        mut.maf1.long <- mut.maf1.long[-3]
        colnames(mut.maf1.long) <- c("samples", "genes", "mutations")
        df <- merge(df, mut.maf1.long, c("samples", "genes"), all.x = TRUE)
        ## annotate coverage information
        print("annotating coverage")
        for (i in 1:length(unique(df$samples))){
            tmp.sample <- unique(df$samples)[i]
            tmp.sample.ver <- master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]
            
            ##oncomap coverage stored separately till we figure out how to display
            ## checks to see how many genes for each sample were not covered in that assay, labels them as such
            if (!tmp.sample.ver){
                ghosts <- df[df$samples == tmp.sample & df$genes %in% not.covered.map, ]
                if (nrow(ghosts) > 0){
                    df[df$samples == tmp.sample & df$genes %in% not.covered.map, ]$mutations <- "nc"
                }else{
                    ## do nothing
                }
            }else{
                ghosts <- df[df$samples == tmp.sample & df$genes %in% not.covered[[tmp.sample.ver]], ]
                if (nrow(ghosts) > 0){
                    df[df$samples == tmp.sample & df$genes %in% not.covered[[tmp.sample.ver]], ]$mutations <- "nc"
                }else{
                    ## do nothing
                }
            }
        }
        
        ## orders rows based on gene most frequently altered genes, removing NAs and uncovered genes
        df_sub <- subset(df, !is.na(df$mutations))
        df_sub <- subset(df, df$mutations != "nc")
        ord <- names(sort(table(df_sub$genes), decreasing = T))
        df.combined <- df
    }
    ## checks if second maf of mutations supplied. If so, duplicate above code to produce second df
    if(!missing(mut.maf2)){
        ## iniatialize df.2
        all.samples <- unique(samples)
        all.genes <- unique(mut.maf2[, input.genes])
        df.2 <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
        colnames(df.2) <- c("samples", "genes")
        
        print("Initializing mut.maf2")
        mut.maf2.long <- melt(mut.maf2[, c(input.samples, input.genes, input.mutations)], id = c(input.samples, input.genes))
        mut.maf2.long <- mut.maf2.long[-3]
        colnames(mut.maf2.long) <- c("samples", "genes", "mutations")
        df.2 <- merge(df.2, mut.maf2.long, c("samples", "genes"), all.x = TRUE)
        ## annotate coverage information
        print("annotating coverage")
        for (i in 1:length(unique(df.2$samples))){
            tmp.sample <- unique(df.2$samples)[i]
            tmp.sample.ver <- master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]
            
            ##oncomap coverage stored separately till we figure out how to display
            ## checks to see how many genes for each sample were not covered in that assay, labels them as such
            if (!tmp.sample.ver){
                ghosts <- df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.map, ]
                if (nrow(ghosts) > 0){
                    df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.map, ]$mutations <- "nc"
                }else{
                    ## do nothing
                }
            }else{
                ghosts <- df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered[[tmp.sample.ver]], ]
                if (nrow(ghosts) > 0){
                    df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered[[tmp.sample.ver]], ]$mutations <- "nc"
                }else{
                    ## do nothing
                }
            }
        }
        df_sub <- subset(df.2, !is.na(df.2$mutations))
        df_sub <- subset(df.2, df.2$mutations != "nc")
        
        ## generate empty row to separate tier 4 from rest
        df.3 <- expand.grid(all.samples, "Tier4")
        df.3$mutations <- "nc"
        colnames(df.3) <- c("samples", "genes", "mutations")
        
        ord.2 <- names(sort(table(df_sub$genes), decreasing = T))
        
        ## checks if first supplied argument or not, either updates or creates tracking variables as appropriate
        if(exists("df.combined")){
            df.combined <- rbind(df.combined, df.3, df.2)
            ord <- c(ord, "Tier4", ord.2)
        }else{
            df.combined <- df.2
            ord <-ord.2
        }
        
    }
    
    ## keep given % of genes based on cutoff if mutations data supplied
    if(!missing(mut.maf1) | !missing(mut.maf2)){
        if (gene.cutoff == 1){
            ## keep all genes
        }else if (gene.cutoff < 1){
            ## keep given percentage of top hits
            idx <- ceiling(length(ord) * gene.cutoff)
            ord <- ord[1:idx]
        }else{
            ## keep given number of top hits
            idx <- min(gene.cutoff, length(ord))
            ord <- ord[1:idx]
        }
    }
    
    ## checks to see if focal.cnv information supplied to comut
    if (!missing(focal.cnv.maf)){
        if (!is.na(focal.cnv.maf)){
            all.samples <- unique(samples)
            ## generates df with all possible cnvs
            all.genes <- unique(focal.cnv.maf$GENE)
            df.cnv <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
            colnames(df.cnv) <- c("samples", "genes")
            
            print("Initializing df.cnv.long")
            colnames(focal.cnv.maf) <- c("samples", "genes", "mutations")
            df.cnv <- merge(df.cnv, focal.cnv.maf, c("samples", "genes"), all.x = TRUE)
            
            # annotes gene coverage
            print("annotating coverage")
            for (i in 1:length(unique(df.cnv$samples))){
                tmp.sample <- unique(df.cnv$samples)[i]
                tmp.sample.ver <- master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]

                ## checks to see how many genes for each sample were not covered in that assay, labels them as such
                ghosts <- df.cnv[df.cnv$samples == tmp.sample & df.cnv$genes %in% not.covered.cnv[[tmp.sample.ver]], ]
                if (nrow(ghosts) > 0){
                    df.cnv[df.cnv$samples == tmp.sample & df.cnv$genes %in% not.covered[[tmp.sample.ver]], ]$mutations <- "nc"
                }else{
                    ## do nothing
                }

            }
            
            df.cnv$genes <- sapply(df.cnv$genes, paste, "-cnv", sep = "", USE.NAMES = F)
            
            ## generate subset with mutations for accurate counting
            df_sub <- df.cnv[!is.na(df.cnv$mutations), ]
            ord.cnv <- names(sort(table(df_sub$genes), decreasing = T))
            
            if(exists("df.combined")){
                df.combined <- rbind(df.combined, df.cnv)
                ord <- c(ord, ord.cnv)
            }else{
                df.combined <- df.cnv
                ord <- ord.cnv
            }
            
        
        }
    }
    
    ## checks to see if broad.cnv information supplied to comut
    if (!missing(broad.cnv.maf)){
        if (!is.na(broad.cnv.maf)){
            ## renames df to enable merging with mutations df, adds to previous dataframe, 
            ## then inserts names at end of dataframe for factor ordering
            
            ## first transpose so samples are columns
            broad.cnv.maf <- t(broad.cnv.maf)
            broad.cnv.maf <- as.data.frame(broad.cnv.maf)
            broad.cnv.maf$samples <- rownames(broad.cnv.maf)
            broad.cnv.maf <- melt(broad.cnv.maf, id = "samples")
            
            colnames(broad.cnv.maf) <- c("samples", "genes", "mutations")

            ## convert 1 and -1 to gain and loss for clarity
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == -1] <- "arm-level loss"
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == 1] <- "arm-level gain"
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == 0] <- NA
            ord.broad.cnv <- names(sort(table(broad.cnv.maf$genes), decreasing = T))
            
            ## add to df
            if(exists("df.combined")){
                df.combined <- rbind(df.combined, broad.cnv.maf)
                ord <- c(ord, ord.broad.cnv)
            }else{
                df.combined <- broad.cnv.maf
                ord <- ord.broad.cnv
            }
        }
    }
    
    ## checks to see if additional phenotype information supplied for comut
    if (!missing(phenotypes)){
        ## transforms to long format so it can be added to df
        df.pheno <- melt(phenotypes, id = colnames(phenotypes)[1])
        
        ## renames df to enable merging with mutations df, adds to previous dataframe, 
        ## then inserts names at end of dataframe for factor ordering
        colnames(df.pheno) <- c("samples", "genes", "mutations")
        df.combined <- rbind(df.combined, df.pheno)
        ord <- c(ord, colnames(phenotypes)[-1])
    }
    
    ## if manual ordering supplied, resorts for factor ordering
    if (!missing(manual.order)){
        ## make sure manual names are actually in df
        if (sum(!(ord %in% df.combined$genes)) > 0){
            stop(paste("Manual order term contains genes not found in data frame: ", ord[!(ord %in% df.combined$genes)]))
        }else {
            ord <- c(manual.order, ord)
            ord <- ord[!duplicated(ord)]
        }   
    }
    
    ## creats factors based on levels, removed genes not present in ord conditions
    df.combined$genes <- factor(df.combined$genes, levels = ord)
    df.combined <- df.combined[order(df.combined$genes), ]
    df.combined <- df.combined[!is.na(df.combined$genes), ]
    
    ## reshape long df to wide format to determine ordering of samples (columns)
    df.wide <- reshape(df.combined, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")
    for (i in 2:ncol(df.wide)){
        ## change NA's and uncovered to filler which will always order last
        df.wide[, i][is.na(df.wide[, i])] <- "z"
        df.wide[, i][df.wide[, i] == "nc"] <- "z"
        ## change mutations to "mutation" so type doesn't factor into cascading order
        df.wide[, i][df.wide[, i] %in% c("frameshift_indel", "missense", "nonsense", "splice_site", "stop_codon", "in_frame_indel", "other", "TSS", 
                                         "damaging mutation", "focal gain", "focal loss", "rearrangemen")] <- "mutation"
    }
    
    ## sorts by single gene if only one significantly mutated, otherwise all columns
    if (ncol(df.wide) == 2){
        df.wide <- df.wide[order(df.wide[, 2]), ]
        
    }else{
        df.wide <- df.wide[do.call("order", df.wide[, -1]), ]
    }
    
    ## takes given order, and sorts samples based on ordering
    ## add on samples with no mutations to factor creation so they don't get NA'd
    missing <- df.combined$sample[!(df.combined$samples %in% df.wide$samples)]
    df.combined$samples <- factor(df.combined$samples, levels = c(df.wide$sample, missing))
    
    ## overwrites ordering if fixed order supplied
    if(!missing(fixed.order)){
        df.combined$samples <- factor(df.combined$samples, levels = fixed.order)
    }
    
    ## final cleanup
    ## switches factor order back for plotting so most frequent is highest on y axis
    df.combined$genes <- factor(df.combined$genes, rev(levels(df.combined$genes)))
    df.combined$mutations[is.na(df.combined$mutations)] <- "wt"
    df.combined <- df.combined[df.combined$samples %in% samples,]
    print("plotting")
    
    ## generate color scheme for plotting
    
    if (!missing(col.vector)){
        print("using supplied color scheme")
        
    }else{
        print("using random color scheme")
        
        ## first set defaults
        default.names <- c("arm-level gain", "arm-level loss", "HA", "2DEL")
        default.colors <- c("red", "blue", "red", "blue")
        names(default.colors) <- default.names
        col.vector <- default.colors
    }

    ## gets arguments that don't have assigned color
    missing.names <- unique(df.combined$mutations)[!(unique(df.combined$mutations) %in% names(col.vector))]
    missing.colors <- distinctColorPalette(length(missing.names))
    names(missing.colors) <- missing.names
    col.vector <- c(col.vector, missing.colors)
    
    ## sets wt to grey, then removes those factors that aren't present
    col.vector[which(names(col.vector) == "wt")] <- "beige"
    col.vector[which(names(col.vector) == "nc")] <- "white"
    col.vector <- col.vector[names(col.vector) %in% unique(df.combined$mutations)]
    colScale <- scale_fill_manual(values = col.vector)
    
    ## generate plot
    mut <- ggplot(df.combined, aes(x=samples, y=genes, height=0.8, width=0.8)) + 
        geom_tile(aes(fill=mutations)) +
        colScale +
        ggtitle(title) +
        theme(
            legend.key = element_rect(fill='NA'),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_blank(),
            legend.position="bottom",
            legend.text = element_text(size=legend.font.size, face="bold"),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            #axis.text.x=element_text(angle = 90, hjust = 1),
            axis.text.x=element_blank(),
            axis.text.y=element_text(colour="Black", size = y.axis.font.size),
            axis.title.x=element_text(face="bold"),
            axis.title.y=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.background=element_blank()
        )
    print(mut) 
    
    ## prints to screen or to filepath
    if(missing(file.path)){
       print(mut) 
    }else{
        print(paste("printing Comut to ", file.path))
        pdf(file.path, width = dimensions[1], height = dimensions[2])
        print(mut)
        dev.off()
    }
    
    if(return.matrix){
        return(df.combined)
    }
    
    if(return.plot){
        return(mut)
    }
}


PlotHistogramCoverage <- function(maf, gene.column, samples, title){
## Takes in a maf file, and generates a histogram plot by percent samples affected, controling for uncovered genes

## Args:
    ## maf: a maf file
    ## gene.column: the column in maf with gene mutation information
    ## samples: a vector with all samples
    
    ## error checking

    if (is.null(maf[[gene.column]])){
        stop("Invalid column name")
    }
    
    if(missing(samples)){
        stop("Missing sample vector")
    }
    
    ## iniatlize dataframe with naive incidence counts, correcting for multiple samples
    genes <- maf[, gene.column]
    tbl <- table(genes)
    df <- data.frame(tbl)
    
    ## Returns the total number of samples that the gene was covered in
    TrueDenominator <- function(gene){
        ## figure out which versions it was in
        idx <- c(gene %in% not.covered[[1]], gene %in% not.covered[[2]], gene %in% not.covered[[3]])
        rev.idx <- !idx
        versions <- 1:3
        
        ## gets the number of samples that were sequenced with the versions containing that gene
        filtered <- subset(master.sheet, SAMPLE_ACCESSION_NBR %in% samples & PANEL_VERSION %in% versions[rev.idx])
        nrow(filtered)
    }
    
    ## creates dataframe with denominator counts for each gene
    denominators <- sapply(unique(genes), TrueDenominator)
    denominators.df <- data.frame(denominators, names(denominators))
    colnames(denominators.df)[2] <- "genes"
    
    ## merges two dataframes, then create column for plotting
    df.merged <- merge(df, denominators.df, by = "genes")
    df.merged$percent <- df.merged$Freq / df.merged$denominators
    df.merged <- df.merged[order(df.merged$percent, decreasing = T), ]
    df.merged$genes <- factor(df.merged$genes, levels = df.merged$genes)
    plot <- ggplot(data = df.merged, aes(x=genes, y = percent)) + geom_bar(stat = "identity") + ggtitle(label =  title) + rameen_theme
    print(plot)
}

CalculateEntropy <- function(tree, correction.method = "none", detailed = NA){
    ## takes in a matrix of samples and features, then determines which column segregregates the data to minimize the entropy
    
    ## Args:
        ## tree: a matrix or data frame. Columns are features of interest, rows are samples
        ## correction.method: one of none, mult, info, or recip. If none, no attempt is made to correct for
            ## for differences in entropy resulting from unbalanced propotion of 1/0. Multiplicative attempts to correct for entropy 
            ## bias by multiplying by (1 / entropy) of the feature of interest. Information corrects by calculating the entropy of the
            ## current row, and subtracting the entropy of each row from that, computing information gain. Reciprocal multiplies the 1/0
        ## count for each column by the reciprocal of the count for the column of interest to balance out effects due purely to count
        ## detailed: optional argument that will produce detailed breakdown of the entropy contribution of each gene to a 
        ## column of interest
    
    ## error checking
    if (ncol(tree) < 3){
        stop("Insufficient number of columns")
    }
    
    if (nrow(tree) < 3){
        stop("Insufficient number of rows")
    }
    
    if (!missing(correction.method)){
        if (!(correction.method %in% c("mult", "info", "recip"))){
            stop("Unknown correction method: must supply either mult, info, or recip")
        }
    }
    
    tree <- apply(tree, 2, as.numeric)
    ## determines if any columns need to be skipped
    sums <- colSums(tree)
    skip.idx <- sums == 0| sums == nrow(tree)
    tree <- tree[, !skip.idx]
    
    ## determines which counts of "i" will be marked
    if (!missing(detailed)){
        detailed.idx <- which(colnames(tree) %in% detailed)
        detailed.counts <- list()
    }else{
        detailed.idx <- c()
    }
    
    ## calculate entropy
    cumulative.entropy <- c()
    for (i in 1:ncol(tree)){
    
        total.entropy <- 0
        detailed.entropy <- c()
        
        ## calculates reciprocal for reciprocal method
        n2.i <- sum(as.numeric(tree[, i]))
        n1.i <- nrow(tree) - n2.i
        s2 <- (n1.i + n2.i) / n2.i
        s1 <- (n1.i + n2.i) / n1.i
        
        for(j in 1:ncol(tree)){
            if (i == j){
                ## skip
            }else{
                ## split data based on identified feature, figure out total on each side
                x <- table(as.numeric(tree[, i]), as.numeric(tree[, j]))
                n1 <- x[1,2]
                n2 <- x[2,2]
                
                if(correction.method == "recip"){
                    n1 <- n1 * s1
                    n2 <- n2 * s2    
                }
                
                p1 <- n1 / (n1 + n2)
                p2 <- n2 / (n1 + n2)
                
                ## sums entropies
                if (!p1 | !p2){
                    entropy <- 0
                }else{
                    entropy <- -(p1*log2(p1) + p2*log2(p2))
                }
                total.entropy <- total.entropy + entropy
                
                if (i %in% detailed.idx){
                    detailed.entropy <- c(detailed.entropy, entropy)
                    names(detailed.entropy)[length(detailed.entropy)] <- colnames(tree)[j]
                }
            }
        }
        ## penalize total based on skewed entropy of factor
        p1 <- n1.i / (n1.i + n2.i)
        p2 <- n2.i / (n1.i + n2.i)
        
        if (correction.method == "mult"){
            total.entropy <- -total.entropy / (p1*log2(p1) + p2*log2(p2))
        }
        
        if (correction.method == "info"){
            total.entropy <- -((ncol(tree) -1) * (p1*log2(p1) + p2*log2(p2))) - total.entropy 
            
        }

        ## additive correction: equivalent to finding information gain in each case
        cumulative.entropy[i] <- total.entropy
        names(cumulative.entropy)[i] <- colnames(tree)[i]
        
        if(i %in% detailed.idx){
            detailed.counts[[length(detailed.counts) + 1]] <- detailed.entropy
            names(detailed.counts)[length(detailed.counts)] <- colnames(tree)[i]
        }
    
    }
    if (!missing(detailed)){
        return(list(cumulative.entropy, detailed.counts, tree))
    }else{
        return(cumulative.entropy)    
    }
    
}

OptimizedHelper <- function(tree, idx.1, idx.2){
    ## functionalizes inner for loop for feed in to apply
    x <- table(as.numeric(tree[, idx.1]), as.numeric(tree[, idx.2]))
    n1 <- x[1,2]
    n2 <- x[2,2]
    
    p1 <- n1 / (n1 + n2)
    p2 <- n2 / (n1 + n2)
    
    ## sums entropies
    if (!p1 | !p2){
        0
    }else{
        -(p1*log2(p1) + p2*log2(p2))
    }
}

OptimizedCalculateEntropy <- function(tree, correction.method = "none", detailed = NA){
    ## trimmed down version for permutations testing
    ## requires tree to be input as numeric object, with full/empty columns removed
   num.cols <- ncol(tree)
 
    ## calculate entropy
    cumulative.entropy <- rep(0, num.cols)
    for (i in 1:num.cols){
        
        total.entropy <- 0
        inner.idx <- 1:num.cols
        inner.idx <- inner.idx[inner.idx != i]
        total.entropy <- sum(sapply(inner.idx, function(x){OptimizedHelper(tree, i, x)}))
        ## adjust based on information gain
        n2.i <- sum(as.numeric(tree[, i]))
        n1.i <- nrow(tree) - n2.i
        p1 <- n1.i / (n1.i + n2.i)
        p2 <- n2.i / (n1.i + n2.i)
        total.entropy <- -((ncol(tree) -1) * (p1*log2(p1) + p2*log2(p2))) - total.entropy 
        cumulative.entropy[i] <- total.entropy
    }
    return(cumulative.entropy)    
}


BranchPlotter <- function(entropy.data, hit.number = 4, incidence = TRUE){
    ## Takes an entropy object and plots the data at each branch point
    
    ## Args:
        ## entropy data: object created by CalculateEntropy
        ## hit.number: the number of features to show detailed information on
        ## incidence: boolean that determines whether incidence plot is generated
    if (!is.list(entropy.data)){
        summary.plot.data <- data.frame(entropy.data, names(entropy.data), stringsAsFactors = FALSE)
    }else{
        summary.plot.data <- data.frame(entropy.data[[1]], names(entropy.data[[1]]), stringsAsFactors = FALSE)
    }
    colnames(summary.plot.data) <- c("information_gain", "feature")
    summary.plot.data$feature <- factor(summary.plot.data$feature, levels = summary.plot.data$feature[order(summary.plot.data$information_gain, decreasing = TRUE)])
    summary.plot.data <- summary.plot.data[summary.plot.data$information_gain > 0, ]
    summary.plot <- ggplot(data = summary.plot.data, aes(x = feature, y = information_gain)) + geom_bar(stat = "identity") + 
        labs(title = "Features ranked by net information gain") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    detailed.plots <- names(entropy.data[[2]])
    return.plots <- list(summary.plot)
    
    
    if (hit.number != 0){
        for (i in 1:length(detailed.plots)){
            subplot.idx <- which(names(entropy.data[[2]]) == detailed.plots[i])
            subplot.data <- data.frame(entropy.data[[2]][[subplot.idx]], names(entropy.data[[2]][[subplot.idx]]), stringsAsFactors = FALSE)
            colnames(subplot.data) <- c("entropy_contribution", "feature")
            subplot.data$feature <- factor(subplot.data$feature, levels = subplot.data$feature[order(subplot.data$entropy_contribution)])
            subplot.data <- subplot.data[subplot.data$entropy_contribution < 0.8, ]
            plot.title <- paste("Entropy by feature for", detailed.plots[i])
            return.plots[[i + 1]] <- ggplot(data = subplot.data, aes(x = feature, y = entropy_contribution)) + geom_bar(stat = "identity") + 
                labs(title = plot.title) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(limits = c(0, 0.8))
        }
    }
    
    if (incidence){
        
        incidence.plot.data <- colSums(entropy.data[[3]])
        incidence.plot.data <- data.frame(names(incidence.plot.data), as.numeric(incidence.plot.data), stringsAsFactors = FALSE)
        colnames(incidence.plot.data) <- c("feature", "count")
        incidence.plot.data <- incidence.plot.data[incidence.plot.data$feature %in% detailed.plots, ]
        incidence.plot.data$feature <- factor(incidence.plot.data$feature, levels = levels(summary.plot.data$feature))
        incidence.plot <- ggplot(data = incidence.plot.data, aes (x = feature, y = count)) + geom_bar(stat = "identity") +
            labs(title = "Number of samples with alteration") + scale_y_continuous(limits = c(0, nrow(entropy.data[[3]])))
        return.plots[[length(return.plots) + 1]] <- incidence.plot
    }
    return(return.plots)
}

PermuteDecisionTree <- function(tree, permutation.number){
    ## Takes a decision tree, permutes the data, then recalculates entropy to determine whether cutoffs are significant
    
    ## Args:
        ## tree: a decision tree
        ## permutation number: number of iterations to permute through
    reps <- sum(colSums(tree) != 0 & colSums(tree) != nrow(tree))
    sample.number <- nrow(tree)
    permuted.values <- rep(0, permutation.number)
    #permuted.values <- rep(0, permutation.number * reps)
    for (i in 1:permutation.number){
        permuted.tree <- apply(tree, 2, function(x){(sample(x, size = sample.number, FALSE))})
        #permuted.values[(reps * (i - 1) + 1):(reps * i)] <- as.numeric(OptimizedCalculateEntropy(permuted.tree))
        permuted.values[i] <- max(as.numeric(OptimizedCalculateEntropy(permuted.tree)))
    }
    
    #sapply(values, function(x){sum(permuted.values >= x)}) / length(permuted.values)
    permuted.values
}





