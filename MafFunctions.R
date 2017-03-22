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
# install.packages("kohonen")
# install.packages("nycflights13")
## Nice packages
library(grDevices)
library(RColorBrewer)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(curl)
library(R.utils)
library(Biobase)
#library(kohonen)
library(nycflights13)

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

PlotComut <- function(maf, samples, input.samples = "Tumor_Sample_Barcode", input.genes = "Hugo_Symbol", input.mutations = "Variant_classification", 
                      gene.cutoff = 1, file.path = NA, unique.muts = NA, phenotypes = NA, title = "Comut Plot"){
## Takes in a maf file, and generates a comut plot of common mutations in each of the samples
## adapted from code at https://benchtobioinformatics.wordpress.com/2015/05/25/how-to-make-a-co-mutation-plot/
    

## Args:C
    ## maf: a maf file
    ## samples: a vector with a list of all samples
    ## input.samples: a vector with all individuals to be plotted (don't take from MAF, could be absent due to no mutations)
    ## input.genes: a column within the maf which specifies genes
    ## input.mutations: a column within the maf describing the type of mutation
    ## gene.cutoff: a threshold for deciding how many genes to include in comut. if beteen 0 and 1, treated as percentage of total. otherwise, number of genes
    ## unique.muts: vector with all the unique mutation classes present in dataset. Useful for consistent colors across multiple comuts
    ## phenotypes: a dataframe with samples as the first column, and subsequent columns taken as additional information to plot in given order
    
    ## Error Checking
    if (is.null(maf[[input.samples]]) |is.null(maf[[input.genes]]) | is.null(maf[[input.mutations]])){
        stop("Invalid column name")
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
    
    ## initailize long dataframe with all possible gene-sample combinations
    all.samples <- unique(samples)
    all.genes <- unique(maf[, input.genes])
    df <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
    colnames(df) <- c("samples", "genes")
    df$mutations <- NA
    
    # Loops through each row in df, determining if given gene is mutated in given sample in maf
    # If so, determines what type of mutation, then writes information to df
    print("Initializing maf")
    for (i in 1:nrow(df)){
        tmp.sample <- df$samples[i]
        tmp.gene <- df$genes[i]
        matches <- maf[maf[, input.samples] == tmp.sample & maf[, input.genes] == tmp.gene, ]
        if (nrow(matches) == 0){
            ## do nothing
        }else{
            df$mutations[i] <- matches[, input.mutations][1]
        }
    }
    df.old <- df
    #samp.old <- unique(df$samples)
    ## annotate coverage information
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
    
    ## keep given % of genes based on cutoff
    if (gene.cutoff == 1){
        ## keep all genes
    }else if (gene.cutoff < 1){
        ## keep given percentage of top hits
        idx <- ceiling(length(ord) * gene.cutoff)
        ord <- ord[1:idx]
    }else{
        ## keep given number of top hits
        ord <- ord[1:idx]
    }
    
    
    ## checks to see if additional phenotype information supplied for comut
    if (!missing(phenotypes)){
        ## transforms to long format so it can be added to df
        df.pheno <- melt(phenotypes, id = colnames(phenotypes)[1])
    }
    
    ##renames df, adds to previous dataframe, then inserts names at end of dataframe for factor ordering
    colnames(df.pheno) <- c("samples", "genes", "mutations")
    df <- rbind(df, df.pheno)
    ord <- c(ord, colnames(phenotypes)[-1])
    
    ## creats factors based on levels, removed genes not present in ord conditions
    df$genes <- factor(df$genes, levels = ord)
    df <- df[order(df$genes), ]
    df <- df[!is.na(df$genes), ]
    
    ## reshape long df to wide format to determine ordering of samples (columns)
    ## first remove non-mutation columns to not intefere with sorting
    df.cleaned <- df[!df$genes %in% colnames(phenotypes)[-1], ]
    df.wide <- reshape(df.cleaned, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")
    for (i in 2:ncol(df.wide)){
        ## change NA's and uncovered to filler which will always order last
        df.wide[, i][is.na(df.wide[, i])] <- "z"
        df.wide[, i][df.wide[, i] == "nc"] <- "z"
    }
    
    ## sorts by single gene if only one significantly mutated, otherwise all columns
    if (ncol(df.wide) == 2){
        df.wide <- df.wide[order(df.wide[, 2]), ]
        
    }else{
        df.wide <- df.wide[do.call(order, df.wide[, -1]), ]
    }
    
    ## takes given order, and sorts samples based on ordering
    ## add on samples with no mutations to factor creation so they don't get NA'd
    missing <- df$sample[!(df$samples %in% df.wide$samples)]
    df$samples <- factor(df$samples, levels = c(df.wide$sample, missing))

    ## switches factor order back for plotting so most frequent is highest on y axis
    df$genes <- factor(df$genes, rev(levels(df$genes)))
    print("plotting")
    
    df$mutations[is.na(df$mutations)] <- "wt"
    
    ## generate color scheme for plotting
    
    if (!missing(unique.muts)){
        ## creates levels based on supplied consistent levels, with wt and nc added on at end, with filler based on number of phenotypes 
        df$mutations <- factor(df$mutations, levels = c(unique.muts, unique(df.pheno$mutations), "nc", "wt"))
        myColors <- c("red", "skyblue", "orange", "yellow", "purple", "gold", "gold", "cyan", brewer.pal(length(levels(df$mutations)) - 8, "Set3"))
        names(myColors) <- levels(df$mutations)
        
        ## sets wt to grey, then removes those factors that aren't present
        myColors[which(names(myColors) == "wt")] <- "beige"
        myColors[which(names(myColors) == "nc")] <- "white"
        myColors <- myColors[levels(df$mutations) %in% unique(df$mutations)]
        colScale <- scale_fill_manual(values = myColors)
    }else{
        df$mutations <- factor(df$mutations, levels = unique(df$mutations))
        myColors <- brewer.pal(length(levels(df$mutations)),"Set3")
        names(myColors) <- levels(df$mutations)
        myColors[which(names(myColors) == "wt")] <- "beige"
        myColors[which(names(myColors) == "nc")] <- "white"
        colScale <- scale_fill_manual(values = myColors)
    }
    
    ## generate plot
    mut <- ggplot(df, aes(x=samples, y=genes, height=0.8, width=0.8)) + 
        geom_tile(aes(fill=mutations)) +
        colScale +
        #scale_colour_discrete(drop=TRUE, values = c("green", "blue", "red", "yellow", "orange", "purple", "gold", "cyan", "grey"), limits = levels(df$mutations)) +
        #scale_fill_manual(name = "mutations", values = c("green", "blue", "red", "yellow", "orange", "purple", "gold", "grey")) +
        xlab("Subject") +
        ggtitle(title) +
        theme(
            legend.key = element_rect(fill='NA'),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_blank(),
            legend.position="bottom",
            legend.text = element_text(size=8, face="bold"),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(angle = 90, hjust = 1),
            axis.text.y=element_text(colour="Black"),
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
    if(is.na(file.path)){
       print(mut) 
    }else{
        print(paste("printing Comut to ", file.path))
        pdf(file.path, width = 14)
        print(mut)
        dev.off()
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


