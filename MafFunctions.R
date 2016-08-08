## Noah Greenwald
## 9/17/15

## Helpful R functions for manipulating MAF output files for data analysis. 
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering/ExacFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//ESPFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//PoNFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//CallCompare.R")

# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("devtools")
# install.packages("R.utils")
# 
## Nice packages
library(grDevices)
library(RColorBrewer)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(curl)
library(R.utils)
library("plyr")

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
