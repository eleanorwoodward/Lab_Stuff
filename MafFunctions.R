## Noah Greenwald
## 9/17/15

## Helpful R functions for manipulating MAF output files for data analysis. 
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/ExacFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/ESPFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/PoNFilter.R")
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/CallCompare.R")

## Nice packages
library(grDevices)
library(RColorBrewer)

FilterMaf <- function(maf, list, column.name, keep = TRUE){
  # Takes a maf, and returns a maf with only those rows who have a value in given column that matches an entry in the list. 
  
  # Args:
  #     maf: a maf
  #     list: list of column values to keep
  #     column.name: column in the maf which will be checked
  #     keep: logical value. When false, discards selected rows
  
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

PlotMaf <- function(maf, column.name, lower.margin = 1, title = ""){
  ## Takes a Maf file and produces a bar plot 
  
  # Args: 
  #     Maf: a maf file in array format. 
  #     column.name: a string with the name of the column that the barplot will be produced from. 
  #     lower.margin: a multiplier that adjusts the marign at the bottom of the graph. 
  #     title: optional, gives the graph a title
  
  ## Error checking
  if (is.null(maf[[column.name]])){
    stop("Invalid column name")
  }
  
  data <- sort(maf[, column.name])
  tbl <- table(data)
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

CountMaf <- function(maf, column, var.list){
  # Counts the number of entries in a given column of the maf that match values in the list.'
  
  # Args:
  #     maf: a maf
  #     column: the column whose entries will be counted
  #     var.list: list of column values that will be counted
  
  # Returns: int with number of occurences
  
  ## Error Checking
  if (is.null(maf[[column]])){
    stop("Invalid column name")
  }
  
  attendance <- var.list %in% maf[[column]]
  if (sum(attendance) < length(attendance)){
    warning("List element ", list[which(!attendance)], " not found")
  }
  
  len <- 0
  for (i in 1:length(var.list)){
    len <- len + length(maf[which(maf[[column]] == var.list[i])])
  }
  return (len)
}

MiniMaf <- function(maf, column.list, keep = TRUE){
  # Generates a paired down version of a maf given a list of columns
  
  # Args:
  #     maf: a maf
  #     column.list: a list of columns in original maf to filter based on
  #     keep: when true, returns maf with only given columns. When false, returns maf without given columns
  
  # Returns: a reduced size maf
  
  ## Error Checking
  attendance <- column.list %in% colnames(maf)
  if (sum(attendance) < length(attendance)){
    stop("List element ", column.list[which(!attendance)], " not found")
  }
  
  
  idx <- names(maf) %in% column.list
  if (keep == TRUE){
    mini <- maf[, idx]
    return (mini)
  }else{
    mini <- maf[, !idx]
    return (mini)
  }
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

PerSampleMaf <- function(maf, column.name, list = maf[[column.name]]){
  ## Takes a maf, and returns a filtered maf that limits occurence to one per sample
  
  ## Args:
  #       maf: a maf file
  #       column.name: name of the column for filtering
  #       list: values that will be kept in filtered maf. Default is all values in column
  
  # Returns: a maf with only those rows which contain an element from list in selected column, 
  # limited to once per sample. 
  filtered <- FilterMaf(maf, list, column.name, TRUE)
  cols <- c("Tumor_Sample_Barcode", column.name)
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

## Takes a folder full of Mafs, and returns one concatenated maf

CombineMaf <- function(directory, file.list = list.files(directory)){
  x <- read.delim(paste(directory, file.list[1], sep = "/"), stringsAsFactors = FALSE, comment.char = "#")
  for (i in 2:length(file.list)){
    temp <- read.delim(paste(directory, file.list[i], sep = "/"), stringsAsFactors = FALSE, comment.char = "#")
    x <- rbind(x, temp)
    }
    x
  }

PairSetFormat <- function(sample.name, cutoff = 3, replace = "M"){
  temp <- substring(sample.name, cutoff)
  temp <- paste(replace, temp, sep = "")
  temp
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

