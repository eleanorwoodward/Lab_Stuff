## Merge PDL-1 data 

Master <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/Master_Sheet_For_R.txt", stringsAsFactors = F)

PDL <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/PDL_for_R.txt", stringsAsFactors = F)

Master[, 6:13] <- "no data"
for (i in 1:nrow(Master)){
    idx <- PDL$mrn %in% Master$mrn[i]
    if (sum(idx) == 1){
        row <- PDL[idx,4:11]
        Master[i, 6:13] <- row
        
    }else if (sum(idx) > 1){
        Master[i, ] <- "Duplicates"
    }
}

write.csv(Master, "C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/R_ouput.csv", row.names = F)

## Correlations analysis
correlations.data <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/Previous versions/Master_Sheet_For_R.txt",
                                stringsAsFactors = F)
correlations.data <- correlations.data[1:222, ]

columns <- colnames(correlations.data)

## Set up variables for use in columns
imaging.features <- c(72:83)
names(imaging.features) <- c("Edema", "Shape", "Intratumoral Heterogeneity", "Multifocality", "Midline shift", "Sinus invasion", 
                             "Bone Invasion", "Hyperostosis", "Cystic", "Necrosis Hemorrhage", "Overt hemorrhage", "Mass Effect")
grades <- c(7, 9)
names(grades) <- c("Molecular Grade", "Grade")

mutation.status <- c(171:179)
names(mutation.status) <- columns[mutation.status]

location <- c(39, 40)
names(location) <- c("Location Class 1", "Location Class 2")

cas <- 180
names(cas) <- "CAS"

mus <- 182
names(mus) <- "MUS"

# 1.	semantic imaging features (columns BT-CE) vs. grade (col G vs. col I) 
#       fisher's exact for individual, log-linear for multiple 

mtrx1 <- matrix(NA, length(imaging.features), 2)
rownames(mtrx1) <- names(imaging.features)
colnames(mtrx1) <- names(grades)
for (i in 1:length(imaging.features)){
    mini <- CleanYourRoom(correlations.data[, c(grades[2], imaging.features[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if (length(dim(tbl)) == 2){
        math.party <- chisq.test(tbl)
        mtrx1[i, 2] <- signif(math.party$p.value, 2)
    }    
}
mtrx1

## multivariate analysis
small <- CleanYourRoom(correlations.data[, c(grades, imaging.features)])
summary(lm(small[, 1] ~ as.character(small[, 3]) + as.character(small[, 4]) + as.character(small[, 5])
           
colnames(small) <- c(names(grades), names(imaging.features))
summary(lm(Grade ~ Edema + Shape, small))
summary(glm(Grade ~ Edema + Shape, data = small, family = binomial(link = "logit")))
                      

> row1 <- c(rep(0, 50), rep(1, 50))
> row2 <- c(rnorm(50, 10, 3), rnorm(50, 14, 3))
> row3 <- c(rep("A", 30), rep("B", 10), rep("C", 10), rep("A", 10), rep("C", 30), rep("B", 10))


# 2.	Mutation status (column Q-Y in Parker sheet individually assessed) vs. grade (col G vs I) - 
#       fisher's exact for individual, log-linear for multiple?

mtrx2 <- matrix(NA, length(mutation.status), 2)
rownames(mtrx2) <- names(mutation.status)
colnames(mtrx2) <- names(grades)
for (i in 1:length(mutation.status)){
    mini <- CleanYourRoom(correlations.data[, c(grades[1], mutation.status[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if (length(dim(tbl)) == 2){
        
        math.party <- chisq.test(tbl)
        mtrx2[i, 1] <- signif(math.party$p.value, 2)
    }
}
mtrx2

# 3.	Mutation status (column Q-Y in Parker sheet) vs. Location (col AM & AN; for classification2 in AN, exclude class 5)
mtrx3 <- matrix(NA, length(mutation.status), 2)
rownames(mtrx3) <- names(mutation.status)
colnames(mtrx3) <- names(location)
for (i in 1:length(mutation.status)){
    ## Filter if location2 used
    #mini <- CleanYourRoom(correlations.data[, c(location[2], mutation.status[i])], 5)
    mini <- CleanYourRoom(correlations.data[, c(location[1], mutation.status[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if (length(dim(tbl)) == 2){
        math.party <- chisq.test(tbl)
        mtrx3[i, 1] <- signif(math.party$p.value, 2)
    }
}
mtrx3


# 4.	Mutation status (individual column Q-Y in Parker sheet) vs. semantic imaging features (columns BT-CE) 
#       Fisher's/Chi for individual, log-linear for multiple

mtrx4 <- matrix(NA, length(mutation.status), length(imaging.features))
rownames(mtrx4) <- names(mutation.status)
colnames(mtrx4) <- names(imaging.features)
for (i in 1:length(imaging.features)){
    
    for (j in 1:length(mutation.status)){
        mini <- CleanYourRoom(correlations.data[, c(imaging.features[i], mutation.status[j])])
        tbl <- table(mini[, 1], mini[, 2])
        if (length(dim(tbl)) == 2){
            math.party <- chisq.test(tbl)
            mtrx4[j, i] <- signif(math.party$p.value, 2)
        }
    }
}
mtrx4

# 5.	CAS score (col CB in Parker sheet) vs. grade (col G vs. I) 
# ANOVA for difference in means between categorical variable
mtrx5 <- matrix(NA, length(grades), 1)
rownames(mtrx5) <- names(grades)
colnames(mtrx5) <- c("CAS")
for (i in 1:length(grades)){
    mini <- CleanYourRoom(correlations.data[, c(cas, grades[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx5[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx5
  
# 6.	CAS score (col CB in Parker sheet) vs. location (col AM & AN, ignore class 5 for col AN) 
#       ANOVA for difference in means between locations
mtrx6 <- matrix(NA, length(location), 1)
rownames(mtrx6) <- names(location)
colnames(mtrx6) <- c("CAS")
for (i in 1:length(location)){
    mini <- CleanYourRoom(correlations.data[, c(cas, location[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx6[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx6
# 7.	CAS score (col CB in Parker sheet) vs. semantic imaging features (columns BT-CE in master sheet) 
mtrx7 <- matrix(NA, length(imaging.features), 1)
rownames(mtrx7) <- names(imaging.features)
colnames(mtrx7) <- c("CAS")
for (i in 1:length(imaging.features)){
    mini <- CleanYourRoom(correlations.data[, c(cas, imaging.features[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx7[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx7

# 	Copy number variation vs. grade / location / semantic imaging features ANOVA for all
##Grade
mtrx8 <- matrix(NA, length(grades), 1)
rownames(mtrx8) <- names(grades)
colnames(mtrx8) <- c("MUS")
for (i in 1:length(grades)){
    mini <- CleanYourRoom(correlations.data[, c(mus, grades[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx8[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx8

## location
mtrx9 <- matrix(NA, length(location), 1)
rownames(mtrx9) <- names(location)
colnames(mtrx9) <- c("MUS")
for (i in 1:length(location)){
    mini <- CleanYourRoom(correlations.data[, c(mus, location[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx9[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx9

## imaging features
mtrx10 <- matrix(NA, length(imaging.features), 1)
rownames(mtrx10) <- names(imaging.features)
colnames(mtrx10) <- c("MUS")
for (i in 1:length(imaging.features)){
    mini <- CleanYourRoom(correlations.data[, c(mus, imaging.features[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx10[i, 1] <- summary(math.party)[[1]][1,5]
    }
}
mtrx10


## Multivariate regressions
summary(lm(correlations.data$MUS ~ correlations.data$KLF4.MUTATION + as.character(correlations.data$grade))
?lm
