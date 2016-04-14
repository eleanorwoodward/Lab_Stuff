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


## format data for analysis
correlations.data <- read.csv("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/Master_Sheet_CSV.csv",
                                stringsAsFactors = F)
correlations.data <- correlations.data[1:222, ]
correlations.data[, 183] <- as.numeric(correlations.data[, 183])
correlations.data[, 184] <- as.numeric(correlations.data[, 184])
correlations.data[, 185] <- as.numeric(correlations.data[, 185])
correlations.data[, 186] <- as.numeric(correlations.data[, 186])
correlations.data[, 187] <- as.numeric(correlations.data[, 187])
correlations.data[, 188] <- as.numeric(correlations.data[, 188])
correlations.data[, c(190:192)] <- NA
colnames(correlations.data)[190:192] <- c("RNA.Score", "IHC.Score", "Combined.Score")


##PD-L1 QC

# Sets median value for columns for 2 TMAs based on different vals
medians.283 <- c(179, .82, NA, 63, 6, NA)
medians.310 <- c(174, .84, 1.6, 162, NA, 85)
for (i in 1:nrow(correlations.data)){
    if (!is.na(correlations.data$TMA[i])){
    ## computes number of metrics pointing towards above or below median results
        if (correlations.data$TMA[i] == 283){
            vals <- correlations.data[i, 183:188]
            list <- vals[1:3] > medians.283[1:3]
            pluses <- sum(list, na.rm = T)
            total <- 3 - sum(is.na(list)) 
            minuses <- total -  pluses
            correlations.data[i, 190] <- round((pluses - minuses) / total, 1)
            list <- vals[4:6] > medians.283[4:6]
            pluses <- sum(list, na.rm = T)
            total <- 3 - sum(is.na(list)) 
            minuses <- total -  pluses
            
            correlations.data[i, 191] <- round((pluses - minuses) / total, 1)
        }else{
            vals <- correlations.data[i, 183:188]
            list <- vals[1:3] > medians.310[1:3]
            pluses <- sum(list, na.rm = T)
            total <- 3 - sum(is.na(list)) 
            minuses <- total -  pluses
            correlations.data[i, 190] <- round((pluses - minuses) / total, 1)
            
            list <- vals[4:6] > medians.310[4:6]
            pluses <- sum(list, na.rm = T)
            total <- 3 - sum(is.na(list)) 
            minuses <- total -  pluses
            correlations.data[i, 191] <- round((pluses - minuses) / total, 1)
        }
    }
}

## Takes RNA and IHC data and creates a combined MIB-1 Score. 

# for (i in 1:nrow(correlations.data)){
#     if (!is.na(correlations.data$TMA[i])){
#         rna <- correlations.data[i, 190]
#         ihc <- correlations.data[i, 191]
#         if (rna == ihc){
#             correlations.data[i, 192] <- rna
#         }else if (ihc == 0 & rna ){
#             
#         }
#     }
# }

## Correlations analysis


columns <- colnames(correlations.data)



## Set up variables for use in columns
imaging.features <- c(72:83)
names(imaging.features) <- c("Edema", "Shape", "Intratumoral Heterogeneity", "Multifocality", "Midline shift", "Sinus invasion", 
                             "Bone Invasion", "Hyperostosis", "Cystic", "Necrosis Hemorrhage", "Overt hemorrhage", "Mass Effect")
imaging.features <- imaging.features[-11]


grades <- c(7, 9)
names(grades) <- c("Molecular Grade", "Grade")

mutation.status <- c(171:179)
names(mutation.status) <- columns[mutation.status]
mutation.status <- mutation.status[-2]

location <- c(39, 40)
names(location) <- c("Location Class 1", "Location Class 2")

cas <- 180
names(cas) <- "CAS"

mus <- 182
names(mus) <- "MUS"

atypical <- c(11, 13,14,16)
names(atypical) <- c("Hypercellularity", "Sheeting", "Prominent Nucleoli", "Necrosis")

pdl <- c(190:192)
names(pdl) <- c("RNA", "IHC", "Combined")

# 1.	semantic imaging features (columns BT-CE) vs. grade (col G vs. col I) 
#       fisher's exact for individual, multinomial logistic regression for multivariate

mtrx1 <- matrix(NA, length(imaging.features), 2)
rownames(mtrx1) <- names(imaging.features)
colnames(mtrx1) <- names(grades)
for (i in 1:length(imaging.features)){
    mini <- CleanYourRoom(correlations.data[, c(grades[1], imaging.features[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if (length(dim(tbl)) == 2){
        math.party <- chisq.test(tbl)
        mtrx1[i, 1] <- signif(math.party$p.value, 2)
    }    
}
mtrx1

## multivariate analysis
small <- CleanYourRoom(correlations.data[, c(grades, imaging.features)])
colnames(small) <- c("molgrade", "grade", names(imaging.features))
model <- multinom(molgrade ~ ., small[,-2])
summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx1a <- matrix(NA, 12, 2)
mult.mtrx1a[,1] <- p[1,-1 ]
mult.mtrx1a[,2] <- p[2,-1 ]
rownames(mult.mtrx1a) <- c(dimnames(summary(model)$coefficients)[[2]][-1])
colnames(mult.mtrx1a) <- c("Grade2", "Grade3")
mult.mtrx1a


model <- multinom(grade ~ ., small[,-1])
summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx1b <- matrix(NA, 12, 2)
mult.mtrx1b[, 1] <- p[1, -1]
mult.mtrx1b[, 2] <- p[2, -1]
rownames(mult.mtrx1b) <- dimnames(summary(model)$coefficients)[[2]][-1]
colnames(mult.mtrx1b) <- c("Grade2", "Grade3")
mult.mtrx1b

# 2.	Mutation status (column Q-Y in Parker sheet individually assessed) vs. grade (col G vs I) - 
#       fisher's exact for individual, multinomial logistic regression for multivariate

mtrx2 <- matrix(NA, length(mutation.status), 2)
rownames(mtrx2) <- names(mutation.status)
colnames(mtrx2) <- names(grades)
for (i in 1:length(mutation.status)){
    mini <- CleanYourRoom(correlations.data[, c(grades[2], mutation.status[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if ((dim(tbl))[2] == 2){
        
        math.party <- chisq.test(tbl)
        mtrx2[i, 2] <- signif(math.party$p.value, 2)
    }
}
mtrx2

##Multivariate
small <- CleanYourRoom(correlations.data[, c(grades, mutation.status)], "0?")
colnames(small) <- c("molgrade", "grade", names(mutation.status))
remove <- c()
for (i in 1:ncol(small)){
    
    if (dim(table(small[, i])) < 2){
        remove <- c(remove, i)
    } 
    
}
small <- small[, -remove]
model <- multinom(molgrade ~ ., small[,-2])
summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx2a <- matrix(NA, ncol(p) - 1, 2)
mult.mtrx2a[, 1] <- p[1, -1]
mult.mtrx2a[, 2] <- p[2, -1]
rownames(mult.mtrx2a) <- c(dimnames(summary(model)$coefficients)[[2]][-1])
colnames(mult.mtrx2a) <- c("Grade2", "Grade3")
mult.mtrx2a

model <- multinom(grade ~ ., small[,-1])
summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx2b <- matrix(NA, 2, ncol(p) - 1)
mult.mtrx2b <- p[, -1]
colnames(mult.mtrx2b) <- c(dimnames(summary(model)$coefficients)[[2]][-1])
rownames(mult.mtrx2b) <- c("Grade2", "Grade3")



# 3.	Mutation status (column Q-Y in Parker sheet) vs. Location (col AM & AN; for classification2 in AN, exclude class 5)
mtrx3 <- matrix(NA, length(mutation.status), 2)
rownames(mtrx3) <- names(mutation.status)
colnames(mtrx3) <- names(location)
for (i in 1:length(mutation.status)){
    ## Filter if location2 used
    mini <- CleanYourRoom(correlations.data[, c(location[1], mutation.status[i])], 5)
    #mini <- CleanYourRoom(correlations.data[, c(location[1], mutation.status[i])])
    tbl <- table(mini[, 1], mini[, 2])
    if (dim(tbl)[2] == 2){
        math.party <- chisq.test(tbl)
        mtrx3[i, 1] <- signif(math.party$p.value, 2)
    }
}
mtrx3



##Multivariate
small <- CleanYourRoom(correlations.data[, c(location, mutation.status)], "0?")
colnames(small) <- c("loc.class.1", "loc.class.2", names(mutation.status))
remove <- c()
for (i in 1:ncol(small)){
    
    if (dim(table(small[, i])) < 2){
        remove <- c(remove, i)
    } 
    
}
small <- small[, -remove]
model <- multinom(loc.class.1 ~ ., small[,-2])
summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx3a <- matrix(NA, 2, ncol(p) - 1)
mult.mtrx3a <- p[, -1]
colnames(mult.mtrx3a) <- c(dimnames(summary(model)$coefficients)[[2]][-1])
rownames(mult.mtrx3a) <- c("Location2", "Location3")


small <- CleanYourRoom(small, 5)
model <- multinom(loc.class.2 ~ ., small[,-1])

summary(model)
#compute p values from standard errors
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- (1- pnorm(abs(z), 0, 1))*2
mult.mtrx3b <- matrix(NA, nrow(p), ncol(p) - 1)
mult.mtrx3b <- p[, -1]
colnames(mult.mtrx3b) <- c(dimnames(summary(model)$coefficients)[[2]][-1])
rownames(mult.mtrx3b) <- c("Location2", "Location3", "Location4")


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
mtrx4[,-c(1:5)]
colnames(mtrx4)[c(3,5)] <- c("Het.Genet", "mid.shift")

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
#       ANOVA for one v one, glm for combined 
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

small <- correlations.data[, c(cas, imaging.features)]
small <- CleanYourRoom(small, "0?")
colnames(small) <- c("CAS", names(imaging.features))
remove <- c()
for (i in 1:ncol(small)){
    
    if (dim(table(small[, i])) < 2){
        remove <- c(remove, i)
    } 
    
}
small <- small[, -remove]
x <- summary(lm(CAS ~ ., data = small))
mult.mtrx7 <- matrix(x$coefficients[-1, 4], 12, 1)
rownames(mult.mtrx7) <- dimnames(x$coefficients)[[1]][-1]
mult.mtrx7

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

## Multivariate
small <- correlations.data[, c(mus, imaging.features)]
small <- CleanYourRoom(small, "0?")
colnames(small) <- c("MUS", names(imaging.features))
remove <- c()
for (i in 1:ncol(small)){
    
    if (dim(table(small[, i])) < 2){
        remove <- c(remove, i)
    } 
    
}
small <- small[, -remove]
x <- summary(lm(MUS ~ ., data = small))
mult.mtrx10 <- matrix(x$coefficients[-1, 4], 12, 1)
rownames(mult.mtrx10) <- dimnames(x$coefficients)[[1]][-1]


## 9 Location vs grade. Fisher's for count data comparison
mtrx11 <- matrix(NA, 2,2)
colnames(mtrx11) <- names(grades)
rownames(mtrx11) <- names(location)
mini <- CleanYourRoom(correlations.data[, c(location[2], grades[1])], "5")
tbl <- table(mini[, 1], mini[, 2])
mtrx11[2,1] <- chisq.test(tbl)$p.value
mtrx11
table(mini$Loclass2, mini$molgrade)



## Atypical features
mtrx12 <- matrix(NA,2,1)
rownames(mtrx12) <- c("MUS", "CAS")
colnames(mtrx12) <- "Atypical Gradient"
mini <- correlations.data[, c(mus, cas, atypical)]
bye.felicia <- function(val){
    if (is.na(val)){
        0
    }else{
        val
    }
}

mini <- apply(mini,2, function(x) sapply(x,bye.felicia))
mini <- cbind(mini, sum(mini[, 3:6]))
mini[, 7] <- rowSums(mini[, 3:6])
mtrx12[1,1] <- summary(lm(mini[, 1] ~ mini[, 7]))$coefficients[2,4]
mtrx12[2,1] <- summary(lm(mini[, 2] ~ mini[, 7]))$coefficients[2,4]

mtrx12


## PD-L1 expression vs Edema

mtrx13 <- matrix(NA, 1,3)
colnames(mtrx13) <- names(pdl)
rownames(mtrx13) <- names(imaging.features)[1]
for (i in 1:2){
    mini <- CleanYourRoom(correlations.data[, c(imaging.features[1], pdl[i])])
    math.party <- fisher.test(table(mini[, 1], mini[, 2]))
    mtrx13[1, i] <-math.party$p.value
}
mtrx13

##PDL-1 and grade
mtrx14 <- matrix(NA, 2,2)
colnames(mtrx14) <- names(pdl[-3])
rownames(mtrx14) <- names(grades)
for (i in 1:2){
    mini <- CleanYourRoom(correlations.data[, c(grades[1], pdl[i])])
    math.party <- fisher.test(table(mini[, 1], mini[, 2]))
    mtrx14[2, i] <-math.party$p.value
}
mtrx14

## PDL-1 and MUS

mtrx15 <- matrix(NA, 1, 2)
rownames(mtrx15) <- names(mus)
colnames(mtrx15) <- names(pdl[-3])
for (i in 1:2){
    mini <- CleanYourRoom(correlations.data[, c(mus, pdl[i])])
    if (length(unique(mini[, 2])) > 1){
        math.party <- aov(mini[, 1] ~ as.character(mini[, 2]))
        mtrx15[1, i] <- summary(math.party)[[1]][1,5]
    }
}
mtrx15


## Generates tables for paper
table(correlations.data[, imaging.features[1]])

## vignette("introduction", package = "dplyr")

correlations.data %>% select(grade, MUS)  %>% group_by(grade) %>% summarise(avgs = mean(MUS, na.rm = T))

correlations.data %>% select(Loclass2, MUS)  %>% group_by(Loclass2) %>% summarise(avgs = mean(MUS, na.rm = T))




















