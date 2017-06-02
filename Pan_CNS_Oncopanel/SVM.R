### SVM analysis for classification of tumor types
## code taken directly from SVM vigenette
## https://cran.r-project.org/web/packages/e1071/vignettes/svmdoc.pdf

install.packages("e1071")
install.packages("rpart")
library(e1071)
library(rpart)
install.packages("mlbench")
data(Glass, package="mlbench")

## split data into a train and test set
index <- 1:nrow(Glass)
testindex <- sample(index, trunc(length(index)/3))
testset <- Glass[testindex,]
trainset <- Glass[-testindex,]
    
## fit the model
svm.model <- svm(Type ~ ., data = trainset, cost = 100, gamma = 1)
svm.pred <- predict(svm.model, testset[,-10])
table(pred = svm.pred, true = testset[,10])
