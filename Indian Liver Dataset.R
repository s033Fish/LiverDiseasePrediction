if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(factoextra)) install.packages("factoextra", repos = "http://cran.us.r-project.org")
if(!require(corrplot)) install.packages("corrplot", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(clusterSim)) install.packages("clusterSim", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(dslabs)) install.packages("dslabs", repos = "http://cran.us.r-project.org")
if(!require(mltools)) install.packages("mltools", repos = "http://cran.us.r-project.org")
if(!require(recosystem)) install.packages("recosystem", repos = "http://cran.us.r-project.org")
if(!require(RCurl)) install.packages("RCurl", repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")
if(!require(gam)) install.packages("gam", repos = "http://cran.us.r-project.org")


#DOWNLOADING DATA
#If you have kaggle-cli installed run this command in the terminal:
#  kaggle datasets download -d uciml/indian-liver-patient-records
#If you don't have kaggle-cli installed run these commands 

#Set your browsing links
dataurl  = "https://drive.google.com/uc?export=download&id=11Ig8V33dHkjekqLIrVVklUBJt8A_ZXJO"

#Downloading the data
td = tempdir()
dat = tempfile(tmpdir=td, fileext=".rds")
download.file(dataurl, dat, mode="wb")
dat = read.csv(dat)
unlink(td)
head(dat)

#Cleaning the Data

#Removed all NAs
dat <- na.omit(dat)

#Changed all of the 2s to 0s
dat$Dataset <- gsub(2,0,dat$Dataset)

#Graphing the distribution of variables
tb <- ggplot(dat, aes(seq_along(Total_Bilirubin), Total_Bilirubin, color=Dataset)) +
  geom_point()

db <- ggplot(dat, aes(seq_along(Direct_Bilirubin), Direct_Bilirubin, color=Dataset)) +
  geom_point()

ap <- ggplot(dat, aes(seq_along(Alkaline_Phosphotase), Alkaline_Phosphotase, color=Dataset)) +
  geom_point()

aa <- ggplot(dat, aes(seq_along(Alamine_Aminotransferase), Alamine_Aminotransferase, color=Dataset)) +
  geom_point()

asa <- ggplot(dat, aes(seq_along(Aspartate_Aminotransferase), Aspartate_Aminotransferase, color=Dataset)) +
  geom_point()

tp <- ggplot(dat, aes(seq_along(Total_Protiens), Total_Protiens, color=Dataset)) +
  geom_point()

a <- ggplot(dat, aes(seq_along(Albumin), Albumin, color=Dataset)) +
  geom_point()

ag <- ggplot(dat, aes(seq_along(Albumin_and_Globulin_Ratio), Albumin_and_Globulin_Ratio , color=Dataset)) +
  geom_point()

#Plotting all of the plots together
ggarrange(tb, db, ap, aa, asa, tp, a, ag,
          ncol = 2, nrow = 4)

#CREATING TEST AND TRAIN SETS
y <- dat$Dataset
set.seed(1)
RNGkind(sample.kind = "Rejection")
test_index <- createDataPartition(y, times = 1, p = 0.1, list = FALSE)
test_set <- dat[test_index, ]
train_set <- dat[-test_index, ]


#GLM Training

#Fitting the model
train_glm <- train(y = train_set[,11], x = train_set[, -11], method = "glm")

#Using the model to make predictions
glm_preds <- predict(train_glm, test_set[, -11])

#Adding the accuracy to a dataframe
acc_table <- data.frame(method = "glm", acc = mean(glm_preds == test_set$Dataset))
acc_table

#Setting the control settings
ctrl <- trainControl(method="repeatedcv",repeats = 3) 


#QDA

#Fitting the model
train_qda <- train(Dataset ~ ., data = train_set, method = "qda", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)

#Using the model to make predictions
qda_preds <- predict(train_qda, test_set[, -11])

#Adding the accuracy to a dataframe
acc_table <- bind_rows(acc_table, data_frame(method="qda", acc = mean(qda_preds == test_set$Dataset)))
acc_table

#KNN

#Fitting the model
train_KNN <- train(Dataset ~ ., data = train_set, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)

#Using the model to make predictions
knn_preds <- predict(train_KNN, test_set)

#Adding the accuracy to the dataframe
acc_table <- bind_rows(acc_table, data_frame(method="knn", acc = mean(knn_preds == test_set$Dataset)))
acc_table

#LDA

#Fitting the model
train_lda <- train(Dataset ~ ., data = train_set, method = "lda", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)

#Using the model to make predictions
lda_preds <- predict(train_lda, test_set[, -11])

#Adding the accuracies to the dataframe
acc_table <- bind_rows(acc_table, data_frame(method="lda", acc = mean(lda_preds == test_set$Dataset)))
acc_table

#GAMLOESS

#Fitting the model
train_gam <- train(Dataset ~ ., data = train_set, method = "gamLoess", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)

#Using the model to make predictions
gam_preds <- predict(train_gam, test_set[, -11])

#Adding the accuracy to the dataframe
acc_table <- bind_rows(acc_table, data_frame(method="gamLoess", acc = mean(gam_preds == test_set$Dataset)))
acc_table

#RF

#Fitting the model
train_rf <- train(Dataset ~ ., data = train_set,method = "rf", tuneGrid = data.frame(mtry = 3), importance = TRUE)

#using the model to make predictions
rf_preds <- predict(train_rf, test_set[, -11])
mean(rf_preds == test_set$Dataset)

#Adding the accuracy to the table
acc_table <- bind_rows(acc_table, data_frame(method="rf", acc = mean(rf_preds == test_set$Dataset)))
acc_table

#ENSEMBLE

#Creating a dataframe with the predictions from the best methods
ensemble <- data.frame(lda_preds, glm_preds, gam_preds, knn_preds, rf_preds)

#Creating a function to find the mean
modefunc <- function(x){
  tabresult <- tabulate(x)
  themode <- which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1) themode <- NA
  return(themode)
}

#Converting the 0s to 2s
ensembleNum<- ifelse(ensemble == 0, 2, 1)

#Finding the modes of all predictions
modes<-apply(ensembleNum, 1, modefunc)

#Converting back to 0s
ensemble$mode <- ifelse(modes == 2, 0, 1)
          
#Adding accuracy to the table
acc_table <- bind_rows(acc_table, data_frame(method="ensemble", acc = mean(ensemble$mode == test_set$Dataset)))
acc_table

