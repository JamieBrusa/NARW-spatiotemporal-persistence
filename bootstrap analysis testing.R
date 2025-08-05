library(tidyverse)
library(here)
library(sf)
library(caret)

set.seed(123)

#Set up a training set (2/3) and a test set (1/3) and match up spatiotemporal sampling units for 3 folds

resurveys <- read.csv("resurveys.csv")

#boot.train.set <- createDataPartition(resurveys$redetects, times = 3, p = 2/3, list = FALSE)
boot.train.set <- createFolds(resurveys$redetects, k = 3)

#Manually iterate over each fold
boot.train.set1 <- as.vector(boot.train.set[[1]])
boot.train.set2 <- as.vector(boot.train.set[[2]])
boot.train.set3 <- as.vector(boot.train.set[[3]])

#Fold 1
boot.train1 <- resurveys[boot.train.set1, ]
boot.test1 <- resurveys[-boot.train.set1, ]

#Match up spatiotemporal sampling units
boot.train1 <- boot.train1[boot.train1$DetNum %in% boot.test1$DetNum,]
boot.test1 <- boot.test1[boot.test1$DetNum %in% boot.train1$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train1$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train1$detnum <- 0
for(i in dumb$detnum){
  boot.train1$detnum[boot.train1$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test1$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test1$detnum <- 0
for(i in dumb2$detnum){
  boot.test1$detnum[boot.test1$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train1$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out1 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train1$detnum)))
for(s in unique(boot.train1$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train1$detnum == s)
    train.samp <- sample(boot.train1$redetects[rows], size = length(rows), replace = TRUE)
    train.out1[i,s] <- mean(train.samp > 0) 
  }
}

train.means1 <- apply(train.out1, 2, mean)

trnm.fac1 <- as.factor(ifelse(train.means1 < 0.5, 0, 1))

iterations <- 10000
test.out1 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test1$detnum)))
for(s in unique(boot.test1$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test1$detnum == s)
    test.samp <- sample(boot.test1$redetects[rows], size = length(rows), replace = TRUE)
    test.out1[i,s] <- mean(test.samp > 0) 
  }
}

test.means1 <- apply(test.out1, 2, mean)
tstm.fac1 <- as.factor(ifelse(test.means1 < 0.5, 0, 1))

CM.boot1 <- confusionMatrix(trnm.fac1, tstm.fac1, positive = "1")


CM.bootAllYears1 <- confusionMatrix(trnm.fac1, tstm.fac1, positive = "1")
boot.AllYears.df1 <- as.data.frame(unlist(CM.bootAllYears1$overall))
colnames(boot.AllYears.df1) <- "All Years"
boot.AllYears.df1 <- as.data.frame(t(boot.AllYears.df1))
boot.AllYears.df1$Year <- "All Years"
boot.AllYears.df1$YearRange <- "Combined"

#Fold 2
boot.train2 <- resurveys[boot.train.set2, ]
boot.test2 <- resurveys[-boot.train.set2, ]

#Match up spatiotemporal sampling units
boot.train2 <- boot.train2[boot.train2$DetNum %in% boot.test2$DetNum,]
boot.test2 <- boot.test2[boot.test2$DetNum %in% boot.train2$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train2$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train2$detnum <- 0
for(i in dumb$detnum){
  boot.train2$detnum[boot.train2$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test2$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test2$detnum <- 0
for(i in dumb2$detnum){
  boot.test2$detnum[boot.test2$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train2$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out2 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train2$detnum)))
for(s in unique(boot.train2$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train2$detnum == s)
    train.samp <- sample(boot.train2$redetects[rows], size = length(rows), replace = TRUE)
    train.out2[i,s] <- mean(train.samp > 0) 
  }
}

train.means2 <- apply(train.out2, 2, mean)

trnm.fac2 <- as.factor(ifelse(train.means2 < 0.5, 0, 1))

iterations <- 10000
test.out2 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test2$detnum)))
for(s in unique(boot.test2$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test2$detnum == s)
    test.samp <- sample(boot.test2$redetects[rows], size = length(rows), replace = TRUE)
    test.out2[i,s] <- mean(test.samp > 0) 
  }
}

test.means2 <- apply(test.out2, 2, mean)
tstm.fac2 <- as.factor(ifelse(test.means2 < 0.5, 0, 1))

CM.boot2 <- confusionMatrix(trnm.fac2, tstm.fac2, positive = "1")


CM.bootAllYears2 <- confusionMatrix(trnm.fac2, tstm.fac2, positive = "1")
boot.AllYears.df2 <- as.data.frame(unlist(CM.bootAllYears2$overall))
colnames(boot.AllYears.df2) <- "All Years"
boot.AllYears.df2 <- as.data.frame(t(boot.AllYears.df2))
boot.AllYears.df2$Year <- "All Years"
boot.AllYears.df2$YearRange <- "Combined"

#Fold 3
boot.train3 <- resurveys[boot.train.set3, ]
boot.test3 <- resurveys[-boot.train.set3, ]

#Match up spatiotemporal sampling units
boot.train3 <- boot.train3[boot.train3$DetNum %in% boot.test3$DetNum,]
boot.test3 <- boot.test3[boot.test3$DetNum %in% boot.train3$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train3$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train3$detnum <- 0
for(i in dumb$detnum){
  boot.train3$detnum[boot.train3$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test3$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test3$detnum <- 0
for(i in dumb2$detnum){
  boot.test3$detnum[boot.test3$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train3$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out3 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train3$detnum)))
for(s in unique(boot.train3$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train3$detnum == s)
    train.samp <- sample(boot.train3$redetects[rows], size = length(rows), replace = TRUE)
    train.out3[i,s] <- mean(train.samp > 0) 
  }
}

train.means3 <- apply(train.out3, 2, mean)

trnm.fac3 <- as.factor(ifelse(train.means3 < 0.5, 0, 1))

iterations <- 10000
test.out3 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test3$detnum)))
for(s in unique(boot.test3$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test3$detnum == s)
    test.samp <- sample(boot.test3$redetects[rows], size = length(rows), replace = TRUE)
    test.out3[i,s] <- mean(test.samp > 0) 
  }
}

test.means3 <- apply(test.out3, 2, mean)
tstm.fac3 <- as.factor(ifelse(test.means3 < 0.5, 0, 1))

CM.boot3 <- confusionMatrix(trnm.fac3, tstm.fac3, positive = "1")


CM.bootAllYears3 <- confusionMatrix(trnm.fac3, tstm.fac3, positive = "1")
boot.AllYears.df3 <- as.data.frame(unlist(CM.bootAllYears3$overall))
colnames(boot.AllYears.df3) <- "All Years"
boot.AllYears.df3 <- as.data.frame(t(boot.AllYears.df3))
boot.AllYears.df3$Year <- "All Years"
boot.AllYears.df3$YearRange <- "Combined"


#Then, we cycle through the years of data collection to leave each year out and use all other years to predict the 
#left-out year. As before, the data are partitioned into a training set (2/3 of original bootstrap output data) and 
#a test set (1/3 of original bootstrap output data), which are balanced to have similar proportions of redetects as 0 
#or 1. A single year is removed from the training set, and all other years except that selected year are removed from 
#the test set. The training and test sets are then filtered to retain only spatiotemporal sampling units present in 
#both training and test sets.

#Set the filter one at a time to remove each year for an independent run
boot.train.set <- createDataPartition(y = resurveys$redetects, p = 2/3, list = FALSE)

boot.train.all <- resurveys[boot.train.set, ]
boot.test.all <- resurveys[-boot.train.set, ]

#Match up spatiotemporal sampling units
boot.train.all <- boot.train.all[boot.train.all$DetNum %in% boot.test.all$DetNum,]
boot.test.all <- boot.test.all[boot.test.all$DetNum %in% boot.train.all$DetNum,]

#Iterate across years
CM.boot <- list()
for(i in 2010:2020){
  boot.train <- boot.train.all %>% filter(Year != i)
  boot.test <- boot.test.all %>% filter(Year == i)
  
  #Match up spatiotemporal sampling units
  boot.train <- boot.train[boot.train$DetNum %in% boot.test$DetNum,]
  boot.test <- boot.test[boot.test$DetNum %in% boot.train$DetNum,]

  #Need to renumber the detnum for the for loop in the resampling approach to function
  detnum <- unique(boot.train$DetNum)
  dumb <- data.frame(detnum = detnum, numbers = 1:length(detnum))
  
  boot.train$detnum <- 0
  for(d in dumb$detnum){
    boot.train$detnum[boot.train$DetNum == d] <- match(d, dumb$detnum)
  }
  
  detnum2 <- unique(boot.test$DetNum)
  dumb2 <- data.frame(detnum2 = detnum2, numbers = 1:length(detnum2))
  
  boot.test$detnum <- 0
  for(d2 in dumb2$detnum){
    boot.test$detnum[boot.test$DetNum == d2] <- match(d2, dumb2$detnum)
  }
  
  detnum <- boot.train$detnum
  
  
  #Run bootstrap resampling
  iterations <- 10000
  train.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train$detnum)))
  for(s in unique(boot.train$detnum)){
    for(t in 1:iterations){
      rows <- which(boot.train$detnum == s)
      train.samp <- sample(boot.train$redetects[rows], size = length(rows), replace = TRUE)
      train.out[t,s] <- mean(train.samp > 0)
    }
  }
  
  train.means <- apply(train.out,2,mean)
  
  trnm.fac <- as.factor(ifelse(train.means < 0.5, 0, 1))
  
  iterations2 <- 10000
  test.out <- matrix(NA, nrow = iterations2, ncol=length(unique(boot.test$detnum)))
  for(s2 in unique(boot.test$detnum)){
    for(t2 in 1:iterations2){
      rows2 <- which(boot.test$detnum == s2)
      test.samp <- sample(boot.test$redetects[rows2], size = length(rows2), replace = TRUE)
      test.out[t2,s2] <- mean(test.samp > 0)
    }
  }
  
  test.means <- apply(test.out,2,mean)
  tstm.fac <- as.factor(ifelse(test.means < 0.5, 0, 1))
  
  
  
  CM.boot[[i]] <- confusionMatrix(trnm.fac, tstm.fac, positive = "1")
  
}

CM.boot.df <-
  as.data.frame(
    do.call(
      rbind,
      lapply(CM.boot,
             FUN = "[[",
             "overall"
      )
    )
  )

CM.boot.df$Year <- as.character(2010:2020)
CM.boot.df$YearRange <- "All"
boot.Years <- rbind(boot.AllYears.df, CM.boot.df)


#Finally, we cycle through the 2014-2020 data to leave each of those years out and use only the previous four years
#to predict the left-out year. As before, the data are partitioned into a training set (2/3 of original bootstrap 
#output data) and a test set (1/3 of original bootstrap output data), which are balanced to have similar proportions 
#of redetects as 0 or 1, prior to extracting specific years. The training and test sets are again filtered to retain
#only spatiotemporal sampling units present in both training and test sets.

boot.train.set <- createDataPartition(resurveys$redetects, p = 2/3, list = FALSE) 

boot.train.all <- resurveys[boot.train.set, ]
boot.test.all <- resurveys[-boot.train.set, ]

#Match up spatiotemporal sampling units
boot.train.all <- boot.train.all[boot.train.all$DetNum %in% boot.test.all$DetNum,]
boot.test.all <- boot.test.all[boot.test.all$DetNum %in% boot.train.all$DetNum,]

#Iterate across years
#For last3 and last5, change years to limit/extend to y-3 or y-5 (and adjust 2014:2020 to 2015:2020)
CM.5yrsboot <- list()
# for(y in 2014:2020){ #for CM.3yrsboot and CM.4yrsboot
  for(y in 2015:2020){ #for CM.5yrsboot
  # boot.train <- boot.train.all %>% filter(Year == (y-1) | Year == (y-2) | Year == (y-3)) #use for last CM.3yrsboot
  # boot.train <- boot.train.all %>% filter(Year == (y-1) | Year == (y-2) | Year == (y-3) | 
  #                                           Year == (y-4)) #use for last CM.4yrsboot
  boot.train <- boot.train.all %>% filter(Year == (y-1) | Year == (y-2) | Year == (y-3) |
                                            Year == (y-4) | Year == (y-5) ) #use for last CM.5yrsboot
  boot.test <- boot.test.all %>% filter(Year == y)
  
  #Match up spatiotemporal sampling units
  boot.train <- boot.train[boot.train$DetNum %in% boot.test$DetNum,]
  boot.test <- boot.test[boot.test$DetNum %in% boot.train$DetNum,]
  
  #Need to renumber the detnum for the for loop in the resampling approach to function
  detnum <- unique(boot.train$DetNum)
  dumb <- data.frame(detnum = detnum, numbers = 1:length(detnum))
  
  boot.train$detnum <- 0
  for(d in dumb$detnum){
    boot.train$detnum[boot.train$DetNum == d] <- match(d, dumb$detnum)
  }
  
  detnum2 <- unique(boot.test$DetNum)
  dumb2 <- data.frame(detnum2 = detnum2, numbers = 1:length(detnum2))
  
  boot.test$detnum <- 0
  for(d2 in dumb2$detnum){
    boot.test$detnum[boot.test$DetNum == d2] <- match(d2, dumb2$detnum)
  }
  
  detnum <- boot.train$detnum
  
  
  #Run bootstrap resampling
  iterations <- 10000
  train.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train$detnum)))
  for(s in unique(boot.train$detnum)){
    for(t in 1:iterations){
      rows <- which(boot.train$detnum == s)
      train.samp <- sample(boot.train$redetects[rows], size = length(rows), replace = TRUE)
      train.out[t,s] <- mean(train.samp > 0)
    }
  }
  
  train.means <- apply(train.out,2,mean)
  
  trnm.fac <- as.factor(ifelse(train.means < 0.5, 0, 1))
  
  iterations2 <- 10000
  test.out <- matrix(NA, nrow = iterations2, ncol=length(unique(boot.test$detnum)))
  for(s2 in unique(boot.test$detnum)){
    for(t2 in 1:iterations2){
      rows2 <- which(boot.test$detnum == s2)
      test.samp <- sample(boot.test$redetects[rows2], size = length(rows2), replace = TRUE)
      test.out[t2,s2] <- mean(test.samp > 0)
    }
  }
  
  test.means <- apply(test.out,2,mean)
  tstm.fac <- as.factor(ifelse(test.means < 0.5, 0, 1))
  
  
  
  CM.5yrsboot[[y]] <- confusionMatrix(trnm.fac, tstm.fac, positive = "1")
  
}

CM.5yrsboot.df <-
  as.data.frame(
    do.call(
      rbind,
      lapply(CM.5yrsboot,
             FUN = "[[",
             "overall"
      )
    )
  )

CM.5yrsboot.df$Year <- as.character(2014:2020) #2015:2020 for Last5 years
CM.5yrsboot.df$YearRange <- "Last5"
boot.Years <- rbind(boot.Years, CM.5yrsboot.df)
boot.Years[1, 8] <- "Combined"


#To try to pinpoint the regional variation, we looked at estimating accuracy for different strata (geographic regions).

filtered.df <- resurveys %>% filter(Stratum == "Nearshore SNE" | Stratum == "Midshore SNE") #Change for different regions

boot.train.set <- createFolds(filtered.df$redetects, k = 3)

#Manually iterate over each fold
boot.train.set1 <- as.vector(boot.train.set[[1]])
boot.train.set2 <- as.vector(boot.train.set[[2]])
boot.train.set3 <- as.vector(boot.train.set[[3]])

#Fold 1
boot.train1 <- resurveys[boot.train.set1, ]
boot.test1 <- resurveys[-boot.train.set1, ]

#Match up spatiotemporal sampling units
boot.train1 <- boot.train1[boot.train1$DetNum %in% boot.test1$DetNum,]
boot.test1 <- boot.test1[boot.test1$DetNum %in% boot.train1$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train1$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train1$detnum <- 0
for(i in dumb$detnum){
  boot.train1$detnum[boot.train1$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test1$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test1$detnum <- 0
for(i in dumb2$detnum){
  boot.test1$detnum[boot.test1$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train1$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out1 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train1$detnum)))
for(s in unique(boot.train1$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train1$detnum == s)
    train.samp <- sample(boot.train1$redetects[rows], size = length(rows), replace = TRUE)
    train.out1[i,s] <- mean(train.samp > 0) 
  }
}

train.means1 <- apply(train.out1, 2, mean)

trnm.fac1 <- as.factor(ifelse(train.means1 < 0.5, 0, 1))

iterations <- 10000
test.out1 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test1$detnum)))
for(s in unique(boot.test1$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test1$detnum == s)
    test.samp <- sample(boot.test1$redetects[rows], size = length(rows), replace = TRUE)
    test.out1[i,s] <- mean(test.samp > 0) 
  }
}

test.means1 <- apply(test.out1, 2, mean)
tstm.fac1 <- as.factor(ifelse(test.means1 < 0.5, 0, 1))

CM.boot1 <- confusionMatrix(trnm.fac1, tstm.fac1, positive = "1")


CM.SNE1 <- confusionMatrix(trnm.fac1, tstm.fac1, positive = "1")
boot.SNE.df1 <- as.data.frame(unlist(CM.SNE1$overall))
colnames(boot.SNE.df1) <- "All Years"
boot.SNE.df1 <- as.data.frame(t(boot.SNE.df1))
boot.SNE.df1$Year <- "All Years"
boot.SNE.df1$Region <- "Southern NE"

#Fold 2
boot.train2 <- resurveys[boot.train.set2, ]
boot.test2 <- resurveys[-boot.train.set2, ]

#Match up spatiotemporal sampling units
boot.train2 <- boot.train2[boot.train2$DetNum %in% boot.test2$DetNum,]
boot.test2 <- boot.test2[boot.test2$DetNum %in% boot.train2$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train2$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train2$detnum <- 0
for(i in dumb$detnum){
  boot.train2$detnum[boot.train2$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test2$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test2$detnum <- 0
for(i in dumb2$detnum){
  boot.test2$detnum[boot.test2$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train2$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out2 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train2$detnum)))
for(s in unique(boot.train2$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train2$detnum == s)
    train.samp <- sample(boot.train2$redetects[rows], size = length(rows), replace = TRUE)
    train.out2[i,s] <- mean(train.samp > 0) 
  }
}

train.means2 <- apply(train.out2, 2, mean)

trnm.fac2 <- as.factor(ifelse(train.means2 < 0.5, 0, 1))

iterations <- 10000
test.out2 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test2$detnum)))
for(s in unique(boot.test2$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test2$detnum == s)
    test.samp <- sample(boot.test2$redetects[rows], size = length(rows), replace = TRUE)
    test.out2[i,s] <- mean(test.samp > 0) 
  }
}

test.means2 <- apply(test.out2, 2, mean)
tstm.fac2 <- as.factor(ifelse(test.means2 < 0.5, 0, 1))

CM.boot2 <- confusionMatrix(trnm.fac2, tstm.fac2, positive = "1")


CM.SNE2 <- confusionMatrix(trnm.fac2, tstm.fac2, positive = "1")
boot.SNE.df2 <- as.data.frame(unlist(CM.SNE2$overall))
colnames(boot.SNE.df2) <- "All Years"
boot.SNE.df2 <- as.data.frame(t(boot.SNE.df2))
boot.SNE.df2$Year <- "All Years"
boot.SNE.df2$Region <- "Southern NE"

#Fold 3
boot.train3 <- resurveys[boot.train.set3, ]
boot.test3 <- resurveys[-boot.train.set3, ]

#Match up spatiotemporal sampling units
boot.train3 <- boot.train3[boot.train3$DetNum %in% boot.test3$DetNum,]
boot.test3 <- boot.test3[boot.test3$DetNum %in% boot.train3$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train3$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train3$detnum <- 0
for(i in dumb$detnum){
  boot.train3$detnum[boot.train3$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test3$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test3$detnum <- 0
for(i in dumb2$detnum){
  boot.test3$detnum[boot.test3$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train3$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out3 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train3$detnum)))
for(s in unique(boot.train3$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train3$detnum == s)
    train.samp <- sample(boot.train3$redetects[rows], size = length(rows), replace = TRUE)
    train.out3[i,s] <- mean(train.samp > 0) 
  }
}

train.means3 <- apply(train.out3, 2, mean)

trnm.fac3 <- as.factor(ifelse(train.means3 < 0.5, 0, 1))

iterations <- 10000
test.out3 <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test3$detnum)))
for(s in unique(boot.test3$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test3$detnum == s)
    test.samp <- sample(boot.test3$redetects[rows], size = length(rows), replace = TRUE)
    test.out3[i,s] <- mean(test.samp > 0) 
  }
}

test.means3 <- apply(test.out3, 2, mean)
tstm.fac3 <- as.factor(ifelse(test.means3 < 0.5, 0, 1))

CM.SNE3 <- confusionMatrix(trnm.fac3, tstm.fac3, positive = "1")


CM.SNE3 <- confusionMatrix(trnm.fac3, tstm.fac3, positive = "1")
boot.SNE.df3 <- as.data.frame(unlist(CM.SNE3$overall))
colnames(boot.SNE.df3) <- "All Years"
boot.SNE.df3 <- as.data.frame(t(boot.SNE.df3))
boot.SNE.df3$Year <- "All Years"
boot.SNE.df3$Region <- "Southern NE"


boot.Regions <- rbind(boot.SE.df, boot.SNE.df, boot.GoM.df, boot.CCB.df, boot.Carolinas.df, 
                      boot.AllYears.df)

#Note, Mid-Atlantic did not have enough data for regional analysis





