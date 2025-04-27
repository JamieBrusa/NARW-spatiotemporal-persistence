library(tidyverse)
library(here)
library(sf)
library(caret)

set.seed(123)

#Set up a training set (2/3) and a test set (1/3) and match up spatiotemporal sampling units

resurveys <- read.csv("resurveys.csv")

boot.train.set <- createDataPartition(resurveys$redetects, p = 2/3, list = FALSE)

boot.train <- resurveys[boot.train.set, ]
boot.test <- resurveys[-boot.train.set, ]

#Match up spatiotemporal sampling units
boot.train <- boot.train[boot.train$DetNum %in% boot.test$DetNum,]
boot.test <- boot.test[boot.test$DetNum %in% boot.train$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train$detnum <- 0
for(i in dumb$detnum){
  boot.train$detnum[boot.train$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test$detnum <- 0
for(i in dumb2$detnum){
  boot.test$detnum[boot.test$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train$detnum


#First, we evaluate the ability of the bootstrap resampling method that generates spatiotemporal persistence 
#estimates from the training data to predict the spatiotemporal persistence estimates generated through the same 
#bootstrap resampling approach using all years of data.

#Run bootstrap resampling
iterations <- 10000
train.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train$detnum)))
for(s in unique(boot.train$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train$detnum == s)
    train.samp <- sample(boot.train$redetects[rows], size = length(rows), replace = TRUE)
    train.out[i,s] <- mean(train.samp > 0) 
  }
}

train.means <- apply(train.out,2,mean)

trnm.fac <- as.factor(ifelse(train.means < 0.5, 0, 1))

iterations <- 10000
test.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test$detnum)))
for(s in unique(boot.test$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test$detnum == s)
    test.samp <- sample(boot.test$redetects[rows], size = length(rows), replace = TRUE)
    test.out[i,s] <- mean(test.samp > 0) 
  }
}

test.means <- apply(test.out,2,mean)
tstm.fac <- as.factor(ifelse(test.means < 0.5, 0, 1))

CM.boot <- confusionMatrix(trnm.fac, tstm.fac, positive = "1")


CM.bootAllYears <- confusionMatrix(trnm.fac, tstm.fac, positive = "1")
boot.AllYears.df <- as.data.frame(unlist(CM.bootAllYears$overall))
colnames(boot.AllYears.df) <- "All Years"
boot.AllYears.df <- as.data.frame(t(boot.AllYears.df))
boot.AllYears.df$Year <- "All Years"
boot.AllYears.df$YearRange <- "Combined"


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

boot.train.set <- createDataPartition(filtered.df$redetects, p = 2/3, list = FALSE)

boot.train <- filtered.df[boot.train.set, ]
boot.test <- filtered.df[-boot.train.set, ]

#Match up spatiotemporal sampling units
boot.train <- boot.train[boot.train$DetNum %in% boot.test$DetNum,]
boot.test <- boot.test[boot.test$DetNum %in% boot.train$DetNum,]

#Need to renumber the detnum for the for loop in the resampling approach to function
detnum <- unique(boot.train$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

boot.train$detnum <- 0
for(i in dumb$detnum){
  boot.train$detnum[boot.train$DetNum == i] <- match(i, dumb$detnum)
}

detnum2 <- unique(boot.test$DetNum)
dumb2 <- data.frame(detnum2 = detnum2,
                    numbers = 1:length(detnum2))

boot.test$detnum <- 0
for(i in dumb2$detnum){
  boot.test$detnum[boot.test$DetNum == i] <- match(i, dumb2$detnum)
}

detnum <- boot.train$detnum

#Run bootstrap resampling
iterations <- 10000
train.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.train$detnum)))
for(s in unique(boot.train$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.train$detnum == s)
    train.samp <- sample(boot.train$redetects[rows], size = length(rows), replace = TRUE)
    train.out[i,s] <- mean(train.samp > 0) 
  }
}

train.means <- apply(train.out,2,mean)

trnm.fac <- as.factor(ifelse(train.means < 0.5, 0, 1))

iterations <- 1000
test.out <- matrix(NA, nrow = iterations, ncol=length(unique(boot.test$detnum)))
for(s in unique(boot.test$detnum)){
  for(i in 1:iterations){
    rows <- which(boot.test$detnum == s)
    test.samp <- sample(boot.test$redetects[rows], size = length(rows), replace = TRUE)
    test.out[i,s] <- mean(test.samp > 0) 
  }
}

test.means <- apply(test.out,2,mean)
tstm.fac <- as.factor(ifelse(test.means < 0.5, 0, 1))


CM.SNE <- confusionMatrix(trnm.fac, tstm.fac, positive = "1")
boot.SNE.df <- as.data.frame(unlist(CM.SNE$overall))
colnames(boot.SNE.df) <- "All Years"
boot.SNE.df <- as.data.frame(t(boot.SNE.df))
boot.SNE.df$Year <- "All Years"
boot.SNE.df$Region <- "Southern NE"


boot.Regions <- rbind(boot.SE.df, boot.SNE.df, boot.GoM.df, boot.CCB.df, boot.Carolinas.df, 
                      boot.AllYears.df)

#Note, Mid-Atlantic did not have enough data for regional analysis





