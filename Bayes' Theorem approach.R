library(tidyverse)
library(here)
library(sf)
library(sp)
library(RColorBrewer)
library(ggmap)
library(stringr)
library(cluster)
library(plotly)
library(raster)


set.seed(123)

#Bring in survey data (JDY stands for Julian day-year)
#Note, coordinates are in epsg: 26718
BTsurvey.pts <- st_read('BTsurvey_pts.shp')


BT.pts <- st_centroid(survey.pts)
BT.r <- st_buffer(BT.pts, 20898)
BT.detections.all <- BT.r %>% filter(Detected == 1)
BT.detections.all$DetNum <- 1:nrow(BT.detections.all)
BT.detections <- BT.detections.all %>% dplyr::select(SVYnum, DetNum)
BT.det.df <- st_drop_geometry(BT.detections)
BT.r <- left_join(BT.r, BT.det.df, by = "SVYnum")



#for unique first detections, sum the whale counts and number of surveys within the next 7 day window for that buffer
#n and t are for the progress bar
n = 0
t = length(unique(na.omit(BT.r$DetNum)))
start <- Sys.time()
resurveys <- list()
probs <- rep(list(list("resurveys" = 0, "redetects" = 0, "prob.second" = 0)), t)
for(i in unique(na.omit(BT.r$DetNum))){
  n <- n+1
  jln <- unique(na.omit(BT.r$Julian[BT.r$DetNum == i]))
  
  d <- as.numeric((jln - 7):(jln + 7))
  d <- case_when(d == 0 ~ 365, d == -1 ~ 364, d == -2 ~ 363, d == -3 ~ 362, d == -4 ~ 361, d == -5 ~ 360, 
                 d == -6 ~ 359, TRUE ~ d)
  d <- case_when(d == 366 ~ 1, d == 367 ~ 2, d == 368 ~ 3, d == 369 ~ 4, d == 370 ~ 5, 
                 d == 371 ~ 6, TRUE ~ d)
  jdq <- BT.detections.all[BT.detections.all$Julian %in% d,] %>%
    distinct(JDY, .keep_all = T)
  
  q <- jdq[st_intersects(jdq, BT.detections[BT.detections$DetNum == i,], sparse = F), ]
  
  qr <- length(unique(na.omit(q$DetNum)))
  pre_q2.l <- rep(list(list("resurveysr" = 0, "redetectsr" = 0)), qr)
  pre_q2 <- data.frame()
  for(r in 1:nrow(q)){
  
    wks <- unique(na.omit(q$JDY[r]))
    wks_mat <- sapply(wks, FUN = function(x){x + c(1:7)}, simplify = T)
  
    pre_qr <- BTsurvey.pts %>% filter(JDY %in% wks_mat)
    pre_q2r <- pre_qr[st_intersects(pre_qr, q[r,], sparse = F), ]
    
    #Calculate number of resurveys 0-7
    if(nrow(pre_q2r) < 1) {resurveysr <- 0
    }else{
      resurveysr <- length(unique(pre_q2r$JDY))} 
    
    
    #Retain correction factor per observation (informative prior)
    if(nrow(pre_q2r) < 1) {pre_q2r <- pre_q2r
    }else{
      num.dets <- pre_q2r[pre_q2r$Detected == 1,]
      zeros <- pre_q2r[pre_q2r$Detected == 0,]
      if(nrow(num.dets) < 1){
        zeros$Correction <- 0
      }else{
        for(z in 1:nrow(zeros)){
        zeros$Correction[z] <- sample(num.dets$Correction, 1)}
      }
      ones <- pre_q2r[pre_q2r$Detected == 1,]
      pre_q2r <- rbind(zeros, ones)}

    #Count the number of days 0-7 where a redetection occurred
    #Informative prior
    #Borrowing information from redetections in the first detection buffer. If there are none, then default to 0.
    if(nrow(pre_q2r) < 1) {redetectsr <- 0
    }else{
      pre_redetectsr <- pre_q2r %>% group_by(JDY) %>% summarise(redetectsr = max(Detected), 
                                                                Correction = sample(Correction, 1))
      zeros.d <- pre_redetectsr[pre_redetectsr$redetectsr == 0,]
      if(nrow(zeros.d) < 1){
        redetectsr <- sum(pre_redetectsr$redetectsr)
        }else{
            for(zdl in 1:length(unique(zeros.d$JDY))){
              zeros.d$redetectsr[zdl] <- zeros.d$Correction[zdl]}

          }
      ones.d <- pre_redetectsr[pre_redetectsr$redetectsr == 1,]
      zeros.d$redetectsr <- ifelse(zeros.d$redetectsr < mean(pre_q2r$Correction), 1, 0)
      new.rd.r <- rbind(zeros.d, ones.d)
      new.rd.r <- new.rd.r %>% group_by(JDY) %>% summarise(redetectsr = max(redetectsr))
      redetectsr <- sum(new.rd.r$redetectsr)}
    
    
    #Not an issue with this dataset, but need to add if statement if nrow(q)=0 could occur (set pre_q2 <- q)
    pre_q2.l[[r]]$resurveysr <- resurveysr
    pre_q2.l[[r]]$redetectsr <- redetectsr
    pre_q2.l[[r]]$JDY <- wks
  }
  pre_q2 <- as.data.frame(do.call(rbind, pre_q2.l))
  
  resurveys <- pre_q2$resurveysr
  redetects <- pre_q2$redetectsr
  

  entries <- nrow(pre_q2)
  probs[[i]]$Julian <- rep(jln, entries)
  probs[[i]]$JDY <- rep(wks, entries)
  probs[[i]]$resurveys <- as.numeric(resurveys)
  probs[[i]]$redetects <- as.numeric(redetects)
  probs[[i]]$prob.second <- (as.numeric(probs[[i]]$redetects) / as.numeric(probs[[i]]$resurveys))
  probs[[i]]$DetNum <- rep(i,entries)
  
  #progress bar
  percent <- n / t * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (n == t){cat('\n')}
  
}
stop <- Sys.time()
#Warnings are just for st_intersection assuming attributes are constant across geometries

probs.df <- data.table::rbindlist(probs)
probs.df$resurveys <- as.numeric(probs.df$resurveys)
probs.df$redetects <- as.numeric(probs.df$redetects)

#Adjust 0 detections according to proportions of corrections for detections in that first detection buffer


#Perform Bayes' Theorem
#pr(2nd det | 1st det) = (pr(1st det | 2nd det) * pr(2nd det)) / pr(1st det)
#pr(1st det | 2nd det) = 1 OR (number redetects in db / total number redetects)
#pr(2nd det) = prob.second  (number redetects in db / number resurveys in db)
#pr(1st det) = 1 OR prob.first (number detections in d / total number detections)

#In this case, pr(first detection) = 1, and pr(detection | redetection) = 1
probs.df$prob.2g1 <- (1 * probs.df$prob.second) / 1

#remove any lines where there were no resurveys
probs.df <- probs.df %>% filter(resurveys != 0)

surv.means.df <- probs.df %>% group_by(DetNum) %>%
  summarise(Julian = mean(Julian), resurveys = mean(resurveys),
            JDY = mean(JDY), prob.2g1 = mean(prob.2g1))




