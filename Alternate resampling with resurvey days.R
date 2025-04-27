library(tidyverse)
library(here)
library(sf)
library(sp)
library(stringr)

set.seed(123)

#Bring in survey data (JDY stands for Julian day-year)
#Note, coordinates are in epsg: 26718
survey.pts <- st_read('survey_pts.shp')

strata.pts <- st_centroid(survey.pts)
strata.r <- st_buffer(strata.pts, 20898)
strata.detections.all <- strata.r %>% filter(Detected == 1)
strata.detections.all$DetNum <- 1:nrow(strata.detections.all)
strata.detections <- strata.detections.all %>% dplyr::select(SVYnum, DetNum)
strata.det.df <- st_drop_geometry(strata.detections)
strata.r <- left_join(strata.r, strata.det.df, by = "SVYnum")


#for unique first detections, sum the whale counts and number of surveys within the next 7 day window for that buffer
#n and t are for the progress bar
n = 0
t = length(unique(na.omit(strata.r$DetNum)))
start <- Sys.time()
resurveys <- list()
samps <- rep(list(list("resurveys" = 0, "redetects" = 0)), t)
for(i in unique(na.omit(strata.r$DetNum))){
  n <- n+1
  jln <- unique(na.omit(strata.r$Julian[strata.r$DetNum == i]))
  
  d <- as.numeric((jln - 7):(jln + 7))
  d <- case_when(d == 0 ~ 365, d == -1 ~ 364, d == -2 ~ 363, d == -3 ~ 362, d == -4 ~ 361, d == -5 ~ 360, 
                 d == -6 ~ 359, TRUE ~ d)
  d <- case_when(d == 366 ~ 1, d == 367 ~ 2, d == 368 ~ 3, d == 369 ~ 4, d == 370 ~ 5, 
                 d == 371 ~ 6, TRUE ~ d)
  jdq <- strata.detections.all[strata.detections.all$Julian %in% d,] %>%
    distinct(DetNum, .keep_all = T)
  
  q <- jdq[st_intersects(jdq, strata.detections[strata.detections$DetNum == i,], sparse = F), ]
  
  qr <- length(unique(na.omit(q$DetNum)))
  pre_q2.l <- rep(list(list("resurveysr" = 0, "redetectsr" = 0)), qr)
  pre_q2 <- data.frame()
  for(r in 1:nrow(q)){
    
    wks <- unique(na.omit(q$JDY[r]))
    wks_mat <- sapply(wks, FUN = function(x){x + c(1:7)}, simplify = T)
    
    pre_qr <- survey.pts %>% filter(JDY %in% wks_mat)
    pre_q2r <- pre_qr[st_intersects(pre_qr, q[r,], sparse = F), ]
    
    #Calculate number of resurveys 0-7
    if(nrow(pre_q2r) < 1) {resurveysr <- 0
    }else{
      resurveysr <- length(unique(pre_q2r$JDY))}
    
    #Count the number of days 0-7 where a redetection occurred
    if(nrow(pre_q2r) < 1) {redetectsr <- 0
    }else{
      pre_redetectsr <- pre_q2r %>% group_by(JDY) %>% summarise(redetectsr = max(Detected)) 
      redetectsr <- sum(pre_redetectsr$redetectsr)}
    
    
    pre_q2.l[[r]]$resurveysr <- resurveysr
    pre_q2.l[[r]]$redetectsr <- redetectsr
    pre_q2.l[[r]]$JDY <- wks
  }
  pre_q2 <- as.data.frame(do.call(rbind, pre_q2.l))
  
  resurveys <- pre_q2$resurveysr
  redetects <- pre_q2$redetectsr
  
  entries <- nrow(pre_q2)
  samps[[i]]$Julian <- rep(jln, entries)
  samps[[i]]$JDY <- rep(wks, entries)
  samps[[i]]$resurveys <- as.numeric(resurveys)
  samps[[i]]$redetects <- as.numeric(redetects)
  samps[[i]]$DetNum <- rep(i,entries)
  
  #progress bar
  percent <- n / t * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (n == t){cat('\n')}
  
}
stop <- Sys.time()
#Warnings are just for st_intersection assuming attributes are constant across geometries

samps.df <- data.table::rbindlist(samps)

#Add a line for each resurvey and include each redetection (instead of summarized as a total)
#This alternative method takes into account variation in the number of resurveys across detection buffers
resurvs <- samps.df %>% filter(resurveys > 0)
resurvs$ID <- 1:nrow(resurvs)

resurvs_new <- as.data.frame(matrix(data = NA, nrow = 0, ncol = length(resurvs), 
                                    dimnames = list(NULL,names(resurvs))))
for(r in 1:nrow(resurvs)){
  n <- resurvs$resurveys[r]
  n2 <- resurvs$redetects[r]
  newrow <- resurvs[r,]
  newrow$redetects <- 0
  newrow2 <- resurvs[r,]
  newrow2$redetects <- n2
  if(n2 < 1){newrows <- as.data.frame(lapply(newrow, rep, n-1))
  }else{
    newrows <- as.data.frame(lapply(newrow, rep, n-n2))
  }
  if(n2 < 1){newrows2 <- NULL
  }else{newrows2 <- as.data.frame(lapply(newrow2, rep, n2-1))
  }
  resurvs_new <- rbind(resurvs_new, resurvs[r,], newrows, newrows2)
}

resurvs_new$redetects <- ifelse(resurvs_new$redetects > 0, 1, 0)

colnames(resurvs_new)[4] <- "JDY_FD" #This field is the JDY from the first detection


#Need to renumber DetNum to eliminate any gaps for the resampling for loop
detnum <- unique(resurvs_new$DetNum)
dumb <- data.frame(detnum = detnum,
                   numbers = 1:length(detnum))

resurvs_new$detnum <- 0
for(i in dumb$detnum){
  resurvs_new$detnum[resurvs_new$DetNum == i] <- match(i, dumb$detnum)
}


##Resampling process##

iterations <- 10000
boot.out <- matrix(NA,nrow = iterations,ncol=length(unique(resurvs_new$detnum)))
for(s in unique(resurvs_new$detnum)){
  for(i in 1:iterations){
    rows <- which(resurvs_new$detnum == s)
    boot.samp <- sample(resurvs_new$redetects[rows], size = length(rows), replace = TRUE)
    boot.out[i,s] <- mean(boot.samp > 0) 
  }
}


#organize output from each analysis
boot.means <- apply(boot.out,2,mean)
boot.lcl <- apply(boot.out,2,function(x) quantile(x,probs=c(0.025)))
boot.ucl <- apply(boot.out,2,function(x) quantile(x,probs=c(0.975)))

#Put output into a dataframe
DetNum <- unique(resurvs_new$DetNum)
boot.out.df <- as.data.frame(cbind(DetNum, boot.means, boot.lcl, boot.ucl))

##Bring in additional information##
FirstDets <- resurvs_new[,3:5]
FirstDets <- dplyr::distinct(FirstDets, .keep_all = TRUE)

boot.out.df <- left_join(FirstDets, boot.out.df, by = "DetNum")



