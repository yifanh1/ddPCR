cloudy <- function (drp, dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE, threshold = NA){
  findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))){  #where bw = is box width, setting the sensitivity of the search
    ###set all vectors to null
    pos.x.max <- NULL ;	pos.y.max <- NULL ;	pos.x.min <- NULL ;	pos.y.min <- NULL
    ###Start of for loop:    we walk down the vector with a window of size "bw"
    for(i in 1:(length(vec)-1)){
      #check if we have reached the end of the vector
      if((i+1+bw)>length(vec)){sup.stop <- length(vec)}else{sup.stop <- i+1+bw}
      #check if we are at beginning of the vector
      if((i-bw) < 1){inf.stop <- 1}else{inf.stop <- i-bw}
      #select window in two parts: values beyond i (superior), and values before i (inferior)
      subset.sup <- vec[(i+1):sup.stop]
      subset.inf <- vec[inf.stop:(i-1)]
      ##############################################################
      #are ALL trailing data smaller than i?
      is.max   <- sum(subset.inf > vec[i]) == 0
      #are ALL leading data smaller than i?
      is.nomin <- sum(subset.sup > vec[i]) == 0
      #are ALL trailing data larger than i?
      no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
      #are ALL leading data larger than i?
      no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
      ##############################################################
      #a maximum is found if  all data before and after i are smaller than i
      if(is.max & is.nomin){
        pos.x.max <- c(pos.x.max, x.coo[i])
        pos.y.max <- c(pos.y.max, vec[i])
      }
      #a maximum is found if  all data before and after i are larger than i
      if(no.max & no.nomin){
        pos.x.min <- c(pos.x.min, x.coo[i])
        pos.y.min <- c(pos.y.min, vec[i])
      }
    }#end of for loop
    ###Output
    return(list("max.X" = pos.x.max, "max.Y" = pos.y.max, "min.X" = pos.x.min, "min.Y" = pos.y.min))}
  require(SuppDists)     #required for the 'moments()' function
  #aim: to change the SD Algr--> get every SD for each peak

  drp <- na.omit(drp)
  fail <- 0 #create failure tracker
  pops <- 0 #create population tracker
  #be sure bandwidth for kernel density is at least 50
  temp <- bw.nrd0(drp)
  if(temp < 50){
    temp <- 50
  }
  #kernel density estimation
  krn <- density(drp, bw = temp)
  
  plot(krn)
  
  krn <- rbind(krn$y, krn$x)
  #we use findpeaks with a high bandwith to single out the populations
  piek <- findpeaks(krn[1, ], bw=20)$max.X
 
  piek <- rbind(piek, krn[1, piek])    #add peak heights
  
  piek <- rbind(piek, krn[2, piek[1,]])#add peak x-locations
  
  
  
  #we also remove all peaks that are smaller than 1% of the max peak height.
  #Outlying fluorescence values will otherwise cause insignificantly small peaks that screw with algorithm robustness
  if(any(piek[2,] < max(piek[2,])/3000)){
    piek <- piek[, -which(piek[2,] < max(piek[2,])/3000)]
    #make sure shit didn't vectorize
    if(is.null(dim(piek))){
      piek <- as.matrix(piek)
    }
  } 

 dal<-NA
  #if there is NA threshold.

  if (is.na(threshold)){
    pops <- length(piek[1,])
    #find dale between peaks!
    pops.1<-pops-1
    for (i in 1:pops.1){
      dal<-rbind(dal,(piek[3,i]+piek[3,i+1])/2)
    }
    
  
  #check the min number
  dal<-dal[-1]

    
  
  print("pops")
  print(pops)

  #get started on the main population matrix
  piik <- piek[3, ]
  print("No threshold entered, automaticly find the following threshold(s)")
  print(dal)
  print("peaks:")
  print(piik)
  
  peaksdivider <- dal
  }
  else{
    if (length(threshold)>=1){
      print("pops")
      pops<-length(threshold)+1
      print(pops)
      
      print("Threshold entered:")
      print(threshold)
      peaksdivider<-threshold
      
      piik <- piek[3, ]
      print("peaks:")
      print(piik)
    }
  }
  
  minnumber<-length(peaksdivider)
  print("min_number")
  print(minnumber)
  #We find the peak base by first estimaging the standard deviation as 1/2 of the width of the peak at 0.6065 percent of the max height
  
  #first we find the y-location (60.65% height) for both peaks
  
  baes <- (dnorm(1)/dnorm(0)) * piek[2, ]
  
  print(peaksdivider)
  BC <-c(NA,NA)

  for (i in 1:minnumber){

    candigroup1 <- head(krn[2, which(krn[2, ] < peaksdivider[i])][order((krn[1, which(krn[2, ] < peaksdivider[i])] - baes[i])^2)] , n = 10) #select 10 possible crossing candidates in the first half of the data (necessary as the data are discrete)
    
    temp <- c(candigroup1[which(candigroup1 < piik[i])][order(candigroup1[which(candigroup1 < piik[i])]-piik[i], decreasing = T)[1]],
            
            candigroup1[which(candigroup1 > piik[i])][order(candigroup1[which(candigroup1 > piik[i])]-piik[i])[1]])#take the ones nearest to the peak (one larger and one smaller)
    
    if(any(is.na(temp))){ # for severely overlapping peaks, only one crossing may be found. If that's the case, we assume symmetry and 'flip' one boundry around the peak centre
      
      
      temp[is.na(temp)] <- piik[i] + piik[i] - temp[!is.na()]
      
    }# END if
    BC<-rbind(BC,temp)
  }
  
  temp <- head(krn[2, which(krn[2, ] > peaksdivider[i])][order((krn[1, which(krn[2, ] > peaksdivider[i])] - baes[i+1])^2)], n = 10) #select 10 possible crossing candidates in the second half of the data
  
  temp <- c(temp[which(temp < piik[i+1])][order(temp[which(temp < piik[i+1])] - piik[i+1], decreasing = T)[1]],
            
            temp[which(temp > piik[i+1])][order(temp[which(temp > piik[i+1])] - piik[i+1])[1]])#take the ones nearest to the peak (one larger and one smaller)
  
  if(any(is.na(temp))){ # for severely overlapping peaks, only one crossing may be found. If that's the case, we assume symmetry and 'flip' one boundry around the peak centre
    
    temp[is.na(temp)] <- piik[i+1] + piik[i+1] - temp[!is.na(temp)]
    
  }# END if
  
  BC <- rbind(BC,temp)
  

  #Expand piik to take all values
  
  piik <- rbind (piik, matrix(c(rep(NA,2*pops)),nrow=2,ncol=pops))
  
  #The higher the kurtosis of the clouds the more sigmas we need to cover 99% of all the droplets, we start with standard values (4) and update below
  
  s<-c(rep(4,pops))
  k<-c(rep(NA,pops))
  #let's populate the population matrix with the population boundaries
  
  for (i in 1:pops){
    piik[2,i]<-piik[1,i]-diff(BC[i+1,])/2*s[i]
    piik[3,i]<-piik[1,i]+diff(BC[i+1,])/2*s[i]
  }

   
  #allright, using this rough delimination of the clouds we can calculate more accurate values
  
  #by iterating the update process three times (i.e. use current estimate to calculate pos and neg,
  
  # then recaculte sd and re-delimite bands). Three iterations is generaly enough to reach
  
  # stability in the first decimal place.
 
  for ( i  in 1 : 5 ) {
 #split 
      drpg1<-drp[which(drp>=piik[2,1]& drp<=piik[3,1])]
      drpg2<-drp[which(drp>=piik[2,2]&drp<=piik[3,2])]
      drpg<-list(drpg1,drpg2)
      if (pops-2>0){drpg3<-drp[which(drp>=piik[2,3]&drp<=piik[3,3])]
      drpg<-list(drpg1,drpg2,drpg3)}
      if (pops-3>0){drpg4<-drp[which(drp>=piik[2,4]&drp<=piik[3,4])]
      drpg<-list(drpg1,drpg2,drpg3,drpg4)}
   
    #calculate kurtosis of droplets
    for(j in 1:pops){
      ppp<-as.vector(unlist(drpg[j], use.names=FALSE))
      k[j]<-abs(moments(ppp)[4])
    
    }

    
    #update s multiplier
    for (j in 1:pops){
      s[j]<-3.8+0.35*log(k[j])+0.045*log(k[j])^2+0.75
    }

    
    #update
    for (j in 1:pops){
      ppp<-as.vector(unlist(drpg[j], use.names=FALSE))
      piik[2,j]<-piik[1,j]-mad(ppp)*s[j]
      piik[3,j]<-piik[1,j]+mad(ppp)*s[j]
      
    }
 #advance the cycle
  }#END of the iterative update
  
 
  #Extra check: severly overlapping populations may cause NA values in piik. NAs are handled badly by the code below
  
  # therefore we replace them by zeroes, this way the error catching below will kick in appropriatly.
  
  if(any(is.na(piik))){
    
    piik[which(is.na(piik))] <- 0
    
    if(!isTRUE(silent)){message("WARNING!: errors in finding population boundaries")
      
      message("this may affect threshold placement")
      
      message("---------")}
    
  }

  ###Calculation of Performance parameters
#  Resolution
  if (pops==2){
    Resol <- 2 * (piik[1, 2] -piik[1, 1]) / ((piik[3, 1] - piik[2, 1]) + (piik[3,2] - piik[2, 2]))
    print('Resol')
    print(Resol)}
    
  else{
    if (pops==1){message("1 population!")}
    else{
    R1<-2 * (piik[1, 2] -piik[1, 1]) / ((piik[3, 1] - piik[2, 1]) + (piik[3,2] - piik[2, 2]))
    R2<-2*(piik[1,3]-piik[1,2])/(piik[3,2]+piik[3,3]-piik[2,2]-piik[2,3])}  
    if (pops==4){
      R3<-2*(piik[1,4]-piik[1,3])/(piik[3,3]+piik[3,4]-piik[2,3]-piik[2,4])
      print('R3')
      print(as.vector(R3))
    }
    print("R1")
    print(as.vector(R1))
    print("R2")
    print(as.vector(R2))}
    
}
  


##################################################################################