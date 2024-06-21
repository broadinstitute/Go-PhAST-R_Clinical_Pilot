
#color scheme
myPalette <- colorRampPalette(c("blue", "white", "red"))(n=299)

#list for data frames
output_list <- list()

#functions
colGeoMeans <- function(myData, na.rm = FALSE) {  ## calculates geometric mean of column of data
  if(na.rm == TRUE){
    myData = log(myData)
    myData[is.na(myData)] = 0
    geoMeanData <- exp( colMeans (myData) )
  }else{
    geoMeanData <- exp( colMeans( log(myData) ) )
  }
  return(geoMeanData)
}

posCtrlCorr <- function(myData) {  ## positive control correction on data (scales by GeoMean)
  # calc GeoMean per lane
  posCtrlGeoMeans <- myData %>%
    filter(grepl("POS", rownames(myData)))  %>%
    colGeoMeans()
  # scale GeoMean by avg GeoMean to get scale factor per column
  colScale <- posCtrlGeoMeans / mean(posCtrlGeoMeans)
  # divide each data point in column by this scale factor
  posCtrlCorrData <- sweep(myData, 2, colScale, FUN = "/")
  return(posCtrlCorrData)
}

negCtrlCorr <- function(myData, round.to.zero = F) {  ## negative control correction on data (subtracts avg)
  # select neg ctrl genes and find col mean
  negCtrlMeans <- myData %>%
    filter(grepl("NEG", rownames(myData)))  %>%
    colMeans(na.rm = T)
  # perform negative control on non-control data
  negCtrlCorrData <- myData %>%
    filter(!( grepl("NEG", rownames(myData)) | grepl("POS", rownames(myData)) )) %>%
    sweep(2, negCtrlMeans)
  if(round.to.zero){
    negCtrlCorrData[is.neg.df(negCtrlCorrData)] = 0
  }
  return(negCtrlCorrData)
}

negCtrl_wReplacement <- function(myData, background.subtraction = F) {  ## negative control correction on data (subtracts avg)
  # select neg ctrl genes and find col mean
  negCtrlMeans <- myData %>%
    filter(grepl("NEG", rownames(myData)))  %>%
    colMeans()
  
  if(background.subtraction){
    # background subtract
    # then all numbers below the background value are rounded to that value
    # meaning that probes have input values of double the background signal or they are rounded up
    negCtrlCorrData <- myData %>%
      filter(!( grepl("NEG", rownames(myData)) | grepl("POS", rownames(myData)) )) %>%
      sweep(2, negCtrlMeans)
    for(x in 1:length(negCtrlMeans)){
      negCtrlCorrData[x][negCtrlCorrData[x]<negCtrlMeans[x]]=negCtrlMeans[x]
    }
  }else{
    # no background subtraction
    # all numbers below the background value are rounded to that value
    negCtrlCorrData <-  myData %>% filter(!( grepl("NEG", rownames(myData)) | grepl("POS", rownames(myData)) ))
    for(x in 1:length(negCtrlMeans)){
      negCtrlCorrData[x][negCtrlCorrData[x]<negCtrlMeans[x]]=negCtrlMeans[x]
    }
  }

  return(negCtrlCorrData)
}

removeFailedProbes <- function(myData, minVal = 50, rm = T, firstRow = F, na.rm = F) { 
  
  # changes NA values already in the data to INF so they are not removed
  if(na.rm){
    myData[is.na.df(myData)] = Inf
  }
  ## removes a probe if ANY value in ANY sample is < minVal, default to 50 but typically use 100 for Ctrls and 40 for Resp
  #if firstRow is true failed probes are only removed from the first row
  if(firstRow){
    myData[1][(myData[1]) < minVal] <- NA
    }else{
    myData[myData < minVal] <- NA
    }
  #myData[,2] = as.numeric(myData[,2])
  
  ## option to remove probes that are below minVal
  if(rm) {
    myData = na.omit(myData)
    if(na.rm){
      myData[is.inf.df(myData)] = NA
    }
    return(myData)
  }
  
  
  if(na.rm){
    myData[is.inf.df(myData)] = NA
  }
  return(myData)
}

calcCtrlCoVs <- function(ctrlProbes) {
  # function that calculates coefficients of variation for control probes across samples, based on ratio of each probe to GeoMean of all probes
  ctrlProbeRatio <- sweep(ctrlProbes,2,colGeoMeans(ctrlProbes),FUN = "/")
  ctrlProbeRatioSDs <- apply(ctrlProbeRatio,1,sd)
  ctrlProbeRatioMeans <- rowMeans(ctrlProbeRatio) 
  ctrlProbeRatioCoVs <- ctrlProbeRatioSDs / ctrlProbeRatioMeans
  return(ctrlProbeRatioCoVs)
}

removeCtrlMaxCoV <- function(ctrlProbesCurrent, probeToRemove = "none") {
  # function to remove the Ctrl probe with the max CoV as a ratio of GeoMean of all probes
  ## input is the current read table of Ctrl probes
  ## output is a list containing (1) updated read table of Ctrl probes; (2) removed probe; (3) list of CoVs
  ctrlProbeRatioCoVs <- calcCtrlCoVs(ctrlProbesCurrent) # calculates the CoV for each probe (by ratio to GeoMean of Ctrls) := Ratio_CoV
  ctrlCoV_iter <- cbind(ctrlProbesCurrent, Ratio_CoV = ctrlProbeRatioCoVs) # attaches this CoV column to Ctrl reads
  #print(ctrlCoV_iter) # display list to console w/ Ratio_CoVs appended to read table of Ctrls
  maxCoV_iter <- max(ctrlProbeRatioCoVs) # value of highest CoV (to be returned in a list from this function)
  #print(c("max CoV this cycle is", maxCoV_iter))
  if(maxCoV_iter > limitCoV) {
    probeToRemove <- names(which.max(ctrlProbeRatioCoVs)) # name of probe with highest CoV, to be removed
  }
  #print(c("probe to remove is:", probeToRemove))
  ctrlProbesNext <- ctrlProbesCurrent[row.names(ctrlProbesCurrent) != probeToRemove,] # removes probe w/ highest Ratio_CoV
  to_return <- list(ctrlProbesNext, probeToRemove, ctrlProbeRatioCoVs)
  return(to_return)
}

# replace probes with a value of 1 so the probe is not discareded
# replaced with the untreated value so no change is registered
subLowValues = function(vctr){
  len = length(vctr)
  vctr[2:len][vctr[2:len] < 2] = vctr[1]
  vctr = unlist(vctr)
  return(vctr)
}

#find infinite values in a dataframe
#likely produced becasue a control probe failed and trying to calculate log fold change produced infinite values
#from https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame
is.inf.df <- function(x) do.call(cbind, lapply(x, is.infinite))
is.na.df <- function(x) do.call(cbind, lapply(x, is.na))
is.zero.df <- function(x) do.call(cbind, lapply(x, function(y){y == 0}))
is.neg.df <- function(x, is_zero_or_less = F) if(is_zero_or_less == F){do.call(cbind, lapply(x, function(y){y < 0}))}else{do.call(cbind, lapply(x, function(y){y <= 0}))}
