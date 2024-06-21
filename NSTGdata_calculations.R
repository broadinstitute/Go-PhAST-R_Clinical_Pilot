setwd("C:/Users/eyoung/OneDrive - The Broad Institute/Documents/R")

#library import
if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

#establish directory
dir = ("G:/.shortcut-targets-by-id/1kcvyeOXWFbS3VaeKPHnjXmRIYRNSmc_q/Bhattacharyya_lab/Eleanor/MGH_BCx_GoPhASTR_pilot")
#probes to be excluded
remove_data = c("KpCip4_R11dn_KPHS_31980_0.95", "KpCip4_R12up_KPHS_40580_0.95", "EcGent4_Rup_b4550_0.92", "EcGent4_Rup_b3686e_0.92", 
                "EcGent4_Rup_b3687e_0.97")
#store all output
sink(paste0(dir, "/calculations_output.txt"))
#file import
known_data_in <- read.csv(paste0(dir, "/", "mel_data/AllMelpaperData_all.csv"), row.names = 1) %>%
  filter(!(tag %in% remove_data)) #read.csv(paste0(dir, "/", "mel_data/AllMelpaperData_2.csv"), row.names = 1)
known_data <- known_data_in %>% filter(use == "use")
dfnames <- list.files(path = paste0(dir, "/","CompositeCSVs"), pattern = "*.csv")
CRE_bckgrnd <- read.csv(paste0(dir, "/","CRE_data/2021.12.13_CRE_bkgd_subt_mean2SD.csv"), row.names = 1, header = FALSE)
#source of functions, colorscheme, and dataframe storage
source(paste0(dir, "/","HeatMapFuncs.R"))

# minimum value for Ctrl probes; any Ctrl probe with any value below this in any sample will be removed from all samples. Should be > neg ctrl probes
minCtrl = 10
# minimum value for Resp probes; any Resp probe with any value below this in any sample will be removed from all samples. Should be > neg ctrl probes
minResp = 0
# note: applying the above cutoffs to raw counts, before calculations
# only apply cutoff to control
untrOnly = F
# initialized value, will be altered in while loop to optimize ctrl probes
maxCoV = 1
# max covariance allowed
limitCoV = 0.25
#threshold for normalized response
minRespUntr = .1
minRespTrt = .1
#antibiotics used
abx <- c("Cip", "Levo", "Gent", "Cz", "Cro", "Atm", "Tzp", "Fep", "Etp", "Mem")
#beta lactams ordered by strength
bls <- c("Cz", "Cro", "Atm", "Tzp", "Fep", "Etp", "Mem")
#species used
species <- c("K. pneumoniae", "E. coli")
#drug classes
classes <- c("Bl", "Fq", "Gent")
#Beta Lactamase dataframe
rows_c <- c("wzi_1" ,"wzy_1" ,"wzi_2" ,"wzy_2" ,"CRE2_KPC_0.95" ,"CRE2_NDM_0.95" ,"CRE2_OXA48_0.95" ,"CRE2_CTXM15_0.95" ,"CRE2_IMP_1_0.95" ,
"CRE2_IMP_3_8_0.95" ,"CRE2_IMP_2_4_0.95" ,"CRE2_IMP_5_0.95" ,"CRE2_IMP_6_0.95" ,"CRE2_IMP_7_0.95" ,"CRE2_VIM_1_0.95" ,"CRE2_VIM_2_3_0.95" ,"OXA10_0.95")

##data frames
# beta lacatamase genes
data_c <- data.frame( matrix(ncol = 0 , nrow = length(rows_c) ) )
row.names(data_c) <- rows_c
# mean abs response for each sample
data_order <- data.frame( matrix(ncol = 3 , nrow = 0  ) )
colnames(data_order) <- c("sample", "r", "class")
# save control genes and note what are thrown out
controls <- data.frame( matrix(ncol = 3 , nrow = 0  ) )
colnames(controls) <- c("probe", "experiment", "pass_or_fail")

for (i in dfnames){
  #read df
  df = read.csv(paste0(dir, "/", "CompositeCSVs/", i), row.names = 1)
  #number of columns
  cols = length(colnames(df))
  #makes empty df for l2fc output
  rspnsPrbs <- row.names(df)[!grepl("NEG_|POS_|_C_|CRE2_|wz|OXA", row.names(df))]
  data_out <- data.frame( matrix(ncol = 0 , nrow = length(rspnsPrbs) ) )
  row.names(data_out) <- rspnsPrbs
  #make empty df for normalized tx (non beta lactamase) probes
  txPrbs <- row.names(df)[!grepl("NEG_|POS_|CRE2_|wz|OXA", row.names(df))]
  data_norm <- data.frame( matrix(ncol = 0 , nrow = length(txPrbs) ) )
  row.names(data_norm) <- txPrbs
  #make empty df for spike-in normalized data- pre control gene normalization
  data_corr <- data.frame( matrix(ncol = 0 , nrow = length(txPrbs) ) )
  row.names(data_corr) <- txPrbs
  # remove negative probes that are (on average) 3x higher than the mean negative probe response
  # remove negative probes that have a variance 3x higher than the average variance
  Neg_Probes = df[(grepl(rownames(df), pattern = "NEG")),] %>% as.matrix()
  Neg_Probes_Mean = Neg_Probes %>% mean(na.rm = TRUE)
  Neg_Probes_Row_Mean = Neg_Probes %>% rowMeans()
  rm_Neg_Probes_Mean = row.names(Neg_Probes)[Neg_Probes_Row_Mean > 3 * Neg_Probes_Mean]
  Neg_Probes_Row_Vars = Neg_Probes %>% matrixStats::rowVars()
  rm_Neg_Probes_Vars = row.names(Neg_Probes)[Neg_Probes_Row_Vars > 3 * mean(Neg_Probes_Row_Vars)]
  cat(paste0(i,
             ":\nNeg control means: ", paste0(round(Neg_Probes_Row_Mean, digits = 1), collapse = ", "), 
             "\nProbes removed by mean: ", rm_Neg_Probes_Mean, 
             "\nNeg control variances: ", paste0(round(Neg_Probes_Row_Vars, digits = 1), collapse = ", "),
             "\nProbes removed by variance: ", rm_Neg_Probes_Vars, "\n"))
  rm_Neg_Probes = unique(c(rm_Neg_Probes_Mean, rm_Neg_Probes_Vars))
  if(length(rm_Neg_Probes) > 0){
    df[rm_Neg_Probes,] = NA
  }
  Neg_Probes = df[(grepl(rownames(df), pattern = "NEG")),] %>% as.matrix()
  #print(i)
  if(grepl("Bl", i)){
    for (j in 1:(cols/8)){
      # partition by experiment
      idx1 = 1 + ((j - 1) * 8)
      idx2 = 8 + ((j - 1) * 8)
      temp = df[,idx1:idx2]
      # starts list to keep track of Ctrl Probes that will be removed
      removedCtrlProbes <- c()
      ##control for low probes
      allRespTrim <- rbind(temp[!rownames(temp) %in% (rspnsPrbs),], removeFailedProbes(temp[rspnsPrbs,], minVal = minResp))
      # pos ctrl correct
      allDataPosCorr <- posCtrlCorr(allRespTrim)
      # neg ctrl correct
      allDataCorr <- negCtrlCorr(allDataPosCorr)
      # split data
      allCtrlCorr <- allDataCorr[grepl("_C_", rownames(allDataCorr)), ]
      allRespCorr <- allDataCorr[grepl("_R", rownames(allDataCorr)), ]
      # pull out CRE probes for later - inelegant b/c "CRE2" wasn't included in a couple of the names, but this works:
      allCRECorr <- allDataCorr[grepl("CRE2|wz|OXA", rownames(allDataCorr)), ]
      # removed failed probes (can change minCtrl and minResp above as needed, with an eye towards neg ctrl values)
      allCtrlTrim <- removeFailedProbes(allCtrlCorr, minCtrl)
      #allRespTrim <- removeFailedProbes(allRespCorr, minResp)
      # copy control list
      ctrlProbes_toIter <- allCtrlTrim
      while(maxCoV > limitCoV) {
        # remove probe w/ max CoV
        ctrlProbeCoV_list <- removeCtrlMaxCoV(ctrlProbes_toIter)
        #print("returned list from loop is:")
        #print(ctrlProbeCoV_list)
        ctrlProbes_toIter <- as.data.frame(ctrlProbeCoV_list[1])
        #print("updated probes list is:")
        #print(ctrlProbes_toIter)
        if(as.character(ctrlProbeCoV_list[2]) != "none")removedCtrlProbes <- c(removedCtrlProbes, as.character(ctrlProbeCoV_list[2]))
        #print("updated list of removed probes:")
        #print(removedCtrlProbes)
        maxCoV <- max(unlist(ctrlProbeCoV_list[3]))
        #print(c("updated max CoV is:", maxCoV))
      }
      #print final list of kept probes:
      #print("final trimmed list of best Ctrl Probes is:")
      #print(rownames(ctrlProbes_toIter))
      #print(c("final max CoV is:", maxCoV))
      optCtrl <- ctrlProbes_toIter
      ctrlNorm <- colGeoMeans(optCtrl)
      allRespNorm <- sweep(allRespCorr, 2, ctrlNorm, FUN = "/")
      # outlier handling
      #if the untreated is below the threshold, it is rounded up
      allRespNorm[,1][allRespNorm[,1] < minRespUntr] = minRespUntr
      #if the treated sample is below the threshold, it is rounded up
      allRespNorm[,2:length(allRespNorm)][apply(allRespNorm[,2:length(allRespNorm)], 2, FUN = function(x) (x < minRespTrt))] = minRespTrt
      # if treated and untreated samples both had to be rounded up, the treated sample is given an NA value
      allRespNorm[,2:length(allRespNorm)][apply(allRespNorm[,2:length(allRespNorm)], 2, FUN = function(x) (x <= minRespTrt) & (allRespNorm[,1] <= minRespUntr))] = NA
      # for fun / for later plots: normalize Ctrl probes by optimized Ctrl probes
      allCtrlNorm <- sweep(allCtrlTrim, 2, ctrlNorm, FUN = "/")
      # save data normalized by nanostring spike ins
      data_corr = cbind(data_corr, allDataCorr[txPrbs,])
      ## but to really visualize these need to then include in fold-changes i.e. add back to allRespNorm
      data_norm = merge(data_norm[txPrbs,], rbind(allCtrlNorm, allRespNorm)[c(rownames(optCtrl), rownames(allRespCorr)),], by = 'row.names', all = TRUE)
      rownames(data_norm) = data_norm$Row.names
      data_norm = data_norm[-1]
      # and normalize CRE probes as well
      allCRENorm <- sweep(allCRECorr, 2, ctrlNorm, FUN = "/")
      # calculate fold induction by comparing to untreated
      allRespNormFold <- allRespNorm[,2:8]/allRespNorm[,1]
      rownames(allRespNormFold) <- rownames(allRespNorm)
      colnames(allRespNormFold) <- colnames(temp)[2:8]
      allRespNormLog <- log2(allRespNormFold)
      #replace infinite values with NaN
      allRespNormLog[is.inf.df(allRespNormLog)] <- NaN
      data_out <- merge(data_out, allRespNormLog, by = "row.names", all = TRUE)
      rownames(data_out) <- unlist(data_out["Row.names"])
      data_out <- subset(data_out, select = -c(1))
      #save normalized CRE data
      data_c <- merge(data_c, allCRENorm, by = "row.names", all = TRUE)
      rownames(data_c) <- unlist(data_c["Row.names"])
      data_c <- subset(data_c, select = -c(1))
      #save which control genes are used and which are not
      failedCtrlProbes = rownames(allCtrlCorr)[!(rownames(allCtrlCorr) %in% rownames(allCtrlTrim))]
      controls_temp <- data.frame(
        controls = c(rownames(ctrlProbes_toIter), removedCtrlProbes, failedCtrlProbes),
        experiment = c(rep( gsub(i, pattern = "[.]csv", replacement = ""), length(rownames(allCtrlCorr)) )),
        pass_or_fail = c( 
          rep("pass", length(rownames(ctrlProbes_toIter))), 
          rep("removed", length(removedCtrlProbes)),
          rep("failed", length(failedCtrlProbes)))
      )
      controls = rbind(controls, controls_temp)
    }
  }else if(grepl("Fq", i)){
    for (j in 1:(cols/3)){
      # partition by experiment
      idx1 = 1 + ((j - 1) * 3)
      idx2 = 3 + ((j - 1) * 3)
      temp = df[,idx1:idx2]
      # starts list to keep track of Ctrl Probes that will be removed
      removedCtrlProbes <- c()
      ##control for low probes
      allRespTrim <- rbind(temp[!rownames(temp) %in% (rspnsPrbs),], removeFailedProbes(temp[rspnsPrbs,], minVal = minResp))
      # pos ctrl correct
      allDataPosCorr <- posCtrlCorr(allRespTrim)
      # neg ctrl correct
      allDataCorr <- negCtrlCorr(allDataPosCorr)
      # split data
      allCtrlCorr <- allDataCorr[grepl("_C_", rownames(allDataCorr)), ]
      allRespCorr <- allDataCorr[grepl("_R", rownames(allDataCorr)), ]
      # removed failed probes (can change minCtrl and minResp above as needed, with an eye towards neg ctrl values)
      allCtrlTrim <- removeFailedProbes(allCtrlCorr, minCtrl)
      #allRespTrim <- removeFailedProbes(allRespCorr, minResp)
      # copy control list
      ctrlProbes_toIter <- allCtrlTrim
      while(maxCoV > limitCoV) {
        # remove probe w/ max CoV
        ctrlProbeCoV_list <- removeCtrlMaxCoV(ctrlProbes_toIter)
        #print("returned list from loop is:")
        #print(ctrlProbeCoV_list)
        ctrlProbes_toIter <- as.data.frame(ctrlProbeCoV_list[1])
        #print("updated probes list is:")
        #print(ctrlProbes_toIter)
        if(as.character(ctrlProbeCoV_list[2]) != "none")removedCtrlProbes <- c(removedCtrlProbes, as.character(ctrlProbeCoV_list[2]))
        #print("updated list of removed probes:")
        #print(removedCtrlProbes)
        maxCoV <- max(unlist(ctrlProbeCoV_list[3]))
        #print(c("updated max CoV is:", maxCoV))
      }
      # print final list of kept probes:
      #print("final trimmed list of best Ctrl Probes is:")
      #print(rownames(ctrlProbes_toIter))
      #print(c("final max CoV is:", maxCoV))
      optCtrl <- ctrlProbes_toIter
      ctrlNorm <- colGeoMeans(optCtrl)
      allRespNorm <- sweep(allRespCorr, 2, ctrlNorm, FUN = "/")
      #outlier handling
      #if the untreated is below the threshold, it is rounded up
      allRespNorm[,1][allRespNorm[,1] < minRespUntr] = minRespUntr
      #if the treated sample is below the threshold, it is rounded up
      allRespNorm[,2:length(allRespNorm)][apply(allRespNorm[,2:length(allRespNorm)], 2, FUN = function(x) (x < minRespTrt))] = minRespTrt
      # if treated and untreated samples both had to be rounded up, the treated sample is given an NA value
      allRespNorm[,2:length(allRespNorm)][apply(allRespNorm[,2:length(allRespNorm)], 2, FUN = function(x) (x <= minRespTrt) & (allRespNorm[,1] <= minRespUntr))] = NA
      # for fun / for later plots: normalize Ctrl probes by optimized Ctrl probes
      allCtrlNorm <- sweep(allCtrlTrim, 2, ctrlNorm, FUN = "/")
      # save data normalized by nanostring spike ins
      data_corr = cbind(data_corr, allDataCorr[txPrbs,])
      ## but to really visualize these need to then include in fold-changes i.e. add back to allRespNorm
      data_norm = merge(data_norm[txPrbs,], rbind(allCtrlNorm, allRespNorm)[c(rownames(optCtrl), rownames(allRespCorr)),], by = 'row.names', all = TRUE)
      rownames(data_norm) = data_norm$Row.names
      data_norm = data_norm[-1]
      # calculate fold induction by comparing to untreated
      allRespNormFold <- allRespNorm[,2:3]/allRespNorm[,1]
      rownames(allRespNormFold) <- rownames(allRespNorm)
      colnames(allRespNormFold) <- colnames(temp)[2:3]
      allRespNormLog <- log2(allRespNormFold)
      #replace infinite values with NaN
      allRespNormLog[is.inf.df(allRespNormLog)] <- NaN
      data_out <- merge(data_out, allRespNormLog, by = "row.names", all = TRUE)
      rownames(data_out) <- unlist(data_out["Row.names"])
      data_out <- subset(data_out, select = -c(1))
      failedCtrlProbes = rownames(allCtrlCorr)[!(rownames(allCtrlCorr) %in% rownames(allCtrlTrim))]
      controls_temp <- data.frame(
        controls = c(rownames(ctrlProbes_toIter), removedCtrlProbes, failedCtrlProbes),
        experiment = c(rep( gsub(i, pattern = "[.]csv", replacement = ""), length(rownames(allCtrlCorr)) )),
        pass_or_fail = c( 
          rep("pass", length(rownames(ctrlProbes_toIter))), 
          rep("removed", length(removedCtrlProbes)),
          rep("failed", length(failedCtrlProbes)))
      )
      controls = rbind(controls, controls_temp)
    }
  }else{
    for (j in 1:(cols/2)){
      idx1 = 1 + ((j - 1) * 2)
      idx2 = 2 + ((j - 1) * 2)
      temp = df[,idx1:idx2]
      # starts list to keep track of Ctrl Probes that will be removed
      removedCtrlProbes <- c()
      ##control for low probes
      allRespTrim <- rbind(temp[!rownames(temp) %in% (rspnsPrbs),], removeFailedProbes(temp[rspnsPrbs,], minVal = minResp))
      # pos ctrl correct
      allDataPosCorr <- posCtrlCorr(allRespTrim)
      # neg ctrl correct
      allDataCorr <- negCtrlCorr(allDataPosCorr)
      # split data
      allCtrlCorr <- allDataCorr[grepl("_C_", rownames(allDataCorr)), ]
      allRespCorr <- allDataCorr[grepl("_R", rownames(allDataCorr)), ]
      allCtrlTrim <- removeFailedProbes(allCtrlCorr, minCtrl)
      #allRespTrim <- removeFailedProbes(allRespCorr, minResp)
      #copy control list
      ctrlProbes_toIter <- allCtrlTrim
      while(maxCoV > limitCoV) {
        # remove probe w/ max CoV
        #print("running while loop") # for troubleshooting
        ctrlProbeCoV_list <- removeCtrlMaxCoV(ctrlProbes_toIter)
        #print("returned list from loop is:")
        #print(ctrlProbeCoV_list)
        ctrlProbes_toIter <- as.data.frame(ctrlProbeCoV_list[1])
        #print("updated probes list is:")
        #print(ctrlProbes_toIter)
        if(as.character(ctrlProbeCoV_list[2]) != "none")removedCtrlProbes <- c(removedCtrlProbes, as.character(ctrlProbeCoV_list[2]))
        #print("updated list of removed probes:")
        #print(removedCtrlProbes)
        maxCoV <- max(unlist(ctrlProbeCoV_list[3]))
        #print(c("updated max CoV is:", maxCoV))
      }
      # print final list of kept probes:
      #print("final trimmed list of best Ctrl Probes is:")
      #print(rownames(ctrlProbes_toIter))
      #print(c("final max CoV is:", maxCoV))
      #allDataCorr[c(rownames(optCtrl), rownames(allRespCorr))]
      optCtrl <- ctrlProbes_toIter
      ctrlNorm <- colGeoMeans(optCtrl)
      allRespNorm <- sweep(allRespCorr, 2, ctrlNorm, FUN = "/")
      #outlier handling
      #if the untreated is below the threshold, it is rounded up
      allRespNorm[,1][allRespNorm[,1] < minRespUntr] = minRespUntr
      #if the treated is below the threshold, it is rounded up
      allRespNorm[,2][allRespNorm[,2] < minRespUntr] = minRespUntr
      # if treated and untreated samples both had to be rounded up, both are given an NA value
      below_background = (allRespNorm[,2] <= minRespTrt) & (allRespNorm[,1] <= minRespUntr)
      allRespNorm[,1][below_background] = NA
      allRespNorm[,2][below_background] = NA
      # for fun / for later plots: normalize Ctrl probes by optimized Ctrl probes
      allCtrlNorm <- sweep(allCtrlTrim, 2, ctrlNorm, FUN = "/")
      ## but to really visualize these need to then include in fold-changes i.e. add back to allRespNorm
      # save data normalized by nanostring spike ins
      data_corr = cbind(data_corr, allDataCorr[txPrbs,])
      # save fully normalized data including Ctrl probes normalized by their own geometric mean
      data_norm = merge(data_norm[txPrbs,], rbind(allCtrlNorm, allRespNorm)[c(rownames(optCtrl), rownames(allRespCorr)),], by = 'row.names', all = TRUE)
      rownames(data_norm) = data_norm$Row.names
      data_norm = data_norm[-1]
      # calculate fold induction by comparing to untreated
      allRespNormFold <- data.frame( allRespNorm[,2]/allRespNorm[,1] )
      rownames(allRespNormFold) <- rownames(allRespNorm)
      colnames(allRespNormFold) <- colnames(temp)[2]
      allRespNormLog <- log2(allRespNormFold)
      #replace infinite values with NaN
      allRespNormLog[is.inf.df(allRespNormLog)] <- NaN
      data_out <- merge(data_out, allRespNormLog, by = "row.names", all = TRUE)
      rownames(data_out) <- unlist(data_out["Row.names"])
      data_out <- subset(data_out, select = -c(1))
      failedCtrlProbes = rownames(allCtrlCorr)[!(rownames(allCtrlCorr) %in% rownames(allCtrlTrim))]
      controls_temp <- data.frame(
        controls = c(rownames(ctrlProbes_toIter), removedCtrlProbes, failedCtrlProbes),
        experiment = c(rep( gsub(i, pattern = "[.]csv", replacement = ""), length(rownames(allCtrlCorr)) )),
        pass_or_fail = c( 
          rep("pass", length(rownames(ctrlProbes_toIter))), 
          rep("removed", length(removedCtrlProbes)),
          rep("failed", length(failedCtrlProbes)))
      )
      controls = rbind(controls, controls_temp)
    }
  }
  
  #remove excluded probes
  data_norm = data_norm[!(row.names(data_norm) %in% remove_data),]
  data_out = data_out[!(row.names(data_out) %in% remove_data),]
  #data order
  data_order = rbind(
    data_order, 
    data.frame( 
      sample = gsub(colnames(data_out), pattern = "_.*" , replacement = "" ),
      r = (colMeans(abs(data_out), na.rm = TRUE)),
      class = rep(gsub( "[.]csv" , "" , gsub( "Ec|Kp" , "" , i )), ncol(data_out)) )
    )
  #append dataframes to list
  # write spike-in corrected data
  write.csv(x = data_corr, 
            file = paste0(dir, "/Background_Threshold_2/CompositeCSVs_corrected/", 
                          gsub(i, pattern = ".csv", replace = ""), "_corrected.csv"))
  # write normalized data and l2fc data
  write.csv(x = data_norm, 
            file = paste0(dir, "/Background_Threshold_2/CompositeCSVs_normalized/", 
                          gsub(i, pattern = ".csv", replace = ""), "_normalized.csv"))
  
  write.csv(x = data_out, 
            file = paste0(dir, "/Background_Threshold_2/CompositeCSVs_l2fc/", 
                          gsub(i, pattern = ".csv", replace = ""), "_l2fc.csv"))
  output_list[[ unlist(strsplit(i, ".csv")) ]] <- data_out
}


##### Beta Lactam Calculations


#in: CRE reads normalized to control probes
#there is only one VIM probe in sumVIM
data_c2 <- data_c %>% 
  filter(!grepl("wzi|wzy|IMP_3_8|IMP_5|VIM_1", row.names(data_c))) 

data_c2 <- rbind(data_c2, colSums(data_c2[grep("IMP", rownames(data_c2)),], na.rm = TRUE))
rownames(data_c2)[grepl("CRE2_VIM_2_3_0.95", rownames(data_c2))] = "KpCRE2_sumVIM_0.95"
rownames(data_c2)[dim(data_c2)[1]] = "KpCRE2_sumIMP_0.95"
data_c2 = data_c2[order(rownames(data_c2)),] %>% 
  filter(!grepl("_IMP_", row.names(data_c2)))
CRE_bckgrnd = CRE_bckgrnd[order(rownames(CRE_bckgrnd)),]
#CRE_bckgrnd = as.data.frame(CRE_bckgrnd[order(rownames(CRE_bckgrnd)), ], row.names = sort(rownames(CRE_bckgrnd)))
data_c3 <- sweep(data_c2, MARGIN = 1,STATS = CRE_bckgrnd, FUN = "-")

data_c4 <- data_c3 %>%
  ##cchange to data_c3 below?
  mutate( tag = row.names(data_c2) ) %>%
  pivot_longer(cols = contains("MGH"), names_to = "sample", values_to = "value") %>%
  separate(sample, c("sample", "drug"), sep = "_" ) %>%
  mutate( #sample = factor(sample, levels = numeric_ordered_samples) ,
          value = ifelse(is.na(value), min(value, na.rm = TRUE), value )) %>%
  group_by(tag, sample) %>%
  dplyr::summarise(value = (mean(value, na.rm = TRUE))) %>%
  mutate(present = ifelse(value > 2, 1, 0)) 
# (value = log(mean(value, na.rm = TRUE)))

maxC4 <- ceiling(max(data_c4$value))

cre_present = data_c4 %>% filter(present == 1) %>% pull(tag) %>% unique()

data_c5 = data_c4 %>% filter(tag %in% cre_present) %>% separate(tag, sep = "_", into = c(NA, "BLase", NA))


##### save data


#save each sample's avg abs l2fc:
write.csv(x = data_order, file = paste0(dir, "/Background_Threshold_2/NSTGsample_AvgAbsL2FC.csv"))
#save each all cre counts normalized by control data (only pulled from beta lactam .csvs like EcBl):
write.csv(x = data_c, file = paste0(dir, "/Background_Threshold_2/NSTGsample_CREnorm.csv"))
#save edited cre count data
write.csv(x = data_c4, file = paste0(dir, "/Background_Threshold_2/NSTGsample_CREnorm_processed.csv"))
#save meaning cre data
write.csv(x = data_c5, file = paste0(dir, "/Background_Threshold_2/NSTGsample_BLase.csv"))
#control probe stats:
controls %>% count(pass_or_fail) %>% arrange(n)

# data_order2 = data_order %>% 
#   filter( !(class %in% c("Fq", "Bl")) ) %>%
#   full_join(sample_renaming_df[c("strains", "titles")], by = c("sample" = "strains")) %>%
#   mutate(title = ifelse(titles == "Ec", "E.coli", "K. pneumoniae"),
#          sample = factor(sample, levels = ordered_samples),
#          r = as.numeric(r),
#          drug = factor(class, levels = c("Amp", bls, "Bl", "Cip", "Levo", "Fq", "Gent")),
#          class = ifelse(drug == "Gent", "Gent",
#                         ifelse(drug %in% c("Levo", "Cip", "Fq"), "Fq", "Bl")))


##### known data data manipulation 


known_data = known_data %>%
  mutate(Strain = ifelse(grepl(time, pattern = "rb", ignore.case = T),
                         time,
                         Strain),
         species = substr(exp, 1, 2),
         drug = str_to_title(drug),
         class = ifelse(drug == "Gent",
                        "Gent",
                        ifelse(drug %in% c("Cip", "Levo"),
                               "Fq",
                               "Beta Lactam"))) %>%
  filter(!(tag %in% remove_data))


##### Calculate SPDs


known_data_spd_data = known_data %>%
  mutate(drug2 = recode(drug,
                        "Cefaz" = "Cz", "Ctx" = "Cro", "Aztr" = "Atm", "Erta" = "Etp", "Mero" = "Mem", "Pip" = "Tzp")) %>%
  #select(!c("drop", "time", "exp"))
  select(!c("time", "exp"))

# calculate the mean of each probe in susceptible strains
known_data_spd_derivS = known_data_spd_data %>%
  #flipping SPD: swap S and R
  filter(S_I_R == "S") %>%
  #filter(S_I_R == "R") %>%
  group_by(species, drug, tag) %>%
  #flipping SPD: R_tag not S_tag
  summarise(s_tag_mean = mean(value))
#summarise(r_tag_mean = mean(value))

# difference between each probe and the probe mean in S (probe displacement)
known_data_spd_displacement = known_data_spd_derivS %>%
  full_join(known_data_spd_data, by = c("species", "drug", "tag")) %>%
  #flipping SPD: R_tag not S_tag
  mutate(tag_disp = value - s_tag_mean)
#mutate(tag_disp = value - r_tag_mean)

#claculate sumSqAxVector
known_data_spd_vector = known_data_spd_displacement %>%
  group_by(species, drug, S_I_R, tag) %>%
  # calculate average probe displacement
  summarise(tag_disp_mean = mean(tag_disp)) %>% 
  #full_join(known_data_spd, by = c("species", "drug", "S_I_R", "tag")) %>% 
  pivot_wider(names_from = S_I_R, values_from = tag_disp_mean, names_prefix = "tag_dis_mean_") %>%
  # calculate axial vector and 
  #flipping SPD: reverse order of R and S to change baseline to R
  mutate(axVector = tag_dis_mean_R - tag_dis_mean_S,
         #mutate(axVector = tag_dis_mean_S - tag_dis_mean_R,
         sqAxVector = axVector^2) %>%  
  ungroup() %>%
  ## calculate sum of the squared axial vector ##
  group_by(species, drug) %>%
  mutate(sumSqAxVector = sum(sqAxVector)) %>%
  ungroup()

#calculate SPD
known_data_spd = known_data_spd_vector %>%
  # full_join(known_data_spd_vector, 
  #         known_data_spd_sum, 
  #         by = c("drug", "species")) %>%
  ## join with displacement data ##
  full_join(known_data_spd_displacement,
            by = c("drug", "species", "tag")) %>%
  mutate(strain_disp = tag_disp * axVector) %>% 
  ## calculte sum of strain displacement ##
  group_by(species, drug, S_I_R, Strain) %>%
  mutate(total_strain_disp = sum(strain_disp)) %>%
  #, abs_total_strain_disp = abs(total_strain_disp)) %>%
  ungroup() %>%
  ## calculate SPD
  mutate(spd = total_strain_disp * abs(total_strain_disp) /(sumSqAxVector^2)) %>%
  ## reorder to follow order of operation ##
  #flipping SPD: changed name of _tag_mean
  select(c('species', 'drug', 'drug2', 'tag', 'S_I_R', 'Strain', 'value', 's_tag_mean', #'r_tag_mean',#'s_tag_mean', 
           'tag_disp', 'tag_dis_mean_R', 'tag_dis_mean_S' , 'axVector', 'sqAxVector', 
           'sumSqAxVector', 'strain_disp', 'total_strain_disp', 'spd'))

single_drug_data = grep(names(output_list), value = T, pattern = "Fq|Bl", invert = T)

#known_data_spd %>% group_by(species, drug, S_I_R, Strain, spd) %>% summarise(n_probes = n()) %>% view()

known_data_spd %>% write.csv(paste0(dir, "/Background_Threshold_2/SPD_calcs_for_mel_data_ReferenceStrains.csv"))

#flipping SPD: changed name of _tag_mean
known_data_spd_subset = known_data_spd %>% select(species, drug, drug2, tag, s_tag_mean, axVector, sqAxVector, sumSqAxVector) %>% distinct()
#known_data_spd_subset = known_data_spd %>% select(species, drug, drug2, tag, r_tag_mean, axVector, sqAxVector, sumSqAxVector) %>% distinct()


sample_data_spd = as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(sample_data_spd) = c("species", "drug", "class", "sample", "spd", "rm_probes")

## calculation changes every time a probe is removed

for (i in single_drug_data){
  temp_data = output_list[[i]]
  s = substr(i, 1, 2)
  d = substr(i, 3, nchar(i))
  c = ifelse(d == "Fep", "Cro", d) # ifelse(d == "Levo", "Cip", d))
  temp_subset = known_data_spd_subset %>%
    filter(species == s & drug2 == c & tag %in% row.names(temp_data))
  sumSqAxVector = temp_subset$sumSqAxVector %>% unique()
  temp_data = temp_data[temp_subset$tag, ]
  n_samples = ncol(temp_data)
  #flipping SPD: changed name of _tag_mean
  temp_disp = sweep(temp_data, 1, temp_subset$s_tag_mean)
  #temp_disp = sweep(temp_data, 1, temp_subset$r_tag_mean)
  temp_strain_disp <- sweep(temp_disp, 1, temp_subset$axVector, FUN = "*")
  temp_sum_strain_disp <- colSums(temp_strain_disp, na.rm = T) # just ignore NA values
  ### check this if things dont work # i was right
  
  #separate into sample data with and without NAs
  #for no na values use precalculated sumSqAxVector
  temp_spd = 
    temp_sum_strain_disp[colSums(is.na(temp_data)) == 0] * 
    abs(temp_sum_strain_disp[colSums(is.na(temp_data)) == 0]) / 
    (sumSqAxVector^2)
  
  temp_df = data.frame(
    samples = names(temp_spd),
    spd = temp_spd,
    rm_probes = "")
  
  for (j in names(temp_sum_strain_disp[colSums(is.na(temp_data)) != 0])){
    sum_strain_disp = temp_sum_strain_disp[j]
    rm_probes = rownames(temp_data)[is.na(temp_data[,j])]
    sumSqAxVector_temp = temp_subset %>% filter(!(tag %in% rm_probes)) %>% pull(sqAxVector) %>% sum()
    spd = sum_strain_disp * abs(sum_strain_disp) / (sumSqAxVector_temp^2)
    temp_df[j,] = c(j, spd, paste0(rm_probes, collapse = ", "))
    
  }
  
  temp_df = temp_df %>% mutate(
    species = s,
    drug = d,
    class = c
  )
  
  sample_data_spd = rbind(sample_data_spd, temp_df)
}

##### calculate SPD for all of Mel's data

all_mel_data_SPD = known_data_spd %>% 
  filter(!(tag %in% remove_data)) %>%
  select(c("species", "drug", "tag", "sumSqAxVector", "s_tag_mean", "axVector")) %>%
  distinct() %>% 
  mutate(drug = str_to_lower(drug)) %>% #view()
  right_join(
    mutate(known_data_in, 
           species = (substr(exp, start = 1, stop = 2)),
           drug = str_to_lower(drug)), 
    by = c("species", "drug", "tag")) %>% 
  #exclude the drug ampicillin
  filter(drug != "amp") %>% 
  mutate(tag_disp = value - s_tag_mean,
         strain_disp = tag_disp * axVector) %>%
  group_by(species, drug, Strain, S_I_R, MIC, exp, use) %>%
  mutate(sum_strain_disp = sum(strain_disp)) %>%
  ungroup() %>%
  mutate(spd = sum_strain_disp * abs(sum_strain_disp) / (sumSqAxVector^2))

all_mel_data_SPD %>% write.csv(paste0(dir, "/Background_Threshold_2/all_mel_data_SPD.csv"))

##### include CRE data

sample_data_spd = sample_data_spd %>%
  mutate(spd = as.numeric(spd),
         drug = factor((drug), levels = abx),
         sample = gsub(samples, pattern = "_.*", replacement = "")) %>%
  full_join(subset(data_c5, present == 1), by = "sample") #%>% view()

export_spd = sample_data_spd %>%
  mutate(sample = gsub(samples, pattern = "[A-z,_]", replacement = "")) %>%
  select(samples, species, sample, drug, spd, rm_probes) %>%
  arrange(drug) %>%
  write.csv(file = paste0(dir, "/Background_Threshold_2/", "MGH_AST_spd.csv"), row.names = F)

sink()
