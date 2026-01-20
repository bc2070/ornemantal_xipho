library(ggplot2)
arrSimDirs <- Sys.glob("sims/drift*")
mMinTractLen <- 0;
mMaxTractLen <- 1e9;

mMinTractLen <- 0;
mMaxTractLen <- 50000;

mMinTractLen <- 50000;
mMaxTractLen <- 100000;

mMinTractLen <- 100000;
mMaxTractLen <- 1000000;

mMinTractLen <- 1000000;
mMaxTractLen <- 10000000;


datAcc <- NULL;
for(sDIR in arrSimDirs) {
  nDrift <- as.numeric(gsub("drift", "", basename(sDIR)));
  datThisAcc <- read.table(paste0(sDIR,"/results_summary_simulated_hybrids_reads_gen50_prop_par1_0.1"), header=F, sep="\t")
  datThisAcc$V7[datThisAcc$V7=="par2par1"] <- "par1par2"
  datThisAcc$Accuracy <- datThisAcc$V8 / (datThisAcc$V8 + datThisAcc$V9)
  datThisAcc$tractsize <- datThisAcc$V3 - datThisAcc$V2 + 1
  datThisAcc$AIMsInTract <- datThisAcc$V8 + datThisAcc$V9
  
  datThisAcc <- datThisAcc[datThisAcc$tractsize>=mMinTractLen & datThisAcc$tractsize <= mMaxTractLen, ];
  #plot(log10(datThisAcc$tractsize), datThisAcc$Accuracy)
  #plot(datThisAcc$AIMsInTract, datThisAcc$Accuracy, xlim=c(0,1000))
  
  datAgg <- aggregate(datThisAcc[, 8:9], by=list(anctype=datThisAcc$V7), FUN=sum)
  datAgg$Accuracy <- datAgg$V8 / (datAgg$V8 + datAgg$V9)
  
  #get conf. int. by resampling
  arrAncTypes <- unique(datAgg$anctype);
  for(sAncType in arrAncTypes) {
    datAccAncType <- datThisAcc[datThisAcc$V7 == sAncType, 8:9];
    arrResampledAcc <- c();
    for(nRep in 1:1000) {
      datAccAncTypeResampled <- datAccAncType[sample(1:nrow(datAccAncType), size = nrow(datAccAncType), replace = T ),];
      nCorrect <- sum(datAccAncTypeResampled$V8);
      nErr <- sum(datAccAncTypeResampled$V9);
      arrResampledAcc <- c(arrResampledAcc, nCorrect / (nCorrect + nErr) )
    }
    
    arrCI <- quantile(arrResampledAcc, c(0.025, 0.975), na.rm=T)
    datAgg[datAgg$anctype==sAncType, "CI_Low"] <- arrCI[1];
    datAgg[datAgg$anctype==sAncType, "CI_High"] <- arrCI[2];
  }
  datAgg$Drift <- nDrift
  datAcc <- rbind(datAcc, datAgg);
}


ggplot(datAcc, aes(x = Drift, y = Accuracy, group = anctype, color=anctype)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = CI_Low, 
                    ymax = CI_High), 
                width = 0.02, linewidth = 0.4) +
  labs(x = "Drift from reference parents", y = "Accuracy")

ggsave(paste0("accu_tractmin_",mMinTractLen,"_tractmax_",mMaxTractLen,".pdf"), width = 4, height=2, units="in")
