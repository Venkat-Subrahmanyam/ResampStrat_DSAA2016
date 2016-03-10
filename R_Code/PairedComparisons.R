library(performanceEstimation)

#HYPOTHESIS 1

pc <- pairedComparisons(exp,baseline="mc.lm",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.svm",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.mars",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.rf",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test
################################

#HYPOTHESIS 2

pc <- pairedComparisons(exp,baseline="mc.lm_UNDERB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.svm_UNDERB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.mars_UNDERB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.rf_UNDERB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.lm_SMOTEB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.svm_SMOTEB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.mars_SMOTEB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test

pc <- pairedComparisons(exp,baseline="mc.rf_SMOTEB",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test
################################

#HYPOTHESIS 3

pc <- pairedComparisons(exp,baseline="mc.arima",maxs=rep(TRUE,3), p.value=0.05) 
pc$F1$WilcoxonSignedRank.test
################################
