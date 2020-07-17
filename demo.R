load("counts_clean_subsample.RData")
load("phenotypic_data.RData")

library(tidyverse)
library(e1071)
library(DESeq2)
library(limma)
library(pracma)

#### Define SVM function
svm.train = function(readcounts, traindata, testdata = NA, referencelevel = "queen", kerneltype = "radial", crossfold = 5, vstCheck = T){
  
  svm.counts.test=NA
  
  # normalise data
  svm.counts = readcounts
  # perform DESeq's variance stabilizing tranformation, which is preferable to logging for gene expression data
  if(vstCheck){svm.counts.vst = vst(svm.counts)}else{svm.counts.vst = varianceStabilizingTransformation(svm.counts)}
  # normalize counts between samples 
  svm.counts.vst.quantiles = normalizeBetweenArrays(svm.counts.vst, method = "quantile")
  # scale counts and remove zero-variance features
  svm.counts.vst.quantiles.scale = t(scale(t(svm.counts.vst.quantiles)))
  svm.counts.vst.quantiles.scale = na.omit(svm.counts.vst.quantiles.scale)
  
  # Divide transcriptomic data into training set (queens and workers from control) and test set (individuals from treatment) 
  svm.counts.train = svm.counts.vst.quantiles.scale[,which(colnames(svm.counts.vst.quantiles.scale) %in% traindata$ID)]
  if(length(testdata)>1){
    svm.counts.test = svm.counts.vst.quantiles.scale[,which((colnames(svm.counts.vst.quantiles.scale) %in% testdata$ID))]
  }
  
  # Perform a grid search to optimise SVM parameters
  svm.counts.tuneResult = tune("svm", 
                               train.x = t(svm.counts.train),
                               train.y = as.numeric(traindata$Role == referencelevel),
                               probability = TRUE, 
                               scale = FALSE,
                               kernel = kerneltype,
                               tunecontrol = tune.control(sampling = "cross",
                                                          cross = crossfold),
                               ranges = list(gamma = 10^(-7:-5),
                                             cost = 2^(3:5))
  )
  
  # Final classifier
  svm.counts.classifier = svm.counts.tuneResult$best.model
  
  svm.counts.prediction = NULL
  if(length(testdata)>1){
    # Make predictions for the test data, if test data were provided.
    svm.counts.prediction = predict(svm.counts.classifier,
                                    t(svm.counts.test),
                                    type = "class", 
                                    probability = TRUE)
  }
  
  #output prediction for test data and cross-validation error for training data
  svm.result = list("prediction" = svm.counts.prediction,
                    "validation_error" = signif(svm.counts.tuneResult$best.performance,4),
                    "traincounts" = svm.counts.train,
                    "testcounts" = svm.counts.test)
  
  #return results
  return(svm.result)
}

#### Perform initial classification
# Divide data into training (control) set and test (queen removal) set
svm.data = phenotypic_data[,c("ID","Role")]
svm.data.train = subset(svm.data, Role %in% c("queen","worker_ctrl"))
svm.data.test = subset(svm.data, !(Role %in% c("queen","worker_ctrl")))

# apply svm to entire set of genes
svm.full = svm.train(counts_clean_subsample,
                     svm.data.train, 
                     svm.data.test, 
                     crossfold = 3,
                     vstCheck = F)

print(paste0("Root mean cross-validation error rate for full model: ",svm.full$validation_error))

#### Perform feature selection
# create copy of training data that we can subject to repeated trimming while preserving original frame
svm.counts.train.iterate = svm.full$traincounts
#record original number of features
nfeatures = nrow(svm.counts.train.iterate)
#target number of features 
nfeatures_target = 100
traindata = svm.data.train
#instantiate data frame to hold data on the error of each model
iterations = data.frame(feature = character(),
                        error_before_removal = numeric())
#iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  
  error = c()
  
  #run repeatedly to account for stochasticity in cross-validation
  for(i in 1:20){
    
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult = tune("svm", 
                                 train.x = t(svm.counts.train.iterate), 
                                 train.y =  as.numeric(traindata$Role == "queen"),
                                 probability = TRUE, 
                                 scale = FALSE,
                                 kernel = "radial", 
                                 tunecontrol = tune.control(sampling = "cross", 
                                                            cross = 3),
                                 ranges = list(gamma = 10^(-5:-7), cost = 2^(4:6)))
    #record error
    error = c(error, svm.counts.tuneResult$best.performance)
  }
  #sample classifier
  svm.counts.classifier = svm.counts.tuneResult$best.model
  #return mean error value
  error = signif(mean(error),4)
  #extract feature weights
  weights = (t(svm.counts.classifier$coefs) %*% svm.counts.classifier$SV)
  #calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature = colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  #remove lowest-weight feature from data frame
  svm.counts.train.iterate = subset(svm.counts.train.iterate, !(rownames(svm.counts.train.iterate) %in% c(weakfeature)))
  #in a dataframe, store removed feature name and error value before removing that feature
  iterations = rbind(iterations, tibble(feature = weakfeature,
                                        error_before_removal = error))
  #tick down
  nfeatures = (nfeatures-1)
  #output every 100 runs to track progress
  if((nfeatures/100)%%1==0){print(paste0("Features remaining: ",nfeatures))}
}
iterLength = 1:nrow(iterations)
# take moving average to smooth out variation
moving_avg = movavg(iterations$error_before_removal, 100, "s") 

# plot data to ensure we have the expected 'hockeystick' shape 
# note that this will be truncated for the subsetted dataset provided for demoing!
hockeyData = data.frame(num = iterLength, error = moving_avg)
hockeyData_plot = hockeyData
hockeyData_plot$num = abs(iterLength - (max(iterLength)+1))
plot(hockeyData_plot$num,hockeyData_plot$error, xlim = rev(c(0, 1000)))
# get minimum of this curve to find the point at which the error window is at its minimum
optimal_removal = which(moving_avg == min(moving_avg));
# list the features to be removed from the original set of genes
features_to_remove = iterations$feature[1:optimal_removal]
# new dataframe with less-useful features removed
counts_clean_subsample = subset(counts_clean_subsample,
                                           !(rownames(counts_clean_subsample) %in% features_to_remove))
# re-perform support vector classification using the new, optimally caste-separating set of features
svm.optimal = svm.train(counts_clean_subsample, 
                        referencelevel = "queen",
                        svm.data.train, 
                        svm.data.test,
                        crossfold = 3,
                        vstCheck = F)

print(paste0("Number of genes included in optimised model: ", nrow(counts_clean_subsample)))
print(paste0("Root mean cross-validation error rate for optimised model: ", svm.optimal$validation_error))
     