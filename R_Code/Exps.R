library(lubridate)
library(xts)
library(performanceEstimation)
library(uba)
library(UBL) # only for distance functions
library(DMwR)

# data set generated with target variable always on the last column.
# This is mandatory for smote versions to work correctly.
# Function to create the embed
create.data<-function(ts,embed){
  t<-index(ts)[-(1:(embed-1))]
  e<-embed(ts,embed)[,embed:1]
  colnames(e)<-paste('V',1:embed,sep='')
  d <-xts(e,t)
  as.data.frame(d)
}

#EVALUATION STATISTICS (Utility-Based Regression Framework)
eval.stats <- function(form,train,test,preds,ph,ls) {
  trues <- resp(form,test)
  trainY <- resp(form,train)
  
  prec <- util(preds,trues,ph,ls,util.control(umetric="P",event.thr=0.9))
  rec  <- util(preds,trues,ph,ls,util.control(umetric="R",event.thr=0.9))
  F05  <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=0.5,event.thr=0.9))
  F1   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=1,event.thr=0.9))
  F2   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=2,event.thr=0.9))
  
  mad=mean(abs(trues-preds))
  mse=mean((trues-preds)^2)
  mape= mean((abs(trues-preds)/trues))*100
  rmse= sqrt(mean((trues-preds)^2))
  mae_phi= mean(phi(trues,phi.parms=ph)*(abs(trues-preds)))
  mape_phi= mean(phi(trues,phi.parms=ph)*(abs(trues-preds)/trues))*100
  mse_phi= mean(phi(trues,phi.parms=ph)*(trues-preds)^2)
  rmse_phi= sqrt(mean(phi(trues,phi.parms=ph)*(trues-preds)^2))
  prec=prec
  rec=rec
  F05=F05
  F1=F1
  F2=F2
  
  c(
    prec=prec,rec=rec,F1=F1
  )
  
}

#EVALUATION FOR ARIMA (Utility-Based Regression Framework)
eval.stats_ARIMA <- function(train,test,preds,ph,ls) {
  trues <- test
  
  prec <- util(preds,trues,ph,ls,util.control(umetric="P",event.thr=0.9))
  rec  <- util(preds,trues,ph,ls,util.control(umetric="R",event.thr=0.9))
  F05  <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=0.5,event.thr=0.9))
  F1   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=1,event.thr=0.9))
  F2   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=2,event.thr=0.9))
  
  mad=mean(abs(trues-preds))
  mse=mean((trues-preds)^2)
  mape= mean((abs(trues-preds)/trues))*100
  rmse= sqrt(mean((trues-preds)^2))
  mae_phi= mean(phi(trues,phi.parms=ph)*(abs(trues-preds)))
  mape_phi= mean(phi(trues,phi.parms=ph)*(abs(trues-preds)/trues))*100
  mse_phi= mean(phi(trues,phi.parms=ph)*(trues-preds)^2)
  rmse_phi= sqrt(mean(phi(trues,phi.parms=ph)*(trues-preds)^2))
  prec=prec
  rec=rec
  F05=F05
  F1=F1
  F2=F2
  
  c(
    prec=prec,rec=rec,F1=F1
  )
  
}

#WORKFLOWS
mc.lm <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.lm_UNDERB <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)

  res
}

mc.lm_UNDERT <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.lm_UNDERTPhi <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.lm_SMOTEB <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.lm_SMOTET <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.lm_SMOTETPhi <- function(form,train,test,...) {
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- lm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_UNDERB <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_UNDERT <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_UNDERTPhi <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_SMOTEB <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_SMOTET <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm_SMOTETPhi <- function(form,train,test,...) {
  require(e1071)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_UNDERB <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_UNDERT <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_UNDERTPhi <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_SMOTEB <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_SMOTET <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars_SMOTETPhi <- function(form,train,test,...) {
  require(earth)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- earth(form,train)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}


mc.rf <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_UNDERB <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_UNDERT <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_UNDERTPhi <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_SMOTEB <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_SMOTET <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf_SMOTETPhi <- function(form,train,test,...) {
  require(randomForest)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- randomForest(form,train,...)
  p <- predict(m,test)
  res <- WFoutput(rownames(test),responseValues(form,test),p)
  eval <- eval.stats(form,train,test,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.arima <- function(form,train,test,...) {
  require(forecast)
  ph <- phi.control(train[,ncol(train)], method="extremes")
  ls <- loss.control(train[,ncol(train)])

  trainY <- resp(form,train)
  trues <- resp(form,test)

  m <- auto.arima(trainY)
  data <- c(trainY,trues)

  p <- fitted(Arima(data,model=m))[(length(trainY)+1):length(data)]
  res <- WFoutput(rownames(test),trues,p)
  eval <- eval.stats_ARIMA(trainY,trues,p,ph,ls)
  workflowInformation(res) <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}


#RESAMPLING FUNCTIONS


# ===================================================
# Performs a random undersampling strategy for regression problems.
# Basically randomly removes a percentage of cases of the "class(es)"
# (bumps below a relevance threshold) selected by the user. 
# Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each "class".
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.5) 
#   alg.myUnd <- randUnderRegress(a7~., clean.algae, C.perc=C.perc)
#   alg.Bal <- randUnderRegress(a7~., clean.algae, C.perc="balance")
#   alg.Ext <- randUnderRegress(a7~., clean.algae, C.perc="extreme")
# 
# P. Branco, Jan 2016
# ---------------------------------------------------
randUnderRegressB <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance", repl=FALSE)
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is relevance determined automatically (default) with uba package or provided by the user
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the under-sampling percentage/s to apply to all/each
  #       "class" obtained with the relevance threshold. This percentage represents the 
  #       percentage of examples that is maintained in each "class". Examples 
  #       are randomly removed in each "class". Moreover, different percentages may 
  #       be provided for each "class". Alternatively, it may be "balance" or "extreme",
  #       cases where the under-sampling percentages are automatically estimated.
  # repl is it allowed to perform sampling with replacement

{
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  #   if(any(is.na(data))){
  #     stop("The data set provided contains NA values!")
  #   }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
#  y <- resp(form,data)
  y <- data[,tgt]
  attr(y,"names") <- rownames(data)
  s.y <- sort(y)

  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]
  
  imp <- sapply(obs.ind, function(x)mean(phi(x,pc)))
  
  und <- which(imp<thr.rel)
  ove <- which(imp>thr.rel)
  
  newdata <- NULL
  for(j in 1:length(ove)){
  newdata <- rbind(newdata,data[names(obs.ind[[ove[j]]]),]) # start with the examples from the minority "classes"
  }
  
  # set the undersampling percentages
  if(is.list(C.perc)){
    if(length(und) > 1 & length(C.perc)==1){ # the same under-sampling percentage is applied to all the "classes" 
      C.perc <- rep(C.perc[1],length(und))
    } else if(length(und)>length(C.perc) & length(C.perc)>1){
      stop("The number of under-sampling percentages must be equal to the number of bumps below the threshold defined!")      
    }else if(length(und)< length(C.perc)){
      stop("the number of under-sampling percentages must be at most the number of bumps below the threshold defined!")
    }# each class has its predefined over-sampling percentage 
  }else if(C.perc == "balance"){
    B <- sum(sapply(obs.ind[ove],length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und],length),5))
  }else if(C.perc== "extreme"){
    Bove <- sum(sapply(obs.ind[ove],length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und],length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length),5))
  }
  
  for(j in 1:length(und)){
    sel <- c()
    if(length(names(obs.ind[[und[j]]]))==1) {
      sel <- names(obs.ind[[und[j]]])
    } else {
      if((C.perc[[j]]*length(obs.ind[[und[j]]]))>length(obs.ind[[und[j]]])) {
        sel <- names(obs.ind[[und[j]]])
      } else {
        sel <- sample(names(obs.ind[[und[j]]]),C.perc[[j]]*length(obs.ind[[und[j]]]), replace=repl)
      }
    }
    newdata <- rbind(newdata, data[sel,])
  }
  
  newdata
}



# ===================================================
# Performs a random undersampling strategy for regression problems.
# Basically randomly removes a percentage of cases of the "class(es)"
# (bumps below a relevance threshold) selected by the user. 
# Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each "class".
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.5) 
#   alg.myUnd <- randUnderRegress(a7~., clean.algae, C.perc=C.perc)
#   alg.Bal <- randUnderRegress(a7~., clean.algae, C.perc="balance")
#   alg.Ext <- randUnderRegress(a7~., clean.algae, C.perc="extreme")
# 
# P. Branco, Jan 2016
# ---------------------------------------------------
randUnderRegressT <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance", repl=FALSE)
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is relevance determined automatically (default) with uba package or provided by the user
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the under-sampling percentage/s to apply to all/each
  #       "class" obtained with the relevance threshold. This percentage represents the 
  #       percentage of examples that is maintained in each "class". Examples 
  #       are randomly removed in each "class". Moreover, different percentages may 
  #       be provided for each "class". Alternatively, it may be "balance" or "extreme",
  #       cases where the under-sampling percentages are automatically estimated.
  # repl is it allowed to perform sampling with replacement

{
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  #   if(any(is.na(data))){
  #     stop("The data set provided contains NA values!")
  #   }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
#  y <- resp(form,data)
  y <- data[,tgt]
  attr(y,"names") <- rownames(data)
  s.y <- sort(y)

  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]
  
  imp <- sapply(obs.ind, function(x)mean(phi(x,pc)))
  
  und <- which(imp<thr.rel)
  ove <- which(imp>thr.rel)
  
  newdata <- NULL
  for(j in 1:length(ove)){
  newdata <- rbind(newdata,data[names(obs.ind[[ove[j]]]),]) # start with the examples from the minority "classes"
  }
  
  # set the undersampling percentages
  if(is.list(C.perc)){
    if(length(und) > 1 & length(C.perc)==1){ # the same under-sampling percentage is applied to all the "classes" 
      C.perc <- rep(C.perc[1],length(und))
    } else if(length(und)>length(C.perc) & length(C.perc)>1){
      stop("The number of under-sampling percentages must be equal to the number of bumps below the threshold defined!")      
    }else if(length(und)< length(C.perc)){
      stop("the number of under-sampling percentages must be at most the number of bumps below the threshold defined!")
    }# each class has its predefined over-sampling percentage 
  }else if(C.perc == "balance"){
    B <- sum(sapply(obs.ind[ove],length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und],length),5))
  }else if(C.perc== "extreme"){
    Bove <- sum(sapply(obs.ind[ove],length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und],length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length),5))
  }
  
  for(j in 1:length(und)){
    chdata <- data[names(obs.ind[[und[j]]]),] # select examples in the bump to undersample
    sdata <- chdata[order(names(obs.ind[[und[j]]])),]
    r <- nrow(sdata)
#    timeRange <- length(unique(as.POSIXct(rownames(sdata))))
    probs <- c()
    for(k in 1:r){
      probs <- c(probs,k/r)
    }

    sel <- c()
    if(length(rownames(sdata))==1) {
      sel <- rownames(sdata)
    } else {
        if((C.perc[[j]]*r)>r) {
          sel <- rownames(sdata)
        } else {
          sel <- sample(rownames(sdata),C.perc[[j]]*r, replace=repl, prob=probs)
        }
      
    }
    
    newdata <- rbind(newdata, data[sel,])
  }
  
  newdata <- newdata[order(rownames(newdata)),]
  newdata
}



# ===================================================
# Performs a random undersampling strategy for regression problems.
# Basically randomly removes a percentage of cases of the "class(es)"
# (bumps below a relevance threshold) selected by the user. 
# Alternatively, it can either balance all the 
# existing classes or it can "smoothly invert" the frequency
# of the examples in each "class".
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.5) 
#   alg.myUnd <- randUnderRegress(a7~., clean.algae, C.perc=C.perc)
#   alg.Bal <- randUnderRegress(a7~., clean.algae, C.perc="balance")
#   alg.Ext <- randUnderRegress(a7~., clean.algae, C.perc="extreme")
# 
# P. Branco, Jan 2016
# ---------------------------------------------------
randUnderRegressTPhi <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance", repl=FALSE)
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is relevance determined automatically (default) with uba package or provided by the user
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the under-sampling percentage/s to apply to all/each
  #       "class" obtained with the relevance threshold. This percentage represents the 
  #       percentage of examples that is maintained in each "class". Examples 
  #       are randomly removed in each "class". Moreover, different percentages may 
  #       be provided for each "class". Alternatively, it may be "balance" or "extreme",
  #       cases where the under-sampling percentages are automatically estimated.
  # repl is it allowed to perform sampling with replacement

{
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  #   if(any(is.na(data))){
  #     stop("The data set provided contains NA values!")
  #   }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
#  y <- resp(form,data)
  y <- data[,tgt]
  attr(y,"names") <- rownames(data)
  s.y <- sort(y)

  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]

  imp <- sapply(obs.ind, function(x)mean(phi(x,pc)))
  
  und <- which(imp<thr.rel)
  ove <- which(imp>thr.rel)
  
  newdata <- NULL
  for(j in 1:length(ove)){
  newdata <- rbind(newdata,data[names(obs.ind[[ove[j]]]),]) # start with the examples from the minority "classes"
  }
  
  # set the undersampling percentages
  if(is.list(C.perc)){
    if(length(und) > 1 & length(C.perc)==1){ # the same under-sampling percentage is applied to all the "classes" 
      C.perc <- rep(C.perc[1],length(und))
    } else if(length(und)>length(C.perc) & length(C.perc)>1){
      stop("The number of under-sampling percentages must be equal to the number of bumps below the threshold defined!")      
    }else if(length(und)< length(C.perc)){
      stop("the number of under-sampling percentages must be at most the number of bumps below the threshold defined!")
    }# each class has its predefined over-sampling percentage 
  }else if(C.perc == "balance"){
    B <- sum(sapply(obs.ind[ove],length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und],length),5))
  }else if(C.perc== "extreme"){
    Bove <- sum(sapply(obs.ind[ove],length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und],length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length),5))
  }
  
  for(j in 1:length(und)){
    chdata <- data[names(obs.ind[[und[j]]]),] # select examples in the bump to undersample
    sdata <- chdata[order(names(obs.ind[[und[j]]])),]
    r <- nrow(sdata)
    #    timeRange <- length(unique(as.POSIXct(rownames(sdata))))
    probs <- c()
    for(k in 1:r){
      probs <- c(probs,k/r)
    }
    s.rel <- phi(sdata[,tgt],pc)
    probsU <- probs*s.rel 
    
    sel <- c()
    if(length(rownames(sdata))==1) {
      sel <- rownames(sdata)
    } else {
      if((C.perc[[j]]*r)>r) {
        sel <- rownames(sdata)
      } else {
        sel <- sample(rownames(sdata),C.perc[[j]]*r, replace=repl, prob=probsU)
      }
      
    }
    
    newdata <- rbind(newdata, data[sel,])
  }
  
  newdata <- newdata[order(rownames(newdata)),]
  newdata 

}


## ===================================================
## Creating a SMOTE training sample for regression problems
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.1, 8) 
#   mysmote.alg <- smoteRegress(a7~., clean.algae, C.perc=C.perc)
#   smoteBal.alg <- smoteRegress(a7~., clean.algae, C.perc="balance")
#   smoteExt.alg <- smoteRegress(a7~., clean.algae, C.perc="extreme")
# 
#   ir<- iris[-c(95:130),]
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.5,2.5))
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.2,4), thr.rel=0.8)
#   smoteBalan.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="balance")
#   smoteExtre.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="extreme")
# 
#   rel <- matrix(0,ncol=3,nrow=0)
#   rel <- rbind(rel,c(2,1,0))
#   rel <- rbind(rel,c(3,0,0))
#   rel <- rbind(rel,c(4,1,0))
#
#   sP.ir <- smoteRegress(Sepal.Width~., ir, rel =rel, C.perc=list(4,0.5,4))
# 
# P. Branco, Jan 2016
# ---------------------------------------------------
smoteRegressB <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance",
                         k=5, repl=FALSE, dist="Euclidean", p=2)
  
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is the relevance determined automatically (default: "auto") with uba 
  #       package or provided by the user through a matrix. See examples.
  # thr.rel is the relevance threshold above which a case is considered
  #       as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the threshold
  #         are increased by this percentage. The undersampling percentage means 
  #         that the normal cases (cases below the threshold) are undersampled by 
  #         this percentage. Alternatively it may be "balance" or "extreme",
  #         cases where the sampling percentages are automatically estimated.
  # k is the number of neighbours to consider as the pool from where
  #               the new generated examples are generated
  # repl is it allowed to perform sampling with replacement
  # dist is the distance measure to be used (defaults to "Euclidean")
  # p is a parameter used when a p-norm is computed
{
 
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  
  if(any(is.na(data))){
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
  if (tgt < ncol(data)) {
    orig.order <- colnames(data)
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  if(is.na(thr.rel)){
    stop("Future work!")
  }
  

  y <- resp(form,data)
  s.y <- sort(y)
  
  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]
  
  newdata <- data.frame()
  
  if(is.list(C.perc)){
    if(length(C.perc)!= nbump) stop("The percentages provided must be the same length as the number of bumps!")
  }else if(C.perc=="balance"){ # estimate the percentages of over/under sampling
    B <- round(nrow(data)/nbump,0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if(C.perc == "extreme"){
    B <- round(nrow(data)/nbump,0)
    rescale <- nbump*B/sum(B^2/sapply(obs.ind,length))
    obj <- round((B^2/sapply(obs.ind, length))*rescale,2)
    C.perc <- round(obj/sapply(obs.ind, length),1)
  }
  
  for(i in 1:nbump){

    if(length(names(obs.ind[[i]]))==1) {
      newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
    } else {

      if(C.perc[[i]]==1){
        newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
      }else if(C.perc[[i]]>1){

        if(length(names(obs.ind[[i]]))<=k) {
          #newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
          newExs <- smote.exsRegressB(data[names(obs.ind[[i]]),],
                                   ncol(data),
                                   C.perc[[i]],
                                   (length(names(obs.ind[[i]]))-1),
                                   dist,
                                   p)
        # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        } else {
          newExs <- smote.exsRegressB(data[names(obs.ind[[i]]),],
                                   ncol(data),
                                   C.perc[[i]],
                                   k,
                                   dist,
                                   p)
        # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        }

        
        
      }else if(C.perc[[i]]<1){

        sel.maj <- sample(1:length(obs.ind[[i]]),
                          as.integer(C.perc[[i]]*length(obs.ind[[i]])),
                          replace=repl)
        newdata <- rbind(newdata, data[names(obs.ind[[i]][sel.maj]),])
        
      }

    }
  }
  
  if (tgt < ncol(data)) {
    newdata <- newdata[,cols]
    data <- data[,cols]
  }
  
  newdata
}



# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
#
# P. Branco, Jan 2016
# ---------------------------------------------------
smote.exsRegressB <- function(data, tgt, N, k, dist, p)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours
  # dist is the distance function used for the neighours computation
  # p is an integer used when a "p-norm" distance is selected
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(data) set of generate
  # examples with rare values on the target
{
  
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2])
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]
  
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  

  ranges <- rep(1,nC)
  if(length(nomatr)){
    for(x in (1:nC)[-c(nomatr)]) ranges[x] <- max(T[,x]) - min(T[,x])
  } else{
    for(x in (1:nC)) ranges[x] <- max(T[,x]) - min(T[,x])
  }

  kNNs <-neighbours(tgt, data, dist, p, k)
    
  nexs <-  as.integer(N-1) # nr of examples to generate for each rare case
  extra <- as.integer(nT*(N-1-nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  new <- matrix(nrow=nexs*nT+extra,ncol=nC)    # the new cases
 
  if(nexs){
    for(i in 1:nT) {
    
        
      for(n in 1:nexs) {
        # select randomly one of the k NNs
        neig <- sample(1:k,1)
      
        # the attribute values of the generated case
        difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
        new[(i-1)*nexs+n,-tgt] <- T[i,-tgt]+runif(1)*difs
        for(a in nomatr) # nominal attributes are randomly selected among the existing values of seed and the selected neighbour 
          new[(i-1)*nexs+n,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
        
        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for(x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i,x] - new[(i-1)*nexs+n,x])/ranges[x]
          d2 <- abs(T[kNNs[i,neig],x] - new[(i-1)*nexs+n,x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
          d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        new[(i-1)*nexs+n,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
        
      }
    }
  }
  if(extra){
    count<-1
    for (i in idx){
    
      # select randomly one of the k NNs
      neig <- sample(1:k,1) 
      
      # the attribute values of the generated case
      difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
      new[nexs*nT+count,-tgt] <- T[i,-tgt]+runif(1)*difs
      for(a in nomatr)
        new[nexs*nT+count,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
      
      
      # now the target value (weighted (by inverse distance) average)
      d1 <- d2 <- 0
      for(x in (1:nC)[-c(nomatr,tgt)]) {
        d1 <- abs(T[i,x] - new[nexs*nT+count,x])/ranges[x]
        d2 <- abs(T[kNNs[i,neig],x] - new[nexs*nT+count,x])/ranges[x]
      }
      if (length(nomatr)) {
        d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
        d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
      }
      # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
      new[nexs*nT+count,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
    
      count <- count+1
    }
  }

  newCases <- data.frame(new)

  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))
  
  colnames(newCases) <- colnames(data)
  newCases
  
}



## ===================================================
## Creating a SMOTE training sample for regression problems
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.1, 8) 
#   mysmote.alg <- smoteRegress(a7~., clean.algae, C.perc=C.perc)
#   smoteBal.alg <- smoteRegress(a7~., clean.algae, C.perc="balance")
#   smoteExt.alg <- smoteRegress(a7~., clean.algae, C.perc="extreme")
# 
#   ir<- iris[-c(95:130),]
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.5,2.5))
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.2,4), thr.rel=0.8)
#   smoteBalan.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="balance")
#   smoteExtre.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="extreme")
# 
#   rel <- matrix(0,ncol=3,nrow=0)
#   rel <- rbind(rel,c(2,1,0))
#   rel <- rbind(rel,c(3,0,0))
#   rel <- rbind(rel,c(4,1,0))
#
#   sP.ir <- smoteRegress(Sepal.Width~., ir, rel =rel, C.perc=list(4,0.5,4))
# 
# P. Branco, jan 2016
# ---------------------------------------------------
smoteRegressT <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance",
                         k=5, repl=FALSE, dist="Euclidean", p=2)
  
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is the relevance determined automatically (default: "auto") with uba 
  #       package or provided by the user through a matrix. See examples.
  # thr.rel is the relevance threshold above which a case is considered
  #       as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the threshold
  #         are increased by this percentage. The undersampling percentage means 
  #         that the normal cases (cases below the threshold) are undersampled by 
  #         this percentage. Alternatively it may be "balance" or "extreme",
  #         cases where the sampling percentages are automatically estimated.
  # k is the number of neighbours to consider as the pool from where
  #               the new generated examples are generated
  # repl is it allowed to perform sampling with replacement
  # dist is the distance measure to be used (defaults to "Euclidean")
  # p is a parameter used when a p-norm is computed
{
 
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  
  if(any(is.na(data))){
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
  if (tgt < ncol(data)) {
    orig.order <- colnames(data)
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  if(is.na(thr.rel)){
    stop("Future work!")
  }
  

  y <- resp(form,data)
  s.y <- sort(y)
  
  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]

  
  newdata <- data.frame()
  
  if(is.list(C.perc)){
    if(length(C.perc)!= nbump) stop("The percentages provided must be the same length as the number of bumps!")
  }else if(C.perc=="balance"){ # estimate the percentages of over/under sampling
    B <- round(nrow(data)/nbump,0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if(C.perc == "extreme"){
    B <- round(nrow(data)/nbump,0)
    rescale <- nbump*B/sum(B^2/sapply(obs.ind,length))
    obj <- round((B^2/sapply(obs.ind, length))*rescale,2)
    C.perc <- round(obj/sapply(obs.ind, length),1)
  }
  
  for(i in 1:nbump){

    if(length(names(obs.ind[[i]]))==1) {
      newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
    } else {
      if(C.perc[[i]]==1){
      newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
      }else if(C.perc[[i]]>1){

        if(length(names(obs.ind[[i]]))<=k) {
          #newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
          newExs <- smote.exsRegressT(data[names(obs.ind[[i]]),],
                                     ncol(data),
                                     C.perc[[i]],
                                     (length(names(obs.ind[[i]]))-1),
                                     dist,
                                     p)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        } else {
          newExs <- smote.exsRegressT(data[names(obs.ind[[i]]),],
                                     ncol(data),
                                     C.perc[[i]],
                                     k,
                                     dist,
                                     p)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        }
        
      }else if(C.perc[[i]]<1){
        
          chdata <- data[names(obs.ind[[i]]),] # select examples in the bump to undersample
          sdata <- chdata[order(names(obs.ind[[i]])),]
          r <- nrow(sdata)
          #    timeRange <- length(unique(as.POSIXct(rownames(sdata))))
          probs <- c()
          for(pr in 1:r){
            probs <- c(probs,pr/r)
          }
          sel <- sample(rownames(sdata),C.perc[[i]]*r, replace=repl, prob=probs)
          newdata <- rbind(newdata, data[sel,])
        
      }
    }

    
  }
  
  if (tgt < ncol(data)) {
    newdata <- newdata[,cols]
    data <- data[,cols]
  }
  
  newdata
}



# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
#
# P. Branco, Jan 2016
# ---------------------------------------------------
smote.exsRegressT <- function(data, tgt, N, k, dist, p)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours
  # dist is the distance function used for the neighours computation
  # p is an integer used when a "p-norm" distance is selected
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(data) set of generate
  # examples with rare values on the target
{
  
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2])
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]
  
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  

  ranges <- rep(1,nC)
  if(length(nomatr)){
    for(x in (1:nC)[-c(nomatr)]) ranges[x] <- max(T[,x]) - min(T[,x])
  } else{
    for(x in (1:nC)) ranges[x] <- max(T[,x]) - min(T[,x])
  }
  
  # order the data by time so that the most recent examples have the higher line numbers
  data <- data[order(rownames(data)),]
  kNNs <-neighbours(tgt, data, dist, p, k)
    
  nexs <-  as.integer(N-1) # nr of examples to generate for each rare case
  extra <- as.integer(nT*(N-1-nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  new <- matrix(nrow=nexs*nT+extra,ncol=nC)    # the new cases
 
  if(nexs){
    for(i in 1:nT) {
    
        
      for(n in 1:nexs) {
        # select randomly one of the k NNs
        #neig <- sample(1:k,1)
        # select the most recent neighbour
        neig <- match(max(kNNs[i,]), kNNs[i,])
      
        # the attribute values of the generated case
        difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
        new[(i-1)*nexs+n,-tgt] <- T[i,-tgt]+runif(1)*difs
        for(a in nomatr) # nominal attributes are randomly selected among the existing values of seed and the selected neighbour 
          new[(i-1)*nexs+n,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
        
        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for(x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i,x] - new[(i-1)*nexs+n,x])/ranges[x]
          d2 <- abs(T[kNNs[i,neig],x] - new[(i-1)*nexs+n,x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
          d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        new[(i-1)*nexs+n,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
        
      }
    }
  }
  if(extra){
    count<-1
    for (i in idx){
    
#       # select randomly one of the k NNs
#       neig <- sample(1:k,1) 
      # select the most recent neighbour
      neig <- match(max(kNNs[i,]), kNNs[i,])
      # the attribute values of the generated case
      difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
      new[nexs*nT+count,-tgt] <- T[i,-tgt]+runif(1)*difs
      for(a in nomatr)
        new[nexs*nT+count,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
      
      
      # now the target value (weighted (by inverse distance) average)
      d1 <- d2 <- 0
      for(x in (1:nC)[-c(nomatr,tgt)]) {
        d1 <- abs(T[i,x] - new[nexs*nT+count,x])/ranges[x]
        d2 <- abs(T[kNNs[i,neig],x] - new[nexs*nT+count,x])/ranges[x]
      }
      if (length(nomatr)) {
        d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
        d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
      }
      # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
      new[nexs*nT+count,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
    
      count <- count+1
    }
  }

  newCases <- data.frame(new)

  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))
  
  colnames(newCases) <- colnames(data)
  newCases
  
}



## ===================================================
## Creating a SMOTE training sample for regression problems
# 
# Examples:
#   library(DMwR)
#   data(algae)
#   clean.algae <- algae[complete.cases(algae),]
#   C.perc=list(0.1, 8) 
#   mysmote.alg <- smoteRegress(a7~., clean.algae, C.perc=C.perc)
#   smoteBal.alg <- smoteRegress(a7~., clean.algae, C.perc="balance")
#   smoteExt.alg <- smoteRegress(a7~., clean.algae, C.perc="extreme")
# 
#   ir<- iris[-c(95:130),]
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.5,2.5))
#   mysmote.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc=list(0.2,4), thr.rel=0.8)
#   smoteBalan.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="balance")
#   smoteExtre.iris <- smoteRegress(Sepal.Width~., ir, dist="HEOM", C.perc="extreme")
# 
#   rel <- matrix(0,ncol=3,nrow=0)
#   rel <- rbind(rel,c(2,1,0))
#   rel <- rbind(rel,c(3,0,0))
#   rel <- rbind(rel,c(4,1,0))
#
#   sP.ir <- smoteRegress(Sepal.Width~., ir, rel =rel, C.perc=list(4,0.5,4))
# 
# P. Branco, jan 2016
# ---------------------------------------------------
smoteRegressTPhi <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance",
                         k=5, repl=FALSE, dist="Euclidean", p=2)
  
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # rel is the relevance determined automatically (default: "auto") with uba 
  #       package or provided by the user through a matrix. See examples.
  # thr.rel is the relevance threshold above which a case is considered
  #       as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the threshold
  #         are increased by this percentage. The undersampling percentage means 
  #         that the normal cases (cases below the threshold) are undersampled by 
  #         this percentage. Alternatively it may be "balance" or "extreme",
  #         cases where the sampling percentages are automatically estimated.
  # k is the number of neighbours to consider as the pool from where
  #               the new generated examples are generated
  # repl is it allowed to perform sampling with replacement
  # dist is the distance measure to be used (defaults to "Euclidean")
  # p is a parameter used when a p-norm is computed
{
 
#  require(uba, quietly=TRUE)
  suppressWarnings(suppressPackageStartupMessages(library('uba')))
  
  if(any(is.na(data))){
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  
  if (tgt < ncol(data)) {
    orig.order <- colnames(data)
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  if(is.na(thr.rel)){
    stop("Future work!")
  }
  

  y <- resp(form,data)
  s.y <- sort(y)
  
  if (is.matrix(rel)){ 
    pc <- phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){ 
    pc <- rel
  }else if(rel=="auto"){
    pc <- phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")
  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes
  
  # collect the indexes in each "class"
  nbumps <- length(bumps) +1 
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]

  
  newdata <- data.frame()
  
  if(is.list(C.perc)){
    if(length(C.perc)!= nbump) stop("The percentages provided must be the same length as the number of bumps!")
  }else if(C.perc=="balance"){ # estimate the percentages of over/under sampling
    B <- round(nrow(data)/nbump,0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if(C.perc == "extreme"){
    B <- round(nrow(data)/nbump,0)
    rescale <- nbump*B/sum(B^2/sapply(obs.ind,length))
    obj <- round((B^2/sapply(obs.ind, length))*rescale,2)
    C.perc <- round(obj/sapply(obs.ind, length),1)
  }
  
  for(i in 1:nbump){

    if(length(names(obs.ind[[i]]))==1) {
      newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
    } else {
      if(C.perc[[i]]==1){
        newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
      }else if(C.perc[[i]]>1){

        if(length(names(obs.ind[[i]]))<=k) {
          #newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
          newExs <- smote.exsRegressTPhi(data[names(obs.ind[[i]]),],
                                     ncol(data),
                                     C.perc[[i]],
                                     (length(names(obs.ind[[i]]))-1),
                                     dist,
                                     p,
                                     pc)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        } else {
          newExs <- smote.exsRegressTPhi(data[names(obs.ind[[i]]),],
                                     ncol(data),
                                     C.perc[[i]],
                                     k,
                                     dist,
                                     p,
                                     pc)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        }
        
      }else if(C.perc[[i]]<1){
        
        chdata <- data[names(obs.ind[[i]]),] # select examples in the bump to undersample
        sdata <- chdata[order(names(obs.ind[[i]])),]
        r <- nrow(sdata)
        #    timeRange <- length(unique(as.POSIXct(rownames(sdata))))
        probs <- c()
        for(pr in 1:r){
          probs <- c(probs,pr/r)
        }
        s.rel <- phi(sdata[,tgt],pc)
        probsU <- probs*s.rel 
        
        sel <- sample(rownames(sdata),C.perc[[i]]*r, replace=repl, prob=probsU)
        newdata <- rbind(newdata, data[sel,])
        
      }
    }

    
  }
  
  if (tgt < ncol(data)) {
    newdata <- newdata[,cols]
    data <- data[,cols]
  }
  
  newdata
}



# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
#
# P. Branco, Jan 2016
# ---------------------------------------------------
smote.exsRegressTPhi <- function(data, tgt, N, k, dist, p, pc)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours
  # dist is the distance function used for the neighours computation
  # p is an integer used when a "p-norm" distance is selected
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(data) set of generate
  # examples with rare values on the target
{
  
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2])
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]
  
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  

  ranges <- rep(1,nC)
  if(length(nomatr)){
    for(x in (1:nC)[-c(nomatr)]) ranges[x] <- max(T[,x]) - min(T[,x])
  } else{
    for(x in (1:nC)) ranges[x] <- max(T[,x]) - min(T[,x])
  }

  # order the data by time so that the most recent examples have the higher line numbers
  data <- data[order(rownames(data)),]
  
  kNNs <-neighbours(tgt, data, dist, p, k)
    
  nexs <-  as.integer(N-1) # nr of examples to generate for each rare case
  extra <- as.integer(nT*(N-1-nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  new <- matrix(nrow=nexs*nT+extra,ncol=nC)    # the new cases
 
  if(nexs){
    for(i in 1:nT) {
      for(n in 1:nexs) {
        # select the nearest neighbour more recent and with higher phi
        y.rel <- phi(data[kNNs[i,],tgt], pc)
        pos.eval <- (kNNs[i,]/max(kNNs[i,]))*y.rel
        neig <- match(max(pos.eval), pos.eval)
        # the attribute values of the generated case
        difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
        new[(i-1)*nexs+n,-tgt] <- T[i,-tgt]+runif(1)*difs
        for(a in nomatr) # nominal attributes are randomly selected among the existing values of seed and the selected neighbour 
          new[(i-1)*nexs+n,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
        
        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for(x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i,x] - new[(i-1)*nexs+n,x])/ranges[x]
          d2 <- abs(T[kNNs[i,neig],x] - new[(i-1)*nexs+n,x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
          d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        new[(i-1)*nexs+n,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
        
      }
    }
  }
  if(extra){
    count<-1
    for (i in idx){
    
      # select the nearest neighbour more recent and with higher phi
      y.rel <- phi(data[kNNs[i,],tgt], pc)
      pos.eval <- (kNNs[i,]/max(kNNs[i,]))*y.rel
      neig <- match(max(pos.eval), pos.eval)
      
      # the attribute values of the generated case
      difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
      new[nexs*nT+count,-tgt] <- T[i,-tgt]+runif(1)*difs
      for(a in nomatr)
        new[nexs*nT+count,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]
      
      
      # now the target value (weighted (by inverse distance) average)
      d1 <- d2 <- 0
      for(x in (1:nC)[-c(nomatr,tgt)]) {
        d1 <- abs(T[i,x] - new[nexs*nT+count,x])/ranges[x]
        d2 <- abs(T[kNNs[i,neig],x] - new[nexs*nT+count,x])/ranges[x]
      }
      if (length(nomatr)) {
        d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
        d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
      }
      # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
      new[nexs*nT+count,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)
    
      count <- count+1
    }
  }

  newCases <- data.frame(new)

  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))
  
  colnames(newCases) <- colnames(data)
  newCases
  
}


#DEFINITION OF VARIABLES FOR RESAMPLING

ds <- create.data(new_d, 10) #Create the embed
form <- as.formula(V10 ~ .) #Define the formula for performanceEstimation

#ds <- knnImputation(ds) #Imputation of NA values (DMwR package)

exp <- performanceEstimation(PredTask(form,ds),
                             c(Workflow("mc.lm"),
                               Workflow("mc.lm_UNDERB"),
                               Workflow("mc.lm_UNDERT"),
                               Workflow("mc.lm_UNDERTPhi"),
                               Workflow("mc.lm_SMOTEB"),
                               Workflow("mc.lm_SMOTET"),
                               Workflow("mc.lm_SMOTETPhi"),
                               Workflow("mc.svm"),
                               Workflow("mc.svm_UNDERB"),
                               Workflow("mc.svm_UNDERT"),
                               Workflow("mc.svm_UNDERTPhi"),
                               Workflow("mc.svm_SMOTEB"),
                               Workflow("mc.svm_SMOTET"),
                               Workflow("mc.svm_SMOTETPhi"),
                               Workflow("mc.mars"),
                               Workflow("mc.mars_UNDERB"),
                               Workflow("mc.mars_UNDERT"),
                               Workflow("mc.mars_UNDERTPhi"),
                               Workflow("mc.mars_SMOTEB"),
                               Workflow("mc.mars_SMOTET"),
                               Workflow("mc.mars_SMOTETPhi"),
                               Workflow("mc.rf"),
                               Workflow("mc.rf_UNDERB"),
                               Workflow("mc.rf_UNDERT"),
                               Workflow("mc.rf_UNDERTPhi"),
                               Workflow("mc.rf_SMOTEB"),
                               Workflow("mc.rf_SMOTET"),
                               Workflow("mc.rf_SMOTETPhi"),
                               Workflow("mc.arima")),
                             EstimationTask("totTime",method=MonteCarlo(nReps=50,szTrain=.5,szTest=.25))
)

