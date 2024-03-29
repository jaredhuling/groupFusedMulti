#' Cross Validation for the groupFusedMulti function
#'
#' @param x input matrix or SparseMatrix of dimension nobs x nvars. Each row is an observation,
#' each column corresponds to a covariate
#' @param y numeric response vector of length nobs
#' @param outcome.groups a vector of length equal to the number of columns in the outcome matrix \code{y}
#' @param lambda A user-specified sequence of lambda values. Left unspecified, the a sequence of lambda values is
#' automatically computed, ranging uniformly on the log scale over the relevant range of lambda values.
#' @param lambda.fused A user-specified sequence of fused lasso tuning parameters.
#' @param type.measure One of \code{c("mse","deviance","class","auc","mae","brier")} indicating measure to evaluate for cross-validation. The default is \code{type.measure = "deviance"}, 
#' which uses squared-error for gaussian models (a.k.a \code{type.measure = "mse"} there), deviance for logistic
#' regression. \code{type.measure = "class"} applies to binomial only. \code{type.measure = "auc"} is for two-class logistic 
#' regression only. \code{type.measure = "mse"} or \code{type.measure = "mae"} (mean absolute error) can be used by all models;
#' they measure the deviation from the fitted mean to the response. \code{type.measure = "brier"} is for models with 
#' \code{family = "coxph"} and will compute the Brier score.
#' @param nfolds number of folds for cross-validation. default is 10. 3 is smallest value allowed. 
#' @param foldid an optional vector of values between 1 and nfold specifying which fold each observation belongs to.
#' @param grouped Like in \pkg{glmnet}, this is an experimental argument, with default \code{TRUE}, and can be ignored by most users. 
#' For all models, this refers to computing nfolds separate statistics, and then using their mean and estimated standard 
#' error to describe the CV curve. If \code{grouped = FALSE}, an error matrix is built up at the observation level from the 
#' predictions from the \code{nfold} fits, and then summarized (does not apply to \code{type.measure = "auc"}). 
#' @param keep If \code{keep = TRUE}, a prevalidated list of arrasy is returned containing fitted values for each observation 
#' and each value of lambda for each model. This means these fits are computed with this observation and the rest of its
#' fold omitted. The folid vector is also returned. Default is \code{keep = FALSE}
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand, such as \pkg{doMC}.
#' @param ... parameters to be passed to groupFusedMulti
#' @importFrom stats predict
#' @importFrom stats stepfun
#' @importFrom stats weighted.mean
#' @import foreach
#' @return An object with S3 class "cv.groupFusedMulti"
#'
#'
#' @import Rcpp
#'
#' @export
#' @examples
#' 
#' set.seed(123)
#'
#' dat.sim <- gen_sparse_multivar_data(nvars = 15L,
#'                    noutcomes = 8L,
#'                    nobs = 100L,
#'                    nobs.test = 100L,
#'                    num.nonzero.vars = 10,
#'                    outcome.groups = rbind(c(1,1,1,2,2,2,2,2),
#'                                           c(1,1,1,2,2,3,3,3)))
#'
#' x        <- dat.sim$x
#' x.test   <- dat.sim$x.test
#' y        <- dat.sim$y
#' y.test   <- dat.sim$y.test
#' beta     <- dat.sim$beta
#' 
#' \donttest{
#'
#' outcome_groups <- rbind(c(1,1,1,2,2,2,2,2),
#'                         c(1,1,1,2,2,3,3,3))
#'                         
#' fit.adapt <- cv.groupFusedMulti(x, y,
#'                                 nlambda        = 50,
#'                                 lambda.fused = c(0.000005, 0.00001, 0.000025, 0.00005, 0.0001),
#'                                 outcome.groups = outcome_groups,
#'                                 adaptive.lasso = TRUE, adaptive.fused = TRUE,
#'                                 gamma          = 0.25,
#'                                 nfolds         = 5)
#' 
#' est.coefs <- predict(fit.adapt, type = "coef")
#' colnames(beta) <- colnames(est.coefs)
#' 
#' round(est.coefs, 3)
#' beta
#' 
#' preds.a <- predict(fit.adapt, x.test, type = 'response')
#' 
#' ## rmse for each outcome
#' sqrt(colMeans((y.test - preds.a) ^ 2))
#' 
#' ## avg rmse
#' mean(sqrt(colMeans((y.test - preds.a) ^ 2)))
#'                    
#' }
cv.groupFusedMulti <- function(x, y,
                               outcome.groups,
                               lambda       = NULL,
                               lambda.fused = NULL,
                               type.measure = c("mse","deviance","class","auc","mae","brier"),
                               nfolds       = 10,
                               foldid,
                               grouped      = TRUE,
                               keep         = FALSE,
                               parallel     = FALSE,
                               ...)
{
  if(missing(type.measure)) type.measure <- "default"
  else type.measure <- match.arg(type.measure)
  if(!is.null(lambda)&&length(lambda)<2) stop("Need more than one value of lambda for cv.groupFusedMulti")
  N <- nrow(x)
  #if(missing(weights))weights=rep(1.0,N)else weights=as.double(weights)
  ###Fit the model once to get dimensions etc of output
  ###Next we construct a call, that could recreate a groupFusedMulti object - tricky
  ### This if for predict, exact=TRUE
  groupFusedMulti.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure","nfolds","foldid","grouped","keep"), names(groupFusedMulti.call), FALSE)

  if(any(which)) groupFusedMulti.call <- groupFusedMulti.call[-which]
  
  if (inherits(y, "data.frame")) y <- as.matrix(y)

  groupFusedMulti.call[[1]] <- as.name("groupFusedMulti")
  groupFusedMulti.object <- groupFusedMulti(x, y, outcome.groups, 
                                            lambda = lambda,
                                            lambda.fused = lambda.fused, 
                                            ...)
  groupFusedMulti.object$call <- groupFusedMulti.call
  

  if(missing(foldid)) 
  {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else 
  {
    nfolds <- max(foldid)
  }

  if(nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  ###First try and do it using foreach if parallel is TRUE
  if (parallel) 
  {
    outlist = foreach (i=seq(nfolds), .packages=c("groupFusedMulti")) %dopar% {
      which <- foldid==i
      if(is.matrix(y)) y_sub <- y[!which,] else y_sub <- y[!which]
      #if(is.offset)offset_sub=as.matrix(offset)[!which,]
      #else offset_sub=NULL
      groupFusedMulti(x[!which,,drop=FALSE], 
                      y_sub, 
                      outcome.groups, 
                      lambda       = lambda,
                      lambda.fused = lambda.fused,
                      ...)
    }
  } else
  {
    for(i in seq(nfolds))
    {
      which <- foldid==i
      if(is.matrix(y)) y_sub=y[!which,] else y_sub=y[!which]
      #if(is.offset)offset_sub=as.matrix(offset)[!which,]
      #else offset_sub=NULL
      #cat("fitting fold", i, "\n")
      outlist[[i]] <- groupFusedMulti(x[!which,,drop=FALSE], 
                                      y_sub, 
                                      outcome.groups, 
                                      lambda       = lambda,
                                      lambda.fused = lambda.fused, 
                                      ...)
    }
  }
  ###What to do depends on the type.measure and the model fit
  fun     <- paste("cv", class(groupFusedMulti.object)[[2]], sep=".")
  cvstuff <- do.call(fun,
                     list(outlist,
                          groupFusedMulti.object$lambda,
                          groupFusedMulti.object$lambda.fused,
                          x,
                          y,
                          outcome.groups,
                          foldid,
                          type.measure,
                          grouped,
                          keep))
  cvm     <- cvstuff$cvm
  cvsd    <- cvstuff$cvsd
  cvname  <- cvstuff$name

  if (is.list(cvm))
  {
    names(cvm) <- names(cvsd) <- groupFusedMulti.object$lambda.fused
    out=list(lambda        = groupFusedMulti.object$lambda,
             lambda.fused  = groupFusedMulti.object$lambda.fused,
             cvm           = cvm,
             cvsd          = cvsd,
             cvup          = lapply(1:length(cvm), function(i) cvm[[i]] + cvsd[[i]]),
             cvlo          = lapply(1:length(cvm), function(i) cvm[[i]] - cvsd[[i]]),
             cvraw         = cvstuff$cvraw,
             name          = cvname,
             groupFusedMulti.fit = groupFusedMulti.object)
    if(keep) out <- c(out,list(fit.preval = cvstuff$fit.preval,
                               foldid     = foldid))
    
    if (type.measure == "auc")
    {
      lam.fused.min.idx <- unname(which.max(unlist(sapply(cvm, max))))
      
      if (is.matrix(groupFusedMulti.object$lambda))
      {
        lamin <- c(getmin(groupFusedMulti.object$lambda[,lam.fused.min.idx], -cvm[[lam.fused.min.idx]], cvsd[[lam.fused.min.idx]]),
                   lambda.fused.min = groupFusedMulti.object$lambda.fused[lam.fused.min.idx],
                   which.min.fused = lam.fused.min.idx)
      } else
      {
        lamin <- c(getmin(groupFusedMulti.object$lambda, -cvm[[lam.fused.min.idx]], cvsd[[lam.fused.min.idx]]),
                   lambda.fused.min = groupFusedMulti.object$lambda.fused[lam.fused.min.idx],
                   which.min.fused = lam.fused.min.idx)
      }
    } else
    {
      lam.fused.min.idx <- unname(which.min(unlist(sapply(cvm, min))))
      if (is.matrix(groupFusedMulti.object$lambda))
      {
        lamin <- c(getmin(groupFusedMulti.object$lambda[,lam.fused.min.idx], cvm[[lam.fused.min.idx]], cvsd[[lam.fused.min.idx]]),
                   lambda.fused.min = groupFusedMulti.object$lambda.fused[lam.fused.min.idx],
                   which.min.fused = lam.fused.min.idx)
      } else
      {
        lamin <- c(getmin(groupFusedMulti.object$lambda, cvm[[lam.fused.min.idx]], cvsd[[lam.fused.min.idx]]),
                   lambda.fused.min = groupFusedMulti.object$lambda.fused[lam.fused.min.idx],
                   which.min.fused = lam.fused.min.idx)
      }
    }
    
  } else
  {
    out=list(lambda        = groupFusedMulti.object$lambda,
             cvm           = cvm,
             cvsd          = cvsd,
             cvup          = cvm + cvsd,
             cvlo          = cvm - cvsd,
             name          = cvname,
             groupFusedMulti.fit = groupFusedMulti.object)
    if(keep) out <- c(out,list(fit.preval=cvstuff$fit.preval,foldid=foldid))
    lamin = if(type.measure=="auc") getmin(groupFusedMulti.object$lambda,-cvm,cvsd)
    else getmin(groupFusedMulti.object$lambda,cvm,cvsd)
  }

  obj <- c(out,as.list(lamin))
  class(obj) <- "cv.groupFusedMulti"
  obj
}

# cv.venncoxph <- function(outlist,lambda,lambda.fused,x,y,groups,foldid,type.measure,grouped,keep=FALSE)
# {
#     typenames <- c(deviance="Coxph Deviance")
#     if(type.measure == "default") type.measure <- "deviance"
#     
#     if(!match(type.measure,c("deviance","brier"),FALSE))
#     {
#         warning("Only 'deviance' or 'brier' available for coxph model; as default 'deviance' will be used")
#         type.measure="deviance"
#     }
# 
#     N      <- nrow(y)
#     nfolds <- max(foldid)
#     if( (N/nfolds < 10) && type.measure == "auc")
#     {
#         warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
#         type.measure="deviance"
#     }
#     if( (N/nfolds <3) && grouped)
#     {
#         warning("Option grouped=FALSE enforced in cv.groupFusedMulti, since < 3 observations per fold",call.=FALSE)
#         grouped <- FALSE
#     }
# 
# 
#     #if(!is.null(offset)){
#     #  is.offset=TRUE
#     #  offset=drop(offset)
#     #}else is.offset=FALSE
#     devmat <- matrix(NA, nfolds, length(lambda))
#     nlams  <- double(nfolds)
# 
#     logpl <- function(link, yyy)
#     {
#         t <- yyy[,1]
#         c <- yyy[,2]
# 
#         index.t <- order(t)
# 
#         link[index.t,] <- apply(link[index.t,], 2, function (tt)
#         {
#             if (max(tt) > 100)
#             {
#                 kkk <- floor(max(tt) / 100)
#                 tt  <- tt - 100*kkk
#             } else 
#             {
#                 tt <- tt
#             }
#         })
# 
#         link.order.exp <- exp(link[index.t,])
#         c.order        <- c[index.t]
#         #       r <- apply(exp(link),2,function(tt){rev(cumsum(rev(tt[index.t])))})
#         logL <- apply(link.order.exp,2,function(tt){r <- rev(cumsum(rev(tt)))
#         sum(log(tt[c.order==1])-log(r[c.order==1]))})
#     }
# 
#     brier <- function(survive,yytest,yytrain,cenFun,cenTime)
#     {
#         index.test  <- order(yytest[,1])
#         ytest       <- yytest[index.test]
#         survive     <- survive[index.test,]
# 
#         index.train <- order(yytrain[,1])
#         ytrain      <- yytrain[index.train]
# 
#         nobs        <- dim(survive)[1]
#         brierScore  <- array(0,c(nobs,1))
#         sub         <- 100 * nobs
# 
#         for (i in 1:nobs)
#         {
#             if (length(survive[i,])<length(ytrain[,1]))
#             {
#                 survm <- c(rep(1,times = length(ytrain[,1])-length(survive[i,])), survive[i,])
#             } else 
#             {
#                 survm <- survive[i,]
#             }
#             surFun  <- stepfun(ytrain[,1],c(1,survm), f=0)
#             proFun1 <- stepfun(ytest[i,1],c(1,0), f=0)
#             proFun2 <- stepfun(ytest[i,1],c(0,1), f=0)
#             knots   <- sort(unique(c(ytrain[,1],ytest[,1],cenTime)))
# 
#             Fun1 <- function(t) { (surFun(t))^2*proFun2(t)}
#             Fun2 <- function(t) { (1-surFun(t))^2*proFun1(t)/cenFun(t)}
#             if (ytest[i,2] == 1)
#             {
#                 brierScore[i] <- brierScore[i] + integrate.step(Fun1,ytest[i,1],max(ytest[,1]),knots) / cenFun(ytest[i,1]) + integrate.step(Fun2,0,ytest[i,1], knots)
#             } else 
#             {
#                 brierScore[i] <- brierScore[i] + integrate.step(Fun2,0,ytest[i,1], knots)
#             }
#     #         for (j in 1:(dim(ytrain)[1]-1)){
#     #             if (ytest[i,1]>ytrain[j,1] && ytest[]){
#     #                 brierScore[i] <- brierScore[i] + (1-survive[i,j])^2*(min(ytrain[j+1,1],ytest[i,1])-ytrain[j,1]) / cenFun(ytrain[j,1])
#     #             }
#     #             if (j == dim(ytrain)[1]-1 & ytest[i,1]>ytrain[j+1,1]){
#     #                 brierScore[i] <- brierScore[i] + (1-survive[i,j+1])^2*(ytest[i,1]-ytrain[j+1,1]) / cenFun(ytrain[j,1]-10^(-6))
#     #             }
#     #         }
#     #         print(6)
#     #
#     #         for (j in 2:(dim(ytrain)[1])){
#     #             if (ytest[i,1]<ytrain[j,1] && ytest[i,2]==1){
#     #                 brierScore[i] <- brierScore[i] + (survive[i,j-1])^2*(ytrain[j,1]-max(ytrain[j-1,1],ytest[i,1])) / cenFun(ytest[i,1])
#     #             }
#     #         }
#     #         print(7)
#     #
#         }
#         brierScore
#     }
# 
#     integrate.step <- function(Fun,lower,upper,knots)
#     {
#         index.lower <- which(knots==lower)
#         index.upper <- which(knots==upper)-1
#         res <- 0
#         for(i in index.lower:index.upper)
#         {
#             res <- res + Fun(knots[i+1] - 0.01) * (knots[i+1]-knots[i])
#         }
#         res
#     }
# 
# 
#     for(i in seq(nfolds))
#     {
#         which  <- foldid==i
#         fitobj <- outlist[[i]]
#         #if(is.offset)off_sub=offset[which]
#         nlami  <- length(outlist[[i]]$lam)
# 
#         if(type.measure == "deviance")
#         {
#             linkout = predict(fitobj, newx = x[!which,,drop=FALSE], group.mat = groups[!which,,drop=FALSE], type="link")
#             #print(anyNA(linkout))
#             link = predict(fitobj, newx = x, group.mat = groups, type="link")
#             #print(anyNA(link))
#             devmat[i,seq(nlami)] = -logpl(link,y) + logpl(linkout,y[!which,,drop=FALSE])
#         } else 
#         {
#             surv = predict(fitobj, newx = x[which,], group.mat = groups[which,,drop=FALSE], type="survival")
#             cenSurv = Surv(time = y[which,1], event = 1-y[which,2])
#             cenSurvFun = c(1,1, survfit(cenSurv~1)$surv)
#             cenTime = c(0, survfit(cenSurv~1)$time)
#             cenFun = stepfun(cenTime, cenSurvFun, f=0)
# 
#             for (j in 1:nlami)
#             {
#                 devmat[i,j] <- mean(brier(surv[[j]],y[which,],y[!which,],cenFun,cenTime))
#                 #print(devmat[i,j])
#                 #print(surv[[j]][6,100])
#             }
#         }
#         nlams[i] <- nlami
#     }
# 
#     #print(nfolds)
#     ###If auc we behave differently
#         ##extract weights and normalize to sum to 1
#         #ywt=apply(y,1,sum)
#         #y=y/ywt
#         #weights=weights*ywt
#         cvraw=devmat
#     #   if(grouped){
#     #        cvob=cvcompute(cvraw,rep(1, nrow(y)),foldid,nlams)
#     #        cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
#     #    }
#     #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
#     cvm  <- apply(cvraw, 2, mean)
#     #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=rep(1, nrow(y)),na.rm=TRUE)/(N-1))
#     cvsd <- apply(cvraw, 2, sd)
#     out <- list(cvm  = cvm,
#                 cvsd = cvsd,
#                 name = typenames[type.measure])
#     if(keep)out$fit.preval <- devmat
#     out
# 
# }

# cv.vennbinomial <- function(outlist,
#                             lambda,
#                             lambda.fused,
#                             x,
#                             y,
#                             groups,
#                             foldid,
#                             type.measure,
#                             grouped,
#                             keep = FALSE)
# {
#   typenames  <- c(mse = "Mean-Squared Error",
#                   mae = "Mean Absolute Error",
#                   deviance = "Binomial Deviance",
#                   auc = "AUC",
#                   class = "Misclassification Error")
#   if(type.measure == "default") type.measure="deviance"
#   if(!match(type.measure, c("mse", "mae", "deviance", "auc", "class"), FALSE))
#   {
#     warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
#     type.measure <- "deviance"
#   }
#   
#   print("lam")
#   print(lambda)
#   print("lamfused")
#   print(lambda.fused)
# 
#   ###These are hard coded in the Fortran, so we do that here too
#   prob_min <- 1e-5
#   prob_max <- 1-prob_min
#   ###Turn y into a matrix
#   nc <- dim(y)
#   if (is.null(nc)) 
#   {
#     y    <- as.factor(y)
#     ntab <- table(y)
#     nc   <- as.integer(length(ntab))
#     y    <- diag(nc)[as.numeric(y), ]
#   }
#   N      <- nrow(y)
#   nfolds <- max(foldid)
#   if( (N/nfolds <10) && type.measure == "auc")
#   {
#     warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
#     type.measure <- "deviance"
#   }
#   if( (N/nfolds <3)&&grouped)
#   {
#     warning("Option grouped=FALSE enforced in cv.groupFusedMulti, since < 3 observations per fold",call.=FALSE)
#     grouped <- FALSE
#   }
# 
# 
#   #if(!is.null(offset)){
#   #  is.offset=TRUE
#   #  offset=drop(offset)
#   #}else is.offset=FALSE
#   if (!is.null(lambda.fused))
#   {
#     predmat <- vector(mode = "list", length = length(lambda.fused))
#     predlist <- vector(mode = "list", length = nfolds)
#     
#     for(i in seq(nfolds))
#     {
#       which    <- foldid == i
#       fitobj   <- outlist[[i]]
#       #fitobj$offset=FALSE
#       predlist[[i]]    <- predict(fitobj, x[which,,drop=FALSE], 
#                                   group.mat = groups[which,,drop=FALSE], type = "response")
#     }
#     
#     for (lf in 1:length(lambda.fused))
#     {
#       if (is.matrix(lambda))
#       {
#         predmat.tmp <- matrix(NA,NROW(y),NROW(lambda))
#       } else
#       {
#         predmat.tmp <- matrix(NA,NROW(y),length(lambda))
#       }
#       nlams   <- double(nfolds)
#       for(i in seq(nfolds))
#       {
#         which    <- foldid == i
#         if (is.matrix(outlist[[i]]$lambda))
#         {
#           nlami    <- NROW(outlist[[i]]$lambda)
#         } else
#         {
#           nlami    <- length(outlist[[i]]$lambda)
#         }
#         predmat.tmp[which, seq(nlami)] <- predlist[[i]][[lf]]
#         nlams[i] <- nlami
#       }
#       predmat[[lf]] <- predmat.tmp
#     }
#   } else
#   {
#     if (is.matrix(lambda))
#     {
#       predmat <- matrix(NA, NROW(y), NROW(lambda))
#     } else
#     {
#       predmat <- matrix(NA, NROW(y), length(lambda))
#     }
#     nlams   <- double(nfolds)
#     for(i in seq(nfolds))
#     {
#       which    <- foldid == i
#       fitobj   <- outlist[[i]]
#       #if(is.offset)off_sub=offset[which]
#       preds    <- predict(fitobj, newx = x[which,,drop=FALSE], group.mat = groups[which,,drop=FALSE], type="response")
#       if (is.matrix(outlist[[i]]$lambda))
#       {
#         nlami    <- length(outlist[[i]]$lambda)
#       } else
#       {
#         nlami    <- length(outlist[[i]]$lambda)
#       }
#       predmat[which,seq(nlami)] <- preds
#       nlams[i] <- nlami
#     }
#   }
#   
#   if (is.list(predmat))
#   {
#     cvraw <- cvm <- cvsd <- N <- cvob <- good <- weights <- vector(mode = "list", length = length(predmat))
#     
#     if (type.measure == "auc")
#     {
#       for (lf in 1:length(predmat))
#       {
#         if (is.matrix(lambda))
#         {
#           cvraw[[lf]] <- matrix(NA, nfolds, NROW(lambda))
#           good[[lf]]  <- matrix(0, nfolds,  NROW(lambda))
#         } else
#         {
#           cvraw[[lf]] <- matrix(NA, nfolds, length(lambda))
#           good[[lf]]  <- matrix(0, nfolds,  length(lambda))
#         }
#         for(i in seq(nfolds))
#         {
#           good[[lf]][i, seq(nlams[i])] <- 1
#           which <- foldid == i
#           for(j in seq(nlams[i]))
#           {
#             #cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
#             cvraw[[lf]][i,j] <- auc.mat(y[which,],predmat[[lf]][which,j], rep(1, sum(which)))
#           }
#         }
#         N[[lf]] <- apply(good[[lf]], 2, sum)
#       }
#     } else
#     {
#       for (lf in 1:length(predmat))
#       {
#         N[[lf]] <- length(y) - apply(is.na(predmat[[lf]]), 2, sum)
#         
#         
#         cvraw[[lf]] <- switch(type.measure,
#                               "mse"      = (y[,1] - (1 - predmat[[lf]])) ^ 2 +(y[,2] - predmat[[lf]]) ^ 2,
#                               "mae"      = abs(y[,1] - (1 - predmat[[lf]])) + abs(y[,2] - predmat[[lf]]),
#                               "deviance" = {
#                                 predmata=pmin(pmax(predmat[[lf]],prob_min),prob_max)
#                                 lp=y[,1]*log(1-predmata)+y[,2]*log(predmata)
#                                 ly=log(y)
#                                 ly[y==0]=0
#                                 ly=drop((y*ly)%*%c(1,1))
#                                 2*(ly-lp)
#                               },
#                               "class" = y[,1]*(predmat[[lf]] > .5) + y[,2] * (predmat[[lf]] <= .5)
#         )
#         if(grouped)
#         {
#           cvob[[lf]]    <- cvcompute(cvraw[[lf]], rep(1, NROW(y)), foldid, nlams)
#           cvraw[[lf]]   <- cvob[[lf]]$cvraw
#           weights[[lf]] <- cvob[[lf]]$weights
#           N[[lf]]       <- cvob[[lf]]$N
#         }
#       }
#     }
#     
#     for (lf in 1:length(predmat))
#     {
#       cvm[[lf]]  <- apply(cvraw[[lf]], 2, mean, w = weights[[lf]], na.rm = TRUE)
#       #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=rep(1, nrow(y)),na.rm=TRUE)/(N-1))
#       cvsd[[lf]] <- sqrt(apply(scale(cvraw[[lf]], cvm[[lf]], FALSE) ^ 2, 2, mean, na.rm = TRUE) / (N[[lf]] - 1))
#     }
#     
#   } else
#   {
#     
#     ###If auc we behave differently
#     if(type.measure == "auc")
#     {
#       if (is.matrix(lambda))
#       {
#         cvraw <- matrix(NA,nfolds, NROW(lambda))
#         good  <- matrix(0,nfolds,  NROW(lambda))
#       } else
#       {
#         cvraw <- matrix(NA,nfolds, length(lambda))
#         good  <- matrix(0,nfolds,  length(lambda))
#       }
#       for(i in seq(nfolds))
#       {
#         good[i, seq(nlams[i])] <- 1
#         which <- foldid == i
#         for(j in seq(nlams[i]))
#         {
#           #cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
#           cvraw[i,j] <- auc.mat(y[which,],predmat[which,j], rep(1, sum(which)))
#         }
#       }
#       N <- apply(good, 2, sum)
#       #weights=tapply(weights,foldid,sum)
#     } else
#     {
#       ##extract weights and normalize to sum to 1
#       #ywt=apply(y,1,sum)
#       #y=y/ywt
#       #weights=weights*ywt
#   
#       N     <- nrow(y) - apply(is.na(predmat), 2, sum)
#       cvraw <- switch(type.measure,
#                       "mse"     = (y[,1] - (1 - predmat)) ^ 2 +(y[,2] - predmat) ^ 2,
#                       "mae"     = abs(y[,1] - (1 - predmat)) + abs(y[,2] - predmat),
#                       "deviance"= {
#                           predmat=pmin(pmax(predmat,prob_min),prob_max)
#                           lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
#                           ly=log(y)
#                           ly[y==0]=0
#                           ly=drop((y*ly)%*%c(1,1))
#                           2*(ly-lp)
#                       },
#                       "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
#       )
#       if(grouped)
#       {
#         cvob    <- cvcompute(cvraw,rep(1, nrow(y)), foldid, nlams)
#         cvraw   <- cvob$cvraw
#         weights <- cvob$weights
#         N       <- cvob$N
#       }
#     }
#     #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
#     cvm  <- apply(cvraw, 2, mean, w = weights, na.rm = TRUE)
#     #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=rep(1, nrow(y)),na.rm=TRUE)/(N-1))
#     cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE) ^ 2, 2, mean, na.rm = TRUE) / (N - 1))
#     
#   }
#   
#   out  <- list(cvm  = cvm,
#                cvsd = cvsd,
#                name = typenames[type.measure])
#   if(keep)out$fit.preval <- predmat
#   out
# 
# }

cv.gfmgaussian <- function(outlist,lambda,lambda.fused,x,y,groups,foldid,type.measure,grouped,keep=FALSE){
  typenames <- c(deviance = "Mean-Squared Error",
                 mse      = "Mean-Squared Error",
                 mae      = "Mean Absolute Error")
  if(type.measure == "default") type.measure <- "mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE))
  {
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure <- "mse"
  }
  #if(!is.null(offset))y=y-drop(offset)
  
  nfolds  <- max(foldid)
  
  noutcomes <- NCOL(y)

  if (!is.null(lambda.fused))
  {
    predmat <- vector(mode = "list", length = length(lambda.fused))
    predlist <- vector(mode = "list", length = nfolds)
    
    for(i in seq(nfolds))
    {
      which    <- foldid == i
      fitobj   <- outlist[[i]]
      #fitobj$offset=FALSE
      predlist[[i]]    <- predict(fitobj, x[which,,drop=FALSE], type = "response")
    }
    
    for (lf in 1:length(lambda.fused))
    {
      if (is.matrix(lambda))
      {
        predmat.tmp <- array(NA, dim = c(NROW(y), NCOL(y), NROW(lambda)))
      } else
      {
        predmat.tmp <- array(NA, dim = c(NROW(y), NCOL(y), length(lambda)))
      }
      nlams   <- double(nfolds)
      for(i in seq(nfolds))
      {
        which    <- foldid == i
        if (is.matrix(outlist[[i]]$lambda))
        {
          nlami    <- NROW(outlist[[i]]$lambda)
        } else
        {
          nlami    <- length(outlist[[i]]$lambda)
        }
        predmat.tmp[which, , seq(nlami)] <- predlist[[i]][[lf]]
        nlams[i] <- nlami
      }
      predmat[[lf]] <- predmat.tmp
    }
  } else
  {
    if (is.matrix(lambda))
    {
      predmat <- array(NA, dim = c(NROW(y), NCOL(y), NROW(lambda)))
    } else
    {
      predmat <- array(NA, dim = c(NROW(y), NCOL(y), length(lambda)))
    }
    nfolds  <- max(foldid)
    nlams   <- double(nfolds)
    for(i in seq(nfolds))
    {
      which    <- foldid == i
      fitobj   <- outlist[[i]]
      #fitobj$offset=FALSE
      preds    <- predict(fitobj, x[which,,drop=FALSE], type = "response")
      if (is.matrix(outlist[[i]]$lambda))
      {
        nlami    <- NROW(outlist[[i]]$lambda)
      } else
      {
        nlami    <- length(outlist[[i]]$lambda)
      }
      predmat[which, , seq(nlami)] <- preds
      nlams[i] <- nlami
    }
  }

  
  if (is.list(predmat))
  {
    cvraw <- cvm <- cvsd <- N <- cvob <- vector(mode = "list", length = length(predmat))
    
    for (lf in 1:length(predmat))
    {
      
      N[[lf]] <- NROW(y)# - apply(is.na(predmat[[lf]]), 2:3, sum)
      cvraw[[lf]] <- predmat[[lf]]
      
      for (l in 1:length(lambda))
      {
        cvraw[[lf]][,,l] <- switch(type.measure,
                              "mse"      =(y - predmat[[lf]][,,l]) ^ 2,
                              "deviance" =(y - predmat[[lf]][,,l]) ^ 2,
                              "mae"      =abs(y - predmat[[lf]][,,l])
        )
      }
      
      if( (NROW(y) / nfolds <3) && grouped)
      {
        warning("Option grouped=FALSE enforced in cv.groupFusedMulti, since < 3 observations per fold",call.=FALSE)
        grouped <- FALSE
      }
      if(grouped)
      {
        cvraw_tmp <- array(NA, dim = c(nfolds, noutcomes, length(lambda)))
        for (o in 1:noutcomes)
        {
          cvob_tmp    <- cvcompute(cvraw[[lf]][,o,], rep(1, NROW(y)), foldid, nlams)
          cvraw_tmp[,o,] <- cvob_tmp$cvraw
          #weights    <- cvob[[lf]]$weights
          N[[lf]]     <- cvob_tmp$N
        }
        cvraw[[lf]] <- cvraw_tmp
      }
      
      
      #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
      cvm[[lf]]  <- apply(cvraw[[lf]], 3, mean, na.rm = TRUE)
      #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
      
      cvscaled <- cvraw[[lf]]
      for (o in 1:noutcomes)
      {
        cvscaled[,o,] <- scale(cvraw[[lf]][,o,], cvm[[lf]], FALSE) ^ 2
      }
      
      cvsd[[lf]] <- sqrt(apply(cvscaled, 3, mean, na.rm = TRUE)/(N[[lf]] - 1))
    }
    
  } else
  {
    N <- NROW(y)# - apply(is.na(predmat), 2, sum)
    
    cvraw <- predmat
    for (l in 1:length(lambda))
    {
      cvraw[,,l] <- switch(type.measure,
                                 "mse"      =(y - predmat[,,l]) ^ 2,
                                 "deviance" =(y - predmat[,,l]) ^ 2,
                                 "mae"      =abs(y - predmat[,,l])
      )
    }
    
    if( (NROW(y)/nfolds <3) && grouped)
    {
      warning("Option grouped=FALSE enforced in cv.groupFusedMulti, since < 3 observations per fold",call.=FALSE)
      grouped <- FALSE
    }
    if(grouped)
    {
      for (o in 1:noutcomes)
      {
        cvob_tmp    <- cvcompute(cvraw[,o,], rep(1, NROW(y)), foldid, nlams)
        cvraw[,o,]  <- cvob_tmp$cvraw
        #weights       <- cvob[[lf]]$weights
        N       <- cvob_tmp$N
      }
    }
    
    #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
    cvm  <- apply(cvraw, 3, mean, na.rm = TRUE)
    #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE) ^ 2, 3, mean, na.rm = TRUE)/(N-1))
  }
  
  
  out  <- list(cvm  = cvm,
               cvsd = cvsd,
               cvraw = cvraw,
               name = typenames[type.measure])
  if(keep) out$fit.preval <- predmat
  out
}



