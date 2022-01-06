
#' Prediction for Group Fused Multi Lasso
#'
#' @param object fitted groupFusedMulti object
#' @param newx new matrix for predictions
#' @param s lambda value for the predictions. defaults to all values computed in the groupFusedMulti object
#' @param sf fused lasso tuning parameter value for the predictions. defaults to all values computed in the groupFusedMulti object
#' @param type type of predictions to be made.
#' @param ... parameters to be passed to groupFusedMulti
#' @return predictions or coefficients
#'
#'
#' @import Rcpp
#' @import methods
#' @import Matrix
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom stats glm
#' @importFrom stats coef
#' @importFrom stats binomial
#'
#' @rdname predict
#' @export
predict.groupFusedMulti <- function(object, newx = NULL,
                                    s = NULL, sf = NULL,
                                    type = c("link", "response", "coefficients", "nonzero", "class", "nvars"),
                                    ...)
{
    type=match.arg(type)
    if(is.null(newx)){
        if(!match(type,c("coefficients","nonzero","nvars"),FALSE))stop("You must supply a value for 'newx'")
    }
    
    #if(exact&&(!is.null(s))){
    #  ###we augment the lambda sequence with the new values, if they are different,and refit the model using update
    #  lambda=object$lambda
    #  which=match(s,lambda,FALSE)
    #  if(!all(which>0)){
    #    lambda=unique(rev(sort(c(s,lambda))))
    #    object=update(object,lambda=lambda)
    #  }
    #}
    
    
    n.obs <- NROW(newx)
    noutcomes <- object$noutcomes
    
    if (!is.null(s) & !is.matrix(s) & is.matrix(object$lambda))
    {
        if (is.null(sf))
        {
            s <- matrix(rep(s, length(object$lambda.fused)), ncol = length(object$lambda.fused))
        } else
        {
            s <- matrix(rep(s, length(sf)), ncol = length(sf))
        }
    }
    
    is_fused <- object$is_fused

    
    ###Check on newx
    if(inherits(newx, "sparseMatrix")) newx = as(newx,"dgCMatrix")
    
    
    if (is_fused)
    {
        nlamfused <- object$nlambda.fused
        a0 <- nbeta <- vector(mode = "list", length = nlamfused)

        # for (lf in 1:nlamfused)
        # {
        #     a0[[lf]] = t(as.matrix(object$intercept[[lf]]))
        #     rownames(a0[[lf]]) = "(Intercept)"
        # }
        
        for (lf in 1:nlamfused)
        {
            nbeta[[lf]] <- object$beta[[lf]]
            dimnames(nbeta[[lf]]) <- list(NULL, NULL, NULL)
        }
        
        if(!is.null(s))
        {
            if (is.null(sf))
            {
                lambda  = object$lambda
                if (!is.matrix(object$lambda))
                {
                    lamlist = lambdaInterp(lambda, s)
                }
                for (lf in 1:nlamfused)
                {
                    if (is.matrix(object$lambda))
                    {
                        lamlist = lambdaInterp(lambda[,lf], s[,lf])
                    }
                    nbeta[[lf]]   = nbeta[[lf]][,,lamlist$left,drop=FALSE] * lamlist$frac + nbeta[[lf]][,,lamlist$right,drop=FALSE] * (1-lamlist$frac)
                    #a0[[lf]]      = a0[[lf]][,lamlist$left, drop=FALSE]    * lamlist$frac + a0[[lf]][,lamlist$right, drop=FALSE]    * (1-lamlist$frac)
                    
                    if (is.matrix(object$lambda))
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s[,lf])))
                        } else
                        {
                            dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda[,lf])))
                        }
                    } else
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s)))
                        } else
                        {
                            dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda)))
                        }
                    }
                }
                if (is.matrix(s))
                {
                    nlam <- NROW(s)
                } else
                {
                    nlam <- length(s)
                }
                
            } else
            {
                ## here we need to do bilinear interpolation (which is the same as 
                ## linear interpolation on each lambda separately)
                lambda       <- object$lambda
                lambda.fused <- object$lambda.fused
                
                
                lamfusedlist <- lambdaInterp(lambda.fused, sf)
                
                nbeta.out <- a0.out <- vector(mode = "list", length = length(sf))
                
                for (f in 1:length(sf))
                {
                    
                    lam.fused.idx.nearest <- which.min(abs(lambda.fused - sf[f]))
                    
                    if (is.matrix(lambda))
                    {
                        lambda      <- lambda[,lam.fused.idx.nearest]
                        lamlist     <- lambdaInterp(lambda, s[,f])
                    } else
                    {
                        lamlist     <- lambdaInterp(lambda, s)
                    }
                    
                    
                    
                    
                    nbeta.left  <- nbeta[[lamfusedlist$left[f]]]
                    nbeta.right <- nbeta[[lamfusedlist$right[f]]]
                    
                    nbeta.interp.lamfused <- nbeta.left * lamfusedlist$frac[f] + nbeta.right * (1 - lamfusedlist$frac[f])
                    #a0.interp.lamfused <- a0[[lamfusedlist$left[f]]] * lamfusedlist$frac[f] + 
                    #    a0[[lamfusedlist$right[f]]] * (1 - lamfusedlist$frac[f])
                    
                    nbeta.out[[f]]   <- nbeta.interp.lamfused[,,lamlist$left,drop=FALSE] * lamlist$frac + 
                        nbeta.interp.lamfused[,,lamlist$right,drop=FALSE] * (1 - lamlist$frac)
                    #a0.out[[f]]      <- a0.interp.lamfused[,lamlist$left, drop=FALSE] * lamlist$frac + 
                    #    a0.interp.lamfused[,lamlist$right, drop=FALSE] * (1 - lamlist$frac)
                    
                    if (is.matrix(object$lambda))
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s[,f])))
                        } else
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda[,f])))
                        }
                    } else
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s)))
                        } else
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda)))
                        }
                    }
                }
                
                nbeta <- nbeta.out
                a0    <- a0.out
                
                if (is.matrix(s))
                {
                    nlam <- NROW(s)
                } else
                {
                    nlam <- length(s)
                }
            }
            
        } else 
        {
            
            if (!is.null(sf))
            {
                lambda.fused <- object$lambda.fused
                lamfusedlist <- lambdaInterp(lambda.fused, sf)
                
                nbeta.out <- a0.out <- vector(mode = "list", length = length(sf))
                
                for (f in 1:length(sf))
                {
                    nbeta.left  <- nbeta[[lamfusedlist$left[f]]]
                    nbeta.right <- nbeta[[lamfusedlist$right[f]]]
                    
                    nbeta.out[[f]] <- nbeta.left * lamfusedlist$frac[f] + nbeta.right * (1 - lamfusedlist$frac[f])
                    #a0.out[[f]] <- a0[[lamfusedlist$left[f]]] * lamfusedlist$frac[f] + 
                    #    a0[[lamfusedlist$right[f]]] * (1 - lamfusedlist$frac[f])
                    
                    if (is.matrix(object$lambda))
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s[,f])))
                        } else
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda[,f])))
                        }
                    } else
                    {
                        if (!is.null(s))
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s)))
                        } else
                        {
                            dimnames(nbeta.out[[f]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lambda)))
                        }
                    }
                }
                
                nbeta <- nbeta.out
                a0    <- a0.out
            } else
            {
                for (lf in 1:nlamfused)
                {
                    if (!is.null(s))
                    {
                        dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = s)))
                    } else
                    {
                        dimnames(nbeta[[lf]]) <- list(object$var.names, object$outcome.names, paste(seq(along = object$lam)))
                    }
                }
            }
            
            if (is.matrix(object$lambda))
            {
                nlam <- NROW(object$lambda)
            } else
            {
                nlam <- length(object$lambda)
            }
        }
        
    } else ## not fused, below:
    {

        
        nbeta <- object$beta
        
        dimnames(nbeta) = list(NULL, NULL, NULL)
        if(!is.null(s))
        {
            lambda  = object$lambda
            lamlist = lambdaInterp(lambda,s)
            nbeta   = nbeta[,,lamlist$left,drop=FALSE] * lamlist$frac + nbeta[,,lamlist$right,drop=FALSE] * (1-lamlist$frac)
            #a0      = a0[,lamlist$left, drop=FALSE] * lamlist$frac + a0[,lamlist$right, drop=FALSE] * (1-lamlist$frac)
            dimnames(nbeta) = list(object$var.names, object$outcome.names, paste(seq(along = s)))
            nlam <- length(s)
        } else 
        {
            dimnames(nbeta) = list(object$var.names, object$outcome.names, paste(seq(along = object$lam)))
            nlam <- length(object$lambda)
        }
        
    }
    
    if (type == "coefficients")
    {
        return(drop(nbeta))
    }
    if (type == "nonzero")
    {
        stop("nonzero not properly working yet")
        nonzeroCoeff(nbeta[,,,drop=FALSE],bystep = TRUE)
        #colnames(nbeta)[1] <- "(Intercept)"
        return(nbeta)
    }
    
    ## and NA variables are due to missingness in a strata. set to zero
    if (is.list(nbeta))
    {
        for (lf in 1:length(nbeta))
        {
            nbeta[[lf]][is.na(nbeta[[lf]])] <- 0
        }
    } else
    {
        nbeta[is.na(nbeta)] <- 0
    }
    
    if (type == "nvars")
    {
        ## returns number of nonzero coefficients
        ## for each strata for each lambda
        
        if (is.list(nbeta))
        {
            nvarssel <- vector(mode = "list", length = length(nbeta))
            for (lf in 1:length(nbeta))
            {
                nvarssel[[lf]] <- apply(nbeta[[lf]], c(1, 3), function(xx) sum(xx != 0))
            }
            
            return(nvarssel)
        } else
        {
            return(apply(nbeta, c(1, 3), function(xx) sum(xx != 0)))
        }
    }
    
    
    
    
    if (is.list(nbeta))
    {
        nfit <- vector(mode = "list", length = length(nbeta))
        
        for (lf in 1:length(nbeta))
        {
            nfit.tmp <- array(NA, dim = c(n.obs, noutcomes, nlam))
   
            for (l in 1:nlam)
            {
                nfit.tmp[,,l] <- cbind(1, newx) %*% nbeta[[lf]][,,l]
            }
                
            nfit[[lf]] <- nfit.tmp
        }
    } else
    {
        nfit <- array(NA, dim = c(n.obs, noutcomes, nlam))
        
        for (l in 1:nlam)
        {
            nfit.tmp[,,l] <- cbind(1, newx) %*% nbeta[,,l]
        }
    }
    
    
    
    #if(object$offset){
    #  if(missing(offset))stop("No offset provided for prediction, yet used in fit of hierSharedLasso",call.=FALSE)
    #  if(is.matrix(offset)&&dim(offset)[[2]]==2)offset=offset[,2]
    #  nfit=nfit+array(offset,dim=dim(nfit))
    #}
    
    
    if (is.list(nfit))
    {
        if (object$family == "gaussian") 
        {
            return(nfit)
        }
        
        ret_list <- vector(mode = "list", length = length(nfit))
        
        if (object$family == "binomial") 
        {
            for (lf in 1:length(nfit))
            {
                ret_list[[lf]] <- switch(type,
                                         response={
                                             pp=exp(-nfit[[lf]])
                                             1/(1+pp)
                                         },
                                         class={
                                             cnum = ifelse(nfit[[lf]] > 0, 2, 1)
                                             clet = object$classnames[cnum]
                                             if(is.matrix(cnum)) clet = array(clet,dim(cnum),dimnames(cnum))
                                             clet
                                         },
                                         nfit[[lf]])
            }
            return(ret_list)
        }
        
    } else
    {
        if (object$family == "gaussian") 
        {
            return(nfit)
        } else if (object$family == "binomial") 
        {
            return(switch(type,
                          response={
                              pp=exp(-nfit)
                              1/(1+pp)
                          },
                          class={
                              cnum=ifelse(nfit>0,2,1)
                              clet=object$classnames[cnum]
                              if(is.matrix(cnum))clet=array(clet,dim(cnum),dimnames(cnum))
                              clet
                          },
                          nfit
            )
            )
        }
    }
}

#' Prediction for Cross Validation Hierarchical Lasso Object
#'
#' @param object fitted cv.groupFusedMulti object
#' @param newx new matrix for predictions
#' @param s lambda value for the predictions. defaults to minimum 
#' @param ... parameters to be passed to predict.groupFusedMulti
#' @return predictions or coefficients
#'
#'
#' @import Rcpp
#' @method predict cv.groupFusedMulti
#' @export
predict.cv.groupFusedMulti <- function(object, newx = NULL, s = c("lambda.min"), ...)
{
    if(is.numeric(s)) lambda <- s
    else
        if(is.character(s))
        {
            s      <- match.arg(s)
            lambda <- object[[s]]
            lambda.fused <- object[["lambda.fused.min"]]
        }
    else stop("Invalid form for s")
    if (object$groupFusedMulti.fit$is_fused)
    {
        return(    drop(predict(object$groupFusedMulti.fit, newx, 
                                s = lambda, 
                                sf = lambda.fused, ...)[[1]]))
    } else
    {
        return(    drop(predict(object$groupFusedMulti.fit, newx, 
                                s = lambda, ...)))
    }
}
