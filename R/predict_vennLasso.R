
#' Prediction for Hierarchical Shared Lasso
#'
#' @param object fitted vennLasso object
#' @param newx new matrix for predictions
#' @param group.mat A matrix of the group memberships for now. Ignore the rest:
#' A list of length equal to the number of groups containing vectors of integers
#' indicating the variable IDs for each group. For example, groups=list(c(1,2), c(2,3), c(3,4,5)) specifies
#' that Group 1 contains variables 1 and 2, Group 2 contains variables 2 and 3, and Group 3 contains
#' variables 3, 4, and 5. Can also be a matrix of 0s and 1s with the number of columns equal to the
#' number of groups and the number of rows equal to the number of variables. A value of 1 in row i and
#' column j indicates that variable i is in group j and 0 indicates that variable i is not in group j.
#' @param s lambda value for the predictions. defaults to all values computed in the vennLasso object
#' @param sf fused lasso tuning parameter value for the predictions. defaults to all values computed in the vennLasso object
#' @param use.refit Should the refitted beta estimates be used for prediction? Defaults to FALSE. If TRUE
#' then the beta estimates from the model refit on just the selected covariates are used
#' @param type type of predictions to be made. \code{type = "median"} is for the median survival time and 
#' \code{type = "survival"} is for the predicted hazard function
#' @param ... parameters to be passed to vennLasso
#' @return predictions or coefficients
#'
#'
#' @import Rcpp
#' @import methods
#' @import Matrix
#' @import survival
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom stats glm
#' @importFrom stats coef
#' @importFrom stats binomial
#'
#' @rdname predict
#' @export
predict.vennLasso <- function(object, newx,
                              group.mat, s = NULL, sf = NULL,
                              use.refit = FALSE,
                              type = c("link", "response", "coefficients", "nonzero", "class", "nvars", "median", "survival"),
                              ...)
{
  type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero","nvars"),FALSE))stop("You must supply a value for 'newx'")
  }
  if( (type == "median" || type == "survival") & object$family != "coxph") stop("Median option is only available for cox model")
  #if(exact&&(!is.null(s))){
  #  ###we augment the lambda sequence with the new values, if they are different,and refit the model using update
  #  lambda=object$lambda
  #  which=match(s,lambda,FALSE)
  #  if(!all(which>0)){
  #    lambda=unique(rev(sort(c(s,lambda))))
  #    object=update(object,lambda=lambda)
  #  }
  #}
  
  
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
  
  if (use.refit & is_fused)
  {
    stop("refit not ready for fused vennLasso")
  }

  if (all(is.na(object$beta.refit)) & use.refit)
  {
      use.refit <- FALSE
      warning("no refitted estimates were computed for this model, defaulting to non-refit coefficients")
  }
  
  
  # only compute if newx needed
  if(!match(type,c("coefficients","nonzero","nvars"),FALSE))
  {
      n.obs <- nrow(newx)
      stopifnot(nrow(group.mat) == n.obs)
      M <- object$n.combinations
      combin.mat <- object$condition.combinations

      data.indices <- vector(mode = "list", length = M)
      if (!is.null(n.obs))
      {
          group.vec <- apply(group.mat, 1, FUN = function(x) {paste(x,collapse = "")})
      } else 
      {
          group.vec <- paste(group.mat,collapse = "")
          n.obs <- 1
      }

      ind.2.remove <- NULL
      for (c in 1:M) 
      {
        data.indices[[c]] <- which(group.vec == paste(combin.mat[c,], collapse = ""))
        if (length(data.indices[[c]]) == 0) 
        {
          ind.2.remove <- c(c, ind.2.remove)
        }
      }
      Mnew <- M
      # remove combinations that have no observations
      if (!is.null(ind.2.remove)) 
      {
        data.indices[ind.2.remove] <- NULL
        Mnew <- length(data.indices)
        #combin.mat <- combin.mat[-ind.2.remove,]
      }

      ###Check on newx
      if(inherits(newx, "sparseMatrix")) newx = as(newx,"dgCMatrix")

  }



  if (is_fused)
  {
    nlamfused <- object$nlambda.fused
    a0 <- nbeta <- vector(mode = "list", length = nlamfused)
    if (use.refit)
    {
      for (lf in 1:nlamfused)
      {
        a0[[lf]] = t(as.matrix(object$intercept.refit[[lf]]))
        rownames(a0[[lf]]) = "(Intercept)"
      }
    } else
    {
      for (lf in 1:nlamfused)
      {
        a0[[lf]] = t(as.matrix(object$intercept[[lf]]))
        rownames(a0[[lf]]) = "(Intercept)"
      }
    }
    
    for (lf in 1:nlamfused)
    {
      if (use.refit)
      {
        nbeta[[lf]] <- object$beta.refit[[lf]]
      } else
      {
        nbeta[[lf]] <- object$beta[[lf]]
      }
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
          a0[[lf]]      = a0[[lf]][,lamlist$left, drop=FALSE]    * lamlist$frac + a0[[lf]][,lamlist$right, drop=FALSE]    * (1-lamlist$frac)
          
          if (is.matrix(object$lambda))
          {
            if (!is.null(s))
            {
              dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = s[,lf])))
            } else
            {
              dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda[,lf])))
            }
          } else
          {
            if (!is.null(s))
            {
              dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = s)))
            } else
            {
              dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda)))
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
            lambda <- lambda[,lam.fused.idx.nearest]
          }
          
          lamlist     <- lambdaInterp(lambda, s[,f])
          
          
          nbeta.left  <- nbeta[[lamfusedlist$left[f]]]
          nbeta.right <- nbeta[[lamfusedlist$right[f]]]
          
          nbeta.interp.lamfused <- nbeta.left * lamfusedlist$frac[f] + nbeta.right * (1 - lamfusedlist$frac[f])
          a0.interp.lamfused <- a0[[lamfusedlist$left[f]]] * lamfusedlist$frac[f] + 
            a0[[lamfusedlist$right[f]]] * (1 - lamfusedlist$frac[f])
          
          nbeta.out[[f]]   <- nbeta.interp.lamfused[,,lamlist$left,drop=FALSE] * lamlist$frac + 
            nbeta.interp.lamfused[,,lamlist$right,drop=FALSE] * (1 - lamlist$frac)
          a0.out[[f]]      <- a0.interp.lamfused[,lamlist$left, drop=FALSE] * lamlist$frac + 
            a0.interp.lamfused[,lamlist$right, drop=FALSE] * (1 - lamlist$frac)
          
          if (is.matrix(object$lambda))
          {
            if (!is.null(s))
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = s[,f])))
            } else
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda[,f])))
            }
          } else
          {
            if (!is.null(s))
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = s)))
            } else
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda)))
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
          a0.out[[f]] <- a0[[lamfusedlist$left[f]]] * lamfusedlist$frac[f] + 
            a0[[lamfusedlist$right[f]]] * (1 - lamfusedlist$frac[f])
          
          if (is.matrix(object$lambda))
          {
            if (!is.null(s))
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = s[,f])))
            } else
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda[,f])))
            }
          } else
          {
            if (!is.null(s))
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = s)))
            } else
            {
              dimnames(nbeta.out[[f]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lambda)))
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
            dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = s)))
          } else
          {
            dimnames(nbeta[[lf]]) <- list(object$combin.names, object$var.names, paste(seq(along = object$lam)))
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
    if (use.refit)
    {
      a0 = t(as.matrix(object$intercept.refit))
    } else
    {
      if (is_fused)
      {
        a0 = t(as.matrix(object$intercept))
      } else
      {
        a0 = t(as.matrix(object$intercept))
      }
    }
    
    rownames(a0) = "(Intercept)"
    
    if (use.refit)
    {
      nbeta <- object$beta.refit
    } else
    {
      nbeta <- object$beta
    }
    
    dimnames(nbeta) = list(NULL, NULL, NULL)
    if(!is.null(s))
    {
      lambda  = object$lambda
      lamlist = lambdaInterp(lambda,s)
      nbeta   = nbeta[,,lamlist$left,drop=FALSE] * lamlist$frac + nbeta[,,lamlist$right,drop=FALSE] * (1-lamlist$frac)
      a0      = a0[,lamlist$left, drop=FALSE] * lamlist$frac + a0[,lamlist$right, drop=FALSE] * (1-lamlist$frac)
      dimnames(nbeta) = list(object$combin.names, object$var.names, paste(seq(along = s)))
      nlam <- length(s)
    } else 
    {
      dimnames(nbeta) = list(object$combin.names, object$var.names, paste(seq(along = object$lam)))
      nlam <- length(object$lambda)
    }
    
  }
  
  
  if (type == "coefficients")
  {
      return(nbeta)
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
      nfit.tmp <- matrix(NA, ncol = nlam, nrow = n.obs)
      for (c in 1:Mnew)
      {
        # make sure observations with each combination of
        # conditions are predicted with the correct coefficients.
        # data.indices cbind(rep(a0, nrow(object$beta)), object$beta)
        
        
        if (!object$one.intercept & (object$family != "coxph"))
        {
          if (!is.null(data.indices[[c]]) & length(data.indices[[c]]) > 0)
          {
            if (n.obs == 1)
            {
              nfit.tmp[data.indices[[c]],] <- c(1, newx) %*%
                matrix(nbeta[[lf]][c,,],ncol=dim(nbeta[[lf]])[3])
            } else
            {
              nfit.tmp[data.indices[[c]],] <- cbind2(1, newx[data.indices[[c]],,drop=FALSE]) %*%
                matrix(nbeta[[lf]][c,,],ncol=dim(nbeta[[lf]])[3])
            }
          }
        } else
        {
          if (!is.null(data.indices[[c]]) & length(data.indices[[c]]) > 0)
          {
            if (n.obs == 1)
            {
              nfit.tmp[data.indices[[c]],] <- c(1, newx) %*%
                rbind(as.vector(a0[[lf]]), matrix(nbeta[[lf]][c,,],ncol=dim(nbeta[[lf]])[3]))
            } else
            {
              nfit.tmp[data.indices[[c]],] <- cbind2(1, newx[data.indices[[c]],,drop=FALSE]) %*%
                rbind(as.vector(a0[[lf]]), matrix(nbeta[[lf]][c,,],ncol=dim(nbeta[[lf]])[3]))
            }
          }
        }
      }
      nfit[[lf]] <- nfit.tmp
    }
  } else
  {
    nfit <- matrix(NA, ncol = nlam, nrow = n.obs)
    for (c in 1:Mnew)
    {
      # make sure observations with each combination of
      # conditions are predicted with the correct coefficients.
      # data.indices cbind(rep(a0, nrow(object$beta)), object$beta)
      
      
      if (!is.null(data.indices[[c]]) & length(data.indices[[c]]) > 0)
      {
        if (!object$one.intercept & (object$family != "coxph"))
        {
          if (n.obs == 1)
          {
            nfit[data.indices[[c]],] <- c(1, newx) %*%
              matrix(nbeta[c,,],ncol=dim(nbeta)[3])
          } else
          {
            nfit[data.indices[[c]],] <- cbind2(1, newx[data.indices[[c]],,drop=FALSE]) %*%
              matrix(nbeta[c,,],ncol=dim(nbeta)[3])
          }
        } else
        {
          if (n.obs == 1)
          {
            nfit[data.indices[[c]],] <- c(1, newx) %*%
              rbind(as.vector(a0), matrix(nbeta[c,,],ncol=dim(nbeta)[3]))
          } else
          {
            nfit[data.indices[[c]],] <- cbind2(1, newx[data.indices[[c]],,drop=FALSE]) %*%
              rbind(as.vector(a0), matrix(nbeta[c,,],ncol=dim(nbeta)[3]))
          }
        }
      } ## end check that indices have any values
      
      
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
    } else if (object$family == "coxph") 
    {
      
      for (lf in 1:length(nfit))
      {
        eta <- apply(nfit[[lf]], 2, function(t){ pmin(t,50) })
        if (type =='link') return(drop(eta))
        if (type =='response') return(drop(exp(eta)))
        
        if (!is.null(s)) 
        {
          if (is.null(sf))
          {
            W <- object$W[[lf]][,lamlist$left,drop = FALSE] * lamlist$frac + 
              object$W[[lf]][,lamlist$right,drop = FALSE] * (1 - lamlist$frac)
          } else
          {
            W.left   <- object$W[[lamfusedlist$left[lf]]]
            W.right  <- object$W[[lamfusedlist$right[lf]]]
            W.interp <- W.left * lamfusedlist$frac[f] + W.right * (1 - lamfusedlist$frac[f])
            
            W <- lamlist$frac * W.interp[,lamlist$left,drop=FALSE] + 
              (1 - lamlist$frac) * W.interp[,lamlist$right,drop=FALSE]
          }
        } else 
        {
          if (is.null(sf))
          {
            W <- object$W[[lf]][,,drop = FALSE]
          } else
          {
            W.left   <- object$W[[lamfusedlist$left[lf]]]
            W.right  <- object$W[[lamfusedlist$right[lf]]]
            W        <- W.left * lamfusedlist$frac[f] + W.right * (1 - lamfusedlist$frac[f])
          }
        }
        #if (type == 'survival' & ncol(W) > 1) stop('Can only return type="survival" for a single lambda value')
        #if (type == 'survival') val <- vector('list', length(eta))
        if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
        hazard <- NULL
        if (type == 'survival' | type == 'median')
        {
          for (j in 1:ncol(eta)) 
          {
            # Estimate baseline hazard
            w <- W[,j]
            r <- rev(cumsum(rev(w)))
            base_hazard <- ifelse(object$status, 1/r, 0)
            tmp <- exp(-cumsum(base_hazard))
            #         temp <- (1-as.matrix(exp(eta[,j]))%*%t(as.matrix(base_hazard))) * ((1-as.matrix(exp(eta[,j]))%*%t(as.matrix(base_hazard))) > 0)
            #         hazard[[j]] <- t(apply(temp,1,cumprod))
            hazard[[j]] <- matrix(0, length(eta[,1]), length(tmp))
            for (i in 1:length(eta[,j]))
            {
              hazard[[j]][i,] <- tmp^(exp(eta[i,j]))
              if(sum(is.infinite(hazard[[j]]))>0) print("Find infinity")
              if(sum(is.na(hazard[[j]]))>0) print("Find NA")
              if(hazard[[j]][i,1]==0) print("eta is too big")
            }
            #print(w[1:4])
            x <- c(0, object$time)
            if (type == 'median') 
            {
              for (i in 1:nrow(eta)) 
              {
                #if (type == 'survival') val[[i]] <- approxfun(x, S, method='constant')
                
                if (any(hazard[[j]][i,] < 0.5)) 
                {
                  val[i,j] <- x[min(which(hazard[[j]][i,] < .5))]
                } else 
                {
                  val[i,j] <- max(x)
                }
              }
            }
          }
        }
        if (type == 'survival') {
          val <- hazard
        }
        if (type == 'median') val <- drop(val)
        ret_list[[lf]] <- val
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
    } else if (object$family == "coxph") 
    {
      eta <- apply(nfit,2,function(t){pmin(t,50)})
      if (type =='link') return(drop(eta))
      if (type =='response') return(drop(exp(eta)))
      
      if (!is.null(s)) {
        W <- lamlist$frac*object$W[,lamlist$left,drop=FALSE] + (1-lamlist$frac)*object$W[,lamlist$right,drop=FALSE]
      } else {
        W <- object$W[,,drop=FALSE]
      }
      #if (type == 'survival' & ncol(W) > 1) stop('Can only return type="survival" for a single lambda value')
      #if (type == 'survival') val <- vector('list', length(eta))
      if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
      hazard <- NULL
      if (type == 'survival' | type == 'median'){
        for (j in 1:ncol(eta)) {
          # Estimate baseline hazard
          w <- W[,j]
          r <- rev(cumsum(rev(w)))
          base_hazard <- ifelse(object$status, 1/r, 0)
          tmp <- exp(-cumsum(base_hazard))
          #         temp <- (1-as.matrix(exp(eta[,j]))%*%t(as.matrix(base_hazard))) * ((1-as.matrix(exp(eta[,j]))%*%t(as.matrix(base_hazard))) > 0)
          #         hazard[[j]] <- t(apply(temp,1,cumprod))
          hazard[[j]] <- matrix(0, length(eta[,1]), length(tmp))
          for (i in 1:length(eta[,j])){
            hazard[[j]][i,] <- tmp^(exp(eta[i,j]))
            if(sum(is.infinite(hazard[[j]]))>0) print("Find infinity")
            if(sum(is.na(hazard[[j]]))>0) print("Find NA")
            if(hazard[[j]][i,1]==0) print("eta is too big")
          }
          #print(w[1:4])
          x <- c(0, object$time)
          if (type == 'median') {
            for (i in 1:nrow(eta)) {
              #if (type == 'survival') val[[i]] <- approxfun(x, S, method='constant')
              
              if (any(hazard[[j]][i,] < 0.5)) {
                val[i,j] <- x[min(which(hazard[[j]][i,] < .5))]
              } else {
                val[i,j] <- max(x)
              }
            }
          }
        }
      }
      if (type == 'survival') {
        val <- hazard
      }
      if (type == 'median') val <- drop(val)
      val
    }
  }
}

#' Prediction for Cross Validation Hierarchical Lasso Object
#'
#' @param object fitted cv.vennLasso object
#' @param newx new matrix for predictions
#' @param group.mat A matrix of the group memberships for now. Ignore the rest:
#' A list of length equal to the number of groups containing vectors of integers
#' indicating the variable IDs for each group. For example, groups=list(c(1,2), c(2,3), c(3,4,5)) specifies
#' that Group 1 contains variables 1 and 2, Group 2 contains variables 2 and 3, and Group 3 contains
#' variables 3, 4, and 5. Can also be a matrix of 0s and 1s with the number of columns equal to the
#' number of groups and the number of rows equal to the number of variables. A value of 1 in row i and
#' column j indicates that variable i is in group j and 0 indicates that variable i is not in group j.
#' @param s lambda value for the predictions. defaults to all values computed in the vennLasso object
#' @param use.refit Should the refitted beta estimates be used for prediction? Defaults to FALSE. If TRUE
#' then the beta estimates from the model refit on just the selected covariates are used
#' @param ... parameters to be passed to predict.vennLasso
#' @return predictions or coefficients
#'
#'
#' @import Rcpp
#' @method predict cv.vennLasso
#' @export
predict.cv.vennLasso <- function(object, newx, group.mat, s = c("lambda.min"), use.refit = FALSE, ...)
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
    if (object$vennLasso.fit$is_fused)
    {
      return(    predict(object$vennLasso.fit, newx, 
                         group.mat, s = lambda, 
                         sf = lambda.fused,
                         use.refit = use.refit, ...)[[1]])
    } else
    {
      return(    predict(object$vennLasso.fit, newx, 
              group.mat, s = lambda, 
              use.refit = use.refit, ...))
    }
}
