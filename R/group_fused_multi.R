


#' Fitting groupFusedMulti 
#'
#' @param x input matrix of dimension nobs by nvars. Each row is an observation,
#' each column corresponds to a covariate
#' @param y numeric response matrix with \code{nobs} rows and \code{noutcomes} number of columns (one outcome per column)
#' @param outcome.groups a vector of length equal to the number of columns in the outcome matrix \code{y}
#' @param family \code{"gaussian"} for least squares problems
#' @param nlambda The number of lambda values. Default is 100.
#' @param nlambda.fused The number of lambda values for the fused lasso. Default is 0, which means that the fused
#' lasso is not used.
#' @param lambda A user-specified sequence of lambda values. Left unspecified, the a sequence of lambda values is
#' automatically computed, ranging uniformly on the log scale over the relevant range of lambda values.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of \code{lambda.max}, the (data derived) entry
#' value (i.e. the smallest value for which all parameter estimates are zero). The default
#' depends on the sample size \code{nobs} relative to the number of variables \code{nvars}. If
#' \code{nobs > nvars}, the default is 0.0001, close to zero. If \code{nobs < nvars}, the default
#' is 0.01. A very small value of \code{lambda.min.ratio} can lead to a saturated fit
#' when \code{nobs < nvars}.
#' @param lambda.fused tuning parameter for fused lasso penalty
#' @param use.alpha.param Should the 'alpha' parameterization for the fused lasso be used? alpha=1 will correspond to no fused lasso and alpha=0 will be all fused lasso
#' @param fuse.none.group Should the coefficients of the none subpopulation be fused with coefficients for subpopulations with just one condition?
#' @param lambda.fused.min minimum value used for the fused lasso tuning parameter
#' @param lambda.fused.max.ratio scalar value that is the ratio of the largest to smallest values
#' for the fused lasso tuning parameter
#' @param penalty.factor vector of weights to be multiplied to the tuning parameter for the
#' group lasso penalty. A vector of length equal to the number of groups
#' @param group.weights A vector of values representing multiplicative factors by which each group's penalty is to
#' be multiplied. Often, this is a function (such as the square root) of the number of predictors in each group.
#' The default is to use the square root of group size for the group selection methods.
#' @param adaptive.lasso Flag indicating whether or not to use adaptive lasso weights. If set to \code{TRUE} and
#' \code{group.weights} is unspecified, then this will override \code{group.weights}. If a vector is supplied to group.weights,
#' then the \code{adaptive.lasso} weights will be multiplied by the \code{group.weights} vector
#' @param adaptive.fused Flag indicating whether or not to use adaptive fused lasso weights. 
#' @param standardize Should the data be standardized? Defaults to \code{FALSE}.
#' @param intercept Should an intercept be fit? Defaults to \code{TRUE}
#' @param one.intercept Should a single intercept be fit for all outcome instead of one
#' for each outcome? Defaults to \code{FALSE}.
#' @param gamma power to raise the MLE estimated weights by for the adaptive lasso. defaults to 1
#' @param rho ADMM parameter. must be a strictly positive value. By default, an appropriate value is automatically chosen
#' @param dynamic.rho \code{TRUE}/\code{FALSE} indicating whether or not the rho value should be updated throughout the course of the ADMM iterations
#' @param maxit integer. Maximum number of ADMM iterations. Default is 500.
#' @param abs.tol absolute convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param rel.tol relative convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param irls.maxit integer. Maximum number of IRLS iterations. Only used if \code{family != "gaussian"}. Default is 100.
#' @param irls.tol convergence tolerance for IRLS iterations. Only used if \code{family != "gaussian"}. Default is 10^{-5}.
#' @param model.matrix logical flag. Should the design matrix used be returned?
#' @param ... not used
#' @return An object with S3 class "groupFusedMulti"
#'
#'
#' @import Rcpp
#'
#' @export
#' @examples
#' library(Matrix)
#' 
groupFusedMulti <- function(x, 
                            y,
                            outcome.groups   = NULL,
                            family           = c("gaussian", "binomial"),
                            nlambda          = 100L,
                            nlambda.fused    = 0,
                            lambda           = NULL,
                            lambda.min.ratio = NULL,
                            lambda.fused     = NULL,
                            use.alpha.param  = TRUE, 
                            fuse.none.group  = TRUE,
                            lambda.fused.min       = 1e-5,
                            lambda.fused.max.ratio = 1e-5,
                            penalty.factor   = NULL,
                            group.weights    = NULL,
                            adaptive.lasso   = FALSE,
                            adaptive.fused   = FALSE,
                            gamma            = 1,
                            standardize      = FALSE,
                            intercept        = TRUE,
                            one.intercept    = FALSE,
                            rho              = NULL,
                            dynamic.rho      = TRUE,
                            maxit            = 500L,
                            abs.tol          = 1e-5,
                            rel.tol          = 1e-5,
                            irls.tol         = 1e-5,
                            irls.maxit       = 100L,
                            model.matrix     = FALSE,
                            ...)
{
    y <- as.matrix(y)
    
    
    family <- match.arg(family)
    this.call = match.call()
    
    if (family != "gaussian") stop("non-gaussian families not available yet")
    
    
    onames <- colnames(y)
    vnames <- colnames(x)
    
    if (is.null(onames))
    {
        onames <- paste0("Outcome_", 1:NCOL(y))
    }
    if (is.null(vnames))
    {
        vnames <- paste0("X_", 1:NCOL(x))
    }
    
    if (intercept)
    {
        vnames <- c("(Intercept)", vnames)
    }
    
    nvars     <- ncol(x)
    nobs      <- nrow(x)
    leny      <- NROW(y)
    noutcomes <- NCOL(y)
    ngroups   <- length(unique(outcome.groups))
    
    if (noutcomes <= 1) stop("outcome y must be a matrix with more than 1 column")
    
    stopifnot(leny == nobs)
    
    nlambda.fused <- as.integer(nlambda.fused[1])
    
    
    if (nlambda.fused > 0 & is.null(lambda.fused))
    {
        
        lambda.fused.max.ratio <- as.double(lambda.fused.max.ratio[1])
        lambda.fused.min       <- as.double(lambda.fused.min[1])
        
        if (lambda.fused.max.ratio <= 0) stop("'lambda.fused.max.ratio' must be strictly positive")
        if (lambda.fused.min <= 0)       stop("'lambda.fused.min' must be strictly positive")
        
        if (use.alpha.param)
        {
            lambda.fused <- seq(1e-5, 1-1e-5, length.out = nlambda.fused)
        } else
        {
            lambda.fused <- exp(seq(log(lambda.fused.min), 
                                    log(lambda.fused.min / lambda.fused.max.ratio), 
                                    length.out = nlambda.fused))
        }
        
    }
    
    ## save lambda as a matrix with ncol = length(lambda.fused)
    ## if both lambda.fused and lambda are provided but lambda is a vector
    # if (!is.null(lambda.fused) & !is.null(lambda))
    # {
    #     if (!is.matrix(lambda))
    #     {
    #         lambda <- matrix(rep(lambda, length(lambda.fused)), ncol = length(lambda.fused))
    #     }
    # }
    
    
    ## we only care about one.intercept if intercept = TRUE
    one.intercept <- !(!one.intercept & intercept)
    
    
    ######################################################################
    ##
    ##                          set up design matrix
    ##
    ######################################################################
    

    
    col_sd_positive <- apply(x,2,sd)>0
    xs <- x[,col_sd_positive]
    
    (pfull <- ncol(x) + intercept)
    
    (p <- ncol(xs) + intercept)
    blocks <- ncol(y)
    
    if (intercept)
    {
        x_bdiag <- Matrix::bdiag(rep(list(cbind(1, xs)), ncol(y)))
    } else
    {
        x_bdiag <- Matrix::bdiag(rep(list(xs), ncol(y)))
    }
    
    y_big <- as.vector(as.matrix(y))
    
    ######################################################################
    ##
    ##                          set up penalty matrices
    ##
    ######################################################################
    
    
    
    
    ## terms for individual selection
    groups_lasso <- lapply(1:(p * blocks), function(x) x)
    
    ## terms for selection of entire var across all models
    groups_all <- lapply((1):p, function(x) c(x + p * (0:(blocks-1))))
    
    
    penalty.factor.out_groups <- NULL
    
    out_groups <- NULL
    if (!is.null(outcome.groups))
    {
        if (is.matrix(outcome.groups))
        {
            if (ncol(outcome.groups) != noutcomes) stop("'outcome.groups' must have number of columns equal to the number of columns in 'y'")
            
            
            for (r in 1:nrow(outcome.groups))
            {
                unique_outcome_groups <- unique(outcome.groups[r,])
                
                ## terms for selection of var across outcome group-specific terms
                for (g in 1:length(unique_outcome_groups))
                {
                    gr_idx <- which(outcome.groups[r,] == unique_outcome_groups[g])
                    
                    new_terms  <- lapply((1):p, function(x) c(x + p * (0:(blocks-1)))[gr_idx] )
                    out_groups <- c(out_groups, new_terms)
                    
                    new_penfact <- rep(1, length(new_terms))
                    if (intercept) new_penfact[1] <- 0
                    penalty.factor.out_groups <- c(penalty.factor.out_groups, new_penfact)
                }
            }
            
        } else
        {
            
            if (length(outcome.groups) != noutcomes) stop("'outcome.groups' must have length equal to the number of columns in 'y'")
            
            unique_outcome_groups <- unique(outcome.groups)
            
            ## terms for selection of var across outcome group-specific terms
            for (g in 1:length(unique_outcome_groups))
            {
                gr_idx <- which(outcome.groups == unique_outcome_groups[g])
                
                new_terms  <- lapply((1):p, function(x) c(x + p * (0:(blocks-1)))[gr_idx] )
                out_groups <- c(out_groups, new_terms)
                
                new_penfact <- rep(1, length(new_terms))
                if (intercept) new_penfact[1] <- 0
                penalty.factor.out_groups <- c(penalty.factor.out_groups, new_penfact)
            }
        }
        
        ## remove groups with just one outcome
        keep_idx_1 <- sapply(out_groups, length) > 1
        penalty.factor.out_groups <- penalty.factor.out_groups[keep_idx_1]
        out_groups                <- out_groups[keep_idx_1]
        
        ## remove duplicates
        dupl_idx <- duplicated(lapply(out_groups, sort))
        penalty.factor.out_groups <- penalty.factor.out_groups[!dupl_idx]
        out_groups                <- out_groups[!dupl_idx]
    }
    
    penalty.factor.lasso <- rep(1, length(groups_lasso))
    
    
    if (intercept) penalty.factor.lasso[(0:(blocks-1)) * (p)+1] <- 0
    
    penalty.factor.groups_all <- rep(1, length(groups_all))
    
    if (intercept) penalty.factor.groups_all[1] <- 0
    
    if (ngroups > 1)
    {
        penalty.factor <- c(penalty.factor.lasso, penalty.factor.out_groups, penalty.factor.groups_all)
        
        groups <- c(groups_lasso, out_groups, groups_all)
    } else
    {
        penalty.factor <- c(penalty.factor.lasso, penalty.factor.out_groups)
        
        groups <- c(groups_lasso, out_groups)
    }
    
    
    
    fit_fused <- FALSE
    
    if (is.matrix(outcome.groups))
    {
        outcome.groups.fuse <- outcome.groups[nrow(outcome.groups),]
    } else
    {
        outcome.groups.fuse <- outcome.groups
    }
    
    group_sizes   <- table(outcome.groups.fuse)
    unique.groups <- names(group_sizes)
    
    groups.to.fuse <- unique.groups[group_sizes>1]
    group_sizes_fuse <- group_sizes[group_sizes>1]
    
    if (!is.null(lambda.fused))
    {
        fit_fused <- TRUE
        ## each column is a variable, each row is a penalty
        fused_mat <- as(sparseMatrix(i=1,j=1,dims = c(p * ( sum(choose(group_sizes_fuse, 2)) ), p * blocks)), "dgCMatrix")
        
        ct <- 0
        for (j in 1:p)
        {
            ## all pairs of coefficients for each variable within each outcome group (not all possible fusings/pairs)
            for (gg in 1:length(groups.to.fuse))
            {
                cur_o_idx  <- which(outcome.groups.fuse == groups.to.fuse[gg])
                cur_blocks <- group_sizes_fuse[gg]
                for (bb in 1:(cur_blocks-1))
                {
                    for (b in (bb+1):(cur_blocks))
                    {
                        ct <- ct + 1
                        fused_mat[ct,j + p * (cur_o_idx[bb] - 1)] <- 1
                        fused_mat[ct,j + p * (cur_o_idx[b] - 1)]  <- -1
                    }
                }
            }
        }
    } else
    {
        fused_mat <- NULL
        lambda.fused <- 0
    }
    
    
    
    
    
    
    ######################################################################
    ##
    ##                          compute model
    ##
    ######################################################################

    fit <- groupFusedMulti::oglasso(x                    = x_bdiag,
                                    y                    = y_big,
                                    group                = groups,
                                    fused                = fused_mat,
                                    family               = family,
                                    nlambda              = nlambda,
                                    lambda               = lambda,
                                    lambda.min.ratio     = lambda.min.ratio,
                                    lambda.fused         = lambda.fused,
                                    use.alpha.param      = use.alpha.param,
                                    penalty.factor       = penalty.factor,
                                    penalty.factor.fused = NULL,
                                    group.weights        = group.weights,
                                    adaptive.lasso       = adaptive.lasso,
                                    adaptive.fused       = adaptive.fused,
                                    gamma                = gamma,
                                    standardize          = standardize,
                                    intercept            = one.intercept,
                                    dynamic.rho          = dynamic.rho,
                                    maxit                = maxit,
                                    abs.tol              = abs.tol,
                                    rel.tol              = rel.tol,
                                    irls.tol             = irls.tol,
                                    irls.maxit           = irls.maxit,
                                    ...)
    
    nlambda <- fit$nlambda
    
    class2 <- switch(family,
                     "gaussian" = "gfmgaussian",
                     "binomial" = "gfmbinomial", "coxph" = "gfmcoxph")
    
    beta.tmp   <- fit$beta

    
    ######################################################################
    ##
    ##                    place betas in a matrix
    ##                   dim = (n vars) x (n outcomes) x (n lambdas)
    ##
    ######################################################################
    
    
    nlambda_fused <- length(lambda.fused)
    
    dim_adjust <- pfull != p
    
    if (fit_fused)
    {
        beta.array <- vector(mode = "list", length = length(lambda.fused))
        for (lf in 1:length(lambda.fused))
        {
            
            
            
            ## place the coefficient results
            ## in an array instead of a long vector
            
            beta_all      <- beta.tmp[[lf]][-1,]
            intercept_all <- beta.tmp[[lf]][1,]
            
            if (dim_adjust)
            {
                beta_all_new <- matrix(0, nrow = pfull * blocks, ncol = nlambda)
                for (b in 1:blocks)
                {
                    ## adding the true for the outcome-specific intercept position
                    beta_all_new[(((b-1)*pfull + 1):(b * pfull))[c(TRUE,col_sd_positive)],] <- beta_all[((b-1)*p + 1):(b * p),]
                }

                
                beta_all <- beta_all_new
                p_use <- pfull
            } else
            {
                p_use <- p
            }
            
            
            beta.array.tmp <- array(0, dim = c(p_use, noutcomes, nlambda))
            
            for (i in 1:nlambda)
            {
                if (intercept)
                {
                    for (b in 1:blocks)
                    {
                        beta.array.tmp[,b,i] <- beta_all[((b-1)*p_use + 1):(b * p_use),i]
                        beta.array.tmp[,b,i] <- beta.array.tmp[,b,i] + intercept_all[i]
                    }
                } else
                {
                    for (b in 1:blocks)
                    {
                        beta.array.tmp[-1,b,i] <- beta_all[((b-1)*p_use + 1):(b * p_use),i]
                    }
                }
            }
            
            if (is.matrix(fit$lambda))
            {
                dimnames(beta.array.tmp) <- list(vnames, onames, paste(seq(along = fit$lambda[,lf])))
            } else
            {
                dimnames(beta.array.tmp) <- list(vnames, onames, paste(seq(along = fit$lambda)))
            }
            
            beta.array[[lf]] <- beta.array.tmp
        }
        
        names(beta.array) <- paste0("lam_fuse", seq(along = lambda.fused))
    } else
    {
        
        
        
        ## place the coefficient results
        ## in an array instead of a long vector
        
        beta_all      <- beta.tmp[-1,]
        intercept_all <- beta.tmp[1,]
        
        if (dim_adjust)
        {
            beta_all_new <- matrix(0, nrow = pfull * blocks, ncol = nlambda)
            for (b in 1:blocks)
            {
                ## adding the true for the outcome-specific intercept position
                beta_all_new[(((b-1)*pfull + 1):(b * pfull))[c(TRUE, col_sd_positive)],] <- beta_all[((b-1)*p + 1):(b * p),]
            }
            
            
            beta_all <- beta_all_new
            p_use <- pfull
        } else
        {
            p_use <- p
        }
        
        beta.array <- array(0, dim = c(p_use, noutcomes, nlambda))
        
        
        for (i in 1:nlambda)
        {
            if (intercept)
            {
                for (b in 1:blocks)
                {
                    beta.array[,b,i] <- beta_all[((b-1)*p_use + 1):(b * p_use),i]
                    beta.array.tmp[,b,i] <- beta.array.tmp[,b,i] + intercept_all[i]
                }
            } else
            {
                for (b in 1:blocks)
                {
                    beta.array[-1,b,i] <- beta_all[((b-1)*p_use + 1):(b * p_use),i]
                }
            }
        }
        
        dimnames(beta.array) <- list(vnames, onames, paste(seq(along = fit$lambda)))
        
    }
    
    
    
    fit$beta                   <- beta.array
    fit$nobs                   <- nobs
    fit$nvars                  <- nvars
    fit$noutcomes              <- noutcomes
    fit$intercept.fit          <- intercept
    fit$standardized           <- standardize
    fit$is_fused               <- fit_fused
    
    
    fit$one.intercept          <- one.intercept
    fit$use.alpha.param        <- use.alpha.param
    
    fit$var.names      <- vnames
    fit$outcome.names  <- onames
    fit$outcome.groups <- outcome.groups
    
    class(fit) <- c("groupFusedMulti", class2)
    fit
}



