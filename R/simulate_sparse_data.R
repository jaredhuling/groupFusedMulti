#' function to generate data with hierarchical sparsity
#'
#' @param nvars number of variables
#' @param noutcomes number of outcomes
#' @param nobs number of observations per outcomes to simulate
#' @param nobs.test number of independent test observations per outcomes to simulate
#' @param outcome.groups groups for outcomes -- either a vector of length \code{noutcomes} of indices of what groups each outcome belongs to
#' or a matrix with \code{noutcomes} columns, with each row as a different grouping of the outcomes -- this allows for potentially overlapping
#' groups of the outcomes
#' @param hier.sparsity.prob probability that a group has its coefficients set to zero for any variable
#' @param individ.sparsity.prob probability that a specific coefficient is set to zero
#' @param group.fused.prob probability that all coefficients in a group are set to be equal to each other
#' @param num.nonzero.vars number variables that will NOT have zero effect across all outcomes
#' @param family family for the response variable; either \code{"gaussian"}, \code{"binomial"}, or \code{"ordinal"}
#' @param x.distribution distribution of covariates \code{x}: either \code{"gaussian"} or \code{"binary"}
#' @param sd standard devation for gaussian simulations
#' @param beta a matrix of true beta values. If given, then no beta will be created and data will be simulated from the given beta
#' @param ordinal.link.cutoffs sorted (ascending) vector of cutoff values for latent variable model when response family is ordinal. ordinal
#' response will have number of levels equal to the length of \code{ordinal.link.cutoffs} plus 1
#' @param x.binary.prob if \code{x.distribution} is binary, this is the probability that each entry in the covariates is 1 (versus 0)
#' @param x.rho scalar, AR(1) correlation parameter for covariate distribution
#' @param y.rho scalar, AR(1) correlation parameter for outcomes
#' @importFrom stats rnorm
#' @importFrom stats var
#' @importFrom stats approx
#' @importFrom stats runif
#' @import MASS
#' @export
#' @examples
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
#' \dontrun{
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
#'
#'
gen_sparse_multivar_data <- function(nvars = 10L,
                                     noutcomes = 8L,
                                     nobs = 100L,
                                     nobs.test = 100L,
                                     outcome.groups = rbind(c(1,1,1,2,2,2,2,2),
                                                            c(1,1,1,2,2,3,3,3)),
                                     hier.sparsity.prob = 0.1,
                                     individ.sparsity.prob = 0.25,
                                     group.fused.prob = 0.5,
                                     num.nonzero.vars = min(nvars, 25L),
                                     family = c("gaussian", "binomial", "ordinal"),
                                     x.distribution = c("gaussian", "binary"),
                                     sd  = 1,
                                     beta = NULL,
                                     ordinal.link.cutoffs = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                                     x.binary.prob = 0.2,
                                     x.rho = 0.5,
                                     y.rho = 0.5)
{

    family <- match.arg(family)
    x.distribution <- match.arg(x.distribution)

    #stopifnot(x.rho[1] >= 0)
    stopifnot(y.rho[1] >= 0)
    
    stopifnot(num.nonzero.vars <= nvars)
    
    
    
    sig_y <- y.rho ^ abs(outer(1:noutcomes, 1:noutcomes, FUN = "-"))
    
    ## training eps
    eps <- MASS::mvrnorm(n = nobs, mu = numeric(noutcomes), Sigma = (sd^2) * sig_y)
    
    ## test eps
    eps.test <- MASS::mvrnorm(n = nobs.test, mu = numeric(noutcomes), Sigma = (sd^2) * sig_y)
    
    if (x.distribution == "gaussian")
    {
      sig_X <- x.rho ^ abs(outer(1:nvars, 1:nvars, FUN = "-"))
      
      ## training x
      x   <- MASS::mvrnorm(n = nobs, mu = numeric(nvars), Sigma = sig_X)
      ## test x
      x.test   <- MASS::mvrnorm(n = nobs.test, mu = numeric(nvars), Sigma = sig_X)
      
    } else if (x.distribution == "binary") 
    {
      ## training x
      x <- matrix(rbinom(nobs * nvars, 1, prob = x.binary.prob), ncol = nvars)
      ## test x
      x.test <- matrix(rbinom(nobs.test * nvars, 1, prob = x.binary.prob), ncol = nvars)
    }
    
    
    
    n_nonzero <- num.nonzero.vars
    beta <- matrix(0, nrow = nvars, ncol = noutcomes)
    beta[1:n_nonzero,] <- sample(c(1,0.5, 0.25, 0.125, -1, -0.5, -0.125), size = n_nonzero * noutcomes, replace = TRUE)
    
    
    
    ## make some effects same across outcomes
    for (j in 1:n_nonzero)
    {
      ## possibly fuse all
      fuse_group <- rbinom(1, size = 1, prob = group.fused.prob)
      if (fuse_group == 1)
      {
        beta[j,] <- beta[j,][1]
      }
    }
    
    if (!is.matrix(outcome.groups))
    {
      outcome_groups <- unique(outcome.groups)
      n_groups <- length(outcome_groups)
      
      
      ## make groups sparse
      for (j in 1:n_nonzero)
      {
        for (g in 1:n_groups)
        {
            beta[j,outcome.groups == outcome_groups[g]] <- beta[j,outcome.groups == outcome_groups[g]] * 
              rbinom(1, size = 1, prob = 1-hier.sparsity.prob)
            
            ## make all coefs in group equal with certain probability
            fuse_group <- rbinom(1, size = 1, prob = group.fused.prob)
            if (fuse_group == 1)
            {
              beta[j,outcome.groups == outcome_groups[g]] <- beta[j,outcome.groups == outcome_groups[g]]
            }
        }
      }
    } else
    {
      for (r in 1:nrow(outcome.groups))
      {
        outcome_groups <- unique(outcome.groups[r,])
        n_groups <- length(outcome_groups)
        
        ## make groups sparse
        for (j in 1:n_nonzero)
        {
          for (g in 1:n_groups)
          {
            beta[j,outcome.groups[r,] == outcome_groups[g]] <- beta[j,outcome.groups[r,] == outcome_groups[g]] * 
              rbinom(1, size = 1, prob = 1-hier.sparsity.prob)
            
            ## make all coefs in group equal with certain probability
            fuse_group <- rbinom(1, size = 1, prob = group.fused.prob)
            if (fuse_group == 1)
            {
              beta[j,outcome.groups[r,] == outcome_groups[g]] <- beta[j,outcome.groups[r,] == outcome_groups[g]][1]
            }
          }
        }
      }
      

    }
    
    ## individual coefficient sparsity-inducing matrix
    sparsity_individ_mask <- matrix(rbinom(nvars * noutcomes, size = 1, 
                                           prob = 1 - individ.sparsity.prob), 
                                    nrow = nvars, ncol = noutcomes)
    
    beta <- beta * sparsity_individ_mask
    
    if (x.distribution == "binary") 
    {
      beta <- beta * 4
    }
    
    
    xbeta <- x %*% beta
    x.testbeta <- x.test %*% beta
    
    
    if (family == "gaussian")
    {
      y <- xbeta + eps
      y.test <- x.testbeta + eps.test
    } else if (family == "binomial")
    {
      y <- apply(xbeta, 2, function(xx) rbinom(NROW(xx), 1, 1/(1+exp(-xx))))
      y.test <- apply(x.testbeta, 2, function(xx) rbinom(NROW(xx), 1, 1/(1+exp(-xx)))) 
    } else if (family == "ordinal")
    {
      y_star <- xbeta + eps
      y.test_star <- x.testbeta + eps.test
      
      y <- apply(y_star, 2, function(yy) {
        yy <- cut(yy, breaks = c(-Inf, ordinal.link.cutoffs, Inf), labels = 1:(length(ordinal.link.cutoffs) + 1))
        yy <- as.numeric(levels(yy)[yy])
      })
      
      y.test <- apply(y.test_star, 2, function(yy) {
        yy <- cut(yy, breaks = c(-Inf, ordinal.link.cutoffs, Inf), labels = 1:(length(ordinal.link.cutoffs) + 1))
        yy <- as.numeric(levels(yy)[yy])
      })
    }
    
    
    list(beta = beta, 
         x = x, y = y, 
         x.test = x.test, 
         y.test = y.test)
}


