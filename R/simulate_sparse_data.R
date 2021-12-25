#' function to generate data with hierarchical sparsity
#'
#' @param nvars number of variables
#' @param noutcomes number of outcomes
#' @param nobs number of observations per outcomes to simulate
#' @param nobs.test number of independent test observations per outcomes to simulate
#' @param outcome.grouping groups for outcomes -- either a vector of length \code{noutcomes} of indices of what groups each outcome belongs to
#' or a matrix with \code{noutcomes} columns, with each row as a different grouping of the outcomes -- this allows for potentially overlapping
#' groups of the outcomes
#' @param effect.size.similarity How should the effect sizes of a variable across the outcomes be related? \code{"sign"} for
#' the signs to match, \code{"mean"} for the coefficients to be the same plus some small uniform noise, or \code{"none"} for there
#' to be no relationship at all
#' @param hier.sparsity.prob probability that a group has its coefficients set to zero for any variable
#' @param prop.zero.vars proportion of all variables that will be zero across all outcomes
#' @param effect.size.max maximum magnitude of the true effect sizes
#' @param family family for the response variable
#' @param sd standard devation for gaussian simulations
#' @param beta a matrix of true beta values. If given, then no beta will be created and data will be simulated from the given beta
#' @param x.covar scalar, pairwise covariance term for covariates
#' @param y.covar scalar, pairwise covariance term for outcomes
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
#'                    prop.zero.vars = 0.25,
#'                    outcome.grouping = rbind(c(1,1,1,2,2,2,2,2),
#'                                             c(1,1,1,2,2,3,3,3),
#'                                             c(1,1,2,3,3,4,4,5),
#'                                             c(1:8)))
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
#'                         c(1,1,1,2,2,3,3,3),
#'                         c(1,1,2,3,3,4,4,5))
#'                         
#' fit.adapt <- cv.groupFusedMulti(x, y,
#'                                 nlambda        = 50,
#'                                 lambda.fused = c(0.000005, 0.00001, 0.000025, 0.00005, 0.0001),
#'                                 outcome.groups = outcome_groups,
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
                                     outcome.grouping = rbind(c(1,1,1,2,2,2,2,2),
                                                              c(1,1,1,2,2,3,3,3),
                                                              c(1,1,2,3,3,4,4,5),
                                                              c(1:8)),
                                     effect.size.similarity = c("sign", "mean", "none"),
                                     hier.sparsity.prob = 0.1,
                                     group.fused.prob = 0.25,
                                     prop.zero.vars = 0.5,
                                     family = c("gaussian", "binomial"),
                                     sd  = 1,
                                     snr = NULL,
                                     beta = NULL,
                                     tau = 10,
                                     x.covar = 0,
                                     y.covar = 0.5
)
{

    effect.size.similarity <- match.arg(effect.size.similarity)
    family <- match.arg(family)

    stopifnot(x.covar[1] >= 0)
    stopifnot(y.covar[1] >= 0)
    
    sig_X <- matrix(x.covar, nrow = nvars, ncol = nvars)
    diag(sig_X) <- 1
    
    sig_y <- matrix(y.covar, nrow = noutcomes, ncol = noutcomes)
    diag(sig_y) <- 1
    
    ## training x/eps
    x   <- MASS::mvrnorm(n = nobs, mu = numeric(nvars), Sigma = sig_X)
    eps <- MASS::mvrnorm(n = nobs, mu = numeric(noutcomes), Sigma = (sd^2) * sig_y)
    
    ## test x/eps
    x.test   <- MASS::mvrnorm(n = nobs.test, mu = numeric(nvars), Sigma = sig_X)
    eps.test <- MASS::mvrnorm(n = nobs.test, mu = numeric(noutcomes), Sigma = (sd^2) * sig_y)
    
    n_nonzero <- floor((1-prop.zero.vars) * nvars)
    beta <- matrix(0, nrow = nvars, ncol = noutcomes)
    beta[1:n_nonzero,] <- sample(c(1,0.5, 0.25, -1, -0.5, -0.25), size = n_nonzero * noutcomes, replace = TRUE)
    
    
    if (!is.matrix(outcome.grouping))
    {
      outcome_groups <- unique(outcome.grouping)
      n_groups <- length(outcome_groups)
      
      ## make groups sparse
      for (j in 1:n_nonzero)
      {
        for (g in 1:n_groups)
        {
            beta[j,outcome.grouping == outcome_groups[g]] <- beta[j,outcome.grouping == outcome_groups[g]] * 
              rbinom(1, size = 1, prob = 1-hier.sparsity.prob)
            
            ## make all coefs in group equal with certain probability
            fuse_group <- rbinom(1, size = 1, prob = group.fused.prob)
            if (fuse_group == 1)
            {
              beta[j,outcome.grouping == outcome_groups[g]] <- beta[j,outcome.grouping == outcome_groups[g]]
            }
        }
      }
    } else
    {
      for (r in 1:nrow(outcome.grouping))
      {
        outcome_groups <- unique(outcome.grouping[r,])
        n_groups <- length(outcome_groups)
        
        ## make groups sparse
        for (j in 1:n_nonzero)
        {
          for (g in 1:n_groups)
          {
            beta[j,outcome.grouping[r,] == outcome_groups[g]] <- beta[j,outcome.grouping[r,] == outcome_groups[g]] * 
              rbinom(1, size = 1, prob = 1-hier.sparsity.prob)
            
            ## make all coefs in group equal with certain probability
            fuse_group <- rbinom(1, size = 1, prob = group.fused.prob)
            if (fuse_group == 1)
            {
              beta[j,outcome.grouping[r,] == outcome_groups[g]] <- beta[j,outcome.grouping[r,] == outcome_groups[g]][1]
            }
          }
        }
      }
    }
    
    
    xbeta <- x %*% beta
    x.testbeta <- x.test %*% beta
    
    y <- xbeta + eps
    y.test <- x.testbeta + eps.test
    
    
    list(beta = beta, 
         x = x, y = y, 
         x.test = x.test, y.test = y.test)
}


