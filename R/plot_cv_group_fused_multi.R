
#' Plot method for cv.groupFusedMulti fitted objects
#'
#' @param sign.lambda Either plot against log(lambda) (default) or its negative if \code{sign.lambda = -1}.
#' @param best.fused if \code{TRUE}, make a CV plot only for the best \code{lambda.fused} value. This includes
#' standard errors in the plot. If \code{FALSE}, plot CV for each value of \code{lambda.fused}. In this case,
#' no standard errors will be plotted.
#' @rdname plot
#' @method plot cv.groupFusedMulti
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics title
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics par
#' @export
#' @examples
#' set.seed(123)
#' 
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 25,
#'                              nobs = 100, 
#'                              hier.sparsity.param = 0.5,
#'                              prop.zero.vars = 0.5,
#'                              effect.size.max = 0.25,
#'                              family = "gaussian")
#'
#' x        <- dat.sim$x
#' x.test   <- dat.sim$x.test
#' y        <- dat.sim$y
#' y.test   <- dat.sim$y.test
#' grp      <- dat.sim$group.ind
#' grp.test <- dat.sim$group.ind.test
#'
#' fit.adapt <- cv.groupFusedMulti(x, y,
#'                           grp,
#'                           adaptive.lasso = TRUE,
#'                           nlambda        = 25,
#'                           nfolds         = 4)
#'                                      
#' plot(fit.adapt) 
#' 
plot.cv.groupFusedMulti <- function(x, sign.lambda = 1, best.fused = FALSE, ...) 
{
    # compute total number of selected variables for each
    # tuning parameter 
    
    if (!is.null(x$lambda.fused))
    {
        if (best.fused)
        {
            cvm   <- x$cvm[[x$which.min.fused]]
            cvup  <- x$cvup[[x$which.min.fused]]
            cvlo  <- x$cvlo[[x$which.min.fused]]
            nzero <- apply(x$groupFusedMulti.fit$beta[[x$which.min.fused]][,-1,], 3, function(bb) sum(bb != 0))
            
            xlab <- expression(log(lambda))
            if(sign.lambda<0)xlab <- paste("-", xlab, sep="")
            if (is.matrix(x$lambda))
            {
                plot.args = list(x    = sign.lambda * log(x$lambda[,x$which.min.fused]),
                                 y    = cvm,
                                 ylim = range(cvup, cvlo),
                                 xlab = xlab,
                                 ylab = x$name,
                                 type = "n")
            } else
            {
                plot.args = list(x    = sign.lambda * log(x$lambda),
                                 y    = cvm,
                                 ylim = range(cvup, cvlo),
                                 xlab = xlab,
                                 ylab = x$name,
                                 type = "n")
            }
            new.args <- list(...)
            if(length(new.args)) plot.args[names(new.args)] <- new.args
            do.call("plot", plot.args)
            
            
            if (is.matrix(x$lambda))
            {
                error.bars(sign.lambda * log(x$lambda[,x$which.min.fused]), 
                           cvup, 
                           cvlo, width = 0.005)
                points(sign.lambda*log(x$lambda[,x$which.min.fused]), cvm, pch=20, col="dodgerblue")
                axis(side   = 3,
                     at     = sign.lambda * log(x$lambda[,x$which.min.fused]),
                     labels = paste(nzero), tick=FALSE, line=0, ...)
            } else
            {
                error.bars(sign.lambda * log(x$lambda), 
                           cvup, 
                           cvlo, width = 0.005)
                points(sign.lambda*log(x$lambda), cvm, pch=20, col="dodgerblue")
                axis(side   = 3,
                     at     = sign.lambda * log(x$lambda),
                     labels = paste(nzero), tick=FALSE, line=0, ...)
            }
            
            
            
            
            abline(v = sign.lambda * log(x$lambda.min), lty=2, lwd = 2, col = "firebrick3")
            abline(v = sign.lambda * log(x$lambda.1se), lty=2, lwd = 2, col = "firebrick1")
            title(x$name, line = 2.5, ...)
        } else
        {
            old.par <- par(no.readonly = TRUE)
            
            layout(t(1:2), widths=c(6,1))
            
            # Set margins and turn all axis labels horizontally (with `las=1`)
            par(mar=rep(.5, 4), oma = rep(3, 4), las=1)
            
            cols <- vpal(length(x$cvm))
            
            ylim <- range(c(unlist(x$cvm)))
            
            maxes <- sapply(x$cvm, max)
            
            rem_idx <- NULL
            for (m in 1:length(maxes))
            {
                if (maxes[m] > 5 * max(maxes[-m]))
                {
                    rem_idx <- c(rem_idx, m)
                }
            }
            
            if (!is.null(rem_idx))
            {
                ylim <- range(c(unlist(x$cvm[-rem_idx])))
            }
            
            xlab <- expression(log(lambda))
            
           
            
            
            
            if(sign.lambda<0)xlab <- paste("-", xlab, sep="")
            if (is.matrix(x$lambda))
            {
                if (x$groupFusedMulti.fit$use.alpha.param)
                {
                    lamplot <- x$lambda[,ncol(x$lambda)] * x$lambda.fused[length(x$lambda.fused)]
                } else
                {
                    lamplot <- x$lambda[,1]
                }
                plot.args = list(x    = sign.lambda * log(lamplot),
                                 y    = x$cvm[[1]],
                                 ylim = ylim,
                                 xlab = xlab,
                                 ylab = x$name,
                                 col  = cols[1],
                                 type = "n")
            } else
            {
                lamplot <- x$lambda
                plot.args = list(x    = sign.lambda * log(lamplot),
                                 y    = x$cvm[[1]],
                                 ylim = ylim,
                                 xlab = xlab,
                                 ylab = x$name,
                                 col  = cols[1],
                                 type = "n")
            }
            new.args <- list(...)
            if(length(new.args)) plot.args[names(new.args)] <- new.args
            do.call("plot", plot.args)
            
            for (lf in 2:length(x$cvm))
            {
                lines(x    = sign.lambda * log(lamplot),
                      y    = x$cvm[[lf]],
                      col  = cols[lf])
            }
            
            # legend("bottomright", paste0("lf: ", round(x$lambda.fused, 4)), 
            #        col = cols, lty = 1, cex = 0.5,
            #        inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
            
            if (x$groupFusedMulti.fit$use.alpha.param)
            {
                image(1, log(x$lambda.fused), t(seq_along(x$lambda.fused)), col=cols, axes=FALSE,
                      main = expression(lambda~fused))
                axis(4, labels = round(x$lambda.fused, 5), at = log(x$lambda.fused))
            } else
            {
                image(1, log(x$lambda.fused), t(seq_along(x$lambda.fused)), col=cols, axes=FALSE,
                      main = expression(lambda~fused))
                axis(4, labels = round(x$lambda.fused, 5), at = log(x$lambda.fused))
            }
            
            par(old.par)
        }
    } else
    {
        nzero <- apply(x$groupFusedMulti.fit$beta[,-1,], 3, function(bb) sum(bb != 0))
        
        xlab <- expression(log(lambda))
        if(sign.lambda<0)xlab <- paste("-", xlab, sep="")
        plot.args = list(x    = sign.lambda * log(x$lambda),
                         y    = x$cvm,
                         ylim = range(x$cvup, x$cvlo),
                         xlab = xlab,
                         ylab = x$name,
                         type = "n")
        new.args <- list(...)
        if(length(new.args)) plot.args[names(new.args)] <- new.args
        do.call("plot", plot.args)
        error.bars(sign.lambda * log(x$lambda), 
                   x$cvup, 
                   x$cvlo, width = 0.005)
        points(sign.lambda*log(x$lambda), x$cvm, pch=20, col="dodgerblue")
        axis(side   = 3,
             at     = sign.lambda * log(x$lambda),
             labels = paste(nzero), tick=FALSE, line=0, ...)
        abline(v = sign.lambda * log(x$lambda.min), lty=2, lwd = 2, col = "firebrick3")
        abline(v = sign.lambda * log(x$lambda.1se), lty=2, lwd = 2, col = "firebrick1")
        title(x$name, line = 2.5, ...)
    }
    
}
