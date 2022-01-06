
#' Plot method for cv.groupFusedMulti fitted objects
#'
#' @param sign.lambda Either plot against log(lambda) (default) or its negative if \code{sign.lambda = -1}.
#' @param best.fused if \code{TRUE}, make a CV plot only for the best \code{lambda.fused} value. This includes
#' standard errors in the plot. If \code{FALSE}, plot CV for each value of \code{lambda.fused}. In this case,
#' no standard errors will be plotted.
#' @param plot.method either \code{'lines'} in which case a different line of the CV error versus the tuning parameter will be plotted
#' for each fused lasso tuning parameter, or \code{'heatmap'} in which case a 2 by 2 heatmap of the errors will be plotted
#' with each tuning parameter on one axis
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
#' 
#' ## examples for plot cv.groupFusedMulti
#' 
#' set.seed(123)
#'
#' dat.sim <- gen_sparse_multivar_data(nvars = 5L,
#'                    noutcomes = 8L,
#'                    nobs = 500L,
#'                    nobs.test = 500L,
#'                    num.nonzero.vars = 5,
#'                    outcome.groups = rbind(c(1,1,1,2,2,2,2,2),
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
#' outcome_groups <- rbind(c(1,1,1,2,2,2,2,2),
#'                         c(1,1,1,2,2,3,3,3),
#'                         c(1,1,2,3,3,4,4,5))
#'                         
#' fit.adapt <- cv.groupFusedMulti(x, y,
#'                                 nlambda        = 25,
#'                                 lambda.fused = c(0.00001, 0.0001, 0.001),
#'                                 outcome.groups = outcome_groups,
#'                                 gamma          = 0.25,
#'                                 nfolds         = 3)
#'                                 
#' plot(fit.adapt)
#' 
#' plot(fit.adapt, plot.method = "heatmap")
#' 
plot.cv.groupFusedMulti <- function(x, sign.lambda = 1, best.fused = FALSE, plot.method = c("lines", "heatmap", "min"), ...) 
{
    # compute total number of selected variables for each
    # tuning parameter 
    
    plot.method <- match.arg(plot.method)
    
    if (!is.null(x$lambda.fused))
    {
        if (plot.method == "lines")
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
                    image(1, x$lambda.fused, t(seq_along(x$lambda.fused)), col=cols, axes=FALSE,
                          main = expression(alpha))
                    axis(4, labels = round(x$lambda.fused, 5), at = (x$lambda.fused))
                } else
                {
                    image(1, log(x$lambda.fused), t(seq_along(x$lambda.fused)), col=cols, axes=FALSE,
                          main = expression(lambda~fused))
                    axis(4, labels = round(x$lambda.fused, 5), at = log(x$lambda.fused))
                }
                
                par(old.par)
            }
        } else if (plot.method == "heatmap")
        {
            cv_mat <- do.call(rbind, x$cvm)
            rn <- round(x$lambda.fused, 5) #gsub("[^0-9\\.]", "", rownames(mp_cv$cv_error))
            cn <- round(x$lambda, 4) #gsub("[^0-9\\.]", "", colnames(mp_cv$cv_error))
            
            # topo.colors(250)
            if (x$groupFusedMulti.fit$use.alpha.param)
            {
                xlab <- expression(lambda)
                ylab <- expression(alpha)
            } else
            {
                xlab <- expression(lambda[1])
                ylab <- expression(lambda[2])
            }
            image(as(cv_mat, "Matrix"), col.regions = vcols, colorkey = TRUE,
                  xlab = xlab, ylab = ylab, scales = list(y = list(labels = rn, at = 1:length(rn)),
                                                          x = list(labels = cn, at = 1:length(cn),
                                                                   rot = 45)),
                  ylab.right = "CV MSE",
                  sub = NULL,
                  ...)
        } else if (plot.method == "min")
        {
            
            if (x$groupFusedMulti.fit$use.alpha.param)
            {
                ylab <- expression("Min CV Error across"~~lambda)
                xlab <- expression(log(alpha))
            } else
            {
                ylab <- expression("Min CV Error across"~~lambda[1])
                xlab <- expression(log(lambda[2]))
            }
            
            min_cves <- sapply(x$cvm, min)
            which_min_cves <- sapply(x$cvm, which.min)
            min_cves_lo <- mapply(function(i,cvl) cvl[i], which_min_cves, x$cvlo)
            min_cves_hi <- mapply(function(i,cvl) cvl[i], which_min_cves, x$cvup)
            
            
            
            plot.args = list(x    = log(x$lambda.fused),
                             y    = min_cves,
                             ylim = range(c(min_cves_lo, min_cves_hi)),
                             xlab = xlab,
                             ylab = ylab,
                             type = "c", col = "black", xaxt = "n")
            
            new.args <- list(...)
            if(length(new.args)) plot.args[names(new.args)] <- new.args
            
            do.call("plot", plot.args)
            
            ## col is scales::alpha("firebrick3", 0.4) = "#CD262666"
            abline(v = log(x$lambda.fused.min), lty = 2, lwd = 1, col = "#CD262666")
            
            axis(1, at = log(x$lambda.fused), labels = x$lambda.fused)
            
            error.bars(log(x$lambda.fused), 
                       min_cves_hi, 
                       min_cves_lo)
            
            points(log(x$lambda.fused), min_cves, pch = 19, col="dodgerblue")
            
            
            
            # plot(min_cves, type = "b", ylim = range(c(min_cves_lo, min_cves_hi)),
            #      ylab = ylab, xlab = xlab)
            #    
            # error.bars(1:length(which_min_cves), min_cves_hi, min_cves_lo)
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
