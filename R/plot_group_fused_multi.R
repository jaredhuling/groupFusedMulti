
#' Plotting method for groupFusedMulti fitted objects
#' 
#' @param x fitted \code{groupFusedMulti} or \code{cv.groupFusedMulti} model object
#' @param plot.type Either \code{"all_variables"}, where all the coefficients for all variables are plotted for a chosen 
#' outcome (given by \code{which.outcome}]) are plotted as a function of the tuning parameter, or \code{"all_outcomes"},
#'  where the coefficients for a given variable 
#' (chosen by \code{which.variable}) are plotted for all outcomes as a function of the tuning parameter
#' @param lam.fused.idx which lambda to use? specify integer between 1 and the number of tuning parameters used for fused lasso
#' @param which.outcome which outcome's coefficients should be plotted? 
#' @param which.variable which outcome's coefficients should be plotted? 
#' @param which.outcome.group if there are multiple outcome groupings, which of the groups should be used? Only used if multiple outcome groupings are present
#' @param xvar What is on the X-axis. "norm" plots against the L1-norm of the coefficients, "lambda" against the log-lambda sequence, and "dev"
#' against the percent deviance explained.
#' @param xlab character value supplied for x-axis label
#' @param ylab character value supplied for y-axis label
#' @param main text for plot title
#' @param labsize label size
#' @param plot.only.nonzero plot path only for values of the tuning parameter for which not all coefficients are zero?
#' @param vline.val value for the x-axis for which to include a vertical line. Defaults to no vertical line. Can be used to track tuning parameter
#' value that was chosen by cross-validation
#' @param base.size base font size for any plots constructed with \code{ggplot2} (i.e. only used when \code{plot.type='all_outcomes'}). Defaults
#' to 12
#' @param ... other graphical parameters for the plot
#' @rdname plot
#' @importFrom graphics matplot mtext strwidth
#' @importFrom grDevices rainbow
#' @importFrom stats na.omit rbinom reshape sd
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid align_plots
#' @importFrom ggplot2 ggplot geom_line ggtitle theme_bw element_text coord_cartesian
#' @importFrom ggplot2 geom_vline scale_x_continuous scale_y_continuous margin
#' @export
#' @examples
#' library(Matrix)
#' 
#' set.seed(123)
#' n.obs <- 200
#' n.vars <- 50
#'
#' true.beta.mat <- array(NA, dim = c(3, n.vars))
#' true.beta.mat[1,] <- c(-0.5, -1, 0, 0, 2, rep(0, n.vars - 5))
#' true.beta.mat[2,] <- c(0.5, 0.5, -0.5, -0.5, 1, -1, rep(0, n.vars - 6))
#' true.beta.mat[3,] <- c(0, 0, 1, 1, -1, rep(0, n.vars - 5))
#' rownames(true.beta.mat) <- c("1,0", "1,1", "0,1")
#' true.beta <- as.vector(t(true.beta.mat))
#'
#' x.sub1 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub2 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub3 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#' x <- as.matrix(rbind(x.sub1, x.sub2, x.sub3))
#'
#' conditions <- as.matrix(cbind(c(rep(1, 2 * n.obs), rep(0, n.obs)),
#'                               c(rep(0, n.obs), rep(1, 2 * n.obs))))
#'
#' y <- rnorm(n.obs * 3, sd = 3) + drop(as.matrix(bdiag(x.sub1, x.sub2, x.sub3)) %*% true.beta)
#'
#' fit <- vennLasso(x = x, y = y, groups = conditions)
#'
#' layout(matrix(1:3, ncol = 3))
#' plot(fit, which.subpop = 1)
#' plot(fit, which.subpop = 2)
#' plot(fit, which.subpop = 3)
#'
plot.groupFusedMulti <- function(x, 
                                 plot.type = c("all_variables", "all_outcomes"),
                                 lam.fused.idx = NULL,
                                 which.outcome = 1,
                                 which.variable = 1,
                                 which.outcome.group = 1,
                                 xvar = c("loglambda", "norm", "lambda", "dev"),
                                 xlab = iname, ylab = "Coefficients",
                                 main = NULL,
                                 labsize = 1,
                                 plot.only.nonzero = TRUE,
                                 vline.val = NULL,
                                 base.size = 12,
                                 ...)
{
    plot_type <- match.arg(plot.type)
    if (is.list(x$beta))
    {
        noutcomes <- dim(x$beta[[1]])[2]
        nvars <- dim(x$beta[[1]])[1]
        
        if (is.null(lam.fused.idx))
        {
            stop("if fused lasso was used, must specify which fused tuning parameter to use via 'lam.fused.idx'")
        }
        
        if (length(lam.fused.idx) > length(x$beta))
        {
            stop("'lam.fused.idx' provided was larger than the number of fused lasso tuning parameters fit")
        }
        
    } else
    {
        noutcomes <- dim(x$beta)[2]
        nvars <- dim(x$beta)[1]
    }
    if (any(which.outcome > noutcomes))
    {
        err.txt <- paste0("Outcome ", which.outcome, " specified, but only ",
                          noutcomes,
                          " outcomes present in the model.")
        stop(err.txt)
    }
    if (any(which.variable > nvars))
    {
        err.txt <- paste0("Variable ", which.variable, " specified, but only ",
                          nvars,
                          " variables present in the model.")
        stop(err.txt)
    }
    xvar <- match.arg(xvar)
    
    
    if (is.matrix(x$outcome.groups))
    {
        if (which.outcome.group < 1 | which.outcome.group > nrow(x$outcome.groups))
        {
            stop(paste0("Invalid value for 'which.outcome.group', must be between 1 and ", nrow(x$outcome.groups), "."))
        }
    }
    
    
    Outcome_Groups <- NULL
    
    if (plot_type == "all_variables")
    {
        if (is.list(x$beta))
        {
            nbeta <- as.matrix(x$beta[[lam.fused.idx]][-1,which.outcome,])
        } else
        {
            nbeta <- as.matrix(x$beta[-1,which.outcome,])
        }
    } else if (plot_type == "all_outcomes")
    {
        if (is.list(x$beta))
        {
            nbeta <- as.matrix(x$beta[[lam.fused.idx]][-1,,][which.variable,,])
        } else
        {
            nbeta <- as.matrix(x$beta[-1,,][which.variable,,])
        }
    }
    
    ## missing means 0 (ie no variation in that subpopulation for that variable)
    nbeta[is.na(nbeta)] <- 0
    
    if (plot.only.nonzero)
    {
        which_nz <- colSums(abs(nbeta)) >= 1e-6
        
        if (any(which_nz))
        {
            which_nz[1] <- FALSE
            
            if (any(!which_nz))
            {
                whichnz <- which(!which_nz)
                which_nz[whichnz[length(whichnz)]] <- TRUE
            }
            
            nbeta <- nbeta[,which_nz]
        }
    }
    
    
    switch(xvar,
           "norm" = {
               index <- colSums(abs(nbeta), na.rm = TRUE)
               iname <- "L1 Norm"
               xlim <- range(index)
           },
           "lambda" = {
               index <- if(is.matrix(x$lambda)){x$lambda[,lam.fused.idx]} else {x$lambda}
               if (plot.only.nonzero) index <- index[which_nz]
               iname <- expression(lambda)
               xlim <- rev(range(index))
           },
           "loglambda" = {
               index <- if(is.matrix(x$lambda)){log(x$lambda[,lam.fused.idx])} else {log(x$lambda)}
               if (plot.only.nonzero) index <- index[which_nz]
               iname <- expression(log(lambda))
               xlim <- rev(range(index))
           },
           "dev" = {
               index = x$sumSquare
               if (plot.only.nonzero) index <- index[which_nz]
               iname = "Sum of Squares"
               xlim <- range(index)
           }
    )
    
    
    if (plot_type == "all_variables")
    {
        if (is.list(x$beta))
        {
            rn <- dimnames(x$beta[[1]])[[2]]
            if (!is.null(rn))
            {
                main.txt <- paste0(rn[which.outcome], ", LamFused=", round(x$lambda.fused[[lam.fused.idx]], 5))
            } else
            {
                main.txt <- paste0("LamFused=", round(x$lambda.fused[[lam.fused.idx]], 5))
            }
        } else
        {
            rn <- dimnames(x$beta)[[2]]
            if (!is.null(rn))
            {
                main.txt <- rn[which.outcome]
            } else
            {
                main.txt <- ""
            }
        }
    } else if (plot_type == "all_outcomes")
    {
        if (is.list(x$beta))
        {
            rn <- dimnames(x$beta[[1]])[[1]]
            if (!is.null(rn))
            {
                lamval   <- round(x$lambda.fused[[lam.fused.idx]], 5)
                main.txt <- paste0(rn[which.variable], ", LamFused=", lamval)
                vnm <- rn[which.variable + 1]
                main.txt <- bquote(.(vnm)~", "~ lambda[Fused]==.(lamval))
            } else
            {
                lamval   <- round(x$lambda.fused[[lam.fused.idx]], 5)
                main.txt <- paste0("LamFused=", round(x$lambda.fused[[lam.fused.idx]], 5))
                main.txt <- bquote(lambda[Fused]==.(lamval))
            }
        } else
        {
            rn <- dimnames(x$beta)[[1]]
            if (!is.null(rn))
            {
                main.txt <- rn[which.variable]
            } else
            {
                main.txt <- ""
            }
        }
    }
    
    
    if (is.null(main))
    {
        main_txt <- main.txt
    } else
    {
        main_txt <- main
    }
    
    
    if (is.matrix(x$outcome.groups))
    {
        outcome.groups <- x$outcome.groups[which.outcome.group,]
    } else
    {
        outcome.groups <- x$outcome.groups
    }
    
    unique.outcome.groups <- unique(outcome.groups)
    
    if (length(unique.outcome.groups) > 1 & plot_type == "all_outcomes")
    {
        cols <- rainbow(length(unique.outcome.groups))
        colseq <- character(length(outcome.groups))
        
        for (og in 1:length(unique.outcome.groups))
        {
            colseq[outcome.groups == unique.outcome.groups[og]] <- cols[og]
        }
    } else
    {
        cols <- rainbow(nrow(nbeta))
        
        ## create sequence that grabs one of ROYGBIV and repeats with
        ## an increment up the rainbow spectrum with each step from 1:7 on ROYGBIV
        n.cols <- 7L
        scramble.seq <- rep(((1:n.cols) - 1) * (length(cols) %/% (n.cols)) + 1, length(cols) %/% n.cols)[1:length(cols)] + 
            (((0:(length(cols)-1)) %/% n.cols))
        
        scramble.seq[is.na(scramble.seq)] <- which(!(1:length(cols) %in% scramble.seq))
        colseq <- cols[scramble.seq]
    }
    
    if (plot_type == "all_variables")
    {
        # Adjust the margins to make sure the labels fit
        labwidth <- ifelse(labsize > 0, max(strwidth(rownames(nbeta), "inches", labsize)), 0)
        margins <- par("mai")
        par("mai" = c(margins[1:3], max(margins[4], labwidth*1.55)))
        
        matplot(index, t(nbeta), lty = 1,
                xlab = xlab, ylab = "",
                xlim = xlim, main = main_txt,
                col = colseq,
                type = 'l', ...)
        
        
        if (is.null(ylab)) 
        {
            mtext(expression(hat(beta)), side = 2, cex = par("cex"), line = 3, las = 1)
        } else 
        {
            mtext(ylab, side = 2, cex = par("cex"), line = 3)
            ylab = ""
        }
        
        if (!is.null(vline.val))
        {
            abline(v = vline.val, lty = "dashed")
        }
        
        #atdf <- pretty(index, n = 10L)
        #plotnz <- approx(x = index, y = x$nzero[[which.model]], xout = atdf, rule = 2, method = "constant", f = approx.f)$y
        #axis(side=3, at = atdf, labels = plotnz, tick=FALSE, line=0, ...)
        
        #title(main, line = 2.5, ...)
        
        
        dotlist <- list(...)
        
        if ("cex.axis" %in% names(dotlist))
        {
            dotlist <- dotlist[-which(names(dotlist) == "cex.axis")]
        }
        
        if ( labsize > 0 && !is.null(rownames(nbeta)) ) 
        {
            for (i in 1:nrow(nbeta)) {
                j <- i
                do.call(axis, c(list(side = 4, at = nbeta[i, ncol(nbeta)], labels = rownames(nbeta)[i],
                                     las = 1, cex.axis = labsize, col.axis = colseq[i],
                                     cex.axis = 0.75,
                                     lty = 1, col = colseq[i]), dotlist))
            }
            # axis(4, at = nbeta[, ncol(nbeta)], labels = rownames(nbeta),
            #      las = 1, cex.axis=labsize, #col.axis = colseq,
            #      lty = 1, col = colseq,  ...)
        }
        par("mai"=margins)
    } else if (plot_type == "all_outcomes")
    {
        coefsplot <- as.data.frame(nbeta)
        coefsplot$Outcome <- rownames(coefsplot)
        coefsplot$Outcome_Groups <- as.factor(outcome.groups)
        lamcols <- head(colnames(coefsplot),-2)
        coefsplot_tall <- reshape(as.data.frame(coefsplot), direction = "long", varying = list(lamcols))#, times = as.numeric(lamcols))
        colnames(coefsplot_tall)[4] <- "Beta"
        
        coefsplot_tall$time <- as.numeric(coefsplot_tall$time)
        index
        
        coefsplot_tall$lambda <- index[coefsplot_tall$time]
        
        coefsplot_label_dat <- coefsplot_tall[coefsplot_tall$time == max(coefsplot_tall$time),]
        
        cfplt <- ggplot(data = coefsplot_tall, 
                        aes(x = lambda, y = Beta, group = Outcome, color = Outcome_Groups)) +
            geom_line() +
            xlab(xlab) +
            ggtitle(main_txt) + 
            ylab(expression(hat(beta))) + 
            theme_bw(base_size = base.size) +
            theme(axis.title.y = element_text(angle = 0, vjust=0.5),
                  legend.position = "bottom")
        
        if (xvar %in% c("lambda", "loglambda"))
        {
            cfplt <- cfplt + coord_cartesian(xlim = rev(range(coefsplot_tall$lambda)))
        }
        
        if (!is.null(vline.val))
        {
            cfplt <- cfplt + geom_vline(xintercept = vline.val, linetype = "dashed")
        }
        
        
        # make the axis plot
        axispl <- ggplot(coefsplot_label_dat,
                         aes(x = 0, y = Beta, label = Outcome, color = Outcome_Groups)) +
            geom_text_repel(min.segment.length = grid::unit(0, "pt"),
                            max.overlaps = nrow(coefsplot_label_dat)*4,
                            force = 3, force_pull = 0.5,
                            size = 0.8*11/.pt  ## ggplot2 theme_grey() axis text
            ) +
            coord_cartesian(xlim = c(0,0.001)) +
            scale_x_continuous(limits = c(0, 0.05), expand = c(0, 0),
                               breaks = NULL, labels = NULL, name = NULL) +
            scale_y_continuous(limits = range(coefsplot_tall$Beta), #expand = c(0, 0),
                               breaks = NULL, labels = NULL, name = NULL) +
            theme(panel.background = element_blank(),
                  legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "pt"))
        
        return(plot_grid(plotlist = align_plots(cfplt, axispl, align = "h", axis = "tb"), nrow = 1, rel_widths = c(5, 1.25)))
    }
        
    
}
