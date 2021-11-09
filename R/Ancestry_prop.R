##### estimate the anscrtry proportion and plot them
#### Credit to Zheng Xiuwen

AdmixProp <- function(x, groups, bound=FALSE)
{
  # check
  eigobj=x$eig

  # 'sample.id' and 'eigenvect' should exist
  stopifnot(!is.null(x$zscores$sample.id))
  stopifnot(is.matrix(x$eig$vectors))

  stopifnot(is.list(groups))
  stopifnot(length(groups) > 1)
  if (length(groups) > (ncol(eigobj$vectors)+1))
  {
    stop("`eig' should have more eigenvectors than ",
         "what is specified in `groups'.")
  }

  grlist <- NULL
  for (i in 1:length(groups))
  {
    if (!is.vector(groups[[i]]) & !is.factor(groups[[i]]))
    {
      stop(
        "`groups' should be a list of sample IDs ",
        "with respect to multiple groups."
      )
    }
    if (any(!(groups[[i]] %in% x$zscores$sample.id)))
    {
      stop(sprintf(
        "`groups[[%d]]' includes sample(s) not existing ",
        "in `sample.id'.", i))
    }

    if (any(groups[[i]] %in% grlist))
      warning("There are some overlapping between group sample IDs.")
    grlist <- c(grlist, groups[[i]])
  }

  stopifnot(is.logical(bound) & is.vector(bound))
  stopifnot(length(bound) == 1)

  # calculate ...

  E <- eigobj$vectors[, 1:(length(groups)-1)]
  if (is.vector(E)) E <- matrix(E, ncol=1)
  mat <- NULL
  for (i in 1:length(groups))
  {
    k <- match(groups[[i]], x$zscores$sample.id)
    Ek <- E[k, ]
    if (is.vector(Ek))
      mat <- rbind(mat, mean(Ek))
    else
      mat <- rbind(mat, colMeans(Ek))
  }

  # check
  if (any(is.na(mat)))
    stop("The eigenvectors should not have missing value!")

  T.P <- mat[length(groups), ]
  T.R <- solve(mat[-length(groups), ] -
                 matrix(T.P, nrow=length(T.P), ncol=length(T.P), byrow=TRUE))

  new.p <- (E - matrix(T.P, nrow=dim(E)[1], ncol=length(T.P),
                       byrow=TRUE)) %*% T.R
  new.p <- cbind(new.p, 1 - rowSums(new.p))
  colnames(new.p) <- names(groups)
  rownames(new.p) <- x$zscores$sample.id

  # whether bounded
  if (bound)
  {
    new.p[new.p < 0] <- 0
    r <- 1.0 / rowSums(new.p)
    new.p <- new.p * r
  }

  new.p
}



PlotAdmix<- function(propmat, group=NULL, col=NULL, multiplot=TRUE,xlab="Individuals", ylab="Ancestry Proportion",
                     showgrp=TRUE, shownum=TRUE, ylim=TRUE, na.rm=TRUE)
{
  # check
  stopifnot(is.numeric(propmat), is.matrix(propmat))
  stopifnot(is.null(group) | is.vector(group) | is.factor(group))
  if (!is.null(group))
    stopifnot(nrow(propmat) == length(group))
  stopifnot(is.null(col) | is.vector(col))
  stopifnot(is.logical(multiplot), length(multiplot)==1L)
  stopifnot(is.logical(showgrp), length(showgrp)==1L)
  stopifnot(is.logical(shownum), length(shownum)==1L)
  stopifnot(is.logical(ylim) | is.numeric(ylim))
  if (is.numeric(ylim))
    stopifnot(length(ylim) == 2L)
  stopifnot(is.logical(na.rm), length(na.rm)==1L)

  if (is.logical(ylim))
  {
    if (isTRUE(ylim))
      ylim <- c(0, 1)
    else
      ylim <- range(propmat, na.rm=TRUE)
  }

  if (!is.null(group))
  {
    if (anyNA(group) & !isTRUE(na.rm))
    {
      group <- as.character(group)
      group[is.na(group)] <- "<NA>"
    }
    grp_name <- levels(factor(group))
    idx <- list()
    for (n in grp_name)
    {
      i <- which(group == n)
      i <- i[order(propmat[i, 1L], decreasing=TRUE)]
      idx <- c(idx, list(i))
    }
    propmat <- propmat[unlist(idx), ]
    grp_len <- lengths(idx, use.names=FALSE)
    xl <- c(0, cumsum(grp_len))
    x <- xl[-1L] - 0.5*grp_len
  }

  if (multiplot)
  {
    opar <- graphics::par(mfrow=c(ncol(propmat), 1L), mar=c(1.25, 5, 1.75, 2),
                oma=c(0, 0, 4, 0))
    on.exit(graphics::par(opar))
    grp <- colnames(propmat)
    if (is.null(grp))
      grp <- paste("group", seq_len(ncol(propmat)))
    ylab <- paste("Prop. of", grp)

    for (i in seq_len(ncol(propmat)))
    {
      barplot(unname(propmat[, i]), space=0, border=NA, ylab=ylab[i],
              ylim=ylim, col=col)
      lines(c(1, nrow(propmat)), c(0, 0))
      lines(c(1, nrow(propmat)), c(1, 1))
      if (!is.null(group))
      {
        abline(v=xl, col="blue")
        if (showgrp)
          text(x, 0.5, labels=grp_name, srt=30)
        if ((i == 1L) & shownum)
        {
          axis(1, c(0, x, nrow(propmat)),
               c("", as.character(lengths(idx)), ""), cex.axis=0.75)
        }
      }
    }
  } else {
    if (is.null(col)) col <- rainbow(ncol(propmat))
    barplot(t(unname(propmat)), col=col,
            xlab=xlab, ylab=ylab, space=0, border=NA)
    if (!is.null(group))
    {
      abline(v=xl, col="black")
      if (showgrp)
        text(x, 0.5, labels=grp_name, srt=30)
      if (shownum)
      {
        axis(1, c(0, x, nrow(propmat)),
             c("", as.character(lengths(idx)), ""), cex.axis=0.75)
      }
    }
  }
  invisible()
}
