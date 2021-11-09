
##### Zscores-qvals

zscores_qvals=function(x,outlier.method="mahalanobis",estim="pairwiseGK",pval.coret.method="bonferroni"){

  crqvalues<-function(scores,K,outlier.method="mahalanobis",estim="pairwiseGK",pval.coret.method="bonferroni")
  {
    #  normalize <- function(x) {
    #   return ((x - min(x)) / (max(x) - min(x)))
    # }
    # scoresnorm <- apply(scores, 2, normalize)

    ##### now this is updated based on Priv?, F et al, 2020 PCADAPT.
    if (outlier.method == "mahalanobis") {
      res <- rep(NA, nrow(scores))
      if (K == 1) {
        res <- (scores - median(scores))^2
      } else if (K > 1) {
        res <- robust::covRob(scores, distance = TRUE, na.action= na.omit, estim=estim)$dist  ####A matrix with no missing values and at least 2 columns
      }
      gif <- median(res, na.rm = TRUE) /stats::qchisq(0.5, df = K)
      reschi2test <- as.numeric(stats::pchisq(res/gif, df = K, lower.tail = FALSE))
      qval <- qvalue::qvalue(reschi2test)
      padj <- stats::p.adjust(reschi2test,method=pval.coret.method)
    }
    else if (outlier.method == "componentwise") {
      res <- apply(scores, 2, FUN = function(h) {h^2})
      gif <- sapply(1:K, FUN = function(h) {
        median(scores[, h]^2, na.rm = TRUE)/stats::qchisq(0.5, df = 1)
      })
      reschi2test <- NULL
      for (k in 1:K) {
        reschi2test <- cbind(reschi2test, stats::pchisq((res/gif)[, k], df = 1, lower.tail = FALSE))
      }
      qval <- qvalue::qvalue(reschi2test)
      padj <- stats::p.adjust(reschi2test,method=pval.coret.method)
    }
    structure(list(pvals=data.frame(p.values=reschi2test, q.values=qval$qvalues,padj=padj), maha = res, gif = gif, K = K), outlier.method=outlier.method,pval.coret.method=outlier.method, class = "GenomeAdapt")

  }

  snpcorr=as.data.frame(t(x$zscores$snpcorr))
  colnames(snpcorr)=x$zscores$sample.id
  rownames(snpcorr)=x$zscores$snp.id
  pvals=crqvalues(snpcorr,K=ncol(snpcorr),outlier.method=outlier.method,estim=estim,pval.coret.method=pval.coret.method)
  chr=x$chr
  return(list(pvals=pvals,chr=chr,class="zscores_qvals"))
}

### plot qvalues

plotmanhattan=function(x,ylim=c(0,200),xlab="",ylab="-log(p-value)",col = x$chr, pch="*",h=10, lcol="blue",...){
  graphics::plot(-log10(x), ylim=ylim, xlab=xlab,ylab = ylab , col=col, pch="*")
  abline(h=h, col=lcol)
}
