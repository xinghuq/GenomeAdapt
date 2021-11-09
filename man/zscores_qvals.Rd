\name{zscores_qvals}
\alias{zscores_qvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating (q)p-values from multiple factors/scores
%%  ~~function to do ... ~~
}
\description{ Converting the Z-scores to (q)p-values. This function calibrates the p-value/q-values considering multiple scores based on Mahalanobis Distance/componentwise method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zscores_qvals(x, outlier.method = "mahalanobis", estim = "pairwiseGK", pval.coret.method = "bonferroni")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The GenomeAdapt object containing the locus scores
%%     ~~Describe \code{scores} here~~
}

  \item{outlier.method}{ The methods used to detect the outliers, "mahalanobis" or "componentwise", default is mahalanobis.
%%     ~~Describe \code{outlier.method} here~~
}
  \item{estim}{ Method used to estimate Mahalanobis distance.The choices are: "mcd" for the Fast MCD algorithm of Rousseeuw and Van Driessen, "weighted" for the Reweighted MCD, "donostah" for the Donoho-Stahel projection based estimator, "M" for the constrained M estimator provided by Rocke, "pairwiseQC" for the orthogonalized quadrant correlation pairwise estimator, and "pairwiseGK" for the Orthogonalized Gnanadesikan-Kettenring pairwise estimator. The default "auto" selects from "donostah", "mcd", and "pairwiseQC" with the goal of producing a good estimate in a reasonable amount of time.
%%     ~~Describe \code{estim} here~~
}
  \item{pval.coret.method}{ Correction method for p-values,choices are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none".The adjustment methods include the Bonferroni correction ("bonferroni") in which the p-values are multiplied by the number of comparisons. Less conservative corrections are also included by Holm (1979) ("holm"), Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"), Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), and Benjamini & Yekutieli (2001) ("BY"), respectively. A pass-through option ("none") is also included. 
%%     ~~Describe \code{pval.coret.method} here~~
}
}
\details{Calculating (q)p-values from GenomeAdapt to identify the outlier loci
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A dataframe with p-values, adjusted p-values and q-values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. http://www.jstor.org/stable/2346101.

R. A. Maronna and R. H. Zamar (2002) Robust estimates of location and dispersion of high-dimensional datasets. Technometrics 44 (4), 307–317.

Capblancq, T., Luu, K., Blum, M. G., & Bazin, E. (2018). Evaluation of redundancy analysis to identify signatures of local adaptation. Molecular Ecology Resources, 18(6), 1223-1233.
%% ~put references to the literature/web site here ~
}
\examples{
##---- Not running as it takes a while ----
# HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),method="EIGMIX",num.thread = 6, autosome.only=TRUE, remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

# Hapmapqval=zscores_qvals(HapmapScan)

# plot.manhattan(Hapmapqval$pvals$pvals$p.values,col=Hapmanqval$chr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ p-values }% use one of  RShowDoc("KEYWORDS")