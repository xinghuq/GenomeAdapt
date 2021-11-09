\name{plotmanhattan}
\alias{plotmanhattan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Making manhattan plot
%%  ~~function to do ... ~~
}
\description{Ploting a Manhattan plot showing the (q)p-values for each SNP.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

plotmanhattan(x, ylim=c(0,200),xlab="",
ylab="-log(p-value)",col = x$chr, pch="*",h=10, lcol="blue", ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The zscores_qvals object
%%     ~~Describe \code{x} here~~
}
  \item{ylim}{ The y axis range
%%     ~~Describe \code{ylim} here~~
}
  \item{xlab}{The x title
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{The y title
%%     ~~Describe \code{xlab} here~~
}
  \item{col}{The colour used for indicating chromosomes
%%     ~~Describe \code{col} here~~
}
  \item{pch}{The shape of of the data points.
%%     ~~Describe \code{pch} here~~
}
  \item{h}{ The cutoff line
%%     ~~Describe \code{h} here~~
}
  \item{lcol}{The colour of the cutoff line
%%     ~~Describe \code{lcol} here~~
}
  \item{...}{passing to other parameters of plot
%%     ~~Describe \code{\dots} here~~
}
}
\examples{
### using Hapmap data
HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),
method="EIGMIX",num.thread = 4, autosome.only=TRUE,
remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

### Not running, it takes a while to finish this

# Hapmapqval=zscores_qvals(HapmapScan)
## plot
# plotmanhattan(Hapmapqval$pvals$pvals$p.values,col=Hapmapqval$chr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}% use one of  RShowDoc("KEYWORDS")

