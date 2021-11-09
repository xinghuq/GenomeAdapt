\name{PlotAdmix}
\alias{PlotAdmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Ancestry Admixture
%%  ~~function to do ... ~~
}
\description{Plot admixture proportions based on their ancestries.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotAdmix(propmat, group = NULL, col = NULL, multiplot = TRUE, showgrp = TRUE, shownum = TRUE, ylim = TRUE, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{propmat}{a sample-by-ancestry matrix of proportion estimates, returned from AdmixProp()
%%     ~~Describe \code{propmat} here~~
}
  \item{group}{a character vector of a factor according to the samples in propmat
%%     ~~Describe \code{group} here~~
}
  \item{col}{specify colors
%%     ~~Describe \code{col} here~~
}
  \item{multiplot}{single plot or multiple plots
%%     ~~Describe \code{multiplot} here~~
}
  \item{showgrp}{show group names in the plot
%%     ~~Describe \code{showgrp} here~~
}
  \item{shownum}{TRUE: show the number of each group in the figure
%%     ~~Describe \code{shownum} here~~
}
  \item{ylim}{TRUE: y-axis is limited to [0, 1]; FALSE: ylim <- range(propmat); a 2-length numeric vector: ylim used in plot()
%%     ~~Describe \code{ylim} here~~
}
  \item{na.rm}{TRUE: remove the sample(s) according to the missing value(s) in group
%%     ~~Describe \code{na.rm} here~~
}
}
\details{The minor allele frequency and missing rate for each SNP passed in snp.id are calculated over all the samples in sample.id.
%%  ~~ If necessary, more details than the description above ~~
}

\references{
Zheng X, Weir BS. Eigenanalysis on SNP Data with an Interpretation of Identity by Descent. Theoretical Population Biology. 2015 Oct 23. pii: S0040-5809(15)00089-1. doi: 10.1016/j.tpb.2015.09.004.

Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., & Weir, B. S. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24), 3326-3328.
%% ~put references to the literature/web site here ~
}

\examples{
##---- Do genome scan and get the ancestry proportation for individuals ----
# Scan genomes (HapMap)

HapmapScan=GenomeAdapt.gds(genfile = snpgdsExampleFileName(),method="EIGMIX",num.thread = 6, autosome.only=TRUE, remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)


# get population information
genofile <- SNPRelate::snpgdsOpen(snpgdsExampleFileName())
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snpgdsClose(genofile)

# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

### estimate the ancestry proportion

Admixpro=AdmixProp(HapmapScan,groups=groups,bound=TRUE)

PlotAdmix(Admixpro,group=as.factor(pop_code),multiplot = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}% use one of  RShowDoc("KEYWORDS")

