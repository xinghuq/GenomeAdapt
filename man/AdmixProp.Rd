\name{AdmixProp}
\alias{AdmixProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate ancestral proportions from the eigen vectors
%%  ~~function to do ... ~~
}
\description{Estimate ancestral (admixture) proportions based on the eigen-analysis
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AdmixProp(x, groups, bound = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object from GenomeAdapt
%%     ~~Describe \code{x} here~~
}
  \item{groups}{A list of sample IDs, such like groups = list( CEU = c("NA0101", "NA1022", ...), YRI = c("NAxxxx", ...), Asia = c("NA1234", ...))
%%     ~~Describe \code{groups} here~~
}
  \item{bound}{if TRUE, the estimates are bounded so that no component < 0 or > 1, and the sum of proportions is one
%%     ~~Describe \code{bound} here~~
}
}
\details{ Estimate ancestral (admixture) proportions based on the eigen-analysis, results give a dataframe of individual probability from the ancestries.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a dataframe of individual probability from the ancestries.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\keyword{AdmixProp}% use one of  RShowDoc("KEYWORDS")

