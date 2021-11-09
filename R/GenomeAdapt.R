### two functions, the conventional linear coefficient method, and the neural network methods

"GenomeAdapt"=function(genfile,method = "EIGMIX", sample.id=NULL, snp.id=NULL,
                       autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                       num.thread=1L, out.fn=NULL, out.prec=c("double", "single"),
                       out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE, ...){
    UseMethod("DeepGenomeScan")
  }


#  GenomeAdapt=function(bedfile=NULL, vcffile=NULL, genofile=NULL,method="EIGMIX",autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,...)
#   {
#  if (is.null(vcffile)==TRUE & is.null(genofile)==TRUE & is.null(bedfile)==FALSE)
#  {SNPRelate::snpgdsBED2GDS(paste0("bedfile","bed"),paste0("bedfile","fam"), paste0("bedfile","bim"), "inputgenofile.gds")
#    genf=SNPRelate::snpgdsOpen("inputgenofile.gds")}
#  if (is.null(bedfile)==TRUE & is.null(genofile)==TRUE & is.null(vcffile)==FALSE)
 # {SNPRelate::snpgdsVCF2GDS(vcffile, "inputgenofile.gds", method="biallelic.only")
#    genf=SNPRelate::snpgdsOpen("inputgenofile.gds")}
#  if (is.null(bedfile)==TRUE & is.null(vcffile)==TRUE & is.null(genofile)==FALSE)
#  {genf=genofile}
 #
#  rv <- SNPRelate::snpgdsGRM(genf, method=method,...) ### "GCTA" - genetic relationship matrix defined in CGTA; "Eigenstrat" - genetic covariance matrix in EIGENSTRAT; "EIGMIX" - two times coancestry matrix defined in Zheng & Weir (2015), "Weighted" - weighted GCTA, as the same as "EIGMIX", "Corr" - Scaled GCTA GRM (dividing each i,j element by the product of the square root of the i,i and j,j elements), "IndivBeta" - two times individual beta estimate relative to the minimum of beta; see details
#  eig <- eigen(rv$grm)  # Eigen-decomposition
#  # close the file
 # chr <- gdsfmt::read.gdsn(index.gdsn(genf, "snp.chromosome"))
#  sampleid=gdsfmt::read.gdsn((index.gdsn(genf,"sample.id")))
#  rownames(eig$vectors)=sampleid
#  zscores <- SNPRelate::snpgdsPCACorr(eig$vectors, genf,...)
#  snpgdsClose(genf)
#  return(list(zscores=zscores,eig=eig,chr=chr))
#}


GenomeAdapt.bed=function(genfile,method="EIGMIX",sample.id=NULL, snp.id=NULL,
                         autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                         num.thread=1L, out.fn=NULL, out.prec=c("double", "single"),
                         out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE,...)
{
SNPRelate::snpgdsBED2GDS(paste0(genfile,"bed"),paste0(genfile,"fam"), paste0(genfile,"bim"), "inputgenofile.gds")
genf=SNPRelate::snpgdsOpen("inputgenofile.gds")
rv <- SNPRelate::snpgdsGRM(genf, sample.id=sample.id, snp.id=snp.id,
                           autosome.only=autosome.only, remove.monosnp=remove.monosnp, maf=maf, missing.rate=missing.rate,
                           method=method, num.thread=num.thread, out.fn=out.fn, out.prec=out.prec,
                           out.compress=out.compress, with.id=with.id, verbose=verbose) ### "GCTA" - genetic relationship matrix defined in CGTA; "Eigenstrat" - genetic covariance matrix in EIGENSTRAT; "EIGMIX" - two times coancestry matrix defined in Zheng & Weir (2015), "Weighted" - weighted GCTA, as the same as "EIGMIX", "Corr" - Scaled GCTA GRM (dividing each i,j element by the product of the square root of the i,i and j,j elements), "IndivBeta" - two times individual beta estimate relative to the minimum of beta; see details
eig <- eigen(rv$grm)  # Eigen-decomposition
# close the file
chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genf, "snp.chromosome"))
sampleid=gdsfmt::read.gdsn((gdsfmt::index.gdsn(genf,"sample.id")))
rownames(eig$vectors)=sampleid
zscores <- SNPRelate::snpgdsPCACorr(eig$vectors, genf,snp.id = rv$snp.id,num.thread=num.thread,...)
SNPRelate::snpgdsClose(genf)
unlink("inputgenofile.gds", force=TRUE)
return(list(zscores=zscores,eig=eig,chr=chr,class = "GenomeAdapt"))
}


GenomeAdapt.vcf=function(genfile,method="EIGMIX",sample.id=NULL, snp.id=NULL,
                         autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                         num.thread=1L, out.fn=NULL, out.prec=c("double", "single"),
                         out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE,...)
{

SNPRelate::snpgdsVCF2GDS(genfile, "inputgenofile.gds", method="biallelic.only")
genf=SNPRelate::snpgdsOpen("inputgenofile.gds")
  rv <- SNPRelate::snpgdsGRM(genf, sample.id=sample.id, snp.id=snp.id,
                             autosome.only=autosome.only, remove.monosnp=remove.monosnp, maf=maf, missing.rate=missing.rate,
                             method=method, num.thread=num.thread, out.fn=out.fn, out.prec=out.prec,
                             out.compress=out.compress, with.id=with.id, verbose=verbose) ### "GCTA" - genetic relationship matrix defined in CGTA; "Eigenstrat" - genetic covariance matrix in EIGENSTRAT; "EIGMIX" - two times coancestry matrix defined in Zheng & Weir (2015), "Weighted" - weighted GCTA, as the same as "EIGMIX", "Corr" - Scaled GCTA GRM (dividing each i,j element by the product of the square root of the i,i and j,j elements), "IndivBeta" - two times individual beta estimate relative to the minimum of beta; see details
  eig <- eigen(rv$grm)  # Eigen-decomposition
  # close the file
  chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genf, "snp.chromosome"))
  sampleid=gdsfmt::read.gdsn((gdsfmt::index.gdsn(genf,"sample.id")))
  rownames(eig$vectors)=sampleid
  zscores <- SNPRelate::snpgdsPCACorr(eig$vectors, genf,snp.id = rv$snp.id,num.thread=num.thread,...)
  SNPRelate::snpgdsClose(genf)
  unlink("inputgenofile.gds", force=TRUE)
  return(list(zscores=zscores,eig=eig,chr=chr,class = "GenomeAdapt"))
}

GenomeAdapt.gds=function(genfile,method="EIGMIX",sample.id=NULL, snp.id=NULL,
                         autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                         num.thread=1L, out.fn=NULL, out.prec=c("double", "single"),
                         out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE,...)
{
  genf=SNPRelate::snpgdsOpen(genfile)
  rv <- SNPRelate::snpgdsGRM(genf, sample.id=sample.id, snp.id=snp.id,
                             autosome.only=autosome.only, remove.monosnp=remove.monosnp, maf=maf, missing.rate=missing.rate,
                             method=method, num.thread=num.thread, out.fn=out.fn, out.prec=out.prec,
                             out.compress=out.compress, with.id=with.id, verbose=verbose) ### "GCTA" - genetic relationship matrix defined in CGTA; "Eigenstrat" - genetic covariance matrix in EIGENSTRAT; "EIGMIX" - two times coancestry matrix defined in Zheng & Weir (2015), "Weighted" - weighted GCTA, as the same as "EIGMIX", "Corr" - Scaled GCTA GRM (dividing each i,j element by the product of the square root of the i,i and j,j elements), "IndivBeta" - two times individual beta estimate relative to the minimum of beta; see details
  eig <- eigen(rv$grm)  # Eigen-decomposition
  # close the file
  chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genf, "snp.chromosome"))
  sampleid=gdsfmt::read.gdsn(gdsfmt::index.gdsn(genf,"sample.id"))
  rownames(eig$vectors)=sampleid
  zscores <- SNPRelate::snpgdsPCACorr(eig$vectors, genf,snp.id = rv$snp.id,num.thread=num.thread,...)
  SNPRelate::snpgdsClose(genf)
  return(list(zscores=zscores,eig=eig,chr=chr,class = "GenomeAdapt"))
}

