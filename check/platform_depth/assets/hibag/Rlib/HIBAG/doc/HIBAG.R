## ----echo=FALSE---------------------------------------------------------------------------------------------
options(width=110)

## -----------------------------------------------------------------------------------------------------------
library(HIBAG)

## ----fig.width=10, fig.height=4, fig.align='center'---------------------------------------------------------
# load the published parameter estimates from European ancestry
#   e.g., filename <- "HumanOmniExpressExome-European-HLA4-hg19.RData"
#   here, we use example data in the package
filename <- system.file("extdata", "ModelList.RData", package="HIBAG")
model.list <- get(load(filename))

# HLA imputation at HLA-A
hla.id <- "A"
model <- hlaModelFromObj(model.list[[hla.id]])
summary(model)
plot(model)    # show the frequency of SNP marker in the model

# SNPs in the model
head(model$snp.id)

head(model$snp.position)

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  #########################################################################
#  # Import your PLINK BED file
#  #
#  yourgeno <- hlaBED2Geno(bed.fn=".bed", fam.fn=".fam", bim.fn=".bim")
#  summary(yourgeno)
#  
#  # best-guess genotypes and all posterior probabilities
#  pred.guess <- hlaPredict(model, yourgeno, type="response+prob")
#  summary(pred.guess)
#  
#  pred.guess$value
#  pred.guess$postprob

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  library(HIBAG)
#  
#  # load HLA types and SNP genotypes in the package
#  data(HLA_Type_Table, package="HIBAG")
#  data(HapMap_CEU_Geno, package="HIBAG")
#  
#  # a list of HLA genotypes
#  # e.g., train.HLA <- hlaAllele(sample.id, H1=c("01:02", "05:01", ...),
#  #                        H2=c("02:01", "03:01", ...), locus="A")
#  #     the HLA type of the first individual is 01:02/02:01,
#  #                 the second is 05:01/03:01, ...
#  hla.id <- "A"
#  train.HLA <- hlaAllele(HLA_Type_Table$sample.id,
#      H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
#      H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
#      locus=hla.id, assembly="hg19")
#  
#  # training genotypes
#  #   import your PLINK BED file,
#  #   e.g., train.geno <- hlaBED2Geno(bed.fn=".bed", fam.fn=".fam", bim.fn=".bim")
#  #   and select the SNPs in the flanking region of 500kb on each side
#  region <- 500    # kb
#  train.geno <- hlaGenoSubsetFlank(HapMap_CEU_Geno, hla.id, region*1000, assembly="hg19")
#  summary(train.geno)
#  
#  # train a HIBAG model
#  set.seed(100)
#  model <- hlaAttrBagging(train.HLA, train.geno, nclassifier=100)

## ----echo=FALSE---------------------------------------------------------------------------------------------
mobj <- get(load(system.file("extdata", "ModelList.RData", package="HIBAG")))
model <- hlaModelFromObj(mobj$A)

## -----------------------------------------------------------------------------------------------------------
summary(model)

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  library(HIBAG)
#  library(parallel)
#  
#  # load HLA types and SNP genotypes in the package
#  data(HLA_Type_Table, package="HIBAG")
#  data(HapMap_CEU_Geno, package="HIBAG")
#  
#  # a list of HLA genotypes
#  # e.g., train.HLA <- hlaAllele(sample.id, H1=c("01:02", "05:01", ...),
#  #                        H2=c("02:01", "03:01", ...), locus="A")
#  #     the HLA type of the first individual is 01:02/02:01,
#  #                 the second is 05:01/03:01, ...
#  hla.id <- "A"
#  train.HLA <- hlaAllele(HLA_Type_Table$sample.id,
#      H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
#      H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
#      locus=hla.id, assembly="hg19")
#  
#  # training genotypes
#  #   import your PLINK BED file,
#  #   e.g., train.geno <- hlaBED2Geno(bed.fn=".bed", fam.fn=".fam", bim.fn=".bim")
#  #   and select the SNPs in the flanking region of 500kb on each side
#  region <- 500    # kb
#  train.geno <- hlaGenoSubsetFlank(HapMap_CEU_Geno, hla.id, region*1000, assembly="hg19")
#  summary(train.geno)
#  
#  
#  # Multithreading
#  cl <- 2    # 2 -- # of threads
#  
#  # Building ...
#  set.seed(1000)
#  hlaParallelAttrBagging(cl, train.HLA, train.geno, nclassifier=100,
#      auto.save="AutoSaveModel.RData")
#  model.obj <- get(load("AutoSaveModel.RData"))
#  model <- hlaModelFromObj(model.obj)
#  summary(model)

## -----------------------------------------------------------------------------------------------------------
library(HIBAG)

# load HLA types and SNP genotypes in the package
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a list of HLA types
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
summary(hlatab$training)
summary(hlatab$validation)

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id,
    HapMap_CEU_Geno$snp.position, hla.id, region*1000, assembly="hg19")
length(snpid)

# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
    samp.sel = match(hlatab$training$value$sample.id,
        HapMap_CEU_Geno$sample.id))
summary(train.geno)

test.geno <- hlaGenoSubset(HapMap_CEU_Geno, samp.sel=match(
    hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  # train a HIBAG model
#  set.seed(100)
#  model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=100)

## ----echo=FALSE---------------------------------------------------------------------------------------------
mobj <- get(load(system.file("extdata", "OutOfBag.RData", package="HIBAG")))
model <- hlaModelFromObj(mobj)

## -----------------------------------------------------------------------------------------------------------
summary(model)

# validation
pred <- hlaPredict(model, test.geno)
summary(pred)

# compare
comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model, call.threshold=0)
comp$overall

## ----fig.align="center", fig.height=4, fig.width=5----------------------------------------------------------
# hierarchical cluster analysis
d <- hlaDistance(model)
p <- hclust(as.dist(d))
plot(p, xlab="HLA alleles")

## ----fig.align="center", fig.height=3.3, fig.width=5--------------------------------------------------------
# violin plot
hlaReportPlot(pred, model=model, fig="matching")

## ----fig.align="center", fig.height=3.3, fig.width=5--------------------------------------------------------
hlaReportPlot(pred, hlatab$validation, fig="call.rate")
hlaReportPlot(pred, hlatab$validation, fig="call.threshold")

## -----------------------------------------------------------------------------------------------------------
# report overall accuracy, per-allele sensitivity, specificity, etc
hlaReport(comp, type="txt")

## ----results='asis'-----------------------------------------------------------------------------------------
# report overall accuracy, per-allele sensitivity, specificity, etc
hlaReport(comp, type="markdown")

## -----------------------------------------------------------------------------------------------------------
# report overall accuracy, per-allele sensitivity, specificity, etc
hlaReport(comp, type="tex", header=FALSE)

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  library(HIBAG)
#  
#  # make a list of HLA types
#  hla.id <- "DQA1"
#  hla <- hlaAllele(HLA_Type_Table$sample.id,
#      H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
#      H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
#      locus = hla.id, assembly = "hg19")
#  
#  # training genotypes
#  region <- 500   # kb
#  snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
#      hla.id, region*1000, assembly="hg19")
#  train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
#      snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
#      samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))
#  
#  set.seed(1000)
#  model <- hlaAttrBagging(hla, train.geno, nclassifier=100)
#  summary(model)
#  
#  # remove unused SNPs and sample IDs from the model
#  mobj <- hlaPublish(model,
#      platform = "Illumina 1M Duo",
#      information = "Training set -- HapMap Phase II",
#      warning = NULL,
#      rm.unused.snp=TRUE, anonymize=TRUE)
#  
#  save(mobj, file="Your_HIBAG_Model.RData")

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  # assume the HIBAG models are stored in R objects: mobj.A, mobj.B, ...
#  
#  ModelList <- list()
#  ModelList[["A"]] <- mobj.A
#  ModelList[["B"]] <- mobj.B
#  ...
#  
#  # save to an R data file
#  save(ModelList, file="HIBAG_Model_List.RData")

## -----------------------------------------------------------------------------------------------------------
sessionInfo()

