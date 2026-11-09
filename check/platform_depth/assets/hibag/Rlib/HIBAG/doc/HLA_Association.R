## ----echo=FALSE---------------------------------------------------------------------------------------------
library(HIBAG)

fn <- system.file("doc", "case_control.txt", package="HIBAG")

## -----------------------------------------------------------------------------------------------------------
dat <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
head(dat)

# make an object for hlaAssocTest
hla <- hlaAllele(dat$sample.id, H1=dat$A, H2=dat$A.1, locus="A", assembly="hg19", prob=dat$prob)
summary(hla)

## -----------------------------------------------------------------------------------------------------------
hlaAssocTest(hla, disease ~ h, data=dat)  # 95% confidence interval (h.2.5%, h.97.5%)

# show details
print(hlaAssocTest(hla, disease ~ h, data=dat, verbose=FALSE))

hlaAssocTest(hla, disease ~ h, data=dat, prob.threshold=0.5)  # regression with a threshold

hlaAssocTest(hla, disease ~ h, data=dat, showOR=TRUE)  # report odd ratio instead of log odd ratio

hlaAssocTest(hla, disease ~ h + pc1, data=dat)  # confounding variable pc1

hlaAssocTest(hla, disease ~ h, data=dat, model="additive")  # use an additive model

hlaAssocTest(hla, trait ~ h, data=dat)  # continuous outcome

## -----------------------------------------------------------------------------------------------------------
aa <- hlaConvSequence(hla, code="P.code.merge")

## -----------------------------------------------------------------------------------------------------------
head(c(aa$value$allele1, aa$value$allele2))

# show cross tabulation at each amino acid position
summary(aa)

# association tests
hlaAssocTest(aa, disease ~ h, data=dat, model="dominant")  # try dominant models

hlaAssocTest(aa, disease ~ h, data=dat, model="dominant", prob.threshold=0.5)  # try dominant models

hlaAssocTest(aa, disease ~ h, data=dat, model="recessive")  # try recessive models

## -----------------------------------------------------------------------------------------------------------
sessionInfo()

