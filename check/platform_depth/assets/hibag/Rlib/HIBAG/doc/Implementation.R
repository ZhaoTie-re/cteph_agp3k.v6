## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------------------------------------
library(HIBAG)
library(ggplot2)

version <- c("v1.24.0\n(baseline)", "v1.24.0\n(POPCNT)", "v1.26.1\n(SSE2)",
    "v1.26.1\n(SSE4&POPCNT)", "v1.26.1\n(AVX)", "v1.26.1\n(AVX2)", "v1.26.1\n(AVX512F)",
    "v1.26.1\n(AVX512BW)")
vs <- factor(version, version)

draw <- function(dat)
{
    ggplot(dat, aes(x=Intrinsics, y=Speedup, fill=Gene)) +
        xlab("Package Version and Intrinsics") + ylab("Speed-up factor") +
        geom_bar(stat="identity", width=0.9, colour="white", position=position_dodge()) +
        geom_hline(yintercept=1, linetype=2, colour="gray33") +
        geom_text(aes(label=Speedup), vjust=-0.55, position=position_dodge(0.9), size=2) +
        scale_fill_brewer(palette="Dark2") + theme_bw()
}

## -----------------------------------------------------------------------------------------------------------
# continue without interrupting
IgnoreError <- function(cmd) tryCatch(cmd, error=function(e) { message("Not support"); invisible() })

## -----------------------------------------------------------------------------------------------------------
IgnoreError(hlaSetKernelTarget("sse4"))
IgnoreError(hlaSetKernelTarget("avx"))
IgnoreError(hlaSetKernelTarget("avx2"))
IgnoreError(hlaSetKernelTarget("avx512f"))
IgnoreError(hlaSetKernelTarget("avx512bw"))

## ----echo=FALSE,fig.align='center',fig.width=9,fig.height=3.2-----------------------------------------------
s <- "HLA-A  HLA-B  HLA-C  HLA-DRB1
1.0	1.0	1.0	1.0
1.7	1.6	1.6	1.6
1.2	1.1	1.0	1.0
2.3	2.2	2.2	2.2
2.7	2.5	2.8	2.6
3.2	2.7	2.8	2.7
3.3	2.8	3.6	2.9
4.2	3.5	4.6	3.9"

b <- read.table(text=s, header=TRUE)
colnames(b) <- gsub(".", "-", colnames(b), fixed=TRUE)
b1 <- data.frame(Intrinsics=rep(vs, ncol(b)), Gene=as.factor(rep(names(b), each=nrow(b))),
    Speedup=unname(unlist(b)))
draw(b1) + ggtitle("Building HIBAG models using ~1,000 samples:")

## ----echo=FALSE,fig.align='center',fig.width=9,fig.height=3.2-----------------------------------------------
s <- "HLA-A	HLA-B	HLA-C	HLA-DRB1	HLA-DQA1	HLA-DQB1
1.0	1.0	1.0	1.0	1.0	1.0
1.6	1.6	1.6	1.6	1.5	1.6
1.1	1.0	1.0	1.0	1.0	1.1
2.2	2.2	2.2	2.2	2.3	2.3
2.6	2.8	2.9	2.8	2.9	2.9
2.7	2.8	2.9	2.9	3.0	3.0
2.9	3.7	3.5	3.4	4.1	3.8
3.5	4.7	4.7	4.6	5.3	5.2"

b <- read.table(text=s, header=TRUE)
colnames(b) <- gsub(".", "-", colnames(b), fixed=TRUE)
b2 <- data.frame(Intrinsics=rep(vs, ncol(b)),
    Gene=factor(rep(names(b), each=nrow(b)), names(b)),
    Speedup=unname(unlist(b)))
draw(b2) + ggtitle("Building HIBAG models using ~5,000 samples:")

## ----echo=FALSE,warning=FALSE,fig.align='center',fig.width=9,fig.height=3.2---------------------------------
s <- "HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQA1 HLA-DQB1
1.0	1.0	1.0	1.0	1.0	1.0
1.5	1.7	1.7	1.7	1.8	1.7
1.2	1.2	1.2	1.1	1.2	1.2
1.9	2.3	2.3	2.3	2.3	2.4
2.2	2.9	2.8	3.0	3.0	2.9
3.3	3.6	3.6	3.6	3.7	3.7
4.1	4.1	4.4	4.3	4.4	4.5
5.4	6.0	6.4	6.5	6.9	7.0"

b <- read.table(text=s, header=TRUE)
colnames(b) <- gsub(".", "-", colnames(b), fixed=TRUE)
b3 <- data.frame(Intrinsics=rep(vs, ncol(b)),
    Gene=factor(rep(names(b), each=nrow(b)), names(b)),
    Speedup=unname(unlist(b)))
draw(b3) + ggtitle("Building HIBAG models using ~10,000 samples:")

## ----echo=FALSE,warning=FALSE,fig.align='center',fig.width=9,fig.height=3.2---------------------------------
s <- "HLAgenes	1thread	8threads	16threads	1thread	8threads	16threads
HLA-A	1.0	7.5	13.8	1.3	9.5	17.5
HLA-B	1.0	7.6	14.6	1.6	12.0	22.6
HLA-C	1.0	7.5	13.9	1.6	11.3	20.1
HLA-DRB1	1.0	7.5	14.2	1.6	11.4	21.3
HLA-DQA1	1.0	7.1	12.2	1.9	11.9	19.0
HLA-DQB1	1.0	7.2	12.9	1.7	11.4	19.0"

b <- read.table(text=s, header=TRUE)
x <- unname(unlist(b[,-1]))
nt <- c(rep(paste("AVX2:", c(1,8,16)), each=nrow(b)), rep(paste("AVX512BW:", c(1,8,16)), each=nrow(b)))
nt <- factor(nt, c("AVX2: 1", "AVX2: 8", "AVX2: 16", "AVX512BW: 1", "AVX512BW: 8", "AVX512BW: 16"))
ss <- levels(nt)
ss[1L] <- "AVX2: 1 (baseline)"
levels(nt) <- ss
dat <- data.frame(Speedup=x, Threads=nt,
    Gene=factor(rep(b$HLAgenes, 6), b$HLAgenes))
ggplot(dat, aes(x=Gene, y=Speedup, fill=Threads)) +
    xlab("HLA genes / HIBAG_v1.26") + ylab("Speed-up factor") +
    geom_bar(stat="identity", width=0.9, colour="white", position=position_dodge()) +
    geom_hline(yintercept=1, linetype=2, colour="gray33") +
    geom_text(aes(label=Speedup), vjust=-0.55, position=position_dodge(0.9), size=2) +
    scale_fill_brewer(palette="Dark2") + theme_bw() +
    ggtitle("Building HIBAG models using ~5,000 samples and AVX2/AVX512BW Intrinsics:")

## -----------------------------------------------------------------------------------------------------------
sessionInfo()

