source("param_conv.R")

conv_sp <- function(par) {
    v <- conv_jds2mspr_sp(par)
    names(v) <- c("Ne Anc S-W-H-I",
                  "t (S-W)-(H-I)",
                  "Ne Anc S-W",
                  "Ne Anc H-I",
                  "t H-I",
                  "Ne H",
                  "Growth rate Ne H",
                  "Ne I",
                  "t S-W",
                  "Ne W",
                  "Ne S",
                  "mig W <- S",
                  "mig S <- W",
                  "mig W <- H",
                  "mig H <- W",
                  "mig I <- H",
                  "mig H <- I",
                  "t start mig W-H",
                  "recomb rate",
                  "Start time growth H",
                  "Start time growth W",
                  "Growth rate Ne W"
                  )
    v
}

conv_hd <- function(par) {
    v <- conv_jds2mspr_hd(par)
    names(v) <- c("Ne Anc S-W-H-I",
                  "t S-(W-H-I)",
                  "Ne Anc H-W",
                  "Ne Anc H-I",
                  "t H-I",
                  "Ne H",
                  "Growth rate Ne H",
                  "Ne I",
                  "t W-H",
                  "Ne W",
                  "Ne S",
                  "mig W <- S",
                  "mig S <- W",
                  "mig W <- H",
                  "mig H <- W",
                  "mig I <- H",
                  "mig H <- I",
                  "t start mig W-S",
                  "recomb rate",
                  "Start time growth H",
                  "Start time growth W",
                  "Growth rate Ne W"
                  )
    v
}

## next are parameters from best fastsimcoal result for model m4

mlpar <- "..."
names(mlpar) <- c("NPOP1", "NPOP2", "NPOP3", "NPOP4", "NANC1", "NANC2", "NANC3", "TMRMG", "MIG01R", "MIG10R", "MIG12R", "MIG21R", "MIG23R", "MIG32R", "GR2", "TDIV3", "TDIV2", "TDIV1", "ENDX1", "MaxEstLhood", "MaxObsLhood")

par <- c(p1=mlpar["NANC3"], p2=mlpar["TDIV3"], p3=mlpar["NANC1"], p4=mlpar["NANC2"], p5=mlpar["TDIV2"], p6=mlpar["NPOP3"],
            p7=mlpar["GR2"], p8=mlpar["NPOP4"], p9=mlpar["TDIV1"], p10=mlpar["NPOP2"], p11=mlpar["NPOP1"],
            p12=mlpar["MIG10R"], p13=mlpar["MIG01R"], p14=mlpar["MIG12R"], p15=mlpar["MIG21R"],
            p16=mlpar["MIG32R"], p17=mlpar["MIG23R"], p18=mlpar["ENDX1"], p19=exp(-2.046097e+01), p20 = mlpar["TMRMG"])

L <- 10000
NS <- par[11]/2
fNS <- 4*NS

fsc_m4_est  <- data.frame(estimate = c(
"Ne Anc S-W-H-I"=mlpar["NANC3"]/2,
      "t S-(W-H-I)"=mlpar["TDIV3"],
      "Ne Anc H-W"=mlpar["NANC1"]/2,
      "Ne Anc H-I"=mlpar["NANC2"]/2,
      "t H-I"=mlpar["TDIV2"],
      "Ne H"=mlpar["NPOP3"]/2,
      "Growth rate Ne H"=mlpar["GR2"],
      "Ne I"=mlpar["NPOP4"]/2,
      "t W-H"=mlpar["TDIV1"],
      "Ne W"=mlpar["NPOP2"]/2,
      "Ne S"=mlpar["NPOP1"]/2,
      "mig W <- S"=mlpar["MIG10R"],
      "mig S <- W"=mlpar["MIG01R"],
      "mig W <- H"=mlpar["MIG12R"],
      "mig H <- W"=mlpar["MIG21R"],
      "mig I <- H"=mlpar["MIG32R"],
      "mig H <- I"=mlpar["MIG23R"],
      "t start mig W-S"=mlpar["TMRMG"],
      "rho/theta"=exp(-2.046097e+01)/3.18e-9
))

###########################################################
## fastsimcoal est model3

ml3 <- "..."
names(ml3) <- c("NPOP1", "NPOP2", "NPOP3", "NPOP4", "NANC1", "NANC2", "NANC3", "TMRMG", "MIG01R", "MIG10R", "MIG12R", "MIG21R", "MIG23R", "MIG32R", "GR2", "TDIV3", "TDIV2", "TDIV1", "ENDX1", "MaxEstLhood", "MaxObsLhood")

m3est <- c("Ne Anc S-W-H-I"=ml3["NANC3"]/2,
      "t (S-W)-(H-I)"=ml3["TDIV3"],
      "Ne Anc S-W"=ml3["NANC1"]/2,
      "Ne Anc H-I"=ml3["NANC2"]/2,
      "t H-I"=ml3["TDIV2"],
      "Ne H"=ml3["NPOP3"]/2,
      "Growth rate Ne H"=ml3["GR2"],
      "Ne I"=ml3["NPOP4"]/2,
      "t S-W"=ml3["TDIV1"],
      "Ne W"=ml3["NPOP2"]/2,
      "Ne S"=ml3["NPOP1"]/2,
      "mig W <- S"=ml3["MIG10R"],
      "mig S <- W"=ml3["MIG01R"],
      "mig W <- H"=ml3["MIG12R"],
      "mig H <- W"=ml3["MIG21R"],
      "mig I <- H"=ml3["MIG32R"],
      "mig H <- I"=ml3["MIG23R"],
      "t start mig C-H"=ml3["TMRMG"],
      "rho/theta"=exp(-2.046097e+01)/3.18e-9)

##########################################################

load("jaatha_estimations/jfrsp_zi_123_folded.RData")
spf <- data.frame(conv_sp(estimatesfrsp$estimate))
f <- system("ls parbs_frsp_folded/jfrsp_parbs_*_123_folded.RData", intern=TRUE)
for(i in 1:length(f)) {
    load(f[i])
    spf[i+1] <- data.frame(conv_sp(bsfrspf$estimate))
}

load("jaatha_estimations/jfrsp_zi_123_folded.RData")
spfi <- data.frame(conv_sp(estimatesfrsp$estimate))
fi <- system("ls parbs_frsp_is_folded/jfrsp_parbs_*_123_is_folded.RData", intern=TRUE)
for(i in 1:length(fi)) {
    load(fi[i])
    spfi[i+1] <- data.frame(conv_sp(bsfrspf$estimate))
}



load("jaatha_estimations/jfrhd_zi_123_folded.RData")
hdf <- data.frame(conv_hd(estimatesfrhd$estimate))
f <- system("ls parbs_frhd_folded/jfrhd_parbs_*_123_folded.RData", intern=TRUE)
for(i in 1:length(f)) {
    load(f[i])
    hdf[i+1] <- data.frame(conv_hd(bsfrhdf$estimate))
}

m4f <- fsc_m4_est
f <- system("ls parbs_mod4_folded/jfrhd_mod4_parbs_*_123_folded.RData", intern=TRUE)
for(i in 1:length(f)) {
    load(f[i])
    m4f[i+1] <- data.frame(conv_hd(bsfrhdf$estimate))
}


plot(c(), xlim=c(0,20000), ylim=c(0,20000))
arrows(spf["t H-I", 1], spf["t S-C", 1], spf["t H-I", 2], spf["t S-C", 2])

arrowplot <- function(d1, d2, p1, p2, main="") {
    m1 <- 1.2 * max(c(unlist(d1[p1,]), unlist(d2[p1,])))
    m2 <- 1.2 * max(c(unlist(d1[p2,]), unlist(d2[p2,])))
    m <- max(m1, m2)
    if(min(m1, m2) > 0.3 * m) {m1=m; m2=m}
    mi1 <- 1.2 * min(c(unlist(d1[p1,]), unlist(d2[p1,])))
    mi2 <- 1.2 * min(c(unlist(d1[p2,]), unlist(d2[p2,])))
    mi <- min(m1, m2)
    if(max(m1, m2) < 0.3 * mi) {mi1=mi; mi2=mi}
    plot(NULL, xlim=c(min(0, mi1), m1), ylim=c(min(0,mi2), m2), xlab=p1, ylab=p2, main=main)
    for(i in 2:ncol(d1)) arrows(d1[p1, 1], d1[p2, 1], d1[p1, i], d1[p2, i])
    for(i in 2:ncol(d2)) arrows(d2[p1, 1], d2[p2, 1], d2[p1, i], d2[p2, i], col="red")
    abline(h=0)
    abline(v=0)
}

par(mfcol=c(3,4))
arrowplot(spf, spu, "t (S-W)-(H-I)", "t H-I", main="From Spain")
arrowplot(spf, spu, "t S-W", "t start mig W-H", main="From Spain")
arrowplot(spf, spu, "Ne Anc S-W-H-I", "Ne Anc S-W", main="From Spain")
arrowplot(spf, spu, "Ne Anc H-I", "Ne Anc S-W", main="From Spain")
arrowplot(spf, spu, "Ne S", "Ne W", main="From Spain")
arrowplot(spf, spu, "Ne H", "Ne I", main="From Spain")
arrowplot(spf, spu, "mig W <- S", "mig S <- W", main="From Spain")
arrowplot(spf, spu, "mig W <- H", "mig H <- W", main="From Spain")
arrowplot(spf, spu, "mig I <- H", "mig H <- I", main="From Spain")
arrowplot(spf, spu, "rho/theta", "Growth rate Ne H", main="From Spain")
## dev.copy2pdf(file="parbs_sp.pdf")

par(mfcol=c(3,4))
arrowplot(hdf, hdu, "t S-(W-H-I)", "t H-I", main="From hooded")
arrowplot(hdf, hdu, "t W-H", "t start mig W-S", main="From hooded")
arrowplot(hdf, hdu, "Ne Anc S-W-H-I", "Ne Anc H-W", main="From hooded")
arrowplot(hdf, hdu, "Ne Anc H-I", "Ne Anc H-W", main="From hooded")
arrowplot(hdf, hdu, "Ne S", "Ne W", main="From hooded")
arrowplot(hdf, hdu, "Ne H", "Ne I", main="From hooded")
arrowplot(hdf, hdu, "mig W <- S", "mig S <- W", main="From hooded")
arrowplot(hdf, hdu, "mig W <- H", "mig H <- W", main="From hooded")
arrowplot(hdf, hdu, "mig I <- H", "mig H <- I", main="From hooded")
arrowplot(hdf, hdu, "rho/theta", "Growth rate Ne H", main="From hooded")
## dev.copy2pdf(file="parbs_hd.pdf")

par(mfcol=c(3,4))
arrowplot(m4f, m4u, "t S-(W-H-I)", "t H-I", main="mod4")
arrowplot(m4f, m4u, "t W-H", "t start mig W-S", main="mod4")
arrowplot(m4f, m4u, "Ne Anc S-W-H-I", "Ne Anc H-W", main="mod4")
arrowplot(m4f, m4u, "Ne Anc H-I", "Ne Anc H-W", main="mod4")
arrowplot(m4f, m4u, "Ne S", "Ne W", main="mod4")
arrowplot(m4f, m4u, "Ne H", "Ne I", main="mod4")
arrowplot(m4f, m4u, "mig W <- S", "mig S <- W", main="mod4")
arrowplot(m4f, m4u, "mig W <- H", "mig H <- W", main="mod4")
arrowplot(m4f, m4u, "mig I <- H", "mig H <- I", main="mod4")
arrowplot(m4f, m4u, "rho/theta", "Growth rate Ne H", main="mod4")
## dev.copy2pdf(file="parbs_m4.pdf")


load("jaatha_estimations/jfrsp_zi_123_folded.RData")
rspf <- data.frame(estimatesfrsp$estimate)
f <- system("ls parbs_frsp_folded/jfrsp_parbs_*_123_folded.RData", intern=TRUE)
for(i in 1:length(f)) {
    load(f[i])
    rspf[i+1] <- data.frame(bsfrspf$estimate)
}

load("jaatha_estimations/jfrsp_zi_123_folded.RData")
rspfi <- data.frame(estimatesfrsp$estimate)
fi <- system("ls parbs_frsp_is_folded/jfrsp_parbs_*_123_is_folded.RData", intern=TRUE)
for(i in 1:length(fi)) {
    load(fi[i])
    rspfi[i+1] <- data.frame(bsfrspf$estimate)
}

load("jaatha_estimations/jfrhd_zi_123_folded.RData")
rhdf <- data.frame(estimatesfrhd$estimate)
f <- system("ls parbs_frhd_folded/jfrhd_parbs_*_123_folded.RData", intern=TRUE)
for(i in 1:length(f)) {
    load(f[i])
    rhdf[i+1] <- data.frame(bsfrhdf$estimate)
}

source("param_conv.R")
cat("Jaatha scaled parameters:", p <- 2*rspf[,1] - apply(rspf[,2:ncol(rspf)], 1, mean), "\n")
cat("msprime scaled parameters:", p <- conv_jds2mspr_hd(p), "\n")
cat("fastsimcoal scaled parameters:", p <- conv_mspr2fsc_hd(p), "\n")
cat("fastsimcoal scaled parameters again, with names:\n")
print(p)

bscorf <- data.frame(conv_sp(rspf[,1]),
                    conv_sp(2*rspf[,1] - apply(rspf[,2:ncol(rspf)], 1, mean)),
                    conv_sp(2*rspf[,1] - apply(rspf[,2:ncol(rspf)], 1, quantile, 0.025)),
                    conv_sp(2*rspf[,1] - apply(rspf[,2:ncol(rspf)], 1, quantile, 0.975)))
colnames(bscorf) <- c("orig", "bscor", "upper", "lower")

bscorfi <- data.frame(conv_sp(rspfi[,1]),
                    conv_sp(2*rspfi[,1] - apply(rspfi[,2:ncol(rspfi)], 1, mean)),
                    conv_sp(2*rspfi[,1] - apply(rspfi[,2:ncol(rspfi)], 1, quantile, 0.025)),
                    conv_sp(2*rspfi[,1] - apply(rspfi[,2:ncol(rspfi)], 1, quantile, 0.975)))
colnames(bscorfi) <- c("orig", "bscor", "upper", "lower")

bscorhf <- data.frame(conv_hd(rhdf[,1]),
                    conv_hd(2*rhdf[,1] - apply(rhdf[,2:ncol(rhdf)], 1, mean)),
                    conv_hd(2*rhdf[,1] - apply(rhdf[,2:ncol(rhdf)], 1, quantile, 0.025)),
                    conv_hd(2*rhdf[,1] - apply(rhdf[,2:ncol(rhdf)], 1, quantile, 0.975)))

colnames(bscorhf) <- c("orig", "bscor", "upper", "lower")

bscorf_raw <- 2*rspf[,1] - apply(rspf[,2:ncol(rspf)], 1, mean)
bscorfi_raw <- 2*rspfi[,1] - apply(rspfi[,2:ncol(rspfi)], 1, mean)
bscorf
bscorfi

bscorff <- signif(data.frame(orig=conv_mspr2fsc_sp(bscorf$orig), bscor=conv_mspr2fsc_sp(bscorf$bscor),
                                  lower=conv_mspr2fsc_sp(bscorf$lower), upper=conv_mspr2fsc_sp(bscorf$upper)),3)



## save(bscorf, bscorfi,bscorf_raw ,bscorfi_raw ,file="jaatha_estimations/jfrsp_zi_123_folded_bs_corrected.RData")

bsplot <- function(d1,  p1, p2, b1, main="", log="", bsci=TRUE, add=FALSE, col="black", xlim=NULL, ylim=NULL) {
    if(bsci) {
        m1 <- max(c(unlist(d1[p1,]),  b1[p1, "upper"]))
        m2 <- max(c(unlist(d1[p2,]),  b1[p2, "upper"]))
    } else {
        m1 <- max(c(unlist(d1[p1,]),  b1[p1, "orig"]))
        m2 <- max(c(unlist(d1[p2,]),  b1[p2, "orig"]))  
    }
    m <- max(m1, m2)
    if(min(m1, m2) > 0.3 * m) {m1=m; m2=m}
    if(bsci) mi <- min(c(unlist(d1[p1,]),  b1[p1, "lower"] ))
    else mi <- min(c(unlist(d1[p1,]),  b1[p1, "orig"] ))
    if(max(m1, m2) < 0.3 * mi) {mi1=mi; mi2=mi}
    if(is.null(xlim)) xlim=c(min(0, mi), m1)
    if(is.null(ylim)) ylim=c(min(0,mi), m2)
    if(add) {
        points(d1[p1, 1], d1[p2, 1], col=col)
    } else {
        plot(d1[p1, 1], d1[p2, 1], xlim=xlim, ylim=ylim, xlab=p1, ylab=p2,
             main=main, log=log, col=col)
    }
    for(i in 2:ncol(d1)) points(d1[p1, i], d1[p2, i], pch=16, cex=0.3, col=col)
    abline(h=0)
    abline(v=0)
    if(bsci) {
        points(b1[p1,2], b1[p2,2], pch=2, col=col)
        lines(rep(b1[p1,2],2), b1[p2,3:4], col=col)
        lines(b1[p1,3:4], rep(b1[p2,2],2), col=col)
    }
}


par(mfcol=c(3,4))
bsplot(spf,   "t (S-W)-(H-I)", "t H-I", bscorf, main="From Spain")
bsplot(spf,   "t S-W", "t start mig W-H", bscorf,  main="From Spain")
bsplot(spf,   "Ne Anc S-W-H-I", "Ne Anc S-W", bscorf, main="From Spain")
bsplot(spf,   "Ne Anc H-I", "Ne Anc S-W", bscorf,  main="From Spain")
bsplot(spf,   "Ne S", "Ne W", bscorf, main="From Spain")
bsplot(spf,   "Ne H", "Ne I", bscorf,  main="From Spain")
bsplot(spf,   "mig W <- S", "mig S <- W", bscorf,  main="From Spain")
bsplot(spf,   "mig W <- H", "mig H <- W", bscorf, main="From Spain")
bsplot(spf,   "mig I <- H", "mig H <- I", bscorf,  main="From Spain")
bsplot(spf,   "recomb rate", "Growth rate Ne H", bscorf,  main="From Spain")
bsplot(spf,   "Growth rate Ne H", "Start time growth H", bscorf,  main="From Spain")
bsplot(spf,   "Growth rate Ne W", "Start time growth W", bscorf,  main="From Spain")
## dev.copy2pdf(file="frsp_bs.pdf")

par(mfcol=c(3,4))
bsplot(spfi,   "t (S-W)-(H-I)", "t H-I", bscorfi, main="From Spain")
bsplot(spfi,   "t S-W", "t start mig W-H", bscorfi,  main="From Spain")
bsplot(spfi,   "Ne Anc S-W-H-I", "Ne Anc S-W", bscorfi, main="From Spain")
bsplot(spfi,   "Ne Anc H-I", "Ne Anc S-W", bscorfi,  main="From Spain")
bsplot(spfi,   "Ne S", "Ne W", bscorfi, main="From Spain")
bsplot(spfi,   "Ne H", "Ne I", bscorfi,  main="From Spain")
bsplot(spfi,   "mig W <- S", "mig S <- W", bscorfi,  main="From Spain")
bsplot(spfi,   "mig W <- H", "mig H <- W", bscorfi, main="From Spain")
bsplot(spfi,   "mig I <- H", "mig H <- I", bscorfi,  main="From Spain")
bsplot(spfi,   "recomb rate", "Growth rate Ne H", bscorfi,  main="From Spain")
bsplot(spfi,   "Growth rate Ne H", "Start time growth H", bscorfi,  main="From Spain")
bsplot(spfi,   "Growth rate Ne W", "Start time growth W", bscorfi,  main="From Spain")
## dev.copy2pdf(file="frsp_is_bs.pdf")

par(mfrow=c(3,4))
bsplot(spf,   "t (S-W)-(H-I)", "t H-I", bscorf, main="A", col="blue")
bsplot(spfi,   "t (S-W)-(H-I)", "t H-I", bscorfi,  add=TRUE)
bsplot(spf,   "t S-W", "t start mig W-H", bscorf,  main="B", col="blue")
bsplot(spfi,   "t S-W", "t start mig W-H", bscorfi,   add=TRUE)
bsplot(spf,   "Ne Anc S-W-H-I", "Ne Anc S-W", bscorf, main="C",  col="blue")
bsplot(spfi,   "Ne Anc S-W-H-I", "Ne Anc S-W", bscorfi,  add=TRUE)
bsplot(spf,   "Ne Anc H-I", "Ne Anc S-W", bscorf, main="D",   col="blue")
bsplot(spfi,   "Ne Anc H-I", "Ne Anc S-W", bscorfi,   add=TRUE)
bsplot(spf,   "Ne S", "Ne W", bscorf, main="E",  col="blue")
bsplot(spfi,   "Ne S", "Ne W", bscorfi,  add=TRUE)
bsplot(spf,   "Ne H", "Ne I", bscorf, main="F",   col="blue")
bsplot(spfi,   "Ne H", "Ne I", bscorfi,   add=TRUE)
bsplot(spf,   "mig W <- S", "mig S <- W", bscorf, main="G",   col="blue")
bsplot(spfi,   "mig W <- S", "mig S <- W", bscorfi,   add=TRUE)
bsplot(spf,   "mig W <- H", "mig H <- W", bscorf, main="H",  col="blue")
bsplot(spfi,   "mig W <- H", "mig H <- W", bscorfi,  add=TRUE)
bsplot(spf,   "mig I <- H", "mig H <- I", bscorf,  main="I",  col="blue", xlim=c(4e-5, 8e-5), ylim=c(8.7e-6, 9e-6))
bsplot(spfi,   "mig I <- H", "mig H <- I", bscorfi,   add=TRUE)
bsplot(spf,   "recomb rate", "Ne Anc S-W-H-I", bscorf, main="J",   col="blue", xlim=c(3e-9, 3.3e-9), ylim=c(28500, 30500))
bsplot(spfi,   "recomb rate", "Ne Anc S-W-H-I", bscorfi,   add=TRUE)
bsplot(spf,   "Growth rate Ne H", "Start time growth H", bscorf, main="K",   col="blue")
bsplot(spfi,   "Growth rate Ne H", "Start time growth H", bscorfi,   add=TRUE)
bsplot(spf,   "Growth rate Ne W", "Start time growth W", bscorf, main="L",   col="blue")
bsplot(spfi,   "Growth rate Ne W", "Start time growth W", bscorfi,   add=TRUE)
## dev.copy2pdf(file="frsp_bs_is_fs.pdf")

### log likelihood-ratios if true model is from-Spain:

## with unfolded spectrum:

s <- system("ls parbs_frsp_unfolded/jfrsp_parbs_*_123_unfolded.RData", intern=TRUE)
h <- system("ls parbs_frhd_unfolded/jfrhd_parbs_*_123_unfolded.RData", intern=TRUE)

lls <- numeric()
llh <- numeric()
for(i in 1:100) {   ### UPDATE NUMBER HERE IF NEEDED
    load(s[i])
    load(h[i])
    lls[i] <-  bsfrsp$loglikelihood
    llh[i] <- bsfrhd$loglikelihood
}

load("jaatha_estimations/jfrsp_zi_123_folded.RData")
estimatesfrsp$loglikelihood
plot(lls, llh)
points(estimatesfrsp$loglikelihood, estimatesfrhd$loglikelihood, col="blue")

### calculation of understandable parameters:

load("jaatha_estimations/jfrsp_zi_123_folded.RData")
p <- conv_sp(estimatesfrsp$estimate)
exp(p["Growth rate Ne W"] * p["Start time growth W"])
exp(p["Growth rate Ne H"] * p["Start time growth H"])
a <- 1-exp(-(p["mig W <- H"]+p["mig W <- S"])*p["t start mig W-H"])
b <- p["mig W <- H"]/(p["mig W <- H"]+p["mig W <- S"])
a*b
