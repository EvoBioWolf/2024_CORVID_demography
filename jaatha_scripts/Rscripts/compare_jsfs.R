
emp <- list() ## will be list of empirical frequeny spectra

### REPLACE PATHS TO FILES WITH PAIRWISE FREQUENCY SPECTRA

## jsfs of S and W
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop1_0.obs")
emp[["SW"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of S and H
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop2_0.obs")
emp[["SH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of S and I 
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_0.obs")
emp[["SI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))
 
## jsfs of W and H    
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop2_1.obs")
emp[["WH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of W and I    
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_1.obs")
emp[["WI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))

## jsfs of H and I     
LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_2.obs")
emp[["HI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))


calc_avg_jsfs <- function(L) {
    ## calculate average frequency spectra of those in the list of files
    Lsum <- list()
    for(n in names(emp)) Lsum[[n]] <- matrix(0, nrow=dim(emp[[n]])[1], ncol=dim(emp[[n]])[2])
    for(f in L) {
        Ljsfs <- list()
        d <- readLines(f)
        i <- which(d=="joint allele frequency spectrum for populations S and W:")
        Ljsfs[["SW"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        i <- which(d=="joint allele frequency spectrum for populations S and H:")
        Ljsfs[["SH"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        i <- which(d=="joint allele frequency spectrum for populations S and I:")
        Ljsfs[["SI"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        i <- which(d=="joint allele frequency spectrum for populations W and H:")
        Ljsfs[["WH"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        i <- which(d=="joint allele frequency spectrum for populations W and I:")
        Ljsfs[["WI"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        i <- which(d=="joint allele frequency spectrum for populations H and I:")
        Ljsfs[["HI"]] <- t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric))
        for(n in names(Ljsfs)) Lsum[[n]] <- Lsum[[n]] + Ljsfs[[n]]
    }
    Lavg <- list()
    for(n in names(Lsum)) Lavg[[n]] <- Lsum[[n]]/length(L)
    return(Lavg)
}

fold_spect <- function(L) {
    ## if L is a matrix, fold the spectrum otherwise assume it's a list and fold its members
    if(is(L,"matrix")) {
        R <- L
        n <- nrow(R)-1
        m <- ncol(R)-1
        ## maximum number of mutations is n+m, so map any (k,l) with k+l > (n+m)/2 to (n-k, m-l),
        ## note that indices (i, j) = (k+1, l+1), so (i,j) with (i,j) > (n+m)/2+2 go to the indices
        ## corresponding to (n-k, m-l)=(n-(i-1),m-(j-1))=(n-i+1,m-j+1), which are (n-i+2,m-j+2)
        for(i in 1:(n+1)) {
            x <- max(1, floor((n+m)/2+3)-i)
            if(x<=m+1) {
                for(j in x:(m+1)) {
                    R[n-i+2,m-j+2] <- R[n-i+2,m-j+2] +  R[i,j]
                    R[i, j] <- 0
                }
            }
        }
    } else {
        R <- list()
        for(n in names(L)) {
            R[[n]] <- fold_spect(L[[n]])
        }
    }
    R
}

compareplot <- function(A, B) {
    par(mfrow=c(2,3))
    for(n in names(A)) {
        C  <-  A[[n]]/(A[[n]]+B[[n]])
        for(i in 1:nrow(C))
            for(j in 1:ncol(C)) {
                if(abs(A[[n]][i,j]-B[[n]][i,j])<1) C[i,j] <- 0.5
            }
        image(C, main=n, breaks=0:11/10-0.05,
              col=rgb(c(0,0.2,0.4,0.6,0.8,rep(1,6)),c(0,0.3,0.6,0.9,1,1,1,0.9,0.6,0.3,0),c(rep(1,6), 0.8,0.6,0.3, 0.2,0)),
              xlab=strsplit(n,"")[[1]][1], ylab=strsplit(n,"")[[1]][2]
              )
    }
}

## simulation results of mod 7x folded
L <- system("ls fs_sim_folded/m7x_simspec_folded_*", intern=TRUE)
m7xavg <- fold_spect(calc_avg_jsfs(L))
compareplot(m7xavg, emp)
##dev.copy2pdf(file="compm7xemp.pdf")

foldsp <- fold_spect(calc_avg_jsfs(system("ls fs_sim_folded/frsp_simspec_folded_*", intern=TRUE)))

compareplot(foldsp, emp)
##dev.copy2pdf(file="compfrspemp.pdf")

foldspbc <- fold_spect(calc_avg_jsfs(system("ls fs_sim_folded_bscor/frsp_simspec_folded_bscor_*", intern=TRUE)))
## foldspbc <- fold_spect(calc_avg_jsfs(system("ls fs_sim_folded/tmp/frspbscor_simspec_folded_*", intern=TRUE)))

compareplot(foldspbc, emp)
##dev.copy2pdf(file="compfrspbsemp.pdf")

foldhd <- fold_spect(calc_avg_jsfs(system("ls fs_sim_folded/frhd_simspec_folded_*", intern=TRUE)))

compareplot(foldhd, emp)
##dev.copy2pdf(file="compfrhdemp.pdf")


foldfscm1 <- fold_spect(calc_avg_jsfs(system("ls fs_sim_folded/fscmod1_simspec_folded_*", intern=TRUE)))

compareplot(foldfscm1, emp)


foldfrspis <- fold_spect(calc_avg_jsfs(system("ls is_sim_folded/frsp_simspec_folded_*", intern=TRUE)))

compareplot(foldfrspis, emp)
##dev.copy2pdf(file="compfrspemp_is.pdf")

## below comparisons to Chyi Yin's spectra based on fastsimcoal simulations


### first check whether we are on the same page regarding empirical spectra

obs <- list() ## will be list of empirical frequeny spectra

## jsfs of S and W
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop1_0.obs")
obs[["SW"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of S and H
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop2_0.obs")
obs[["SH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of S and I 
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_0.obs")
obs[["SI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))
 
## jsfs of W and H    
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop2_1.obs")
obs[["WH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))

## jsfs of W and I    
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_1.obs")
obs[["WI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))

## jsfs of H and I     
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_2.obs")
obs[["HI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))

compareplot(obs, emp) ## YES! all fine


### model 7x simulated with fastsimcoal

m7xf <- list() ## will be list of empirical frequeny spectra

## jsfs of S and W
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop1_0.txt")
m7xf[["SW"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m7xf[["SW"]] <- m7xf[["SW"]]*(sum(obs[["SW"]]) - obs[["SW"]][1,1])/(1-m7xf[["SW"]][1,1])
    
## jsfs of S and H
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop2_0.txt")
m7xf[["SH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m7xf[["SH"]] <- m7xf[["SH"]]*(sum(obs[["SH"]]) - obs[["SH"]][1,1])/(1-m7xf[["SH"]][1,1])

## jsfs of S and I 
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_0.txt")
m7xf[["SI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m7xf[["SI"]] <- m7xf[["SI"]]*(sum(obs[["SI"]]) - obs[["SI"]][1,1])/(1-m7xf[["SI"]][1,1])
 
## jsfs of W and H    
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop2_1.txt")
m7xf[["WH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m7xf[["WH"]] <- m7xf[["WH"]]*(sum(obs[["WH"]]) - obs[["WH"]][1,1])/(1-m7xf[["WH"]][1,1])

## jsfs of W and I    
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_1.txt")
m7xf[["WI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m7xf[["WI"]] <- m7xf[["WI"]]*(sum(obs[["WI"]]) - obs[["WI"]][1,1])/(1-m7xf[["WI"]][1,1])

## jsfs of H and I     
LCS <- readLines("../fastsimcoal/folded/no_chr18/4PopModel7x_jointMAFpop3_2.txt")
m7xf[["HI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m7xf[["HI"]] <- m7xf[["HI"]]*(sum(obs[["HI"]]) - obs[["HI"]][1,1])/(1-m7xf[["HI"]][1,1])

compareplot(m7xf, emp)
##dev.copy2pdf(file="compm7xemp_fsim.pdf")

### model 8x simulated with fastsimcoal


m8xf <- list() ## will be list of empirical frequeny spectra

## jsfs of S and W
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop1_0.txt")
m8xf[["SW"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m8xf[["SW"]] <- m8xf[["SW"]]*(sum(obs[["SW"]]) - obs[["SW"]][1,1])/(1-m8xf[["SW"]][1,1])
    
## jsfs of S and H
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop2_0.txt")
m8xf[["SH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m8xf[["SH"]] <- m8xf[["SH"]]*(sum(obs[["SH"]]) - obs[["SH"]][1,1])/(1-m8xf[["SH"]][1,1])

## jsfs of S and I 
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop3_0.txt")
m8xf[["SI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m8xf[["SI"]] <- m8xf[["SI"]]*(sum(obs[["SI"]]) - obs[["SI"]][1,1])/(1-m8xf[["SI"]][1,1])
 
## jsfs of W and H    
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop2_1.txt")
m8xf[["WH"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:32], "[", 2:32), 1, as.numeric))
m8xf[["WH"]] <- m8xf[["WH"]]*(sum(obs[["WH"]]) - obs[["WH"]][1,1])/(1-m8xf[["WH"]][1,1])

## jsfs of W and I    
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop3_1.txt")
m8xf[["WI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m8xf[["WI"]] <- m8xf[["WI"]]*(sum(obs[["WI"]]) - obs[["WI"]][1,1])/(1-m8xf[["WI"]][1,1])

## jsfs of H and I     
LCS <- readLines("../fastsimcoal/folded/no_chr18/from_hooded/4PopModel8x_jointMAFpop3_2.txt")
m8xf[["HI"]] <- t(apply(sapply(strsplit(LCS, "[\t ]")[2:12], "[", 2:32), 1, as.numeric))
m8xf[["HI"]] <- m8xf[["HI"]]*(sum(obs[["HI"]]) - obs[["HI"]][1,1])/(1-m8xf[["HI"]][1,1])

compareplot(m8xf, emp)
##dev.copy2pdf(file="compm8xemp_fsim.pdf")
