### functions to calculate summary statistics for jaatha analyses

Rcpp::sourceCpp("../src/c++/fgcv.cc")

confuse_direction <- function(J, p) {
    ## for a fraction p of the mutations in the joint site-frequency spectrum J
    ## confuse ancestral and derived
    J2 <- J
    n <- nrow(J)
    m <- ncol(J)
    for(i in 1:n)
        for(j in 1:m) {
            if( 2*i != n | 2*j != m) {
                k <- rbinom(1, J[i, j], p)
                J2[i, j] <- J2[i, j] - k
                J2[n+1-i, m+1-j] <- J2[n+1-i, m+1-j] + k
            }
        }
    J2
}

categorize_fgc <- function(B) {
    ## for single populations use boundaries 0.1, 0.2, 0.3,
    ## for pairs of populations use boundaries 0.005, 0.01, 0.015
    x <- apply(apply(B, 1, "<", rep(c(0.1, 0.005), c(4,3))), 1, sum, na.rm=TRUE)
    y <- apply(apply(B, 1, "<", rep(c(0.2, 0.01), c(4,3))), 1, sum, na.rm=TRUE) - x
    z <- apply(apply(B, 1, "<", rep(c(0.3, 0.015), c(4,3))), 1, sum, na.rm=TRUE) - y
    c(x,y,z,apply(apply(B, 1, ">=", rep(c(0.3, 0.015), c(4,3))), 1, sum, na.rm=TRUE))
}

coarse_jsfs_info_folded <-  function() {
    ind1 <- c("0","1:2","3:(n-3)","(n-2):(n-1)","n")
    ind2 <- c("0","1:2","3:(m-3)","(m-2):(m-1)","m")
    for(i in 1:5) {
        for(j in 1:5) {
            if( i+j < 6 & !(i==1 & j==1)) {
                cat("(",ind1[[i]],",", ind2[[j]],") + (",ind1[[6-i]],",", ind2[[6-j]],")", "\n")
            }
        }
    }
}


coarse_jsfs <- function(J, folded=FALSE, cap=0.0) {
    ## cap is probability that the ancestral state is confused
    n <- nrow(J)
    m <- ncol(J)
    ind1 <- list(1,2:3,4:(n-3),(n-2):(n-1),n)
    ind2 <- list(1,2:3,4:(m-3),(m-2):(m-1),m)
    s <- numeric()
    if(folded) {
        for(i in 1:5) {
            for(j in 1:5) {
                if( i+j < 6 & !(i==1 & j==1)) {
                    s <- c(s, sum(J[ind1[[i]], ind2[[j]]]) + sum(J[ind1[[6-i]], ind2[[6-j]]]))
                }
            }
        }
        s <- c(s, sum(J[ind1[[3]], ind2[[3]]]),
               sum(J[ind1[[2]], ind2[[4]]]) + sum(J[ind1[[4]], ind2[[2]]]),
               sum(J[ind1[[1]], ind2[[5]]]) + sum(J[ind1[[5]], ind2[[1]]]))
    } else {
        if(cap > 0.0) {
            dJ <- matrix(0, nrow=n, ncol=m)
            for(i in 1:n) {
                for(j in 1:m) {
                    if((i!=1 | j!=1) & (i!=n-i+1 | j!=m-j+1) & (i!=n | j!=m)) {
                        d <- rbinom(1, J[i,j], cap)
                        dJ[i, j] <- dJ[i, j] - d
                        dJ[n-i+1, m-j+1] <- dJ[n-i+1, m-j+1] + d
                    }
                }
            }
            J  <- J + dJ
        }
        for(i in 1:5) {
            for(j in 1:5) {
                if( !(i==1 & j==1) & !(i==5 & j==5) ) {
                    s <- c(s, sum(J[ind1[[i]], ind2[[j]]]))
                }
            }
        }
    }
    s
}
