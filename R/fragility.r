#' Compute the normalized fragility row for adjacency matrix A
#' 
#' @param A Numeric. Adjacency Matrix  
#' @param n_search Integer. Number of eigenvalues tried to find the minimum norm vector 
fragilityRowNormalized <- function(A, n_search = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    #return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = n_search + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(n_search)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  maxf <- max(fragNorm)
  fragNorm2 <- (maxf - fragNorm) / maxf

  return(fragNorm2)
}

#' Compute the fragility row for adjacency matrix A
#'
#' @inheritParams fragilityRowNormalized
fragilityRow <- function(A, n_search = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = n_search + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(n_search)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  return(fragNorm)
}


#' Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
#'
#' @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calc_adj_frag}
#' @param elecSoz Integer.  Vector soz electrodes (for good electrodes)
#' 
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#' @export
#'
#' @examples
#' data("pt01Frag")
#' data("pt01Epoch")
#' sozindex<-attr(pt01Epoch,"sozindex")
#' pt01fragstat<-fragStat(frag=pt01Frag, elecSoz=sozindex)
fragStat <- function(frag, elecSoz){
  if (is(frag, "Fragility")) frag <- frag$frag
  if(!inherits(frag,"matrix")) stop("Fragmust be matrix or Fragility object")
  
  steps <- ncol(frag)   # Time window number
  elecSozC<-which(!(seq_len(nrow(frag))%in%elecSoz))
  hmapSoz<-frag[elecSoz,,drop=FALSE]
  hmapSozC<-frag[elecSozC,,drop=FALSE]
  muSoz<-colMeans(hmapSoz)
  muSozC<-colMeans(hmapSozC)
  sdSoz<-apply(hmapSoz,2L,sd)
  sdSozC<-apply(hmapSozC,2L,sd)
  Q<-seq(.1,1,by=.1)
  qmatrix<-rbind(
    apply(hmapSoz,2L,quantile,Q),
    apply(hmapSozC,2L,quantile,Q)
  )
  rowPrefix<-rep(c("SOZ","SOZC"),each=10)
  dimN<-dimnames(qmatrix)
  dimnames(qmatrix)<-list(
    Quantiles= paste0(rowPrefix,dimN[[1L]]),
    Step     =dimN[[2L]]
  )
  FragStat(
    qmatrix=qmatrix,
    cmeansoz=muSoz,
    cmeansozc=muSozC,
    csdsoz=sdSoz,
    csdsozc=sdSozC
      )
}


predictRidge <- function(xt, A) {
    ## the data matrix
    if (nrow(A) == ncol(A) + 1) {
        x <- cbind(1, as.matrix(xt))
    } else {
        x <- as.matrix(xt)
    }
    x %*% A
}