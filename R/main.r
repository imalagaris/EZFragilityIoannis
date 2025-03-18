#' Calculate adjacency matrices and fragility matrix from iEEG recording
#'
#' The function calculates the neural fragility column
#' from an adjacency matrix in each time window
#'
#' @source Recreation of the method described in
#' Li A, Huynh C, Fitzgerald Z, Cajigas I, Brusko D, Jagid J, et al.
#' Neural fragility as an EEG marker of the seizure onset zone.
#'  Nat Neurosci. 2021 Oct;24(10):1465–74
#' (\href{https://pubmed.ncbi.nlm.nih.gov/34354282/}{pubmed}).
#' We have found solutions to fill up missing details in the paper method description
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns, or an Epoch object
#' @param window Integer. The number of time points to use in each window
#' @param step Integer. The number of time points to move the window each time
#' @param lambda Numeric. The lambda value to use in the ridge regression.
#' If NULL, the lambda will be chosen automatically
#' ensuring that ensuring that the adjacent matrix is stable (see details)
#' @param nSearch Integer. Number of minimization to compute the fragility row
#' @param progress Logical. If TRUE, print progress information. If `parallel` is TRUE, this option only support the `doSNOW` backend.
#' @param parallel Logical. If TRUE, use parallel computing. 
#' Users must register a parallel backend with the foreach package
#' 
#'
#' @return A list containing the normalized ieegts,
#' adjacency matrices, fragility, and R^2 values
#'
#' @examples
#' ## A simple example
#' data <- matrix(rnorm(100), nrow = 10)
#' window <- 10
#' step <- 5
#' lambda <- 0.1
#' calcAdjFrag(ieegts = data, window = window,
#' step = step, lambda = lambda, progress = TRUE)
#'
#' ## A more realistic example with parallel computing
#' \dontrun{
#' ## Register a SNOW backend with 4 workers
#' library(parallel)
#' library(doSNOW)
#' cl <- makeCluster(4, type="SOCK")
#' registerDoSNOW(cl)
#' 
#' data("pt01Epoch")
#' window <- 250
#' step <- 125
#' title <- "PT01 seizure 1"
#' calcAdjFrag(ieegts = pt01Epoch, window = window,
#'   step = step, parallel = TRUE, progress = TRUE)
#' 
#' ## stop the parallel backend
#' stopCluster(cl)
#' }
#' 
#'
#'
#' @details
#' 1/ For each time window i, a discrete stable Linear time system
#' (adjacency matrix) is computed named \eqn{A_i}
#' such that
#' \eqn{A_i x(t) = x(t+1)}. The 'lambda' option is the regularization parameter
#' for the ridge regression. 
#' `lambda=NULL`(default) will find a lambda value that ensures 
#' the stability of the estimated \eqn{A_i}.
#'
#' 2/For each stable estimated \eqn{A_i}, the minimum norm perturbation \eqn{\Gamma_{ik}} (k index of the electrodes)
#' for column perturbation is computed.
#' Each column is normalized \eqn{\frac{max(\Gamma_{i})-\Gamma_{ik}}{max(\Gamma_i)}}
#'
#' @export
calcAdjFrag <- function(ieegts, window, step, lambda = NULL, nSearch = 100L, progress = FALSE, parallel = FALSE) {
    ## check the input types
    stopifnot(is.matrix(ieegts) | is(ieegts, "Epoch"))
    stopifnot(isWholeNumber(window))
    stopifnot(isWholeNumber(step))
    stopifnot(step > 0)
    stopifnot(is.null(lambda) | is.numeric(lambda))


    if (is(ieegts, "Epoch")) {
        epoch <- ieegts
        ieegts <- ieegts$matrix
    }else {
        epoch <- Epoch(ieegts)
    }

    ## The input matrix must have at least window rows
    stopifnot(nrow(ieegts) >= window)
    scaling <- 10^floor(log10(max(abs(ieegts))))
    ieegts  <- ieegts / scaling
    # Electrode count and names
    elCnt <- ncol(ieegts)
    elNms <- colnames(ieegts)
    # Number/sequence of steps
    nParts  <- floor((nrow(ieegts) - window) / step) + 1L
    partitions   <- seq_len(nParts)
    # Pre-allocate output
    dm   <- c(elCnt, elCnt, nParts)
    dmn  <- list(Electrode  = elNms, Step = partitions)
    dmnA <- list(Electrode1 = elNms, Electrode2 = elNms, Step = partitions)
    
    # Indices of window at time 0
    i0 <- seq_len(window - 1L)

    # Pre-allocate output
    A    <- array(.0, dim = dm,     dimnames = dmnA)
    R2   <- array(.0, dim = dm[-1], dimnames = dmn)
    f = fR <- R2
    lbd <- rep(0, nParts) |> setNames(partitions)
    
    ## switch between parallel and sequential computing
    if(parallel){
        `%run%` <- foreach::`%dopar%`
    }else{
        `%run%` <- foreach::`%do%`
    }

    ## initialize the progress bar
    if (progress){
        pb <- progress_bar$new(
        format = "Step = :current/:total [:bar] :percent in :elapsed | eta: :eta",
        total = nParts, 
        width = 60)
        
        progress <- function(n){
            pb$tick()
        } 
        opts <- list(progress = progress)
        on.exit(pb$terminate())
    }else{
        opts <- list()
        
        progress <- function(n){
            NULL
        } 
    }

    ## Initial data for data aggregation
    init <- list(A = A, R2 = R2, f = f, lbd = lbd)
    foreach(
        iw = partitions, 
        .combine = .combine, 
        .init = init, 
        .inorder = FALSE,
        .options.snow = opts
        ) %run% {
            progress(iw)
            ## Not necessary, but for clear the R check error
            iw <- get('iw')
            
            si   <- i0 + (iw - 1L) * step
            xt   <- ieegts[si, , drop = FALSE]
            xtp1 <- ieegts[si + 1L, , drop = FALSE]
            
            adjMatrix <- EZFragility:::ridgeSearch(xt, xtp1, lambda)
            R2Column <- EZFragility:::ridgeR2(xt, xtp1, adjMatrix)
            fColumn  <- EZFragility:::fragilityRow(adjMatrix, nSearch)

            list(
                iw = iw, 
                adjMatrix = adjMatrix, 
                R2Column = R2Column, 
                fColumn = fColumn
            )
    } -> results

    ## unpack the results
    A <- results$A
    R2 <- results$R2
    f <- results$f
    lbd <- results$lbd
    
    ## column rank of fragility matrix
    fR <- apply(f, 2, rank) / elCnt

    ## start time point/indices for each partition
    startTimes <- (seq_len(nParts) - 1L) * step + 1L
    if (!is.null(epoch$times)) {
        startTimes <- epoch$times[startTimes]
    }
    
    ## TODO: why the row of frag is the electrode names? not the column of frag?
    ## This does not match the input data
    Fragility(
        ieegts = ieegts,
        adj = A,
        R2 = R2,
        frag = f,
        frag_ranked = fR,
        lambdas = lbd,
        startTimes = startTimes,
        electrodes = epoch$electrodes
    )
}

.combine <- function(results, x){
    iw <- x$iw
    adjMatrix <- x$adjMatrix
    R2Column <- x$R2Column
    fColumn <- x$fColumn

    results$A[,, iw]  <- adjMatrix
    results$R2[, iw]  <- R2Column
    results$f[,  iw]  <- fColumn
    results$lbd[[iw]] <- attr(adjMatrix, "lambda")
    results
}

#' Find Serzure Onset Zone
#' 
#' @param x Fragility object
#' @param method Character. The method to use to find the onset zone.
#' Must be one of 'max', 'min', or "mean"
#' @param proportion Numeric. The proportion of electrodes to consider as the onset zone. 
#' The electrode number will be rounded to the nearest integer.
#' @param ... Additional arguments
#' 
#' @return A vector of electrode names, or indices if the electrode names are NULL
#' @export 
findOnset <- function(x, method = c("mean","max", "min"), proportion = 0.1, ...) {
    method <- match.arg(method)

    stopifnot(is(x, "Fragility"))

    frag <- x$frag
    elCnt <- nrow(frag)
    nSOZ <- round(elCnt * proportion)
    stopifnot(nSOZ > 0 & nSOZ <= elCnt)

    if (method == "max") {
        stat <- apply(frag, 2, max)
    } else if (method == "min") {
        stat <- apply(frag, 2, min)
    } else if (method == "mean") {
        stat <- apply(frag, 2, mean)
    }

    sozIndex <- order(stat, decreasing = TRUE)[seq_len(nSOZ)]
    if(!is.null(frag)){
        sozIndex <- rownames(frag)[sozIndex]
    }

    sozIndex
}