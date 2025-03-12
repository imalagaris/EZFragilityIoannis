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
#' with time points as rows and electrodes names as columns
#' @param window Integer. The number of time points to use in each window
#' @param step Integer. The number of time points to move the window each time
#' @param lambda Numeric. The lambda value to use in the ridge regression.
#' If NULL, the lambda will be chosen automatically
#' ensuring that ensuring that the adjacent matrix is stable (see details)
#' @param nSearch Integer. Number of minimization to compute the fragility row
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
#' step = step, lambda = lambda)
#'
#' ## A more realistic example, but it will take a while to run
#' \dontrun{
#' data("pt01Epochm1sp2s")
#' window <- 250
#' step <- 125
#' lambda <- NULL
#' nSearch <- 100
#' title <- "PT01 seizure 1"
#' resfrag <- calcAdjFrag(ieegts = pt01Epochm1sp2s, window = window,
#'   step = step, lambda = lambda,nSearch=nSearch)
#' }
#'
#'
#' @details
#' 1/ For each time window i, a discrete stable Linear time system
#' (adjacency matrix) is computed named \eqn{A_i}
#' such that
#' \eqn{A_i x(t) = x(t+1)}
#' option Lambda=NULL ensures that the matrix is stable
#'
#' 2/For each stable estimated \eqn{A_i}, the minimum norm perturbation \eqn{\Gamma_{ik}} (k index of the electrodes)
#' for column perturbation is computed.
#' Each column is normalized \eqn{\frac{max(\Gamma_{i})-\Gamma_{ik}}{max(\Gamma_i)}}
#'
#' @export
calcAdjFrag <- function(ieegts, window, step, lambda = NULL, nSearch = 100L) {
    stopifnot(isWholeNumber(window))
    stopifnot(isWholeNumber(step))
    stopifnot(is.null(lambda) | is.numeric(lambda))
    stopifnot(nrow(ieegts) >= window)
    self <- environment()
    ConLogger <- ConsoleLogger(self)
    scaling <- 10^floor(log10(max(ieegts)))
    ieegts  <- ieegts / scaling
    # Electrode count and names
    elCnt <- ncol(ieegts)
    elNms <- colnames(ieegts)
    # Number/sequence of steps
    nsteps  <- floor((nrow(ieegts) - window) / step) + 1L
    STEPS   <- seq_len(nsteps)
    # Pre-allocate output
    dm   <- c(elCnt, elCnt, nsteps)
    dmn  <- list(Electrode  = elNms, Step = STEPS)
    dmnA <- list(Electrode1 = elNms, Electrode2 = elNms, Step = STEPS)
    A    <- array(.0, dim = dm,     dimnames = dmnA)
    R2   <- array(.0, dim = dm[-1], dimnames = dmn)
    f = fR <- R2
    lbd <- rep(0, nsteps) |> setNames(STEPS)
    # Indices of window at time 0
    i0 <- seq_len(window - 1L)
    ConLogger$ProcessInfo()
    for (iw in STEPS) {
        si   <- i0 + (iw - 1L) * step
        xt   <- ieegts[si, ]
        xtp1 <- ieegts[si + 1L, ]
        ConLogger$RidgeStart()
        adjMatrix <- ridgesearchlambdadichomotomy(xt, xtp1, lambda)
        A[,, iw]  <- adjMatrix
        R2[, iw]  <- ridgeR2(xt, xtp1, adjMatrix)
        ConLogger$FragStart()
        f[,  iw]  <- fragilityRow(adjMatrix, nSearch)
        fR[, iw]  <- rank(f[, iw]) / elCnt # ranks should probably be here...
        lbd[[iw]] <- attr(adjMatrix, "lambda")
        ConLogger$StepEnd()
    }
    ConLogger$TotalTime()
    Fragility(
        ieegts = ieegts,
        adj = A,
        R2 = R2,
        frag = f,
        frag_ranked = fR,
        lambdas = lbd
    )
}


# Get number of seconds from a reference timestamp in specified format
getTimeSecs <- \(ref) {
    difT <- difftime(Sys.time(), ref, units = "secs") |> as.double()
    sprintf("%.2f", difT)
}

# Optional utility which prints info while the above is running.
ConsoleLogger <- \(e) {
    start <- Sys.time()
    stepStart = TimeKeeper = RidgeRun = FragRun <-  NULL;
    InfoDash = TabDash <- NULL
    window <- e$window;
    step <- e$step;
    samples <- nrow(e$ieegts)
    self <- environment()
    
    # Save time of Ridge (and step) start and print step number (run before ridge)
    RidgeStart <- \() {
        self$stepStart <- Sys.time()
        self$TimeKeeper <- self$stepStart
        sprintf("%7d ", e$iw) |> cat()
    }
    # Save time of Fragility start and print Ridge runtime (run before fragility)
    FragStart <- \() {
        self$RidgeRun <- getTimeSecs(self$TimeKeeper)
        self$TimeKeeper <- Sys.time()
        sprintf("%11s ", self$RidgeRun) |> cat()
    }
    # Print runtime for Fragility and for the whole step (run at the end of the step)
    StepEnd <- \() {
        sprintf("%11s ", getTimeSecs(self$TimeKeeper)) |> cat()
        sprintf("%7s",   getTimeSecs(self$stepStart))  |> cat("\n")
    }
    # Print the total runtime of the whole process
    TotalTime <- \() {
        total <- Sys.time() - start
        fmt <- "  Total Runtime:  %21.2f %s"
        self$StepDash |> shift(2) |> cat("\n")
        sprintf(fmt, as.double(total), attr(total, "units")) |> cat("\n")
    }
    # Prints process specifications
    initTab <- \(nb = 2) {
        header <- list("Samples", "Window", "Shift", "Steps")
        values <- list(samples, window, step, e$nsteps)
        hFmt <- do.call(sprintf, c(list(" %7s |%7s |%6s |%6s "), header))
        vFmt <- do.call(sprintf, c(list(" %7d  %7d  %6d  %6d "), values))
        self$InfoDash <- paste(rep("-", nchar(hFmt)), collapse = "")
        c(self$InfoDash, hFmt, vFmt) |> shift(nb) |> cat(sep = "\n")
    }
    # Prints the headers of the step time table
    stepTab <- \(nb = 2) {
        header <- c("Step", "Adjacency", "Fragility", "Total")
        hFmt <- do.call(sprintf, c(list("%5s | %9s | %9s | %5s"), header))
        len <- nchar(hFmt)
        self$StepDash <- paste(rep("-", len + 1), collapse = "")
        subTitle <- paste(rep("-", len - 12), collapse = "")
        subtext <- "Runtime (sec)"
        id <- (1 + 0.5 * (nchar(subTitle) - nchar(subtext))) |> ceiling()
        substr(subTitle, id, id + nchar(subtext)) <- subtext
        c(self$StepDash, shift(subTitle, 6), hFmt) |> shift(nb) |> cat(sep = "\n")
    }
    ProcessInfo <- \() { initTab(); stepTab() }
    self
}
