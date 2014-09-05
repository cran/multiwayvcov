###################################################################### 
#' @title Bootstrapped multi-way standard error clustering
#' 
#' @description Return a bootstrapped multi-way cluster-robust variance-covariance matrix
#'
#' @param model The estimated model, usually an \code{lm} or \code{glm} class object
#' @param cluster A vector, \code{matrix}, or \code{data.frame} of cluster variables,
#' where each column is a separate variable.  If the vector \code{1:nrow(data)}
#' is used, the function effectively produces a regular 
#' heteroskedasticity-robust matrix.
#' @param parallel Scalar or list.  If a list, use the list as a list
#' of connected processing cores/clusters.  A scalar indicates no
#' parallelization.  See the parallel package.
#' @param use_white \code{Logical} or \code{NULL}.  See description below.
#' @param force_posdef Logical.  Force the eigenvalues of the variance-covariance
#' matrix to be positive.
#' @param R \code{Integer}. The number of bootstrap replicates; passed directly to \code{boot}.
#' @param debug \code{Logical}.  Print internal values useful for debugging to 
#' the console.
#'
#' @keywords clustering multi-way robust standard errors bootstrap boot block
#' 
#' @details
#' This function implements cluster bootstrapping (also known as the block bootstrap)
#' for variance-covariance matrices, following Cameron, Gelbach, & Miller (CGM) (2008).
#' Usage is generally similar to the \code{cluster.vcov} function in this package, but this
#' function does not support degrees of freedome corrections or leverage adjustments.
#' 
#' In the terminology that CGM (2008) use, this function implements
#' \emph{pairs cluster bootstrap-se}.
#' 
#' Multi-way clustering is handled as described by Petersen (2009) and generalized
#' according to Cameron, Gelbach, & Miller (2011).  This means that cluster.boot
#' estimates a set of variance-covariance matrices \emph{for the variables} separately 
#' and then sums them (subtracting some matrices and adding others).
#' The method described by CGM (2011) estimates a set of variance-covariance matrices
#' \emph{for the residuals} (sometimes referred to as the meat of the sandwich estimator)
#' and sums them appropriately.  Whether you sum the meat matrices and then compute
#' the model's variance-covariance matrix or you compute a series of model matrices
#' and sum those is mathematically irrelevant, but may lead to (very) minor numerical
#' differences.
#' 
#' Unlike the \code{cluster.vcov} function, this function does not depend upon the \code{estfun}
#' function from the \bold{sandwich} package, although it does make use of the \code{vcovHC}
#' function for computing White HC0 variance-covariance matrices. 
#' 
#' Parallelization (if used) is handled by the \bold{boot} package.
#' 
#' @return a \eqn{k} x \eqn{k} variance-covariance matrix of type \code{matrix}
#' 
#' @export
#' @author Nathaniel Graham \email{npgraham1@@gmail.com}
#' 
#' @references
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based improvements 
#' for inference with clustered errors. The Review of Economics and Statistics, 90(3), 414-427.
#' 
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust inference with multiway 
#' clustering. Journal of Business & Economic Statistics, 29(2).
#' 
#' Petersen, M. A. (2009). Estimating standard errors in finance panel data sets: Comparing 
#' approaches. Review of financial studies, 22(1), 435-480.
#' 
#' @importFrom sandwich vcovHC vcovHC.default
#' @importFrom boot boot
#' @importFrom compiler cmpfun
#' 
#' @examples
#' \dontrun{
#' library(lmtest)
#' data(petersen)
#' m1 <- lm(y ~ x, data = petersen)
#' 
#' # Cluster by firm
#' boot_firm <- cluster.boot(m1, petersen$firmid)
#' coeftest(m1, boot_firm)
#' 
#' # Cluster by year
#' boot_year <- cluster.boot(m1, petersen$year)
#' coeftest(m1, boot_year)
#' 
#' # Double cluster by firm and year
#' boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year))
#' coeftest(m1, boot_both)
#' 
#' # Go multicore using the parallel package
#' require(parallel)
#' cl <- makeCluster(4)
#' boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year), parallel = cl)
#' stopCluster(cl)
#' coeftest(m1, boot_both)
#' }
cluster.boot <- function(model, cluster, parallel = FALSE, use_white = NULL, 
                         force_posdef = FALSE, R = 300,
                         debug = FALSE) {
  
  cluster <- as.data.frame(cluster)
  cluster_dims <- ncol(cluster)
  
  # total cluster combinations, 2^D - 1
  tcc <- 2 ** cluster_dims - 1
  # all cluster combinations
  acc <- list()
  for(i in 1:cluster_dims) {
    acc <- append(acc, combn(1:cluster_dims, i, simplify = FALSE))
  }
  if(debug) print(acc)
  
  # We need to subtract matrices with an even number of combinations and add
  # matrices with an odd number of combinations
  vcov_sign <- sapply(acc, function(i) (-1) ** (length(i) + 1))
  
  # Drop the original cluster vars from the combinations list
  acc <- acc[-1:-cluster_dims]
  if(debug) print(acc)
  
  # Handle omitted or excluded observations
  if(!is.null(model$na.action)) {
    if(class(model$na.action) == "exclude") {
      cluster <- cluster[-model$na.action,]
    } else if(class(model$na.action) == "omit") {
      cluster <- cluster[-model$na.action,]
    }
    cluster <- as.data.frame(cluster)  # silly error somewhere
  }
  if(debug) print(class(cluster))
  
  # Make all combinations of cluster dimensions
  if(cluster_dims > 1) {
    for(i in acc) {
      cluster <- cbind(cluster, Reduce(paste0, cluster[,i]))
    }
  }
  
  df <- data.frame(M = integer(tcc),
                   N = integer(tcc),
                   K = integer(tcc))
  
  for(i in 1:tcc) {
    df[i, "M"] <- length(unique(cluster[,i]))
    df[i, "N"] <- length(cluster[,i])
    df[i, "K"] <- model$rank
  }
  
  if(is.null(use_white)) {
    if(cluster_dims > 1 && df[tcc, "M"] == prod(df[-tcc, "M"])) {
      use_white <- TRUE
    } else {
      use_white <- FALSE
    }
  }
  
  if(use_white) {
    df <- df[-tcc,]
    tcc <- tcc - 1
  }
  
  if(debug) {
    print(acc)
    print(paste("Original Cluster Dimensions", cluster_dims))
    print(paste("Theoretical Cluster Combinations", tcc))
    print(paste("Use White VCOV for final matrix?", use_white))
  }
  

  if("model" %in% names(model)) {
    X <- model$model
  } else {
    X <- model.frame(model)
  }
  
  boot.outs <- list()
  
  par_type <- "no"
  if(length(parallel) > 1) {
    clusterExport(parallel, varlist = c("cluster", "model"), envir = environment())
    par_type <- "snow"
  }
  
  dat_loc <- which(names(model$call) == "data")
  args <- model$call[c(-1, -dat_loc)] 
  estimator <- eval(model$call[[1]])
  
  est.func <- cmpfun(function(grp, i, data2, clustvar, arglist, estimator) {
    j <- unlist(lapply(i, function(n) which(n == clustvar)))
    coef(estimator(arglist, data = data2[j,]))
  })
  
  for(i in 1:tcc) {
    boot.outs[[i]] <- boot(unique(cluster[,i]), est.func, R = R,
                           parallel = par_type, cl = parallel,
                           data2 = X, clustvar = cluster[,i],
                           arglist = args, estimator = estimator)
  }
  
  if(debug) {
    print(df)
    print(vcov_sign)
  }
  
  vcov_matrices <- list()
  for(i in 1:tcc) {
    vcov_matrices[[i]] <- vcov_sign[i] * cov(boot.outs[[i]]$t)
  }
  
  if(use_white) {
    i <- i + 1
    vcov_matrices[[i]] <- vcov_sign[i] * vcovHC(model, type = "HC0")
  }
  
  if(debug) {
    print(vcov_matrices)
  }
  
  vcov_matrix <- Reduce('+', vcov_matrices)
  
  if(force_posdef) {
    decomp <- eigen(vcov_matrix, symmetric = TRUE)
    if(debug) print(decomp$values)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    vcov_matrix <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }
  
  return(vcov_matrix)
}