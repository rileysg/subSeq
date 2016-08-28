# this file contains functions that are used in two ways

# a function to estimate proportions of "truely null" hypothesis.
#This should be performed on the results of a single experiment
# PVALS an n vector of pvalues
#METHODS an n vector of methods used to perform differential expression testing
estimatePi0 <- function(pvals, methods){
  #TODO(riley) it would be nice to check that only a single experiment is being analysed
  if ( !is.numeric(pvals))
    stop("PVALS should be supplied as a numeric vector")
  if ( !(is.character(methods) || is.factor(methods)))
    stop("METHODS should be supplied as character or factor")
  if ( length(methods) != length( pvals))
    stop( "there must be equal numbers of pvals and methods, because the methods describe how the pvalue was found")
  ret <- data.frame(pvalue= pvals, method= as.character(methods)) %>%
    group_by(method) %>%
    summarize( pi0=qvalue::qvalue(pvalue, lambda = seq(0.05,0.9, 0.05))$pi0) %>% group_by()
  return( ret)
}

getOracle <- function(methods, counts, design){
  out <- subsample(counts, design, method=methods, replications = 1,
            subsampleDesigner= identityDesign,
            chunker= fixedLengthChunker(length = 1),
            writer = NULL, getQvalues= NULL)
  
  out0 = out %>% group_by(method) %>%
    summarize(pi0=qvalue::qvalue(pvalue, lambda = seq(0.05,0.9, 0.05))$pi0) 
  # out = out %>% inner_join(out0, by = c("method")) %>%
  #   mutate(qvalue=qvalue.fixedpi0(pvalue, pi0s=unique(pi0))) %>% group_by()

  out = out %>% inner_join(out0, by = c("method")) %>% group_by(method) %>%
    mutate(qvalue=qvalue.fixedpi0(pvalue, pi0s= unique(pi0))) %>% group_by()
  }
  

#construct a function that generates q-values from p-values using a using a "gold standard" estimate of pi0,
# the gold standard estimate depends on supplying counts and designs, and estimationg pi0 using the functions
# provided in methods
# VALUES
# a function that takes arguments, pvals and method, and returns q-vals
commonPi0Est <- function(methods, counts, design){
  #estimate a common pi0
  # run subseq without subsampling using all of the methods provided in METHODS
  ora <- getOracle(methods, counts, design)
  function(pvals, method){
    #check that a pi0 estimate is defined for METHOD
    #use only the gold standard estimate corresponding to METHOD
    frm <- data_frame( method, pvals) %>%
      inner_join(distinct(select(ora, method, pi0)), by = "method") %>%
      group_by( method) %>%
      mutate(qvalue=qvalue.fixedpi0(pvals, pi0s=unique(pi0))) %>% group_by() %>%
      select( pi0, qvalue) %>%
      as_data_frame()
  }
}

# a function to do estimation of pi0 and q-values in a subsample-wise manner
exptPi0Est <- function(pvals, method){
  pi0 <- estimatePi0(pvals, method)$pi0
  frm <- data_frame( pi0, pvals) %>%
    mutate(qvalue = qvalue.fixedpi0(pvals, pi0s=unique(pi0))) %>% group_by() %>%
    select( pi0, qvalue) %>%
    as_data_frame()
}

# OLD DGR CODE ----------------------------------------
# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when 
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.
qvalue.filtered1 = function(p) {
  # given a vector of p-values
  q.without1 = qvalue(p[p < 1])
  # when p == 1, the FDR is pi0
  q = rep(q.without1$pi0, length(p))
  q[p < 1] = q.without1$qvalue
  q
}


# NEW AJB CODE ----------------------------------------
#Estimate the common pi0,

# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.

#TODO(riley) what is the purpose of the pfdr flag

qvalue.fixedpi0 <- function (p, fdr.level = NULL, pfdr = FALSE, pi0s=NULL, ...)
{
  # initalize values ---------------------------------
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  # check for valid input  ---------------------------------
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  }
  else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level >
                                   1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  m <- length(p)
  u <- order(p)
  v <- rank(p, ties.method = "max")
  if (pfdr) {
    qvals <- (pi0s * m * p)/(v * (1 - (1 - p)^m))
  }
  else {
    qvals <- (pi0s * m * p)/v
  }
  qvals[u[m]] <- min(qvals[u[m]], 1)
  for (i in (m - 1):1) {
    qvals[u[i]] <- min(qvals[u[i]], qvals[u[i + 1]])
  }
  qvals_out[rm_na] <- qvals
  lfdr <- lfdr(p = p, pi0 = pi0s, ...)
  lfdr_out[rm_na] <- lfdr
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level))
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0s, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out)
  }
  return(qvals_out)
}