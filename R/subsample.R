#' Subsample reads and perform statistical testing on each sample
#'
#' Perform subsampling at multiple proportions on a matrix of count
#' data representing mapped reads across multiple samples in many
#' genes. For each sample, perform some statistical operations.
#'
#' @param counts Matrix of unnormalized counts
#' @param treatments Vector (factor) of experimental treatments corresponding to
#' collumns of counts
#' @param proportions Vector of subsampling proportions in (0, 1]
#' @param bioReplicates Vector specifying number of samples from each treatment used in subsampling
#' @param method One or more methods to be performed at each subsample,
#' such as edgeR or DESeq (see Details)
#' @param replications Number of replications to perform at each depth
#' @param seed An initial seed, which will be stored in the output
#' so that any individual simulation can be reproduced.
#' @param qvalues Whether q-values should be calculated for multiple hypothesis
#' test correction at each subsample.
#' @param env Environment in which to find evaluate additional hander functions
#' that are given by name
#' @param ... Other arguments given to the handler, such as \code{treatment}
#'
#' @return A subsample S3 object, which is a data.table containing
#'
#' \item{pvalue}{A p-value calculated for each gene by the handler}
#' \item{coefficient}{An effect size (usually log fold change) calculated
#' for each gene by the handler}
#' \item{ID}{gene ID}
#' \item{count}{the number of reads to this specific gene in this subsample}
#' \item{depth}{the overall sequencing depth of this subsample}
#' \item{method}{the method used (the name of the handler)}
#'
#' @details Method represents the name of a handler function, which can be
#' custom-written by the user.
#'
#' If a gene has a count of 0 at a particular depth, we set the p-value to 1
#' and the coefficient to 0 to stay consistent between programs. If the gene has
#' a count that is not 0 but the p-value is NA, we set the p-value to 1 but
#' keep the estimated coefficient.
#'
#' @examples
#'
#' data(hammer)
#'
#' hammer.counts = Biobase::exprs(hammer)[, 1:4]
#' hammer.design = Biobase::pData(hammer)[1:4, ]
#' hammer.counts = hammer.counts[rowSums(hammer.counts) >= 5, ]
#' 
#' ss = subsample(counts=        hammer.counts,
#'                treatments =   hammer.design$protocol,
#'                proportions=   c(.01, .1, 1),
#'                bioReplicates= c(2),
#'                method=        c("edgeR", "voomLimma"))
#'
#' @import data.table
#' @importFrom dplyr group_by do mutate filter
#' @import magrittr
#' @importFrom qvalue qvalue lfdr
#' @importFrom Biobase exprs pData
#' @export
subsample <-
  function(counts, treatments, proportions, bioReplicates, method="edgeR", replications=1,
           replacement= FALSE, ballancedproportions= FALSE, seed=NULL,
           qvalues = TRUE, use.common.pi0= TRUE, env=parent.frame(), ...) {
    # error checking ------------------------------
    if (length(proportions) == 0) {
      stop("No proportions to sample")
    }
    if (length( bioReplicates) == 0){
      stop("No numbers of bioReplicates to sample")
    }
    if (any(proportions > 1 | proportions == 0)) {
      stop("Proportions must be in range (0, 1]")
    }
    if ( length( treatments) != dim( counts)[2]){
      stop("Collumns of counts and elements of treatments should correspond
           so they should have the same length")
    }
    if ( any( !( bioReplicates %in% seq(2, dim(counts)[2])))){
      stop( "The number of biological replicates must be greater
            than one and less than the number of samples")
    }
    
    # check that the counts is an unnormalized numeric matrix
    counts = as.matrix(counts)
    error = max(abs(round(counts) - counts))
    if (error > 1e-4 | any(counts < 0)) {
      stop("Counts should be unnormalized integer counts (not, e.g., RPKM)")
    }
    if (is.null(seed)) {
      # come up with a random initial seed (can't use current one, since
      # the point is to save it for later)
      seed = sample(.Machine$integer.max, 1)
    }
    
    # create a list mapping function name to function
    if (is.function(method)) {
      # if given a single function, make that into a list
      methods = list(method)
      names(methods) = deparse(substitute(method))
    } else {
      methods = lapply(method, function(m) {
        # function can be given directly: otherwise, check subSeq package
        # then caller's env
        if (is.function(m)) {
          handler = m
        }
        #         else if (exists(m, mode="function", envir=as.environment("package:subSeq"))) {
        #            handler = get(m, mode="function", envir=as.environment("package:subSeq"))
        #       }
        else if (exists(m, mode="function", envir=env)) {
          handler = get(m, mode="function", envir=env)
        }
        else if (m %in% c("edgeR", "edgeRDispersions", "voomLimma", "DESeq2", "edgeR.glm")) {
          handler = get(m, mode="function")
        }
        else {
          stop(paste("Could not find handler", m))
        }
        handler
      })
      names(methods) = method
    }
    # Specify parameters for resampling------------------------------
    # perform one for each method x proportion x bioReplicates x replication
    params = expand.grid(method=names(methods), proportion=proportions, biological.replicate= bioReplicates, replication=1:replications)
    # If all treatments are sampled at highest number of biological replicates,
    # replications of of full depth might be identical
    all.max.brep <- all( table( treatments) == max( bioReplicates))
    if ( !replacement && !ballancedproportions && all.max.brep){
      params = params %>% filter(!(proportion == 1 & biological.replicate == max( biological.replicate) & replication > 1))
    }
    # apply method to each subsample the specified number of times
    #m.ret = as.data.table(do.call(rbind, lapply(1:nrow(prop.reps),
    perform.ballanced.subsampling <- function(method, proportion, bioReplicates, replication) {
      # resample biological replicates: specify which biological replicates will be used in the subsampling
      inds <- getIndecesFromCatagoricalTreatment( treatments, bioReplicates, replication, seed, replacement)
      treatment <- treatments[inds]
      #Get proportions for each index in inds
      #Calculating colsums once would help efficientcy, but it needs to be done in correct place to be readable
      if( ballancedproportions == TRUE){
        total.counts <- colSums( counts)
        #ind.proportion * total.counts == min( total.counts)
        ind.proportions <- proportion * min( total.counts) / total.counts[ inds]
      } else {
        ind.proportions <- rep( proportion, length.out= length( inds))
      }
      # subsample reads and use handler
      subcounts = generateSubsampledMatrix(counts, inds, ind.proportions, seed, replication)
      id = which(rowSums(subcounts) >= 5)
      subcounts = subcounts[id,] ### Filter counts at zero
      if (length(id) == 0) return("Error: counts too low at subsampling proportion")
      
      handler = methods[[method]]
      ret = handler(subcounts, treatment, ...)
      ## add gene names (ID) and per-gene counts
      ## If there is one handler row per gene, the remaining collumns of ret, "ID" and "counts", can be infered.
      infer.per.gene = dim( ret)[1] == dim( subcounts)[1]
      if ( !any( "ID" == colnames(ret))){
        if (infer.per.gene){
          ret$ID = rownames(subcounts)
        } else {
          stop("if a handler doesn't return one row per gene then it must specify an ID collumn")
        }
      }
      if ( !any( "count" == colnames(ret))){
        if (infer.per.gene){
          ret$count = as.integer(rowSums(subcounts))
        } else {
          ret$count = NA
        }
      }
      ret$depth = sum(subcounts)
      
      # in any cases of no reads, fix coefficient/pvalue to 0/1
      if ("count" %in% colnames(ret)) {
        if ("pvalue" %in% colnames(ret)) {
          ret$pvalue[ret$count == 0 | is.na(ret$pvalue) | ret$pvalue == 1] = NA #new qvalue accepts NA values
        }
        if ("coefficient" %in% colnames(ret)) {
          ret$coefficient[ret$count == 0 | is.infinite(ret$coefficient)] = 0
        }
      }
      ret
    }

    ret = params %>% group_by(method, proportion, biological.replicate, replication) %>%
      do(perform.ballanced.subsampling( .$method, .$proportion, .$biological.replicate, .$replication))

    ## cleanup
    if (qvalues) {
      # calculate q-values
      if( use.common.pi0){ #a switch for experimenting with different pi0 estimation methods
        max.proportion <- max( ret$proportion)
        max.bio.rep <- max( ret$biological.replicate)
        # TODO(riley) is it necesary to pick a particular replicate, or
        #will the estimate of pi0 be more accurate when performed using all of the replicates?
        rep.id <- ret %>% filter(proportion == max.proportion, biological.replicate == max.bio.rep) %>%
          group_by() %>%
          sample_n(1)
        ret0 = ret %>% filter(proportion == max.proportion,
                              biological.replicate == max.bio.rep,
                              replication == rep.id$replication) %>%
          group_by( method) %>%
          summarize(pi0=qvalue::qvalue(pvalue, lambda = seq(0.05,0.9, 0.05))$pi0) %>% group_by()
        ret = ret %>% inner_join(ret0, by = c("method"))
        ret = ret %>% group_by(proportion, biological.replicate, method, replication) %>%
          mutate(qvalue=qvalue.fixedpi0(pvalue, pi0s=unique(pi0))) %>% group_by()
      }
      else {
        ret.pi0 = ret %>% group_by(proportion, biological.replicate, method, replication) %>%
          summarize(pi0=qvalue::qvalue(pvalue, lambda = seq(0.05,0.9, 0.05))$pi0) %>% group_by()
        ret = ret %>% inner_join( ret.pi0, by= c("proportion", "biological.replicate", "method", "replication"))
        ret = ret %>% group_by(proportion, biological.replicate, method, replication) %>%
          mutate(qvalue=qvalue.filtered1(pvalue)) %>% group_by()
      }

    }

    # turn into a subsamples object
    ret = as.data.table(as.data.frame(ret))
    class(ret) = c("subsamples", class(ret))
    attr(ret, "seed") = seed

    return(ret)
  }

# a function that returns a subsample function
# a subsample function should do the following:
# 1. Check for errors in input
# 2. Get parameters to fully specify the subsampling: 
#   * getSamplingParams
#     * INPUTS: ( treatments, bioReplicates, replication, seed, replacement)
#     * OUTPUTS: ( a matrix of subsampling proportions. Rows correspond to subsamplings, collumns correspond to samples)
#                 an additional collumn of the matrix indexes replications with a common design
# 3. Generate subsampled matrices:                    
#   * mySubsample
#     * INPUTS: (counts, proportions, seed, replication)
#     * OUTPUTS: count matrix with gene names
# 4. Perform a test for differential gene expression: should return "coefficient" and "p-value" collumns
#   * handler
#     * INPUTS: counts, design
#     * OUTPUTS: ID, coefficients, p-value, count,
# 5. Get q-values, if necesary:                       should add a column to the result of testing for differential gene expression
#   * qvalues
#     * INPUTS: p-value
#     * OUTPUTS: q-value
# 6. Summarize a subsampling
#   * getSsAttributes
#     * INPUTS:   subsampled counts, design
#     * OUTPUTS:  total reads (depth), number of samples (samples)
# 7. Combine descriptors of a subsampling ( 2,6) with descriptors that are particular to each collumn (4,6)

require( dplyr)
require( digest)
subsampleFactory <- function(getSamplingParams, mySubsample, handlers, qvalues){
  function(counts, design){
    # Error checking ------------------------------
    # TODO(riley): make sure that each function accepts the correct arguments
    if (dim( design)[1] != dim( counts)[2])
      stop("the design matririx should have rows corresponding to samples")
    if ( is.null(rownames( counts))){
      warning("the count matrix does not have IDs specified using rownames. Indeces will be used")
      rownames( counts) <- as.character(1:dim(counts)[1])
    }
    if ( !is.list(handlers)){ 
      stop("handlers should be a list")
    } else {
      if ( is.null( names( handlers)))
        names( handlers) <- deparse(substitute(handlers))
    }
    #Error checking should be the same for all types of subseq functions
    
    # Subsampling Design ------------------------------
    # get a matrix of proportions with rows corresponding to subsamplings and collumns corresponding to samples
    subsampling.design <- getSamplingParams( design = design, counts = counts) #other parameters should include: seed, proportion, ballancing, replication, and replacement
    replications <- subsampling.design %>%
      group_by_( paste0(colnames(subsampling.design), collapse= "," )) %>%
      transmute_( replication = "row_number()")

      mutate_( replication= n())
      mutate_( replication = dense_rank("a")) 
    %>%
      group_by_( col)
    %>%
      do_()
    
      
    # Subsampling  ------------------------------
    result.array <- apply( subsampling.design, 1, function(props){
      # generate a subsampled matrix of counts according to a row of subsampling.design
      subcounts <- mySubsample(counts = counts, props = props)
      s.remaining <- which( props != 0)
      #TODO(riley) this is assuming that the design matrix has rows corresponding to samples
      subcounts.design <- design[s.remaining, , drop = FALSE]
      # return an analysis using each handler that has been provided
      analysis.list <- lapply( names(handlers), function(name){
        # for this to work, the sample orderings of subcounts and subcounts.design have to be consistent
        dex = handlers[[name]](subcounts, subcounts.design)
        dex = cbind( method= rep( name, times= dim(dex)[1]), dex)
        ## subseq no longer infers colnames from gene names or returns counts, the handler should do this itself...
        dex$qvals <- qvalues( dex)
        dex$depth = sum(subcounts)
        return(dex)
      })
      if ( length(analysis.list) > 1)
        analysis <- bind_rows( analysis.list)
      else
        analysis <- analysis.list[[1]]
      analysis
      })
    num.ss <- dim(result.array)[3]
    num.col <- dim( result.array)[2]
    num.
    apermarray( result.array, )
    }
}

#Goal, run subsample using do() for speed. Send a proportions vector of arbitrary length using a data.table containing vectors

factoryOut <- subsampleFactory( getSamplingParams= pairedTuplesDesign,
                                 mySubsample= lookSimpleSubsampler,
                                 handlers= list( lookHandler= lookHandler),
                                 qvalues= lookQval)
factoryOut( counts = lookCounts, design = lookDesign)

#Test that the ballanced proportions handler is working
factoryOut.twoTreatment <- subsampleFactory(getSamplingParams = twoTreatmentAllPermutations,
                                            mySubsample =       lookSimpleSubsampler,
                                            handlers =          list( lookHandler= lookHandler),
                                            qvalues =           lookQval)
factoryOut.twoTreatment( counts = lookCounts, design = lookDesign)

factoryOut.twoTreatment.edgeR <- subsampleFactory(getSamplingParams = twoTreatmentAllPermutations,
                                                  mySubsample =       lookSimpleSubsampler,
                                                  handlers =          list( edgeR= lookedgeR),
                                                  qvalues =           lookQval)
look <- factoryOut.twoTreatment.edgeR( counts = lookCounts, design = lookDesign)

# a function that generates all ballanced designs 
twoTreatmentAllPermutations <- function(design, counts){
  # Error checking------------------------------
  if (!is.data.frame( design))
    stop(" the design matrix should be a data frame")
  if (length( unique(design$treatment)) != 2)
    stop(" the design matrix should have only two treatments")
  
  # Permutations for each treatment ------------------------------
  #use indeces in the original count matrix as labels for the samples
  treatment.inds <- lapply( unique( design$treatment), function(myTreatment){
    which( design$treatment == myTreatment)
  })
  # for each treatment, get all permutations of the associated indeces
  all.combinations <- lapply( treatment.inds, function(inds){
    all.combinations <- expand.grid( lapply( inds, function(x) c(0,1)))
    num.samples <- rowSums( all.combinations)
    min.samples <- 2
    # ballanced design constrains maximum number of samples
    max.samples <- min( table( design$treatment))
    all.combinations[ num.samples >= min.samples & num.samples <= max.samples, ]
  })
  # Pairs of sample permutations from each treatment ------------------------------
  num.combinations <- sapply( all.combinations, function(combs) dim(combs)[1])
  a.rows <- rep(1:num.combinations[1], times = num.combinations[2])
  b.rows <- rep(1:num.combinations[2], each = num.combinations[1])
  # bind rows based on treatment
  out <- bind_cols( all.combinations[[1]][a.rows, ], all.combinations[[2]][b.rows, ])
  #permute columns to correspond to rows of DESIGN
  ii <- sort( unlist(treatment.inds), index.return= TRUE)$ix
  out[ , ii]
}

# a function that specifies subsampling at even proportions across samples
twoTreatmentAllPermutationsBallanced <- function( design, counts){
  unballanced <- twoTreatmentAllPermutations( design)
  total.counts <- colSums( counts)
  #ind.proportion * total.counts == min( total.counts)
  ind.proportions <- proportion * min( total.counts) / total.counts[ inds]
  design %*% diag( ind.proportions)
}

# write some simple subsampling functions as examples
lookTuplesDesign <- function( design, counts){
  treatments <- as.factor(design$treatment)
  #find the treatment with the minimum number of samples and use it to specify pairs
  treatment.list <- lapply( levels( treatments), function( myTreatment) which( treatments == myTreatment))
  min.length <- min( sapply( treatment.list, length))
  out <- apply( array(1:min.length, dim= c(min.length,1)), 1, function(ind){
    zeros <- rep( 0, times= length( treatments))
    one.inds <- sapply( treatment.list[1:2], function(list.el) list.el[ind])
    zeros[ one.inds] <- 1
    zeros
  })
  return( t( out))
}

lookHandler <- function( counts, design){
  treatment <- design
  sums <- mapply( unique( treatment), FUN= function(myTreat){
    treat.inds <- which( treatment == myTreat)
    rowSums( counts[, treat.inds, drop=FALSE])
  })
  data.frame( diff= abs( sums[, 1] - sums[, 2]))
}

lookedgeR <-
  function(counts, design, pair=NULL) {
    # filter counts below some threshold
    MIN.COUNT <- 20
    total.count <- rowSums( counts)
    pass.threshold <- which( total.count > MIN.COUNT)
    count.matrix <- counts[pass.threshold, ]
    id <- rownames( count.matrix)
    total.count <- total.count[ pass.threshold]
    # TODO(riley) check that the number of remaining rows is sufficient to continue...
    treatment = factor(design$treatment)
    if (is.null(pair)) {
      if (length(levels(treatment)) > 2) {
        stop(paste("If here are more than two levels of treatment,",
                   "must be explicitly given pair to compare"))
      }
      pair = levels(treatment)
    }
    
    d = edgeR::DGEList(counts=count.matrix, group=treatment)
    d = edgeR::calcNormFactors(d)
    d = edgeR::estimateCommonDisp(d)
    d = edgeR::estimateTagwiseDisp(d)
    edgeR.result = edgeR::exactTest(d, pair=pair)$table
    ret = data.frame(id= id,
                     coefficient=edgeR.result$logFC,
                     pvalue=edgeR.result$PValue,
                     count= total.count)
    ret
  }

lookSimpleSubsampler <- function( counts, props){
  pos.inds <- which( props != 0)
  pos.props <- props[ pos.inds]
  generateSubsampledMatrix(counts= counts,
                           indeces= pos.inds,
                           proportions= pos.props,
                           seed= 123,
                           replication=1) 
}

lookQval <- function( results){
  length.out <- dim( results)[1]
  data.frame( qval = runif( length.out))
}

lookCounts <- array( 1:16, dim=c(4,4))
lookDesign <- data.frame( treatment= c(1,1,2,2))

look <- pairedTuplesDesign( lookDesign, lookCounts)
look <- lookHandler( lookCounts, lookDesign)
look <- lookSimpleSubsampler( lookCounts, c(1,0,1,0))

# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.
qvalue.filtered1 = function(p) {
  # given a vector of p-values
  q.without1 = qvalue(p)$qvalue
  # when p == 1, the FDR is pi0
  #  q = rep(q.without1$pi0, length(p))
  # q[p < 1] = q.without1$qvalue
  q.without1
}

# when calculating pi0, it is prudent to filter out cases where p-values are exactly 1.
# for example, methods of differential expression tend to give p-values of 1 when
# the read count is 0. This would lead to over-estimating pi0 and thus the resulting
# q-values, even though those genes have no impact on the rest of the genes.

qvalue.fixedpi0 <- function (p, fdr.level = NULL, pfdr = FALSE, pi0s=NULL, ...)
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
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
  lfdr <- qvalue::lfdr(p = p, pi0 = pi0s, ...)
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
