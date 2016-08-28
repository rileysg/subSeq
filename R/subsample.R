
#' Subsample reads and perform statistical testing on each sample
#'
#' Perform subsampling at multiple proportions on a matrix of count
#' data representing mapped reads across multiple samples in many
#' genes. For each sample, perform some statistical operations.
#'
#' @param counts Matrix of unnormalized counts
#' @param proportions Vector of subsampling proportions in (0, 1]
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
#' ss = subsample(hammer.counts, c(.01, .1, 1), treatment=hammer.design$protocol,
#'                  method=c("edgeR", "DESeq2", "voomLimma"))
#'
#' @import data.table
#' @importFrom dplyr group_by do mutate filter bind_rows group_by_ mutate_
#' @importFrom tidyr gather
#' @import magrittr
#' @importFrom qvalue qvalue
#' @importFrom Biobase exprs pData
#' @export
subsample <-
    function(counts, design, method="edgeR", replications = 1,
             subsampleDesigner, chunker, writer,
             seed=NULL, getQvalues = NULL, env=parent.frame(), ...) {
        # a switch that specifies whether subsample should return anything or rely on writer to deal with output
        COLLECT.SUBSAMPLES <- TRUE
        # error checking------------------------------
        # check that the counts is an unnormalized numeric matrix
        counts = as.matrix(counts)
        error = max(abs(round(counts) - counts))
        if (error > 1e-4 | any(counts < 0)) {
            stop("Counts should be unnormalized integer counts (not, e.g., RPKM)")
        }
        
        #ensure that every column of COUNTS has a unique label
        if (is.null( colnames(counts)))
          colnames(counts) <- as.character( 1:dim(counts)[2])
        rep.colnames <- which(table(colnames(counts)) > 1)
        if (max(rep.colnames) > 1)
          stop("Each count column (sample) must have unique name")
        
        if (is.null(seed)) {
            # come up with a random initial seed (can't use current one, since
            # the point is to save it for later)
            # borrowing the method from gabor csardi
            # here: https://github.com/airoldilab/stat221/blob/6e98d122f75c1fe38903380e7912816398ab3ec8/R/utils.r
            # Generate a random seed based on microseconds
            now <- as.vector(as.POSIXct(Sys.time())) / 1000
            seed <- as.integer(abs(now - trunc(now)) * 10^8)
        }

        # create a list mapping function name to function
        if (is.function(method)) {
            # if given a single function, make that into a list
            methods = list(method)
            names(methods) = deparse(substitute(method))
        }
        else {
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
                else if (m %in% c("edgeR", "voomLimma", "DESeq2", "edgeR.glm")) {
                    handler = get(m, mode="function")
                }
                else {
                    stop(paste("Could not find handler", m))
                }
                handler
            })
            if ( !is.null(names(method)))
              names(methods) = names(method)
            else
              names(methods) = method
        }
        
        if (!(is.function(getQvalues) || is.null(getQvalues)))
          stop("the subsample argument getQvalues must be either a function or NULL")
        
        # specify how subsampling should be done ------------------------------
        # get unique combinations of subsampling
        props = subsampleDesigner( design, counts)
        colnames( props) <- paste0("prop.", colnames(props))
        #TODO(riley): check that all subsampling proportions are unique

        # cross subsampling proportions with method and replication
        params = expand.grid.df( props,
                                 data.frame( method = names(methods)),
                                 data.frame( replication = seq_len(replications)))

        #define how subsampling will be done ------------------------------
        perform.subsampling <- function(method, proportion, replication) {
            # generate subcounts, dropping collumns and rows that are subsampled at proportion 0
            subcounts = generateSubsampledMatrix(counts, proportion, seed, replication)
            id = which(rowSums(subcounts) >= 5)
            subcounts = subcounts[id,] ### Filter counts at zero
            if (length(id) == 0) return("Error: counts too low at subsampling proportion")
            handler = methods[[method]]
            #get treatments associated with remaining samples
            ind.original <- which(colnames(counts) %in% colnames(subcounts))
            subtreat <- design$treatment[ind.original]
            #test for differential expression
            ret = handler(subcounts, treatment= subtreat, ...)
            
            # add gene names (ID) and per-gene counts
            # If there is one handler row per gene, the remaining columns of ret, "ID" and "counts", can be infered.
            infer.per.gene = dim( ret)[1] == dim( subcounts)[1]
            if ( !any( "ID" == colnames(ret))){
              if (infer.per.gene){
                ret$ID = rownames(subcounts)
              } else {
                stop("if a handler doesn't return one row per gene then it must specify an ID column")
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
        
        #make a subsample object from a dataframe
        frameToSubsample <- function(frame){
          frame = as.data.table(as.data.frame(frame))
          class(frame) = c("subsamples", class(frame))
          attr(frame, "seed") = seed
          frame
        }
        
        #perform analysis tasks in batch using CHUNKER ------------------------------
        # to do subsampling and writing to the disk in batch, use a list of vectors of indeces
        chunks <- chunker( params)
        ss.out <- lapply( chunks, function(myChunk){
          chunk.sss <- lapply(myChunk, function(ii){
            myParams <- params[ii,]
            ss.props <- as.numeric( myParams[as.character(colnames(props))])
            res <- perform.subsampling( myParams$method,
                                        ss.props,
                                        myParams$replication)
            res <- cbind( myParams, res, row.names= NULL)
            #if a method for getting qvals has been suppiled, use it
            if (is.function(getQvalues)) {
              qvals <- getQvalues(res$pvalue, res$method)
              res <- bind_cols( res, qvals)
            }
            res
            })

          # If a writer is specified, cache subsample objects individually
          if (is.function(writer)){
            lapply(chunk.sss, function(mySS) writer( frameToSubsample(mySS)))
          }
            
          out <- NULL
          if (COLLECT.SUBSAMPLES)
            out <- bind_rows(chunk.sss)
          out
        })
        ss.out <- frameToSubsample( bind_rows(ss.out))
        return(ss.out)
    }
