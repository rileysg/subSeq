#' combine multiple subsamples objects
#' 
#' Given two or more subsamples objects, combine them into one larger object, on
#' which we can perform all the usual analyses and plots.
#' 
#' @param ... Two or more subsamples objects
#' 
#' @details If there are columns in some subsamples objects that are not in others,
#' the missing values will be filled with NA
#' @examples 
#' # see ?subsample to see how ss is generated
#' data(ss)
#' 
#' # combine multiple subsampling objects (in this example they happen to be the same object)
#' ss_new <- combineSubsamples(ss, ss)
#' 
#' @return subSeq object
#' 
#' @export
combineSubsamples <-
    function(...) {
        lst = list(...)

        # confirm that the global seed is the same
        seeds = sapply(lst, function(s) attr(s, "seed"))
        if (length(unique(seeds)) != 1) {
            stop(paste("Cannot combine subsamples objects that have",
                 "different random seeds: methods are not comparable"))
        }
        
        #confirm that all elements have the same class
        classs <- sapply(lst, function(s) class(s))
        if (dim(unique(classs, MARGIN=2))[2] != 1){
          stop(paste("Cannot combine subsamples objects that have",
                     "different classes"))
        }

        all.colnames = Reduce(union, lapply(lst, colnames))
        with.NA = lapply(lst, function(s) {
            ifnot = lapply(all.colnames, function(i) rep(NA, nrow(s)))
            as.data.frame(mget(all.colnames, envir=as.environment(s),
                                ifnotfound=ifnot), stringsAsFactors = FALSE)
        })

        ret = as.data.table(do.call(rbind, with.NA))
        # save the class
        class(ret) = class(lst[[1]])
        # save the seed
        attr(ret, "seed") = attr(lst[[1]], "seed")
        ret
    }
