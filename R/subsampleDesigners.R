#subsample designers are used to specify subsampling proportions and should have the following behavior
# A collumn for each sample in the experiment representing the probability of resampling a given read
# subsamplers should NOT specify replications

# to simplify the interfaces between functions, subsamplings are specified using data.frames

#------------------------------
#Helper functions
#------------------------------

# a function that generates all ballanced designs 
# a more readable implementation could recursively generate "grey codes" then use an
# expand.grid.df function to generate all pairs of sample subsets
#DESIGN refers to experimental design and is used to identify treatments
permute.indicators <- function(design, counts){
  # Error checking------------------------------
  if (!is.data.frame( design))
    stop(" the design matrix should be a data frame")
  if (length( unique(design$treatment)) != 2)
    stop(" the design matrix should have only two treatments")
  if (min(table(design$treatment)) < 2)
    stop("each treatment must have at least two associated samples")
  
  # Permutations for each treatment ------------------------------
  #get indeces associated with each treatment
  treatment.inds <- lapply( unique( design$treatment), function(myTreatment){
    which( design$treatment == myTreatment)
  })

  # ballanced design constrains maximum number of samples
  # also, differential expression methods generally require >= 2 samples
  max.samples <- min( table( design$treatment))
  min.samples <- 2
  # for each treatment, get all permutations of the associated indeces
  # store them based on number of samples included
  all.combinations <- lapply( treatment.inds, function(inds){
    all.combinations <- expand.grid( lapply( inds, function(x) c(0,1)))
    num.samples <- rowSums( all.combinations)
    all.combinations[ num.samples >= min.samples & num.samples <= max.samples, ]
    out <- lapply( 2:max.samples, function(n)
      all.combinations[which(num.samples == n), ])
    out
  })

  # Pairs of sample permutations from each treatment ------------------------------
  # each must have the same total number of samples
  out <- lapply(seq_len(max.samples - 1), function(ii){
    bind.cols <- expand.grid.df(all.combinations[[1]][[ii]], all.combinations[[2]][[ii]])
  })
  out <- do.call(rbind, out)
  #combinations taken for each treatment seprately,
  #want to map back to original order
  new.order <- unlist(treatment.inds)
  restore.order <- sort.int(new.order, index.return = TRUE)$ix
  out <- out[ , restore.order]
  colnames( out) <- colnames(counts)
  out
}

#this function is taken from stack exchange:
# a version of expand.grid for data frames
# http://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

# each row of FRAME is scaled by each of the elements in scales
# if FRAME is (m x n) and SCALES is a vector of length p then the output is dimension (pm x n) 
expand.scale <- function(dsg, scales){
  out <- lapply( scales, function(x) x*dsg)
  bind_rows(out)
}

# a function that takes an unballanced subsampling matrix (DSG) and
# modifies its proportions such that each sample will have the same
# expected number of reads
ballance <- function(dsg, counts){
  col.counts <- colSums( counts)
  #note that ind.proportion * total.counts == min( total.counts)
  normalize.col <- min( col.counts) / col.counts 
  out <- as.matrix(dsg) %*% diag( normalize.col)
  out <- as.data.frame(out)
  colnames(out) <- colnames(dsg)
  out
}

# a function that partitions a design into NUM.SETS distinct subsets,
# each of which is associated to an ID
# The idea is to use this function to distribute subsampling among several nodes
designSplitter <- function(id, design, n.sets){
  #Error checking -----------------------------------------
  if (round(n.sets) != n.sets || n.sets < 1)
    stop("n.sets must be a natural number")
  if (!(id %in% 1:n.sets))
    stop("id must be an integer")
  #Partition the rows of design, then return those associated with ID
  n.row <- dim(design)[1]
  chunks <- recursiveChunker(length = n.row, n.sets)
  design[ chunks[[id]], ,drop = FALSE]
}

#------------------------------
#Designers
#------------------------------

# a function that does subsampling for a design in which 
# all permutations of samples are fullly crossed with the specified proportions
permutationXproportion <- function(design, counts, proportions){
  # get all ballanced subsampling permutations, each row contains indicators of whether a sample is included
  out <- permute.indicators(design, counts)
  # sample all permutations at each subsampling proportion
  out <- expand.scale(out, proportions)
  out
}

# a function that subsamples equal numbers of expected reads
# for all combinations of samples crossed with subsampling within a given proportion
ballancedPermutationXproportion <- function(design, counts, proportions, max.decimals= 9){
  # get all ballanced subsampling permutations, each row contains indicators of whether a sample is included
  out <- permute.indicators(design, counts)
  # sample all permutations at each subsampling proportion
  out <- expand.scale(out, proportions)
  # each included sample should have same expected number of reads
  out <- ballance(out, counts)
  # for convenience in naming files, can limit digits in a proportion
  round( out, max.decimals)
}

#a subsampler that uses all of the reads of the original experiment
identityDesign <- function(design, counts){
  # rows of the design matrix represent samples
  n.cols <- dim(design)[1]
  ret <- rep_len(1, n.cols)
  ret <- data.frame( matrix( ret, nrow = 1))
  colnames( ret) <- colnames( counts)
  ret
}