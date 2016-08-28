# a chunker that either includes a fixed number of indeces in its chunk, or all remaining indeces
# returns a list with indeces corresponding to chunks
fixedLengthChunker <- function(length){
  function(design){
    total <- dim(design)[1]
    num.chunks <- ceiling(total/length)
    start.inds <- seq.int(from=1, by=length, length.out= num.chunks)
    end.inds <- c(tail(start.inds, -1) - 1, total)
    lapply(seq_along(start.inds), function(ii) start.inds[ii]:end.inds[ii])
  }}

recursiveChunker <- function(length, n){
  if (n != round(n) || n < 0)
    stop("number of partitions must be a natural number")
  if( length != round(length) || length < 0)
    stop("can only partition a natural number of elements")
  if (n > length)
    stop("cannot partition n integers into more than n subsets")
  if (length == 0){
    #final.list <- lapply(1:n, function(start) (1 + r * (start -1)):(start * (r - 1)))
    return(NULL)
  } else {
    r <- length / n
    step <- ceiling(r)
    next.len <- length - step
    tail.inds <- list( c((next.len + 1):length))
    head.inds <- recursiveChunker(length = next.len, n = n - 1)
    return( c(head.inds, tail.inds))
  }
}
