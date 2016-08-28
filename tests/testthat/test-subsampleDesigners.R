#generate a synthetic design matrix
synthMat <- function(n) matrix( seq_len(n^2), nrow=n)
synthFrame <- function(n) as.data.frame( matrix( seq_len(n^2), nrow=n),)
synthTreat <- function(n) {
  split <- floor(n/2)
  treatments <- c(rep("a", split), rep("b", n - split))
  data.frame( treatment = treatments)
}
# generate a synthectic count matrix
synthCounts <- function(m, n, seed = NULL){
  TRIALS <- 100
  PROB <- 0.1
  if ( !is.null(seed))
    seed <- set.seed(seed)
  randCounts <- rbinom(m*n, TRIALS, PROB)
  dim.names <- list(c(),letters[1:n])
  data.frame( matrix(randCounts, ncol = n, dimnames = dim.names))
}

# combinatorial constraints on designs --------------------------------------------
# calculate the number of rows with a particular sum
# n.a and n.b represent the number of samples with treatments a, and b respectively
rowCount <- function(n.a, n.b, k) {
  # ballanced design => each row has an even sum
  if (round(k/2) - k/2 != 0) {
    warning("all rows must have an even number of elements")
    out <- 0
  }
  # check constrainsts on number of elements to be chosen
  else if ( k/2 > max(n.a, n.b) | k < 4){
    warning("cannot choose more elements than exist for a given treatment")
    out <- 0
  } else {
    out <- choose(n.a, k/2) * choose(n.b,k/2)
  }
  out
}

context("generating permutations for subsampling")

test_that("permute.indicators can be called",{
  simple.fr <- synthMat(4)
  simple.tr <- synthTreat(4)
  permute.indicators(design = simple.tr, counts = simple.fr)
})

test_that("permute.indicators catches invalid designs",{
  simple.fr <- synthMat(4)
  inval.tr <- data.frame( treatments = rep("a", times = 4))
  expect_error( permute.indicators(inval.tr,simple.fr))
  
  small.frame <- synthMat(3)
  small.treat <- data.frame( treatments = c("sick","sick","well"))
  expect_error( permute.indicators(small.treat, small.frame))
})

test_that("subsamplings from a ballanced design have correct dimensions",{
  iis <- seq(4,20,by=2)
  #assuming both treatments have n samples, the number of rows is given by
  count.rows <- function(n) choose(2*n,n) - n^2 - 1
  total.rows <- sapply(iis, function(ii) count.rows(ii/2))
  designs <- lapply(iis, function(i) permute.indicators(design = synthTreat(i), counts = synthFrame(i)))
  
  # elements should be indicators
  expect_true(all(unique(unlist(designs)) %in% c(0,1)))
  
  # check dimensions
  expect_equal(sapply(designs, nrow), total.rows)
  expect_equal(sapply(designs, ncol), iis)

  # check column sums
  # after fixing a particular element to be 1, there remain n-1 choose k-1 other potential arangements of the two
  # treatments
  # The total number is sum(i=2:n, choose(n-1,i-1) * choose(n,i))
  # == choose(2*n - 1, n) - choose(n,1)
  col.constrained.permutations <- function(n) choose(2*n - 1, n) - n
  cs.expected <- sapply(iis, function(ii) rep(col.constrained.permutations(ii/2), times=ii))
  cs.observed <- sapply(designs, function(dsg) setNames(colSums(dsg), c()))
  expect_equal( cs.observed, cs.expected)
  
  # check row sums
  marginal.counts.observed <- sapply(designs, function(dsg) as.numeric(table(rowSums(dsg))))
  marginal.counts.expected <- lapply(iis, function(total) {
    sapply(seq(4,total, by = 2), function(k) rowCount(total/2, total/2, k))
  })
  expect_equal( marginal.counts.observed, marginal.counts.expected)
})

test_that("margins are correct for randomly selected sets of treatments",{
  #generate randomly permuted sample lables subject to the constraints
  #at least two of each type of sample
  randomTreat <- function(n, prob.a = 0.2){
    if (n < 4)
      stop("there must be at least four samples")
    required.treatments <- c("a","a","b","b")
    random.treatmets <- sample(c("a","b"), size= n - 4, replace = TRUE, prob = c(prob.a, 1-prob.a))
    tr <- c( required.treatments, random.treatmets)
    # derange required and random treatments
    tr <- sample(tr, n, replace= FALSE)
    data.frame(treatment = tr)
  }
  
  # generate some synthetic samples with different numbers of columns and random treatments
  iis <- seq(4,20)
  mats <- lapply(iis, synthMat)
  seed <- .Random.seed
  trts <- lapply(iis, randomTreat)
  dsgs <- lapply(seq_len(length(mats)), function(i) permute.indicators(design = trts[[i]], counts = mats[[i]]))
  n.as <- sapply(trts, function(chars) as.numeric( table(chars)["a"]))
  n.bs <- iis - n.as
  
  #TODO(riley) it would be nice to check this using a simpler expression.
  count.rows.general <- function(n.a, n.b){
    k.max <- max(n.a,n.b)
    ks <- 2:k.max
    sum( sapply(ks, function(k) rowCount(n.a,n.b,2*k)))
  }
  
  # check dimensions
  expect_equal(sapply(dsgs, ncol), iis)
  expect.rows <- sapply(seq_along(n.as), function(ii) count.rows.general(n.as[ii], n.bs[ii]))
  expect_equal(sapply(dsgs, nrow), expect.rows)
  
  #check marginal sums
  rowsum.observed <- sapply(dsgs, function(dsg) as.numeric(table(rowSums(dsg))))
  rowsum.expected <- sapply(seq_along(n.as), function(ii) {
    n.max <- min( n.as[ii], n.bs[ii])
    total.samples <- seq(4, 2 * n.max, by = 2)
    sapply(total.samples, function(k) {
      out <- rowCount(n.as[ii], n.bs[ii], k)
    })
  })
  expect_equal(rowsum.observed, rowsum.expected)
})

context("scaling proportions for subsamplings")

test_that("BALLANCE behaves as expected for a ballanced design",{
  seed <- .Random.seed
  n.col <- 5
  n.row <- 3
  seed <- 1234
  cts <- Reduce(cbind, Map( function(...) synthCounts(n.row,1,seed), 1:n.col))
  colnames(cts) <- letters[1:n.col]
  
  dsg <- data.frame(matrix(rep(1,times=3*n.col), ncol= n.col, dimnames = dimnames(cts)))
  dsg.bal <- ballance(dsg, cts)
  expect_equivalent(dsg, dsg.bal)
})

test_that("BALLANCE behaves as expected for an unballnced design",{
  seed <- .Random.seed
  n.col <- 5
  n.row <- 3
  seed <- 1234
  cts <- Reduce(cbind, Map( function(...) synthCounts(n.row,1,seed), 1:n.col))
  two.scale <- sapply(0:4, function(n) 2^n)
  scale.counts <- as.matrix(cts) %*% diag(two.scale)
  scale.counts <- as.data.frame(scale.counts)
  colnames(scale.counts) <- colnames(cts)
  
  dsg <- data.frame(matrix(rep(1,times=3*n.col), ncol= n.col, dimnames = dimnames(cts)))
  
  expect.scaled <- as.matrix(dsg) %*% diag(1/two.scale)
  expect.scaled <- as.data.frame(expect.scaled)
  colnames(expect.scaled) <- colnames(dsg)
  
  scale.dsg <- ballance(dsg, scale.counts)
  expect_equivalent(scale.dsg, expect.scaled)
})
  #scale counts and expect proportions to be scaled by the inverse


test_that("expand.scale works correctly for a simple example",{
  dsg <- data.frame(a=1,b=1,c=1)
  scales <- seq(0,1,0.1)
  scaled <- expand.scale(dsg, scales)
  expected.scale <- scales %*% as.matrix(dsg)
  expected.scale <- as.data.frame(expected.scale)
  colnames(expected.scale) <- colnames(dsg)
  expect_equivalent(scaled, expected.scale)
})

matchRows <- function(vec, frame){
  compare <- apply(frame, 1, function(row) all(row == vec))
  matches <- which(compare)
  matches
}

list.compare <- function(a, b){
  sapply( seq_along(a), function(ii) a[ii] == b[ii])
}

test_that("designSplitter correctly perofrms a partition for some small designs",{
  #generate a matrix of unique rows, test many possible numbers of rows
  sapply(1:30, function(n){
    mat <- synthMat(n)
    # partition the matrix in all possible ways
    observed.parts <- lapply(1:n, function(n.parts){
      # get each of the partitions and see which rows of the orignal matrix they contain
      part.rows <- sapply(1:n.parts, function(id){
        part <- designSplitter(id, mat, n.parts)
        i.rows <- seq_len( dim(part)[1])
        sapply( i.rows, function(i) matchRows(part[i,], mat))
      })
      #count the number of times each row appears
      inds.counts <- table(unlist(part.rows), deparse.level = 0)
      expected.parts <- table(1:n)
      expect_equal(inds.counts, expected.parts)
    })
  })
})