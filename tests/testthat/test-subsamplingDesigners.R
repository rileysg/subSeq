context("Design subsampling")

test_that("a simple unballanced design produces the correct combinations",{
  lookDesign <- data.frame( treatment= c(1,2,1,2,2))
  ss.design <- twoTreatmentAllPermutations( lookDesign)
  # all rows must sum to 4
  expect_equal( rowSums(ss.design), rep(4, 3))
  # collumns 2,4, and 5 must sum to 2
  #TODO(riley) this currently fails due to missing rownames
  expect_true( colSums(ss.design)[c(2,4,5)] == rep( 2,3))
})


#THE FOLLOWING IS STILL UNDER CONSTRUCTION ------------------------------
lookHandler <- function( counts, design){
  treatment <- design
  sums <- mapply( unique( treatment), FUN= function(myTreat){
    treat.inds <- which( treatment == myTreat)
    rowSums( counts[, treat.inds, drop=FALSE])
  })
  data.frame( diff= abs( sums[, 1] - sums[, 2]))
}

lookCounts <- sapply( means, FUN= function( lambda){
  rpois( num.genes, lambda)})
num.genes <- 10
means <- c( 1,4,1,4,4)

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