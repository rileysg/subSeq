# Setup ------------------------------
# synthetic counts
means <- c( 1,4,1,4,4)
num.genes <- 10
lookCounts <- sapply( means, FUN=function( lambda){
  rpois( num.genes, lambda)})
lookDesign <- data.frame( treatment= c(1,2,1,2,2))
exptId <- 99999
proportions <- list( c(1),c(1),c(1),c(1),c(1))

example.subsample <- 

test_that("saving a subsampling results in the generation of a new file",{
  #specify the file to output, clear space to rewrite if necesary
  toSave
})