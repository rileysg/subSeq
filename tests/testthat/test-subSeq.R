#in addition to subseq, the following dependancies must be loaded:
# TODO:(riley) is there a way to require that these are loaded as part of the subsample function definition
# require(dplyr)
# require(tidyr)
# TODO:(riley) need to show that samples arn't being deranged durring the subsampling process

### custom expectations ###

expect_na = function(x){
  eval(bquote(expect_equal(sum( !is.na(.(x))), 0)))
}
expect_not_na = function(x){
  eval(bquote(expect_equal(sum( is.na(.(x))), 0)))
}

### setup ###

data(hammer)

hammer.counts = hammer@assayData$exprs[, 1:4]
hammer.design = hammer@phenoData@data[1:4, ]


# create a smaller set of counts to perform (faster) subsampling on
counts = hammer.counts[rowSums(hammer.counts) > 5000, ]
proportions = c(.01, .1, 1)
treatment = data.frame( treatment= hammer.design$protocol)

context("Individual handlers")

# for handlers, work with a small set that covers a range of sizes
smallcounts = hammer.counts[rowSums(hammer.counts) <= 100, ]
largecounts = hammer.counts[rowSums(hammer.counts) >= 100, ]

testcounts = rbind(smallcounts[sample(nrow(smallcounts), 20), ],
                   largecounts[sample(nrow(largecounts), 280), ])

test.output = function(output, numgenes=NULL) {
    # run tests on the output of a handler
    expect_that(output, is_a("data.frame"))
    expect_true("pvalue" %in% colnames(output))
    expect_true("coefficient" %in% colnames(output))
    expect_true(all(output$pvalue[!is.na(output$pvalue)] >= 0))
    expect_true(all(output$pvalue[!is.na(output$pvalue)] <= 1))

    if (!is.null(numgenes)) {
        expect_equal(nrow(output), numgenes)
    }
}

test_that("edgeR.glm handler works", {
   mm = model.matrix(~ treatment$treatment)
   test.output(subSeq::edgeR.glm(testcounts, mm=mm), nrow(testcounts))
   confounder = rnorm(ncol(testcounts))
   mm2 = model.matrix(~ treatment$treatment + confounder)
   test.output(subSeq::edgeR.glm(testcounts, mm=mm2), nrow(testcounts))

    # check that the edgeR applied to a random confounder is roughly uniform
    # (and thus that it actually is using that confounder)
   conf.out = subSeq::edgeR.glm(counts, mm=mm2, column=3)
   test.output(conf.out, nrow(counts))
   expect_true(mean(conf.out$pvalue) < .7 & mean(conf.out$pvalue) > .3)
})
test_that("edgeR handler works", {
   test.output(subSeq::edgeR(testcounts, treatment$treatment), nrow(testcounts))
})

test_that("DESeq2 handler works", {
   test.output(DESeq2(testcounts, treatment$treatment), nrow(testcounts))
})

test_that("voomLimma handler works", {
   test.output(voomLimma(testcounts, treatment$treatment), nrow(testcounts))
})

# DEXSeq has problems; leaving it out

# test_that("DEXSeq handler works", {
#     data( "pasillaExons", package="pasilla" )
#     n = 200
#     exon.counts = head(pasillaExons@assayData$counts, n)
#     design = pasillaExons@phenoData@data[c("condition", "type")]
#     geneIDs = head(pasillaExons@featureData@data$geneID, n)
#     exonIDs = head(pasillaExons@featureData@data$exonID, n)
# 
#     ret.exon = DEXSeq(exon.counts, design, geneIDs, exonIDs)
#     test.output(ret.exon, n)
# })
context("SLURM script")
library( "biorep")
library( "S4Vectors")
library( "IRanges")
library( "GenomicRanges")
library( "SummarizedExperiment" )

#specify parameters for subsampling
EXPERIMENT.ID <- 60590
NUM.PROC <- 1
my.id <- 1
OUT.DIR <- file.path(".","testCache")
ids <- c("SRR1555218", "SRR1555193", "SRR1555192", "SRR1555185", "SRR1555210", 
         "SRR1555204", "SRR1555199", "SRR1555206", "SRR1555202", "SRR1555194")
design <- GEOD60590Loader() %>% filter(file.id %in% ids)
load(paste0(DATA.DIR,"E-GEOD-", EXPERIMENT.ID, "-atlasExperimentSummary.Rdata"))
counts <- assays(experimentSummary$rnaseq)$counts
counts <- counts[,design$file.id]
proportions <- c(1)
subsampleDesigner <- function(design, counts){
  full.design <- ballancedPermutationXproportion(design, counts, proportions)
  designSplitter(id= my.id, full.design, NUM.PROC)
}
chunker <- fixedLengthChunker(length = 10)
# setting ORACLE to NULL should cause each subsample to be analysed individually
#TODO(riley) check this
writer <- summarizeAndSave(base.dir = file.path(".","testCache"),
                           expt.id = EXPERIMENT.ID, oracle = NULL)
getQvalues <- exptPi0Est

#perform subsampling
slurm.ss.out <- subsample(counts, design, method=c("edgeR", "voomLimma"), replications = 1,
          subsampleDesigner, chunker, writer, getQvalues = getQvalues)

break()

context("q-values")
#specify default arguments for calls to the subseq function
#most tests involve a small number of deviations from the default
dft                   <- list(counts = counts)
dft$design            <- treatment
dft$method            <- c("edgeR", "voomLimma")
dft$replications      <- 1
dft$subsampleDesigner <- function(design, counts)
                  permutationXproportion(design, counts, proportions = proportions)
dft$chunker           <- fixedLengthChunker(length = 1)
# use previous arguments to sp
dft.oracle <- getOracle(dft$method, dft$counts, dft$design)
dft$writer            <- summarizeAndSave(base.dir = file.path(".","testCache"),
                                  expt.id = "00001", oracle = dft.oracle)
dft$getQvalues        <- commonPi0Est(methods= c(edgeR= edgeR, voomLimma = voomLimma),
                                          counts= counts,
                                          design= treatment)
dft$seed <- 1234


test_that("an oracle can be generated using subsampling",{
  #when an oracle is not specified, the subsampling conditions with
  #highest proportion of reads should be used to specify an oracle for each differential expression testing method
  args.o                  <- dft
  args.o$method            <- c("edgeR")
  args.o$subsampleDesigner <- identityDesign
  args.o$writer            <- list(NULL)
  ss.o = do.call( subsample, args.o)
  #generate the oracel using getOracle directly
  edgeR.o <- getOracle(args.o$method, args.o$counts, args.o$design)
  expect_equivalent(select( ss.o, qvalue), select(edgeR.o, qvalue))
})

context("Subsampling")
ss = do.call( subsample, dft)
ss.summ = summary(ss)

#get fully crossed design with methods and replicatoins
ss.design <- expand.grid.df(dft$subsampleDesigner(dft$design, dft$counts),
               data.frame( method= dft$method),
               data.frame( replication= seq_len(dft$replication)))
proportion.inds <- !(colnames(ss.design) %in% c("method","replication"))
#load the summary and subsample objects,
load.params <- lapply(seq_len(dim(ss.design)[1]), function(ii)
  list( proportions = ss.design[ii, proportion.inds],method = ss.design[ii,"method"], replication = ss.design[ii,"replication"],
        expt.id = "00001", base.dir = file.path(".","testCache")))
ss.reload.list <- lapply( load.params, function(params) do.call( readSs, params))
ss.reload <- do.call( combineSubsamples, ss.reload.list)
ss.summary.reload.list <- lapply( load.params, function(params) do.call( readSsSummary, params))
#perhaps combineSubsamples is not intended to be called on "summary.subsamples" objects
ss.summary.reload <- do.call( combineSubsamples, ss.summary.reload.list)

#TODO(riley) summary.subsamples objects include a FDR.level atribute, should this be included? What does the vignette suggest?
test_that("qvalues of generated for an oracle match those generated for a subsample",{
  single.method <- "edgeR"
  ii <- which(ss.design$SRX020102 == 1 & ss.design$method == single.method)
  edgeR.ora <- dft.oracle %>% filter( method == single.method)
  reload <- ss.reload.list[[ii]]
  expect_equivalent(select(edgeR.ora, qvalue), select(reload, qvalue))
})

test_that("cached subsamples objects are identical after being reloaded", {
  expect_equivalent(ss.reload, ss)
  # reorder the reloaded summaries to match summary(ss)
  # use naming convention for subsampling design parameters
  match.colnames <- grep( "prop\\.|method|replication", colnames(summary(ss)), value = TRUE)
  ind.matches <- join.keys(ss.summary.reload, summary(ss), by= match.colnames)
  #find positions of summary(ss) rows in the reloaded table
  new.order <- match(ind.matches$y, ind.matches$x)
  expect_equivalent(ss.summary.reload[new.order, ], summary(ss))
})
#compare the result of subsampling to the unsubsampled object,
# the proportoins are all the same, so the results should match

test_that("subsamples returns a data table with the right structure", {
    expect_that(ss, is_a("subsamples"))
    expect_that(ss, is_a("data.table"))

    # check that the proportions and depths have the properties we expect
    all.proportions <- select( ss, contains("prop"))
    expect_equal(sort(unique(unlist(all.proportions))), proportions)
    expect_equal(max(ss$depth), sum(counts))
    expect_false(is.null(getSeed(ss)))

    # check the proportions are about what you expect (has some noise)
    proportion.changes = log2(ss$depth / max(ss$depth) / ss$proportion)
    expect_that(all(proportion.changes < .05), is_true())

    # check that the per-gene counts add up to the total depth
    countsums = ss[, list(total=sum(count)), by=c("depth", "method")]
    expect_equal(countsums$depth, countsums$total)

    # check that no replication was performed
    expect_true(all(ss$replication == 1))
})

ss.summ <- summary(ss)
test_that("quality metrics improve with increasing depth", {
    for (m in unique(ss.summ$method)) {
        sm = ss.summ[method == m]
        expect_that(all(diff(sm$pearson) > 0), is_true())
        expect_that(all(diff(sm$percent) > 0), is_true())
        expect_that(all(diff(sm$MSE) < 0), is_true())
    }
})

test_that("summaries can be created with other p-value corrections", {
    ss.summ.BH = summary(ss, p.adjust.method="BH")
    ss.summ.bon = summary(ss, p.adjust.method="bonferroni")
    expect_that(all(ss.summ.BH$significant < ss.summ$significant), is_true())
    expect_that(all(ss.summ.bon$significant < ss.summ.BH$significant), is_true())

    # check that you can't give it a nonsense method
    expect_that(summary(ss, p.adjust.method="nomethod"), throws_error("should be one of"))
})

#Test sampling with replication. Most of the default arguments are reused
rep <- dft
rep$replications <- 2
rep$subsampleDesigner = function(design, counts)
  permutationXproportion(design, counts, proportions = c(0.1))
rep$getQvalues
ss.rep <- do.call(subsample, rep)
test_that("Replications (multiple at each proportion) works", {
    sumrep = summary(ss.rep)
    # (distinct sampling proportions) x (methods) x (replications) = 1 * 2 * 2 = 4
    expect_equal(nrow(sumrep), 4)
    rep2 = sumrep[replication == 2]
    expect_equal(nrow(rep2), 2)
    expect_true(all( unlist( select( rep2, contains("prop."))) == .1))

    # TODO(riley) need to update this test to work with new method for labling proportion collumns
    # confirm there aren't replicates of 1.0
    #expect_false(any(ss.rep$proportion == 1 & ss.rep$replication == 2))

    # confirm that averaging works
    av.rep = summary(ss.rep, average=TRUE)
    expect_equal(nrow(av.rep), 2)
    expect_null(av.rep$replication)
})

test_that("subSeq can handle low counts", {
    low.counts = hammer.counts[rowSums(hammer.counts) < 2000, ]
    low.counts = low.counts[sample(nrow(low.counts), 500), ]

    low.proportions       <- c(.01, .1, 1)
    low                   <- dft
    low$counts            <- low.counts
    low$subsampleDesigner <- function(design, counts)
      permutationXproportion(design, counts, proportions = low.proportions)
    low.oracle            <- getOracle(low$method, low$counts, low$design)
    low$writer            <- summarizeAndSave(base.dir = file.path(".","testCache"),
                                              expt.id = "00001", oracle = low.oracle)
    ss.low = do.call(subsample, low)

    # test that plots still work
    expect_that(plot(summary(ss.low)), is_a("ggplot"))

    # significance might get wonky at this level, but correlations had better line up
    summ.low = summary(ss.low)
    for (m in unique(summ.low$method)) {
        sm.low = summ.low[method == m]
        expect_true(all(unlist(select(sm.low, contains("prop"))) %in% low.proportions))
        expect_true(all(diff(sm.low$pearson) > 0))
        expect_true(all(diff(sm.low$MSE) < 0))
    }

    # there should be no NAs or infinities
  #  expect_not_na(summ.low)
    expect_false(any(apply(summ.low, 2, is.nan)))
    expect_false(any(apply(summ.low, 2, is.infinite)))
})

test_that("Combining subsamples works", {
    seed = getSeed(ss)
    # try three other proportions
    more.proportions = c(.05, .3, .5)
    combine <- dft
    combine$subsampleDesigner <- function(design, counts)
      permutationXproportion(design, counts, proportions = more.proportions)
    combine$seed <- seed

    ss2 = do.call(subsample, combine)

    combined = combineSubsamples(ss, ss2)
    expect_equal(getSeed(ss), getSeed(ss2))
})


context("Custom Handlers")

test_that("Can provide custom error handlers", {
  fake.pvalues = runif(nrow(dft$counts))
  custom <- function(counts, treatment) {
    return(data.frame(pvalue=fake.pvalues, coefficient=-.2))
  }
  # note: substitute() does not always work properly when using do.call: see ?do.call
  # so construct the function call from elements of the default argument list
  ss.custom = subsample(method= custom,
                        writer= saveWithoutSummary(base.dir = file.path(".","testCache"),
                                                   expt.id = "00001"),
                        getQvalues=NULL,
                        counts= dft$counts,
                        design= dft$design,
                        replications = 1,
                        subsampleDesigner= dft$subsampleDesigner,
                        chunker= dft$chunker)

  test.output(custom(dft$counts, dft$treatment), length(fake.pvalues))
  expect_true(all(ss.custom$method == "custom"))
  expect_equal(fake.pvalues, ss.custom[depth == min(ss.custom$depth)]$pvalue)

  # check it can be given as a string as well
  ss.custom2 = subsample(method= "custom",
                         writer= saveWithoutSummary(base.dir = file.path(".","testCache"),
                                                    expt.id = "00001"),
                         getQvalues=NULL,
                         counts= dft$counts,
                         design= dft$design,
                         replications = 1,
                         subsampleDesigner= dft$subsampleDesigner,
                         chunker= dft$chunker)
  expect_true(all(ss.custom2$method == "custom"))
  expect_equal(fake.pvalues, ss.custom2[depth == min(ss.custom2$depth)]$pvalue)
})
test_that("Handlers can have columns that others don't", {
    fake.pvalues = runif(nrow(counts))
    othercols = replicate(3, runif(nrow(counts)))
    custom1 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 3], other=othercols[, 1]))
    }
    custom2 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 2], other=othercols[, 2]))
    }
    custom3 = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=othercols[, 1], other3=othercols[, 3]))
    }
    #TODO(riley) generating oracles from non-default methods currently is broken
    # the problem could be fixed by providing the correct ENV argument to the subsample function
    ss.custom = subsample(method= c("edgeR", "custom1", "custom2", "custom3"),
                          writer= saveWithoutSummary(base.dir = file.path(".","testCache"),
                                                     expt.id = "00001"),
                          getQvalues= commonPi0Est(methods= c("edgeR", "custom1", "custom2", "custom3"),
                                                   counts= dft$counts,
                                                   design= dft$treatment),
                          counts= dft$counts,
                          design= dft$design,
                          replications = 1,
                          subsampleDesigner= dft$subsampleDesigner,
                          chunker= dft$chunker)

    # we expect it to fill in missing values with NAs
    expect_true(all(c("other", "other3") %in% colnames(ss.custom)))
    expect_na(ss.custom[method == "edgeR", other])
    expect_na(ss.custom[method == "custom1", other3])
    expect_na(ss.custom[method == "custom2", other3])
    expect_na(ss.custom[method == "custom3", other])
    for (p in proportions) {
        expect_equal(ss.custom[method == "custom1" & proportion == p, other], othercols[, 1])
        expect_equal(ss.custom[method == "custom2" & proportion == p, other], othercols[, 2])
        expect_equal(ss.custom[method == "custom3" & proportion == p, other3], othercols[, 3])
    }
})
test_that("Handlers don't have to return one row per gene", {
    # some handlers, such as for gene set enrichment, don't necessarily return
    # one row per gene. Check it can return more and less
    for (n in c(nrow(counts) / 2, nrow(counts) * 2)) {
        fake.pvalues = runif(n)
        othercol = runif(n)
        coefs = rnorm(n)
        custom.different = function(counts, treatment) {
            return(data.frame(pvalue=fake.pvalues, coefficient=coefs, other=othercol,
                              ID=as.character(1:n)))
        }

        args.custom <- dft
        args.custom$method <- c("edgeR", "custom.different")
        ss.custom = do.call( subsample, args.custom)
        ss.edgeR = ss.custom[method == "edgeR"]
        ss.custom.different = ss.custom[method == "custom.different"]

        expect_equal(nrow(ss.edgeR), nrow(counts) * length(proportions))
        expect_equal(nrow(ss.custom.different), length(coefs) * length(proportions))

        expect_not_na(ss.edgeR$ID)
        expect_not_na(ss.custom.different$ID)
        expect_not_na(ss.custom.different$other)
        expect_na(ss.edgeR$other)
        expect_na(ss.custom.different$count)

        expect_not_na(ss.custom[method == "custom.different", other])

        # check that it doesn't affect the summary of the data
        custom.summ = summary(ss.custom)
        expect_not_na(custom.summ$pearson)
        expect_not_na(custom.summ$spearman)
    }

    # check that the handler does need to return an ID column in that case
    custom.noID = function(counts, treatment) {
        return(data.frame(pvalue=fake.pvalues, coefficient=coefs, other=othercol))
    }
    args.err <- dft
    args.err$method <- "custom.noID"
    expect_that(do.call(subsample,args.err),
                throws_error("if a handler doesn't return one row per gene then it must specify an ID column"))
})


context("Reproducibility from seeds")

s.edgeR = ss[method == "edgeR"]
s.voom = ss[method == "voomLimma"]

test_that("seeds are reproducible between methods", {
    summ.edgeR = ss.summ[method == "edgeR"]
    summ.voom = ss.summ[method == "voomLimma"]

    expect_equal(s.edgeR$depth, s.voom$depth)
    expect_equal(s.edgeR$count, s.voom$count)
})

test_that("seeds are reproducible between runs", {
    # perform it again with the same seed, and see that it matches the
    # first replication all the way through
    args.seed2 <- dft
    args.seed2$seed <- getSeed(ss)
    ss2 = do.call( subsample, args.seed2)

    expect_equal(ss$count, ss2$count)
    expect_equal(ss$pvalue, ss2$pvalue)
    expect_equal(ss$coefficient, ss2$coefficient)

    # for sanity check, check the same with the summary
    ss.summ2 = summary(ss2)
    expect_equal(ss.summ$significant, ss.summ2$significant)
    expect_equal(ss.summ$MSE, ss.summ2$MSE)
})

test_that("generateSubsampledMatrix retrieves the correct subsampled matrices", {
    p = sapply(select(ss, contains("prop.")), min)
    names(p) <- NULL
    
    seed = getSeed(ss)
    subm = generateSubsampledMatrix(counts, p, seed)
    expect_that(subm, is_a("matrix"))
    expect_equal(dim(subm), dim(counts))
    expect_equal(sum(subm), min(ss$depth))

    expect_equal(rowSums(subm), s.edgeR[prop.SRX020102 == unique(p)]$count, check.names = FALSE)
    expect_equal(rowSums(subm), s.edgeR[prop.SRX020102 == unique(p)]$count, check.names = FALSE)

    # confirm that edgeR on that matrix gives the same results
    subm.results = edgeR(subm, treatment=unlist(dft$design))
    expect_equal(subm.results$pvalue, s.edgeR[prop.SRX020102 == unique(p)]$pvalue)
    expect_equal(subm.results$coefficient, s.edgeR[prop.SRX020102 == unique(p)]$coefficient)

    # confirm summary object also contains correct seed
    summ = summary(ss)
    expect_equal(subm, generateSubsampledMatrix(counts, p, getSeed(summ)))

    # confirm generateSubsampledMatrix works if you explicitly
    # tell it replication is 1
    subm.rep1 = generateSubsampledMatrix(counts, p, seed, replication=1)
    expect_equal(subm, subm.rep1)
})

test_that("Performing multiple replicates is reproducible", {
    args.rep2 <- rep
    args.rep2$seed <- getSeed(ss.rep)
    ss.rep.2 = do.call(subsample, args.rep2)
    expect_equal(ss.rep$pvalue, ss.rep.2$pvalue)
})

context("Plotting")

test_that("plotting is possible without errors", {
    expect_that(plot(ss.summ), is_a("ggplot"))
})


context("Error handling")

test_that("Raises an error on edgeR if there are >2 treatments", {
    new.treatment = c("A", "A", "B", "C")
    # check with multiple handlers
    expect_that(subsample(counts, proportions, treatment=new.treatment,
                          method="edgeR"), throws_error("more than two levels"))
})

test_that("Raises an error if it cannot find the handler", {
    expect_that(subsample(counts, proportions, treatment=treatment,
                          method="nomethod"),
                throws_error("Could not find handler nomethod"))
})

test_that("error messages are thrown when proportions are incorrect", {
    expect_that(subsample(counts, c(), treatment=treatment, method="edgeR"),
                throws_error("No proportions"))
    expect_that(subsample(counts, c(.1, 1, 2), treatment=treatment, method="edgeR"),
                throws_error("Proportions must be in range"))
    expect_that(subsample(counts, c(0, 1), treatment=treatment, method="edgeR"),
                throws_error("Proportions must be in range"))
})

test_that("error message is thrown if counts were normalized", {
    # confirming that there was no normalization
    expect_that(subsample(scale(counts), proportions, treatment=treatment,
                          method="edgeR"), throws_error("unnormalized"))
    sc.counts = scale(counts, center=FALSE)
    expect_that(subsample(sc.counts, proportions, treatment=treatment,
                          method="edgeR"), throws_error("unnormalized"))
})

test_that("combineSubsamples raises an error when combining different seeds", {
    more.proportions = c(.05, .3, .5)
    ss2.arg <- dft
    ss2.arg$subsampleDesigner <- function(design, counts)
      permutationXproportion(design, counts, proportions = more.proportions)
    ss2 = do.call(subsample, ss2.arg)

    expect_false(getSeed(ss) == getSeed(ss2))
    expect_that(combineSubsamples(ss, ss2), throws_error("different.*seed"))
})

# handlers:
# -write baySeq, DEXSeq

# plots:
# -volcano plots
# -per-gene plots
