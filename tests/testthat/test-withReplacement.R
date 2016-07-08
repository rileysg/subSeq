
reproduceError <- function(){

load("/media/share/GradSchool/Storey/riley/storey/bioRepSeq/software/biorep/inst/extdata/E-GEOD-52202-atlasExperimentSummary.Rdata" )
  rSumExp <- experimentSummary$rnaseq
  counts <- assays( rSumExp)$counts
  als.design <- structure(list(AtlasAssayGroup = c("g2", "g2", "g2", "g2", "g2", 
                                   "g2", "g2", "g2", "g1", "g1", "g1", "g1", "g1", "g1", "g1", "g1"
), progenitor_cell_type = c("induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell", "induced pluripotent stem cell", 
                            "induced pluripotent stem cell"), cell_type = c("motor neuron", 
                                                                            "motor neuron", "motor neuron", "motor neuron", "motor neuron", 
                                                                            "motor neuron", "motor neuron", "motor neuron", "motor neuron", 
                                                                            "motor neuron", "motor neuron", "motor neuron", "motor neuron", 
                                                                            "motor neuron", "motor neuron", "motor neuron"), disease = c("normal", 
                                                                                                                                         "normal", "normal", "normal", "normal", "normal", "normal", "normal", 
                                                                                                                                         "ALS", "ALS", "ALS", "ALS", "ALS", "ALS", "ALS", "ALS"), organism = c("Human", 
                                                                                                                                                                                                               "Human", "Human", "Human", "Human", "Human", "Human", "Human", 
                                                                                                                                                                                                               "Human", "Human", "Human", "Human", "Human", "Human", "Human", 
                                                                                                                                                                                                               "Human"), file.id = c("SRR1027591", "SRR1027592", "SRR1027593", 
                                                                                                                                                                                                                                     "SRR1027594", "SRR1027595", "SRR1027596", "SRR1027597", "SRR1027598", 
                                                                                                                                                                                                                                     "SRR1027599", "SRR1027600", "SRR1027601", "SRR1027602", "SRR1027603", 
                                                                                                                                                                                                                                     "SRR1027604", "SRR1027605", "SRR1027606"), biological.replicate = structure(c(5L, 
                                                                                                                                                                                                                                                                                                                   5L, 6L, 6L, 7L, 7L, 8L, 8L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L), .Label = c("C9 28i", 
                                                                                                                                                                                                                                                                                                                                                                                           "C9 29i", "C9 30i", "C9 52i", "Control 00i", "Control 03i", "Control 14i", 
                                                                                                                                                                                                                                                                                                                                                                                           "Control 83i"), class = "factor"), lib.size = c(13432337, 13420175, 
                                                                                                                                                                                                                                                                                                                                                                                                                                           12851479, 12841250, 12037102, 12010712, 12692826, 12676512, 9971404, 
                                                                                                                                                                                                                                                                                                                                                                                                                                           9936337, 11192983, 11169677, 15227880, 15162335, 14508124, 14467152
                                                                                                                                                                                                                                                                                                                                                                                           ), treatment = c("normal", "normal", "normal", "normal", "normal", 
                                                                                                                                                                                                                                                                                                                                                                                                            "normal", "normal", "normal", "ALS", "ALS", "ALS", "ALS", "ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                            "ALS", "ALS", "ALS"), design.note = c("NA", "NA", "NA", "NA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  "NA"), question = c("a mutation that causes ALS", "a mutation that causes ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "a mutation that causes ALS", "a mutation that causes ALS", "a mutation that causes ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "a mutation that causes ALS", "a mutation that causes ALS", "a mutation that causes ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "a mutation that causes ALS", "a mutation that causes ALS", "a mutation that causes ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "a mutation that causes ALS", "a mutation that causes ALS", "a mutation that causes ALS", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "a mutation that causes ALS", "a mutation that causes ALS"), 
sample = c("cell culture", "cell culture", "cell culture", 
           "cell culture", "cell culture", "cell culture", "cell culture", 
           "cell culture", "cell culture", "cell culture", "cell culture", 
           "cell culture", "cell culture", "cell culture", "cell culture", 
           "cell culture")), class = "data.frame", .Names = c("AtlasAssayGroup", 
                                                              "progenitor_cell_type", "cell_type", "disease", "organism", "file.id", 
                                                              "biological.replicate", "lib.size", "treatment", "design.note", 
                                                              "question", "sample"), row.names = c(NA, -16L))



get.pooled.counts <- function( design, counts){
  # find collumns of a matrix associated with a given biological replicate
  getInds <- function( biological.replicate){
    design.inds <- which( design$biological.replicate == biological.replicate)
    design.names <- design$file.id[ design.inds]
    which( dimnames( counts)[[2]] %in% design.names)
  }
  
  add <- function( x) Reduce("+", x, accumulate = FALSE)
  
  #add collums corresponding to a particular biological replicate
  pool <- function( brep){
    brep.inds <- getInds( brep)
    cols <- Map( function(col) counts[, col],
                 brep.inds)
    add( cols)
  }
  #get collumns corresponding to each biological replicate
  out.list <- Map( function(myBrep) pool( myBrep),
                   names( get.pooled.treats( design)))
  data.frame( out.list)
}
get.breps <- function( design){
  as.character( unique( design$biological.replicate))
}

get.pooled.treats <- function( design){
  #get the names of the biological replicates using the same mehtod as get.pooled.counts
  breps <- get.breps( design)
  
  #find a treatment associated with a particular biological replicate
  getTreatment <- function( brep){
    treat.ind <- which( brep == as.character( design$biological.replicate))
    design$treatment[ treat.ind]
  }
  
  all.treats <- Map( getTreatment, breps)
  sapply( all.treats, unique)
}

pooled.count <- get.pooled.counts( als.design, counts)
pooled.treat <- get.pooled.treats( als.design)
bio.reps <- 2:4
proportions <- bio.reps/4
replications <- 1

ss <- subsample( pooled.count, pooled.treat, proportions = proportions, bioReplicates = bio.reps,
                 method = c("edgeRDispersions", "voomLimma"), replications = replications,
                 replacement = TRUE, ballancedproportions = TRUE,
                 use.common.pi0 = TRUE)


}