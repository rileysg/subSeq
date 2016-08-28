#used to specify a common naming convention for saving and loading files
# PROPORTIONS should be specified as a list with elements corresponding to subsamplings that occured
# for a particular sample
masterName <- function(expt.id, proportions, method, replication, obj.type) {
  #specify a threshold for the maximum number of decimal places to be used in names
  MAX.DECIMALS <- 9
  match.ind <- which(obj.type %in% c( "summary.subsamples", "subsamples"))
  if (1 != length(match.ind))
    stop( "naming conventions are defined only for objects of type subsamples or summary.subsamples")
  # check that the input proportions are specified to three decimal places
  if (any( unlist( proportions) - round( unlist( proportions), MAX.DECIMALS) != 0))
    stop( paste0("proportions cannot be specified using more than ", MAX.DECIMALS),
          "fractional digits")
  #convert subsampling proportions from ragged array format to string format
  # separate proportions
  proportion.strings <- sapply( proportions, function(vec) paste( vec, collapse= "-"))
  # separate samples
  proportions.long <- paste( proportion.strings, collapse= "__")
  paste0("E-GEOD-", expt.id, "-", obj.type[match.ind],
         "_proportions-", proportions.long,
         "_method-", method,
         "_replication-", replication, ".rds")
}

#use an md5 hash to limit the number of files per directory
hashDir <- function(filename){
  NUM.SUBDIR <- 3
  LENGTH.SUBDIR <- 2
  hash <- digest::sha1( filename, "md5")
  # TODO(riley) what is the range of possible hash lengths?
  if ( nchar(hash) < 2 * NUM.SUBDIR)
    stop("A subsample does not have a sufficiently long hash to be placed in a directory")
  #directories are specified by pairs of characters, begining with the first character
  start.chars <- 1 + ((0:(NUM.SUBDIR - 1) * (LENGTH.SUBDIR))) # (start index) + (offset)
  substr.pairs <- lapply( start.chars, function(ii) substr(hash, ii, ii + LENGTH.SUBDIR - 1))
  do.call( file.path, substr.pairs)
}

# a function that does the fundamental aspects of saving an object
# that results from a single subsampling
# proportions should be a list of vectors
masterSave <- function(obj, base.dir, expt.id){
  # TODO(riley) should the subsamples object store information about the proportion of each gene in the .RDS file?
  #make sure the directory exists and is formatted correctly
  if ( !dir.exists( base.dir))
    stop("BASE.DIR not found. It should refer to the directory used to store the hashed directory structure")
  
  #Extract design parameters and check they are unique
  #select replications
  obj <- ungroup(obj)
  replication <- obj %>% select( replication) %>% distinct()
  replication <- replication$replication
  if (length(replications) != 1){
    stop("can only save a subsampling with a single replication")
  }
  #select proportions
  proportions <- obj  %>% select( contains("prop.")) %>% distinct() %>% as.data.frame()
  if (dim( proportions)[1] > 1){
    stop("can only save a subsampling with a single combination of subsampling proportions")
  }
  method <- obj  %>% select( method) %>% distinct() %>% as.data.frame()
  method <- method$method
  if (length(method) > 1){
    stop("can only save a subsampling with a single method")
  }
  #save the file
  type.obj <- attr(obj, "class")
  out.file <- masterName( expt.id= expt.id, obj.type = type.obj, method= method, proportions = proportions, replication = replication)
  #use of hashDir avoids large numbers of files in a single directory
  out.dir <- file.path( base.dir, hashDir( out.file))
  ifelse(!dir.exists(out.dir), dir.create(out.dir, recursive=TRUE), FALSE)
  out <- file.path( base.dir, hashDir( out.file), out.file)
  saveRDS(obj, file = out)
}

# a helper function that can be used to summarize a subsampling and
# save both the original subsampling and the summary
summarizeAndSave <- function(base.dir, expt.id, oracle = NULL){
  function(ss){
    masterSave(obj = ss, base.dir = base.dir, expt.id = expt.id)
    sm <- summary( ss, oracle)
    masterSave(obj = sm, base.dir = base.dir, expt.id = expt.id)
  }
}

saveWithoutSummary <- function(base.dir, expt.id){
  function(ss){
    masterSave(obj = ss, base.dir = base.dir, expt.id = expt.id)
  }
}



# a function that does the fundamental aspects of loading objects that involve repeated subsampling
# of reads and biological replicates
masterRead <- function(base.dir, expt.id, method, obj.type, proportions, replication){
  #ensure that the diectory is correctly formatted
  if ( !dir.exists( base.dir))
    stop("BASE.DIR not found. It should refer to the directory used to store the hashed directory structure")
  
  #get file name by knowing the naming convention
  in.file <- masterName( expt.id = expt.id, obj.type = obj.type, method = method, proportions = proportions, replication = replication)
  hash.dir <- hashDir(filename = in.file)
  readRDS( file.path(base.dir, hash.dir, in.file))
}

readSsSummary <- function(proportions, replication, method, expt.id, base.dir){
  masterRead( base.dir= base.dir, expt.id = expt.id, method = method, obj.type = "summary.subsamples", proportions, replication)
}

#given the metadata on a particular subsampling, load the data
readSs <- function(proportions, replication, method, expt.id, base.dir){
  masterRead( base.dir= base.dir, expt.id = expt.id, method = method, obj.type = "subsamples", proportions, replication)
}
