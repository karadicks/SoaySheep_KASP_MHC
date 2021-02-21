
library(reshape)

### THIS IS THE FUNCTION. It requires reshape

ExtractHaplotypes <- function(results.file, make.ref.alleles = T){
  require(reshape)
  
  message("Reading and formatting Beagle results file...")
  
  #~~ read in the results file
  
  beagle.res <- read.table(results.file, skip=10, stringsAsFactors = F)
  
  beagle.refalleles <- beagle.res[,4:5]
  
  #~~ transpose data and retain only the haplotype information
  haplotypes.raw <- data.frame(t(beagle.res[,10:ncol(beagle.res)]))
  
  #~~ add ID information and tidy up col names
  
  #beagletab <- read.table(paste0(input.prefix, ".txt"), header = T, stringsAsFactors = F)
  #beaglemap <- read.table(paste0(input.prefix, ".map"), stringsAsFactors = F)
  
  idvec <- read.table(results.file, skip=9, nrows = 1, stringsAsFactors = F, comment.char = "")
  idvec <- idvec[1, 10:ncol(idvec)]
  idvec <- apply(idvec, 2, function(x) strsplit(x, split = "_")[[1]][2])
  idvec <- as.vector(idvec)
  
  haplotypes.raw <- cbind(idvec, haplotypes.raw)
  names(haplotypes.raw) <- c("ID", as.character(beagle.res[,3]))
  
  #~~ get rid of duplicate ID in each ID string
  
  haplotypes.raw$ID <- as.character(haplotypes.raw$ID)
  haplotypes.raw$ID <- sapply(haplotypes.raw$ID, function(x) strsplit(as.character(x), split = " ")[[1]][1])
  haplotypes.raw$ID <- gsub("X", "", haplotypes.raw$ID)
  haplotypes.raw$ID <- gsub("^\\.", "-", haplotypes.raw$ID)
  
  #~~ convert to character
  
  for(i in 1:ncol(haplotypes.raw)) haplotypes.raw[,i] <- as.character(haplotypes.raw[,i])
  
  #~~ The first two characters (i.e. 0|0) is the phased genotype. 
  #   Retain these and substitute in the reference alleles (if specified).
  
  message("Extracting phased genotypes...")
  
  haplotypes.extract <- data.frame(matrix(NA,
                                          nrow = nrow(haplotypes.raw),
                                          ncol = ncol(haplotypes.raw)))
  
  names(haplotypes.extract) <- names(haplotypes.raw)
  haplotypes.extract$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    haplotypes.extract[,i] <- substr(haplotypes.raw[,i], 1, 3)
    
    if(make.ref.alleles == TRUE){
      haplotypes.extract[,i] <- gsub(0, as.character(beagle.refalleles[i-1,1]), haplotypes.extract[,i])
      haplotypes.extract[,i] <- gsub(1, as.character(beagle.refalleles[i-1,2]), haplotypes.extract[,i])
    }
  }
  
  
  #~~ extract the genotype probabilities
  
  message("Extracting genotype probabilities...")
  
  
  genotype.probs <- data.frame(matrix(NA,
                                      nrow = nrow(haplotypes.raw),
                                      ncol = ncol(haplotypes.raw)))
  names(genotype.probs) <- names(haplotypes.raw)
  genotype.probs$ID <- haplotypes.raw$ID
  
  for(i in 2:ncol(haplotypes.raw)){
    genotype.probs[,i] <- sapply(haplotypes.raw[,i], function(x) strsplit(x, split = ":")[[1]][2])
  }
  
  #~~ create a results table for the haplotypes
  
  message("Identifying haplotypes...")
  
  haplo.results <- data.frame(ID = haplotypes.extract$ID,
                              Haplo1 = NA,
                              Haplo2 = NA)
  
  for(i in 1:nrow(haplotypes.extract)){
    haplo.results$Haplo1[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 1, 1), collapse = ""))
    haplo.results$Haplo2[i] <- apply(haplotypes.extract[i,2:ncol(haplotypes.extract)], 1, function(x) paste(substring(x, 3, 3), collapse = ""))
  }
  
  #~~ melt into one column
  require(reshape)
  haplo.results.1col <- melt(haplo.results, id.vars="ID")
  names(haplo.results.1col) <- c("ID", "HaploID", "Haplo")
  
  for(i in 1:ncol(haplo.results.1col)) haplo.results.1col[,i] <- as.character(haplo.results.1col[,i])
  
  return(list(haplotypes = haplo.results.1col,
              phased.genotypes = haplotypes.extract,
              genotype.probabilities = genotype.probs))
  message("...done.")
}
