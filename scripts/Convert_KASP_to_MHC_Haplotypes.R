# Script to convert the raw KASP genotypes into MHC haplotypes

library(tidyverse)

# Read in KASP data
kasp.data.raw <- read.table("data/KASP_formatted.txt", sep = "\t",
                        stringsAsFactors = F)

kasp.data.raw[kasp.data.raw == "0:0"] <- NA  # Convert missing data (0:0) to NA


# Remove parentage & sex data to a new df
kasp.ind.info <- select(kasp.data.raw, ID, MOTHER, FATHER, DatabaseSex)

kasp.data <- select(kasp.data.raw, -MOTHER, -FATHER, -DatabaseSex)

#~~ remove DQA1 SNPs to a new df
DQA1 <- select(kasp.data, ID, MasterPlate, Well, DQA1_171, DQA1_195)
kasp.data <- subset(kasp.data, select = -c(DQA1_171, DQA1_195))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Process Duplicates                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ First of all, find samples that were rubbish

bad.snps <- apply(kasp.data[,4:ncol(kasp.data)], MARGIN = 2, function(x) length(which(is.na(x)))/length(x))
sort(bad.snps)


kasp.data$Missing <- apply(kasp.data[,4:ncol(kasp.data)], MARGIN = 1, function(x) length(which(is.na(x))))

kasp.data %>% 
  count(Missing)

#~~ remove anything with equal or more than 10 missing genos

kasp.data <- subset(kasp.data, Missing < 10)
kasp.data <- subset(kasp.data, select = -Missing)

kasp.data.beforeDups <- kasp.data


#~~ duplicate IDs

dup.ids <- data.frame(table(kasp.data$ID))
dup.ids <- subset(dup.ids, Freq > 1)
names(dup.ids)[1] <- "ID"


#~~ Try doubles first. Check also if mismatch with parents

double.tab <- subset(dup.ids, Freq == 2) # indivudals occurring twice
double.tab$Mismatches <- NA
double.tab$Matches    <- NA

double.tab$Fam1 <- NA
double.tab$Fam2 <- NA

double.tab$Fam1Count <- NA
double.tab$Fam2Count <- NA

for(i in 1:nrow(double.tab)){
  
  print(paste("Running row", i, "of", nrow(double.tab)))
  x <- subset(kasp.data, ID == double.tab$ID[i]) # subset kasp.data with indiviudals in double.tab but only indivudal i
  
  double.tab$Fam1[i] <- x[1,2] 
  double.tab$Fam2[i] <- x[2,2] # paste first row and MasterPlate (col 2) in Fam1, then paste the second row and MasterPlate in Fam2
  
  double.tab$Fam1Count[i] <- length(na.omit(unlist(c(x[1,4:ncol(x)])))) # number of genotypes occurrance 1
  double.tab$Fam2Count[i] <- length(na.omit(unlist(c(x[2,4:ncol(x)])))) # number of genotypes occurrance 2
  
  na.cols <- unique(c(which(is.na(x[1,])), which(is.na(x[2,])))) # create a list of the columns with NA
  if(length(na.cols) > 0) x <- x[,-na.cols] # remove the column with NA
  x1 <- table(x[1,4:ncol(x)] == x[2,4:ncol(x)]) # create a table which has the counts of FALSE (mismatches) and TRUE (matches)
  
  if("FALSE" %in% names(x1)) double.tab$Mismatches[i] <- x1[["FALSE"]]
  if("TRUE" %in% names(x1)) double.tab$Matches[i] <- x1[["TRUE"]] # add these counts to the double.tab df
  
  rm(x, x1)
}


bad.pairs <- subset(double.tab, Mismatches > 0)

head(bad.pairs)

bad.ids <- subset(kasp.data, ID %in% bad.pairs$ID)


#~~ Do they match mum?

bad.pairs <- plyr::join(bad.pairs, kasp.ind.info)
bad.pairs <- subset(bad.pairs, !is.na(MOTHER))

bad.pairs$MOTHER <- ifelse(bad.pairs$MOTHER %in% kasp.data$ID, bad.pairs$MOTHER, NA)
bad.pairs <- subset(bad.pairs, !is.na(MOTHER))
head(bad.pairs)


for(i in 1:nrow(bad.pairs)){
  
  x <- subset(kasp.data, ID == bad.pairs$ID[i])
  y <- subset(kasp.data, ID == bad.pairs$MOTHER[i])
  
  x <- rbind(x, y)
  x <- data.frame(t(x[,4:ncol(x)]), stringsAsFactors = F)
  x <- na.omit(x)
  
  
  x$Match1 <- NA
  x$Match2 <- NA
  
  
  for(j in 1:nrow(x)){
    
    x$Match1[j] <- ifelse(any(strsplit(x[j,1], split = ":")[[1]] %in% strsplit(x[j,3], split = ":")[[1]]), "yes", "no")
    x$Match2[j] <- ifelse(any(strsplit(x[j,2], split = ":")[[1]] %in% strsplit(x[j,3], split = ":")[[1]]), "yes", "no")
    
  }
  
  bad.pairs$MumMatch1[i] <- length(which(x$Match1 == "yes"))
  bad.pairs$MumMatch2[i] <- length(which(x$Match2 == "yes"))
  
  #~~ Dad
  
  rm(x)
  
  x <- subset(kasp.data, ID == bad.pairs$ID[i])
  y <- subset(kasp.data, ID == bad.pairs$FATHER[i])
  
  x <- rbind(x, y)
  x <- data.frame(t(x[,4:ncol(x)]), stringsAsFactors = F)
  x <- na.omit(x)
  
  
  x$Match1 <- NA
  x$Match2 <- NA
  
  if(ncol(x) == 5){
    
    for(j in 1:nrow(x)){
      
      x$Match1[j] <- ifelse(any(strsplit(x[j,1], split = ":")[[1]] %in% strsplit(x[j,3], split = ":")[[1]]), "yes", "no")
      x$Match2[j] <- ifelse(any(strsplit(x[j,2], split = ":")[[1]] %in% strsplit(x[j,3], split = ":")[[1]]), "yes", "no")
      
    }
    
    bad.pairs$DadMatch1[i] <- length(which(x$Match1 == "yes"))
    bad.pairs$DadMatch2[i] <- length(which(x$Match2 == "yes"))
    
  } else {
    bad.pairs$DadMatch1[i] <- NA
    bad.pairs$DadMatch2[i] <- NA
  }
  
  
}


# Decide what to keep.

double.tab$Keep <- ifelse(double.tab$Fam2Count > double.tab$Fam1Count, 2, 1)

for(i in 1:nrow(bad.pairs)){
  
  sum1 <- sum(bad.pairs$MumMatch1[i], bad.pairs$DadMatch1[i], na.rm = T)
  sum2 <- sum(bad.pairs$MumMatch2[i], bad.pairs$DadMatch2[i], na.rm = T)
  
  double.tab$Keep[which(double.tab$ID == bad.pairs$ID[i])] <- ifelse(sum1 > sum2, 1, ifelse(sum2 > sum1, 2, ifelse(double.tab$Fam2Count > double.tab$Fam1Count, 2, 1)))
  
}




#~~ Triples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

triple.tab <- subset(dup.ids, Freq == 3)

triple.tab$Fam1 <- NA
triple.tab$Fam2 <- NA
triple.tab$Fam3 <- NA

triple.tab$Fam1Count <- NA
triple.tab$Fam2Count <- NA
triple.tab$Fam3Count <- NA

triple.tab$Mismatches12 <- NA
triple.tab$Mismatches13 <- NA
triple.tab$Mismatches23 <- NA


for(i in 1:nrow(triple.tab)){
  
  print(paste("Running row", i, "of", nrow(triple.tab)))
  x <- subset(kasp.data, ID == triple.tab$ID[i])
  
  triple.tab$Fam1[i] <- x[1,2]
  triple.tab$Fam2[i] <- x[2,2]
  triple.tab$Fam3[i] <- x[3,2]
  
  triple.tab$Fam1Count[i] <- length(na.omit(unlist(c(x[1,4:ncol(x)]))))
  triple.tab$Fam2Count[i] <- length(na.omit(unlist(c(x[2,4:ncol(x)]))))
  triple.tab$Fam3Count[i] <- length(na.omit(unlist(c(x[3,4:ncol(x)]))))
  
  na.cols <- unique(c(which(is.na(x[1,])), which(is.na(x[2,])), which(is.na(x[3,]))))
  if(length(na.cols) > 0) x <- x[,-na.cols]
  
  x12 <- table(x[1,4:ncol(x)] == x[2,4:ncol(x)])
  x13 <- table(x[1,4:ncol(x)] == x[3,4:ncol(x)])
  x23 <- table(x[2,4:ncol(x)] == x[3,4:ncol(x)])
  
  
  if(length(x12) > 1) triple.tab$Mismatches12[i] <- x12[["FALSE"]]
  if(length(x13) > 1) triple.tab$Mismatches13[i] <- x13[["FALSE"]]
  if(length(x23) > 1) triple.tab$Mismatches23[i] <- x23[["FALSE"]]
  
  rm(x, x12, x13, x23)
}


triple.tab$Keep <- c(3, 1, 1)


#~~ Golden Sheep

golden <- subset(kasp.data, ID == 7658)

for(i in 4:ncol(golden)) print(table(golden[,i]))

golden$Keep <- 0
golden$Keep[1] <- 1

kasp.data <- kasp.data[-which(kasp.data$ID == 7658 & kasp.data$MasterPlate != "Plate 35"),]

#~~~~~~~~~ DETERMINE LINES TO REMOVE


for(i in 1:nrow(double.tab)){
  
  if(double.tab$Keep[i] == 0){
    kasp.data <- subset(kasp.data, ID != double.tab$ID[i])
  }
  
  
  if(double.tab$Keep[i] == 1){
    kasp.data <- kasp.data[-which(kasp.data$ID == double.tab$ID[i] & kasp.data$MasterPlate == double.tab$Fam2[i]),]
  }
  
  if(double.tab$Keep[i] == 2){
    kasp.data <- kasp.data[-which(kasp.data$ID == double.tab$ID[i] & kasp.data$MasterPlate == double.tab$Fam1[i]),]
  }
  
  if(nrow(kasp.data) == 0) print(i)
  
}


for(i in 1:nrow(triple.tab)){
  
  if(triple.tab$Keep[i] == 0){
    kasp.data <- subset(kasp.data, ID != triple.tab$ID[i])
  }
  
  
  if(triple.tab$Keep[i] == 1){
    kasp.data <- kasp.data[-which(kasp.data$ID == triple.tab$ID[i] & !kasp.data$MasterPlate == triple.tab$Fam1[i]),]
  }
  
  if(triple.tab$Keep[i] == 2){
    kasp.data <- kasp.data[-which(kasp.data$ID == triple.tab$ID[i] & !kasp.data$MasterPlate == triple.tab$Fam2[i]),]
  }
  
  if(triple.tab$Keep[i] == 3){
    kasp.data <- kasp.data[-which(kasp.data$ID == triple.tab$ID[i] & !kasp.data$MasterPlate == triple.tab$Fam3[i]),]
  }
  
  if(nrow(kasp.data) == 0) print(i)
  
}


rm(triple.tab, bad.ids, bad.pairs, dup.ids, bad.snps, i, j, na.cols, sum1, sum2, golden)

table(table(kasp.data$ID))



#### This code has also removed individuals which occur twice on the same plate, rather than keeping the best rows. need to add them back in. 

missingIDs <- setdiff(kasp.data.beforeDups$ID, kasp.data$ID)
missingIDs <- subset(double.tab, ID %in% missingIDs)

missingIDs$Fam1 <- NA
missingIDs$Fam2 <- NA

missingIDs$Well1 <- NA
missingIDs$Well2 <- NA

missingIDs$Fam1Count <- NA
missingIDs$Fam2Count <- NA


for(i in 1:nrow(missingIDs)){
  
  print(paste("Running row", i, "of", nrow(missingIDs)))
  x <- subset(kasp.data.beforeDups, ID == missingIDs$ID[i]) # subset kasp.data with indiviudals in double.tab but only indivudal i
  
  missingIDs$Fam1[i] <- x[1,2] 
  missingIDs$Fam2[i] <- x[2,2] # paste first row and MasterPlate (col 2) in Fam1, then paste the second row and MasterPlate in Fam2
  
  missingIDs$Well1[i] <- x[1,3] 
  missingIDs$Well2[i] <- x[2,3] 
  
  missingIDs$Fam1Count[i] <- length(na.omit(unlist(c(x[1,4:ncol(x)])))) # number of genotypes occurrance 1
  missingIDs$Fam2Count[i] <- length(na.omit(unlist(c(x[2,4:ncol(x)])))) # number of genotypes occurrance 2
  
  na.cols <- unique(c(which(is.na(x[1,])), which(is.na(x[2,])))) # create a list of the columns with NA
  if(length(na.cols) > 0) x <- x[,-na.cols] # remove the column with NA
  x1 <- table(x[1,4:ncol(x)] == x[2,4:ncol(x)]) # create a table which has the counts of FALSE (mismatches) and TRUE (matches)
  
  if("FALSE" %in% names(x1)) missingIDs$Mismatches[i] <- x1[["FALSE"]]
  if("TRUE" %in% names(x1)) missingIDs$Matches[i] <- x1[["TRUE"]] # add these counts to the double.tab df
  
  rm(x, x1)
}


#~~~~~~~~~ DETERMINE LINES TO REMOVE

kasp.data.missing <- subset(kasp.data.beforeDups, ID %in% missingIDs$ID)


for(i in 1:nrow(missingIDs)){
  
  if(missingIDs$Keep[i] == 0){
    kasp.data.missing <- subset(kasp.data.missing, ID != missingIDs$ID[i])
  }
  
  
  if(missingIDs$Keep[i] == 1){
    kasp.data.missing <- kasp.data.missing[-which(kasp.data.missing$ID == missingIDs$ID[i] & kasp.data.missing$Well == missingIDs$Well2[i]),]
  }
  
  if(missingIDs$Keep[i] == 2){
    kasp.data.missing <- kasp.data.missing[-which(kasp.data.missing$ID == missingIDs$ID[i] & kasp.data.missing$Well == missingIDs$Well1[i]),]
  }
  
  if(nrow(kasp.data.missing) == 0) print(i)
  
}


kasp.data <- rbind(kasp.data, kasp.data.missing)
table(table(kasp.data$ID))



rm(double.tab, missingIDs)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. PLINK                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#######################
######################
#######################


#~~ Load info on loci but remove the DQA loci

info <- read.table("data/SNP.Info.txt", header = T, stringsAsFactors = F, sep = "\t")
head(info)
info <- subset(info, SNP.ID != "DQA1_171" & SNP.ID != "DQA1_195")

info$Map <- sapply(info$Position, function(x) strsplit(x, split = "_")[[1]][2])

info$Chromosome <- sapply(info$Position, function(x) strsplit(x, split = "_")[[1]][1])
info$Chromosome <- as.numeric(gsub("Chr", "", info$Chromosome))

info <- arrange(info, Chromosome, Map)

#~~ Create map file by extracting first 4 columns of tped file

mapfile <- info[,c("Chromosome", "SNP.ID", "Map")]

write.table(mapfile, "output/KASP.map", col.names  = F, row.names = F, quote = F)


#~~ Create a PLINK file (ped)

# add in the sex information downloaded from the database
# kasp.ind.info$sex <- kasp.ind.info$DatabaseSex
# kasp.ind.info$sex[which(kasp.ind.info$sex == 1)] <- 0
# kasp.ind.info$sex[which(kasp.ind.info$sex == 2)] <- 1
# 
# kasp.ind.info$sex[which(kasp.ind.info$ID == 3614)] <- 0
# kasp.ind.info$sex[which(kasp.ind.info$ID == 5871)] <- 1
# 
# # head(kasp.ind.info)
# 
# kasp.data <- plyr::join(kasp.data, kasp.ind.info)
# table(kasp.data$sex, useNA = "always")
# 
# 
# kasp.data$sex[is.na(kasp.data$sex)] <- 0
# 
# abeldata.file <- read.table("./data/Plates1-77merged.QC3.KASP.pheno", header = T, stringsAsFactors = F)
# 
# names(abeldata.file) <- c("ID", "blah", "oldsex")
# 
# kasp.data <- join(kasp.data, abeldata.file)
# 
# kasp.data$sex[which(kasp.data$sex != kasp.data$oldsex)] <- kasp.data$oldsex[which(kasp.data$sex != kasp.data$oldsex)]
# 


plink.file <- data.frame(Family = 1,
                         ID = kasp.data$ID,
                         MOTHER = 0,
                         FATHER = 0,
                         sex = 0,
                         trait = 0)


plink.file <- cbind(plink.file, kasp.data[,mapfile$SNP.ID])

# head(plink.file)

#~~ Recode as SNP genotypes. For indels, will make it A/T for absence/presence

for(i in 7:ncol(plink.file)){
  plink.file[,i] <- gsub(":", " ", plink.file[,i])
}

plink.file$B <- gsub("-", "A", plink.file$B)
plink.file$B <- gsub("AGGAA", "T", plink.file$B)

plink.file[is.na(plink.file)] <- "0 0"

write.table(plink.file, "output/KASP.ped", col.names  = F, row.names = F, quote = F)



#~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Run through PLINK     #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

## Make sure plink (plink.exe v ) is present in the 

############## Read the data into plink and run it's QCs
### read the .ped file into plink and make a .bed file
plink.file <- "output/KASP"

system("cmd", input = paste0("plink --file output/KASP --sheep --mind 0.14 --recode --out output/KASP_QC"))
# input = 5480 sheep, out = 5378 sheep, 22 SNPs

#~~ run thru with -het to calculate homozygosity per individual
system("cmd", input = paste0("plink --file output/KASP_QC  --sheep --het --recode --out output/KASP_QC"))



### read the .ped file, extract MHC SNPs by selecting Chr 20 only, and make a vcf file
system("cmd", input = paste0("plink --file output/KASP_QC --sheep --chr 20 --recode vcf --out output/KASP_QC"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. BEAGLE                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#######################
######################
#######################

# Run through Beagle on the command line using the following code: 

# java -Xmx3g -jar beagle.21Oct15.abc.jar gt=output/KASP.vcf out=output/KASP_QC_bgl



# extract the haplotypes using custom function

source("scripts/ParseBeagle.R")

haplo.raw <- ExtractHaplotypes("output/KASP_QC_bgl.vcf")
haplo.res <- cast(haplo.raw$haplotypes, ID ~ HaploID)
haplo.res.raw <- cast(haplo.raw$haplotypes, ID ~ HaploID)


# open the reference SNP to MHC haplotypes
ref.haplo <- read.table("data/Ref_MHC_haplotypes.txt", header=T, stringsAsFactors = F)

haplo.freq <- arrange(as.data.frame(table(haplo.raw$haplotypes$Haplo)), desc(Freq))
names(haplo.freq)[1] <- "Haplo"
# Add in the haplotype name
haplo.freq <- merge(haplo.freq, ref.haplo, all.x=T, by="Haplo")
haplo.freq <- arrange(haplo.freq, HaplotypeID)


haplo.melt <-arrange(haplo.raw$haplotypes, ID)
haplo.melt <- merge(haplo.melt, ref.haplo, all.x=T, by="Haplo")
haplo.melt <- arrange(select(haplo.melt, ID, HaploID, Haplo, HaplotypeID), ID)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Mendelian Errors     Haplotypes      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#pedigree <- read.table("./data/4_Updated_Pedigree_Feb2017.txt", header = T, stringsAsFactors = F)
pedigree <- kasp.ind.info


haplo.res <- cast(haplo.melt, ID ~ HaploID)

haplo.res$MHChaplo <- paste(haplo.res$Haplo1, haplo.res$Haplo2, sep = ":")

# haplo.resNA <- haplo.res[grep("NA", haplo.res$MHChaplo),]
# haplo.res <- haplo.res[-grep("NA", haplo.res$MHChaplo),]


haplo.mend <- haplo.res[c("ID", "MHChaplo")]

haplo.parent <- haplo.mend

ped.melt <- melt(pedigree, id.vars = "ID")
names(ped.melt) <- c("ID", "Parent", "ParentID")

haplo.mend <- plyr::join(haplo.mend, ped.melt)

names(haplo.parent) <- c("ParentID", "Parent.value")

haplo.mend <- plyr::join(haplo.mend, haplo.parent)
haplo.mend <- na.omit(haplo.mend)

#~~ Are the genotypes incompatible?

comp.func <- function(x, y){
  if(any(is.na(x), is.na(y))){
    NA
  } else {
    ifelse(any(strsplit(x, split = ":")[[1]] %in% strsplit(y, split = ":")[[1]]), "match", "mismatch")
  }
}

haplo.mend$ParentMatch <- mapply(comp.func, haplo.mend$MHChaplo, haplo.mend$Parent.value)

table(haplo.mend$ParentMatch, haplo.mend$Parent)

rm(haplo.parent, ped.melt)

## Firstly, which Mendelian errors are due to unknown haplotypes (NAs)

haplo.mend$MismatchReason <- ifelse(haplo.mend$ParentMatch == "match", NA, ifelse(grepl("NA", haplo.mend$MHChaplo), "NewHaplo", ifelse(grepl("NA", haplo.mend$Parent.value), "NewHaplo", "Other")))

haplo.mend.newh <- subset(haplo.mend, MismatchReason == "NewHaplo")
haplo.mend$MismatchReason <- ifelse(haplo.mend$ID %in% haplo.mend.newh$ID, "NewHaplo", haplo.mend$MismatchReason)
haplo.mend.newh <- subset(haplo.mend, MismatchReason == "NewHaplo")

haplo.mend.errors <- subset(haplo.mend, MismatchReason == "Other")

haplo.mend.ok <- subset(haplo.mend, is.na(MismatchReason) )

rownames(haplo.mend.errors)<-NULL
haplo.mend.errors
## reveals 3 individuals with haplotype mismatches. One individual (4597) mismatches both parents. 


#~~ remove the individuals that mismatches parents
haplo.res <- subset(haplo.res, ID != 4597 & ID != 5864 & ID != 8450)


# remove the individuals with novel haplotypes
haplo.res2 <- subset(haplo.res, !is.na(Haplo1) & !is.na(Haplo2))

haplo.res2$Hetz <- ifelse(haplo.res2$Haplo1 == haplo.res2$Haplo2, "Homozygote", "Heterozgyote")










#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Recombinant haplotypes - DQA SNPs    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ DQA SNP genotypes were stored at stage 1 in the dataframe DQA1 (updated IDs but duplicates still present)
head(DQA1)

# need to match the plate/well of an ID to those included in the final mhc data
DQA1$ID_plate_well <- paste(DQA1$ID, DQA1$MasterPlate, DQA1$Well, sep="-")
# create the same in the kasp.data df
kasp.data$ID_plate_well <- paste(kasp.data$ID, kasp.data$MasterPlate, kasp.data$Well, sep="-")

# subset DQA1 to include on those with a match in kasp.data
DQA1 <- subset(DQA1, ID_plate_well %in% kasp.data$ID_plate_well)

# merge the DQA1 data with the MHC haplotypes data. 
kasp.DQA1 <- merge(DQA1, haplo.res2, by="ID", all.y=T) 
kasp.DQA1 <- select(kasp.DQA1, -ID_plate_well)
kasp.DQA1$Haplo1 <- as.character(kasp.DQA1$Haplo1)
kasp.DQA1$Haplo2 <- as.character(kasp.DQA1$Haplo2)



#~ Assign haplotype classes as either H, null or other
kasp.DQA1$HaploA <- ifelse(kasp.DQA1$Haplo1 == "H", "H", ifelse(kasp.DQA1$Haplo1 == "A"| kasp.DQA1$Haplo1 == "G" | kasp.DQA1$Haplo1 == "E", "null", "other"))

kasp.DQA1$HaploB <- ifelse(kasp.DQA1$Haplo2 == "H", "H", ifelse(kasp.DQA1$Haplo2 == "A"| kasp.DQA1$Haplo2 == "G" | kasp.DQA1$Haplo2 == "E", "null", "other"))

kasp.DQA1$MHChaploAB <- paste(kasp.DQA1$HaploA, kasp.DQA1$HaploB, sep = ":")

kasp.DQA1$MHChaploAB <- ifelse(kasp.DQA1$MHChaploAB == "other:null", "null:other", ifelse(kasp.DQA1$MHChaploAB == "H:null", "null:H", ifelse(kasp.DQA1$MHChaploAB == "H:other", "other:H", kasp.DQA1$MHChaploAB)))



#~create a column that defines whether the genotype is expected or not
kasp.DQA1$expected171 <- "unexpected"
kasp.DQA1$expected171 <- ifelse(kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_171 == "G:G", "expected", ifelse(kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171 == "G:C", "expected", ifelse(kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171 == "G:G", "expected", ifelse(kasp.DQA1$MHChaploAB == "null:null" & is.na(kasp.DQA1$DQA1_171), "expected", ifelse(kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171 == "C:C", "expected", ifelse(kasp.DQA1$MHChaploAB == "H:H" & kasp.DQA1$DQA1_171 == "C:C", "expected",  "unexpected"))))))

kasp.DQA1$expected171 <- ifelse((is.na(kasp.DQA1$DQA1_171) & kasp.DQA1$MHChaploAB != "null:null"), "failed", kasp.DQA1$expected171)


#~create a column by which to colour the figures - expected genotype is blue
kasp.DQA1$expected195 <- "unexpected"
kasp.DQA1$expected195 <- ifelse(kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_195 == "A:A", "expected", ifelse(kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_195 == "C:A", "expected", ifelse(kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_195 == "A:A", "expected", ifelse(kasp.DQA1$MHChaploAB == "null:null" & is.na(kasp.DQA1$DQA1_195), "expected", ifelse(kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_195 == "C:C", "expected", ifelse(kasp.DQA1$MHChaploAB == "H:H" & kasp.DQA1$DQA1_195 == "C:C", "expected",  "unexpected"))))))

kasp.DQA1$expected195 <- ifelse((is.na(kasp.DQA1$DQA1_195) & kasp.DQA1$MHChaploAB != "null:null"), "failed", kasp.DQA1$expected195)


#### Create a column that accounts for both DQA1 loci
kasp.DQA1$DQA1_expected <- "unknown"

kasp.DQA1$DQA1_expected <- ifelse(kasp.DQA1$expected171 == "unexpected" | kasp.DQA1$expected195 == "unexpected" , "unexpected", ifelse(kasp.DQA1$expected171 == "failed" | kasp.DQA1$expected195 == "failed" , "failed", ifelse(kasp.DQA1$MHChaplo == "H:H" & kasp.DQA1$DQA1_171 == "G:C", "Recomb/dropout", ifelse(kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171 == "G:G", "Recomb/dropout", ifelse(kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171 == "G:G", "Recomb/dropout","expected")))))

kasp.DQA1$DQA1_171F <- as.character(kasp.DQA1$DQA1_171)
kasp.DQA1$DQA1_195F <- as.character(kasp.DQA1$DQA1_195)

kasp.DQA1[c("DQA1_171F", "DQA1_195F")][is.na(kasp.DQA1[c("DQA1_171F", "DQA1_195F")])] <- "??"

kasp.DQA1$DQA1_171F <- as.factor(kasp.DQA1$DQA1_171F)
kasp.DQA1$DQA1_195F <- as.factor(kasp.DQA1$DQA1_195F)

kasp.DQA1$DQA1_171F <- factor(kasp.DQA1$DQA1_171F, levels= c("C:C", "G:C", "G:G", "??"))
kasp.DQA1$DQA1_195F <-  factor(kasp.DQA1$DQA1_195F, levels= c("C:C", "C:A", "A:A", "??"))

#Expected
kasp.DQA1$DQA1_expected <- 
  ifelse(kasp.DQA1$MHChaploAB == "H:H" & kasp.DQA1$DQA1_171F == "C:C" & kasp.DQA1$DQA1_195F == "C:C", "expected", ifelse(
    kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "C:C" & kasp.DQA1$DQA1_195F == "C:C", "expected", ifelse(
      kasp.DQA1$MHChaploAB == "null:null" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "??", "expected", ifelse(
        kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "expected", ifelse(
          kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171F == "G:C" & kasp.DQA1$DQA1_195F == "C:A", "expected", ifelse(
            kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "expected", "other"))))))

# Expected & Unxpected
kasp.DQA1$DQA1_expected <- 
  ifelse(kasp.DQA1$MHChaploAB == "H:H" & kasp.DQA1$DQA1_171F == "C:C" & kasp.DQA1$DQA1_195F == "C:C", "Expected", ifelse(
    kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "C:C" & kasp.DQA1$DQA1_195F == "C:C", "Expected", ifelse(
      kasp.DQA1$MHChaploAB == "null:null" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "??", "Expected", ifelse(
        kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "Expected", ifelse(
          kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171F == "G:C" & kasp.DQA1$DQA1_195F == "C:A", "Expected", ifelse(
            kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "Expected", ifelse(
              kasp.DQA1$MHChaploAB == "H:H" & kasp.DQA1$DQA1_171F == "G:C" & kasp.DQA1$DQA1_195F == "C:A", "Recomb/dropout", ifelse(
                kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "Recomb/dropout", ifelse(
                  kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "A:A", "Recomb/dropout" , ifelse(
                    kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "C:C" & kasp.DQA1$DQA1_195F == "??", "Fail",  ifelse(
                      kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "??", "Fail", ifelse(
                        kasp.DQA1$MHChaploAB == "null:H" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "C:C", "Fail",ifelse(           
                          kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "??", "Fail",  ifelse(
                            kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "??", "Fail", ifelse(  
                              kasp.DQA1$MHChaploAB == "null:other" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "A:A", "Fail", ifelse( 
                                kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171F == "G:C" & kasp.DQA1$DQA1_195F == "??", "Fail", ifelse(
                                  kasp.DQA1$MHChaploAB == "other:H" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "C:A", "Fail", ifelse(
                                    kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_171F == "G:G" & kasp.DQA1$DQA1_195F == "??", "Fail",  ifelse(
                                      kasp.DQA1$MHChaploAB == "other:other" & kasp.DQA1$DQA1_171F == "??" & kasp.DQA1$DQA1_195F == "A:A", "Fail",
                                      "Unexpected"  )))))))))))))))))))


# Re-call individuals as recombinants if: H:H, null:H or other:H. Then run mendelian check and see which ones seem to pass through a family

kasp.DQA1$NewHaplo1 <- kasp.DQA1$Haplo1
kasp.DQA1$NewHaplo2 <- kasp.DQA1$Haplo2

kasp.DQA1$NewHaplo1 <- ifelse(kasp.DQA1$DQA1_expected == "Recomb/dropout" & kasp.DQA1$NewHaplo1 == "H", "R", kasp.DQA1$NewHaplo1)
kasp.DQA1$NewHaplo2 <- ifelse(kasp.DQA1$DQA1_expected == "Recomb/dropout" & kasp.DQA1$NewHaplo2 == "H", "R", kasp.DQA1$NewHaplo2)

kasp.DQA1$NewMHC <- paste(kasp.DQA1$NewHaplo1, kasp.DQA1$NewHaplo2, sep = ":")

kasp.DQA1$NewMHC <- ifelse(kasp.DQA1$NewMHC == "R:R", "R:H", kasp.DQA1$NewMHC)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mendelian Errors for DQA1                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

haplo.mend <- kasp.DQA1[c("ID", "NewMHC")]

haplo.parent <- haplo.mend

ped.melt <- melt(pedigree, id.vars = "ID")
names(ped.melt) <- c("ID", "Parent", "ParentID")

haplo.mend <- plyr::join(haplo.mend, ped.melt)

names(haplo.parent) <- c("ParentID", "Parent.value")

haplo.mend <- plyr::join(haplo.mend, haplo.parent)
haplo.mend <- na.omit(haplo.mend)

#~~ Are the genotypes incompatible?

comp.func <- function(x, y){
  if(any(is.na(x), is.na(y))){
    
    NA
  } else {
    ifelse(any(strsplit(x, split = ":")[[1]] %in% strsplit(y, split = ":")[[1]]), "match", "mismatch")
  }
}

haplo.mend$ParentMatch <- mapply(comp.func, haplo.mend$NewMHC, haplo.mend$Parent.value)
rm( haplo.parent, ped.melt)

subset(haplo.mend, ParentMatch == "mismatch") # note the only mismatch is between 4179 & his mother. The recombination is known to have originated in 4179



##############
##### convert the recombinants to haplotype B and export dataset
###############


kasp.DQA1$NewNewMHC <- kasp.DQA1$NewMHC

kasp.DQA1$NewNewMHC <- gsub("R", "B", kasp.DQA1$NewNewMHC)

kasp.DQA1$DQA1_status <- ifelse(kasp.DQA1$NewHaplo1 == "R" | kasp.DQA1$NewHaplo2 == "R", "recomb", ifelse(kasp.DQA1$DQA1_expected == "Fail", "Fail", "normal"))

kasp.DQA1$NewNewHaplo1 <- gsub("R", "B", kasp.DQA1$NewHaplo1)
kasp.DQA1$NewNewHaplo2 <- gsub("R", "B", kasp.DQA1$NewHaplo2)

kasp.DQA1.export <- select(kasp.DQA1, ID, NewNewHaplo1, NewNewHaplo2, NewNewMHC, DQA1_status)

names(kasp.DQA1.export)[names(kasp.DQA1.export)=="NewNewHaplo1"] <- "Haplo1"
names(kasp.DQA1.export)[names(kasp.DQA1.export)=="NewNewHaplo2"] <- "Haplo2"
names(kasp.DQA1.export)[names(kasp.DQA1.export)=="NewNewMHC"] <- "MHC_diplo"






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Add on the Sanger only haplotypes    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Sang <- read.table("./data/SangerOnlyHaplos.txt", header=T, stringsAsFactors = F, sep=",")

Sang <- select(Sang, ID, Haplo1, Haplo2, MHC_diplo, DQA1_status)

kasp.DQA1.export_Sang <- rbind(kasp.DQA1.export, Sang)


# Final output file
write.table(kasp.DQA1.export_Sang, "./output/SNP_Sang_MHC_haplos.txt", quote=F, sep="\t", row.names = F)









