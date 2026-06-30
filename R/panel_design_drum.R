# Pedigree inference panel design
# with SNPs called by GATK
# phasing by whatshap
# mh He calculated by CalcHe_mh_vcf.jar
# greedy algorithm to find loci and design primers
# R v4
# 
# edits: 
#  filtering strategy, 
#  and inputs
# 
# for ceres: module load primer3 bowtie2 r/4.3.0 gcc
# for andromeda: module load Bowtie2; module load R/4.1.2-foss-2021b
# sbatch  -t 12:00:00 --cpus-per-task 1 --mem=20G --wrap="module load Bowtie2; module load R/4.1.2-foss-2021b; Rscript panel_design.R"

# order is important here so slice defaults to dplyr::slice
#.libPaths(c("/project/oyster_gs_sim/R_packages/4.3/", .libPaths()))
library(Rsamtools)
library(tidyverse)
source("utils.R")

# prefix
# 0 or 1 for filterNoRecomb
# 0 or 1 for filterAdj
# 75 or 150 for readlength
cmdArgs <- commandArgs(trailingOnly=TRUE)

###
## Begin inputs
###
prefix <- cmdArgs[1] # prefix on output files
filterNoRecomb <- as.logical(as.numeric(cmdArgs[2])) # only loci with Na <= Nsnps + 1
filterAdj <- as.logical(as.numeric(cmdArgs[3])) # only loci without adjacent SNPs
readLength <- as.numeric(cmdArgs[4]) # anticipated read length for amplicon sequencing
# pop names and he calc output files 
inputfiles <- tibble(
  pop = c("WGULF","EGULF")
)

if(readLength == 75){
  inputfiles <- inputfiles %>% mutate(heFile = paste0("He_50_", pop, ".txt"))
} else if(readLength == 150){
  inputfiles <- inputfiles %>% mutate(heFile = paste0("He_125_", pop, ".txt"))
}
summaryFunction <- min # function to summarize He or equivalent across groups
numLoci <- 500 # target number of loci in the panel
refGenome <- "/work/marinegenomics/jmatt1/wgs/mega_soce_cmh/resources/genome.fasta" # reference genome fasta file
bowtie2_ref <- "/work/marinegenomics/jmatt1/wgs/mega_soce_cmh/resources/genome_bt2"
primer3path <- "/home/jmatt1/.conda/envs/microhap/bin/primer3_core"
fwdAdapter <- "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" # note this is illumina read 1 primer
revAdapter <- "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" # and read 2 primer
# chromosome names and lengths
chrTable <- tibble(
  chr = c(
    "CP136866.1",
    "CP136867.1",
    "CP136868.1",
    "CP136869.1",
    "CP136870.1",
    "CP136871.1",
    "CP136872.1",
    "CP136873.1",
    "CP136874.1",
    "CP136875.1",
    "CP136876.1",
    "CP136877.1",
    "CP136878.1",
    "CP136879.1",
    "CP136880.1",
    "CP136881.1",
    "CP136882.1",
    "CP136883.1",
    "CP136884.1",
    "CP136885.1",
    "CP136886.1",
    "CP136887.1",
    "CP136888.1",
    "CP136889.1"
  ),
  l = c(
    31019971,
    35478940,
    33107411,
    26829222,
    27851347,
    32124745,
    28686924,
    27691033,
    27086610,
    22788609,
    15492723,
    28117437,
    30983267,
    25631868,
    29804843,
    28869540,
    31839088,
    29116862,
    32745965,
    30640032,
    26875509,
    30575656,
    35187552,
    22062508
  )
)
###
## End inputs
###

# reading in He values
for(i in 1:nrow(inputfiles)){
  temp <- read_tsv(inputfiles$heFile[i], col_types = "ccddddc") %>% 
    select(-numFailGeno, -numInvalidPhase)
  colnames(temp)[3:5] <- paste0(inputfiles$pop[i], "_", colnames(temp)[3:5])
  if(i == 1){
    allHeData <- temp
  } else {
    allHeData <- allHeData %>% full_join(temp, by = c("Chr", "Pos"))
  }
}

allSnps <- allHeData %>% select(Chr, Pos)


# filtering

# remove loci missing in any population
boolMissData <- (allHeData %>% mutate(across(.cols = matches("_He$"), .fns = is.na)) %>%
                   select(matches("_He$")) %>% as.matrix() %>% rowSums()) == 0
table(boolMissData)

allHeData <- allHeData %>% filter(boolMissData)

# remove loci with fewer than 10 inds in any population
boolFewInds <- (allHeData %>% mutate(across(.cols = matches("_numInds$"), .fns = (function(x) x < 10))) %>%
                  select(matches("_numInds$")) %>% as.matrix() %>% rowSums()) == 0
table(boolFewInds)

allHeData <- allHeData %>% filter(boolFewInds)



# remove loci with adjacent snps
adjSNPS <- sapply(str_split(allHeData$Pos, ","), function(x){
  if(length(x) < 2) return(FALSE)
  x <- as.numeric(x)
  return(any(x[2:length(x)] - x[1:(length(x) - 1)] < 2))
})
table(adjSNPS)

if(filterAdj){
  print("removing adjacent")
  allHeData <- allHeData %>% filter(!adjSNPS)
}


# remove loci where number of alleles is > numSnps + 1
ns <- sapply(str_split(allHeData$Pos, ","), length)

temp <- allHeData %>% select(matches("_AlleleFreq$"))
if(ncol(temp) > 1){
  for(i in 2:ncol(temp)){
    temp[[1]] <- paste(temp[[1]], temp[[i]], sep = ",")
  }
  temp <- temp[[1]]
}
nAlleles <- sapply(str_split(temp, ","), function(x){
  return(n_distinct(gsub(":.+", "", x)))
})
rm(temp)

table(ns)
table(ns <= 7)
table(nAlleles)
table(nAlleles <= (ns + 1))
table(ns <= 7 & nAlleles <= (ns + 1))

if(filterNoRecomb){
  print("removing recomb")
  allHeData <- allHeData %>% filter(nAlleles <= (ns + 1))
}							

rm(ns, nAlleles, boolMissData, boolFewInds, adjSNPS)

print("number of candidates")
print(nrow(allHeData))


expHet <- tibble() # chrom, pos, mean genotyped, min genotyped, numAlleles, each pop He, summary statistic

sumHe <- allHeData %>% select(Chr, Pos, matches("_He$")) %>%
  pivot_longer(3:ncol(.), names_to = "pop", values_to = "He") %>% 
  group_by(Chr, Pos) %>% summarise(summaryHe = summaryFunction(He)) %>%
  ungroup()

# summary(sumHe$summaryHe)
# sum(sumHe$summaryHe > 0.4)
# hist(sumHe$summaryHe)


# now run greedy algorithm

# assign start and end values for each window and sort
sumHe <- sumHe %>%
  mutate(wStart = as.numeric(gsub(",.+$", "", Pos)),
         wEnd = as.numeric(gsub("^.+,", "", Pos))) %>%
  arrange(Chr, wStart)






# distribute loci according to chr length
chrTable <- chrTable %>% mutate(num = round(numLoci * (l / sum(l))))
# account for rounding error
chrTable$num[which.max(chrTable$num)] <- chrTable$num[which.max(chrTable$num)] + numLoci - sum(chrTable$num)

# index reference genome if not already indexed
if(!file.exists(paste0(refGenome, ".fai"))) indexFa(refGenome)

# make a list of all SNPs
allSnps_temp <- tibble()
for(i in 1:nrow(chrTable)){
  allSnps_temp <- allSnps_temp %>% bind_rows(tibble(chr = chrTable$chr[i],
                                                    pos = allSnps %>% filter(Chr == chrTable$chr[i]) %>% pull(Pos) %>%
                                                      str_split(",") %>% unlist %>% as.numeric %>% unique))
}
allSnps <- allSnps_temp
rm(allSnps_temp)


lociReject <- rep(0, 2) # no primers, no uniquely mapping primers
fullPanel <- tibble()
fullPrimers <- tibble()

# for each chromosome
for(i in 1:nrow(chrTable)){
  cat("Chr:", i, "\n")
  
  # calculate initial scores
  chrCands <- sumHe %>% filter(Chr == chrTable$chr[i]) %>% 
    mutate(score = scoreGreedy(start = 1, end = chrTable$l[i], maf = summaryHe, pos = wStart), 
           lastStart = 1, lastEnd = chrTable$l[i])
  
  panel <- tibble()
  while(nrow(panel) < chrTable$num[i]){
    # this is useful if it get's stuck
    # cat("Current number of markers selected for chromosome", i, ":", nrow(panel), "\n")
    if(nrow(chrCands) < 1){
      warning("Ran out of candidates in chr:", chrTable$chr[i], "\n")
      break
    }
    # select locus
    chosen <- chrCands %>% arrange(desc(score)) %>% slice(1)
    # remove from candidates list (don't need to check Chr b/c all the same Chr)
    chrCands <- chrCands %>% filter(Pos != chosen$Pos)
    
    snpPos <- as.numeric(str_split(chosen$Pos, ",")[[1]]) # snps in locus
    
    # not checking HWE of microhaps b/c checked HWE of SNPs and filtered prior
    # HWE section of script deleted
    
    # primer design
    
    # chunk of reference to fetch
    # start is base such that the last snp is at the end of the read
    # end is 200 past the last snp
    targetRange <- c(max(1, snpPos[length(snpPos)] - readLength + 1), 
                     min(chrTable$l[i], snpPos[length(snpPos)] + 200)) 
    
    # get reference
    ref <- toString(scanFa(file = refGenome, 
                           param = GRanges(seqnames = chrTable$chr[i], 
                                           ranges = IRanges(start = targetRange[1], end = targetRange[2]))))
    # avoid primers overlapping any SNPs
    snpsInRegion <- allSnps %>% 
      filter(chr == chrTable$chr[i], pos >= targetRange[1], pos <= targetRange[2]) %>%
      pull(pos)
    ref <- str_split(ref, "")[[1]]
    ref[snpsInRegion - targetRange[1] + 1] <- "N"
    ref <- paste(ref, collapse = "")
    
    # Primer3 input
    # note primer3 default is 0-based
    # product sizes of 60 – 150 with optimum of 80
    # primer lengths of 10-25 with optimum of 20
    # primer tm of 55-60 with optimum of 57
    
    primer3In <- paste0("SEQUENCE_ID=", chosen$Chr, ":", chosen$Pos, "\n",
                        "SEQUENCE_TEMPLATE=", ref, "\n",
                        "SEQUENCE_TARGET=", snpPos[1] - targetRange[1], ",", snpPos[length(snpPos)] - snpPos[1] + 1, "\n", # start and length
                        "PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=10
PRIMER_MAX_SIZE=25
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_PRODUCT_SIZE_RANGE=60-150
PRIMER_PRODUCT_OPT_SIZE=80
PRIMER_MIN_TM=58
PRIMER_MAX_TM=65
PRIMER_OPT_TM=60
P3_FILE_FLAG=0
PRIMER_FIRST_BASE_INDEX=0
=", sep = "")
    # adjust product size ranges if needed
    if(readLength == 150){
      primer3In <- gsub("PRIMER_PRODUCT_SIZE_RANGE=60-150", "PRIMER_PRODUCT_SIZE_RANGE=120-300", primer3In)
      primer3In <- gsub("PRIMER_PRODUCT_OPT_SIZE=80", "PRIMER_PRODUCT_OPT_SIZE=160", primer3In)
    }
    
    # design primer pairs
    # remember: START POSITION OF PRIMERS IN PRIMER3 OUTPUT HERE IS 0-BASED
    primer3out <- system(paste0("echo \"", primer3In, "\" | ", primer3path), intern = TRUE)
    
    # parse primer3 output
    primerTable <- parsePrimer3(primer3out)
    
    # reject locus if no valid primer pairs found
    if(nrow(primerTable) < 1){
      lociReject[1] <- lociReject[1] + 1
      next
    }
    
    primerTable <- primerTable %>% mutate() %>% 
      mutate(across(.cols = c(num, tm, start, len), .fns = as.numeric)) %>%
      arrange(num, orient)
    
    # check for multiple alignments
    # running all primers for this locus at once (vs running one by one as needed)
    #   b/c bowtie2 has a relatively long startup time
    left <- primerTable %>% filter(orient == "left") %>% pull(seq)
    right <- primerTable %>% filter(orient == "right") %>% pull(seq)
    samFull <- system2("bowtie2", args = c("--end-to-end", "-x", bowtie2_ref, "-k", "6", "-c", 
                                           "-1", paste(left, collapse = ","), "-2", paste(right, collapse = ","), 
                                           "-X", "1000", "--no-hd", "--no-mixed", "--no-discordant"), stdout = TRUE,
                       stderr = FALSE)
    
    acceptPrimers <- FALSE
    for(j in 0:((nrow(primerTable)/2)-1)){ # for each primer pair
      tempPrimer <- primerTable %>% filter(num == j) %>% arrange(orient) # left will be row 1, right 2
      # make sure target is in product
      if(any((snpPos - targetRange[1]) <= (tempPrimer$start[1] + tempPrimer$len[1] - 1) | # last base of fwd primer
             (snpPos - targetRange[1]) >= (tempPrimer$start[2] - tempPrimer$len[2] + 1))){ # first base of rev primer
        save.image("onError.rda")
        stop("target snps not flanked by primers. error in code?")
      }
      
      # get bowtie2 output for this primer pair
      sam <- str_split(samFull[grepl(paste0("^", j, "\t"), samFull)], "\t")
      if(all(sapply(sam, function(x) return(x[3])) == "*")){
        warning("0 alignments for primers chosen for ", chosen$Chr, " ", chosen$Pos, ". ",
                "Manual inspection of these primers recommended.")
      } else{
        # now remove those with edit distances of more than 1
        # need to take each pair at a time
        trem <- c()
        for(k in seq(1, length(sam) - 1, 2)){
          temp <- sam[k:(k + 1)]
          if (as.numeric(gsub("NM:i:", "", temp[[1]][grepl("NM:i:", temp[[1]])])) > 1 ||
              as.numeric(gsub("NM:i:", "", temp[[2]][grepl("NM:i:", temp[[2]])])) > 1){
            trem <- c(trem, k, k + 1)
          }
        }
        sam <- sam[!((1:length(sam)) %in% trem)]
      }
      
      # reject if more than one pair of alignments
      if(length(sam) <= 2){
        acceptPrimers <- TRUE
        break
      }
      
    }
    
    if(acceptPrimers){
      # add to panel
      # make a descriptive name for the locus with no difficult characters
      locusName <- paste0(chosen$Chr, "_", chosen$wStart, "_", chosen$wEnd)
      panel <- panel %>% bind_rows(chosen %>% mutate(Locus = locusName))
      fullPrimers <- fullPrimers %>% 
        bind_rows(tempPrimer %>% mutate(genomeStart = targetRange[1] + start,
                                        Chr = chosen$Chr,
                                        Pos = chosen$Pos,
                                        Locus = locusName))
      
      # remove overlapping loci from candidates
      chrCands <- chrCands[!((chrCands$wStart >= chosen$wStart & chrCands$wStart <= chosen$wEnd) | # start in window OR
                               (chrCands$wEnd >= chosen$wStart & chrCands$wEnd <= chosen$wEnd)),] # end in window
      # recalculate scores
      # note that there is only one interval to update b/c the new locus start position can only
      # have been in one interval
      toUpdate <- chrCands %>% filter(lastStart < chosen$wStart & lastEnd > chosen$wStart) %>%
        select(lastStart, lastEnd) %>% distinct() # this is some bs to help vectorize operations for R
      tempBool <- chrCands$lastStart == toUpdate$lastStart & chrCands$lastEnd == toUpdate$lastEnd & chrCands$wStart < chosen$wStart
      chrCands$score[tempBool] <- scoreGreedy(start = toUpdate$lastStart,
                                              end = chosen$wStart,
                                              maf = chrCands$summaryHe[tempBool],
                                              pos = chrCands$wStart[tempBool])
      chrCands$lastEnd[tempBool] <- chosen$wStart
      tempBool <- chrCands$lastStart == toUpdate$lastStart & chrCands$lastEnd == toUpdate$lastEnd & chrCands$wStart > chosen$wStart
      chrCands$score[tempBool] <- scoreGreedy(start = chosen$wStart,
                                              end = toUpdate$lastEnd,
                                              maf = chrCands$summaryHe[tempBool],
                                              pos = chrCands$wStart[tempBool])
      chrCands$lastStart[tempBool] <- chosen$wStart
      
    } else {
      lociReject[2] <- lociReject[2] + 1
    }
  }
  
  fullPanel <- fullPanel %>% bind_rows(panel)
}

cat("Loci failed HWE:", lociReject[3], "\n",
    "Loci failed primer design:", lociReject[1], "\n", 
    "Loci failed uniquely mapping primers: ", lociReject[2], "\n")

write.table(fullPanel, paste0(prefix, "panel.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fullPrimers, paste0(prefix, "primers.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


# build input files for microhapWrap.py


# fullPanel <- read_tsv("ind_minpanel.txt" , col_types = "cc")
# fullPrimers <- read_tsv("ind_minprimers.txt")

# check to see what % of reference alleles are in the microhap alleles
# create genotyping inputs for microhapWrap
# create reference inputs for HiSat2

# define reverse complement of reverse adaptor
revComp_revAdapter <- revComp(revAdapter)

locusCount <- rep(0,2) # number loci w/ reference allele in microhap discovery alleles, total number of loci
posFile <- tibble()
ampRef <- ""
for(i in 1:nrow(fullPanel)){
  # get microhap alleles
  tempC <- allHeData %>% right_join(fullPanel %>% slice(i) %>% select(Chr, Pos), 
                                    by = c("Chr", "Pos")) %>% select(matches("_AlleleFreq$")) %>%
    unlist() %>% paste(collapse = ",")
  
  tempC <- unlist(str_split(tempC, ","))
  # filter alleles to not include very low frequency - NOT DOING THIS
  # tempC <- tempC[sapply(str_split(tempC, ":"), function(x) as.numeric(x[2]) > 0.01)]
  tempC <- unique(gsub(":.+", "", tempC))
  
  # get reference sequence
  snpPos <- as.numeric(str_split(fullPanel$Pos[i], ",")[[1]]) # snps in locus
  ref <- toString(scanFa(file = refGenome, param = GRanges(seqnames = fullPanel$Chr[i],
                                                           ranges = IRanges(start = snpPos[1], end = snpPos[length(snpPos)]))))
  # make characters slicable
  ref <- str_split(ref, "")[[1]]
  ref <- paste(ref[snpPos - snpPos[1] + 1], collapse = "") # extract the snps and make a microhap allele
  
  # determine if reference allele is in microhap alleles
  locusCount <- locusCount + c(ref %in% tempC, 1)
  
  
  # microhapWrap inputs
  
  # ampRef file
  tempPrim <- fullPrimers %>% 
    filter(Locus == fullPanel$Locus[i]) %>%
    arrange(orient)
  curAmpRef <- toString(scanFa(file = refGenome, param = GRanges(seqnames = fullPanel$Chr[i], 
                                                                 ranges = IRanges(start = tempPrim$genomeStart[1], end = tempPrim$genomeStart[2]))))
  curAmpRef <- paste0(curAmpRef, revComp_revAdapter)
  if(nchar(curAmpRef) < readLength){
    warning("Reference sequence for ", fullPanel$Locus[i], " is too short as you will likely read through the adapter.")
  } else {
    curAmpRef <- substr(curAmpRef, 1, readLength)
  }
  
  ampRef <- paste0(ampRef, ">", fullPanel$Locus[i], "\n", curAmpRef, "\n")
  
  # pos file - Locus Pos (1-based) Type ValidAlt
  vAlt <- matrix(unlist(str_split(tempC, "")), ncol = length(tempC))
  vAlt2 <- rep("", nrow(vAlt))
  for(j in 1:nrow(vAlt)){
    x <- unique(vAlt[j,])
    vAlt2[j] <- paste(x[x != substr(ref, j, j)], collapse = ",")
  }
  
  # in rare cases, a SNP is included but all the samples with an alternate allele
  #   are unphased in the microhap, so the alt allele isn't represented in the 
  # microhap output
  # in this cases, it wasn't included in the He calculation either
  # Issue warning to manually add it back in
  if(any(vAlt2 == "")) warning("manually add alt allele for SNP in ", fullPanel$Locus[i])
  
  posFile <- posFile %>% bind_rows(tibble(Locus = fullPanel$Locus[i], 
                                          RefPos = snpPos - tempPrim$genomeStart[1] + 1,
                                          Type = "S", 
                                          ValidAlt = vAlt2))
  
}

cat(ampRef, file = paste0(prefix, "ampRef.fa"), append = FALSE)
write.table(posFile, paste0(prefix, "pos.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
# fwd file for microhapWrap
fullPrimers %>% filter(orient == "left") %>% select(Locus, seq) %>%
  write.table(paste0(prefix, "fwd.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

cat("Total number of loci:", locusCount[2], "\n",
    "Number with reference allele represented: ", locusCount[1], 
    paste0("(", round((locusCount[1] / locusCount[2]) * 100, 2), "%)"), "\n")

warnings()
