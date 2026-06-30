# utility functions
library(tidyverse)

#' scoring function for greedy algorithem
#' 
#' Score function for greedy algorithm
#' from Matukumalli et al. 2009 https://doi.org/10.1371/journal.pone.0005350
#' @param start start of interval
#' @param end end of interval
#' @param maf maf of markers within interval in original algorithm, but can be any statistic
#' @param pos position of marker within interval (order
#'   corresponds to order of `maf`)
scoreGreedy <- function(start, end, maf, pos){
	return(
		maf * (end - start - abs((2 * pos) - end - start))
	)
}

#' parse primer3 stdout output for primer pairs
#' does not assume consistent order from primer3
#' could be sped up considerably by assuming consistent order
#' @param x stdout from primer3 for example from 
#'   calling primer3out <- system(paste0("echo \"", primer3In, "\" | primer3_core"), intern = TRUE)
#' @return a data frame withe primer number, sequence, left vs right, tm, start position, and length
parsePrimer3 <- function(x){
	outTable <- tibble()
	# get number of pairs returned
	pairsReturned <- as.numeric(gsub("PRIMER_PAIR_NUM_RETURNED=", "", x[grepl("PRIMER_PAIR_NUM_RETURNED=", x)]))
	if(pairsReturned > 0){
		for(i in 0:(pairsReturned - 1)){
			
			# seq
			lSeq <- gsub(".+=", "", x[grepl(paste0("PRIMER_LEFT_", i, "_SEQUENCE="), x)])
			rSeq <- gsub(".+=", "", x[grepl(paste0("PRIMER_RIGHT_", i, "_SEQUENCE="), x)])
			
			# tm
			lTM <- gsub(".+=", "", x[grepl(paste0("PRIMER_LEFT_", i, "_TM="), x)])
			rTM <- gsub(".+=", "", x[grepl(paste0("PRIMER_RIGHT_", i, "_TM="), x)])
			
			# pos and len
			lpos_len <- str_split(gsub(".+=", "", x[grepl(paste0("PRIMER_LEFT_", i, "="), x)]), ",")[[1]]
			rpos_len <- str_split(gsub(".+=", "", x[grepl(paste0("PRIMER_RIGHT_", i, "="), x)]), ",")[[1]]
			
			# now bind into a data.frame
			outTable <- outTable %>% bind_rows(
				tibble(num = i, seq = c(lSeq, rSeq), orient = c("left", "right"),
					   tm = c(lTM, rTM), start = c(lpos_len[1], rpos_len[1]), len = c(lpos_len[2], rpos_len[2]))
			)
		}
	}
	return(outTable)
}

#' reverse complement of DNA string
#' This function is very slow and not optimized
#' For computationally intense operations, use the BioStrings package
#' @param x a string
revComp <- function(x){
	rc <- ""
	dict <- data.frame(a = c("A", "C", "G", "T", "[", "]"),
					   b = c("T", "G", "C", "A", "]", "["))
	for(i in nchar(x):1){
		a <- substr(x,i,i)
		if (a %in% dict$a){
			b <- dict$b[match(a, dict$a)]
		} else {
			b <- a
		}
		rc <- paste0(rc, b)
	}
	return(rc)
}

