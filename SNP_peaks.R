#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Get top peaks") %>%
  add_argument("--prefix", help = "ex. <prefix>_chr1.RData") %>%
  add_argument("--p_min", help = "maximum p-value below which to search for peaks",
               default = 5e-4) %>%
  add_argument("--top", help = "number of SNPs peaks to save", default = NULL) %>%
  add_argument("--window", help = "bp distance within which SNPs should be grouped into a peak", default = 15000) %>%
  add_argument("--out_file", help = "file to save SNP peaks", default = "SNP_peaks.rds")

argv <- parse_args(argp)


combine_assoc_by_chr<- function(argv$prefix){

	load(paste0(prefix,"_chr1.RData"))
	assoc_full <- assoc

	for(i in 2:22){
		load(paste0(prefix, "_chr", i, ".RData"))
		assoc_full <- rbind(assoc_full, assoc)
	}

	load(paste0(prefix, "_chrX.RData"))
	assoc_full <- rbind(assoc_full, assoc)

assoc_comb <<- assoc_full
}



get_peaks <- function(Assoc, p_min = 5e-4, top = NULL, window = 15000){
	Assoc_highP <- Assoc[Assoc$Score.pval < p_min,]

	close_SNPs <- Assoc[1,] #Initialize dataframe to hold close-together SNPs

	SNP_peaks <- Assoc_highP[0,] #Initialize dataframe to hold 1 peak from each close_SNPs

	i <- 1

	for(i in 1:(nrow(Assoc_highP)-1)){
		if(Assoc_highP$pos[i+1] < Assoc_highP$pos[i] + window){ #if SNP(i) (already in close SNPs) is within 15kb of the next SNP (i+1)...
			close_SNPs <- rbind(close_SNPs, Assoc_highP[i+1,]) #add SNP(i+1) to close_SNPs dataframe
		} else { #If the next SNP is not close to the previous...
			minP <- close_SNPs[close_SNPs$Score.pval == min(close_SNPs$Score.pval),] #Find the minimum P-value from the previous set of close SNPs (the peak)
			if(nrow(minP) > 1){ #If there's more than 1 minimum (happens in p = 0)...
				minP <- minP[ceiling(nrow(minP)/2),] #Take the middlemost value
			}
			SNP_peaks <- rbind(SNP_peaks, minP)
			close_SNPs <- Assoc_highP[i+1,] #Re-initialize dataframe to hold close-together SNPs
		}
		i <- i+1
	}

	if(!is.null(top) & nrow(SNP_peaks) > top){
		SNP_peaks <- SNP_peaks[order(SNP_peaks$Score.pval),][1:top,] #Take just top 100 (by lowest p val) SNP peaks
	} else {
      	SNP_peaks <- SNP_peaks[order(SNP_peaks$Score.pval),]
	}

Assoc_peaks <<- SNP_peaks
saveRDS(Assoc_peaks, file = arv$out_file)
}
	

get_peaks(assoc_comb, p_min = argv$p_min, top = argv$top,window = argv$window)

