#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
library(dplyr)


# read arguments
argp <- arg_parser("Get top peaks") %>%
  add_argument("--assoc_file", help = "name of file with association results for all chromosomes") %>%
  add_argument("--p_min", help = "maximum p-value below which to search for peaks",
               default = 5e-4) %>%
  add_argument("--window", help = "bp distance within which SNPs should be grouped into a peak", default = 100000) %>%
  add_argument("--out_prefix", help = "file to save SNP peaks", default = "SNP_peaks")
#  add_argument("--p_name", help = "name of column for p-value", default = "Score.pval") %>%
#  add_argument("--pos_name", help = "name of column for position", default = "pos") %>%


argv <- parse_args(argp)



#If association results are in seperate files by chromosome, combine them
assoc_comb <- readRDS(argv$assoc_file)


get_peaks <- function(Assoc, p_min = 5e-4, top = NULL, window = 15000){
	Assoc_highP <- Assoc[Assoc$Score.pval < p_min,]

	close_SNPs <- Assoc[1,] #Initialize dataframe to hold close-together SNPs

	SNP_peaks <- Assoc_highP[0,] #Initialize dataframe to hold 1 peak from each close_SNPs

	i <- 1

	for(i in 1:(nrow(Assoc_highP)-1)){
		if(Assoc_highP$pos[i+1] < Assoc_highP$pos[i] + window){ #if SNP(i) (already in close SNPs) is within --window of the next SNP (i+1)...
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

   SNP_peaks <- SNP_peaks[order(SNP_peaks$Score.pval),]


Assoc_peaks <<- SNP_peaks %>% 
 mutate_if(is.numeric, function(x) floor(x) + signif(x - floor(x), 3))

saveRDS(Assoc_peaks, file = paste0(argv$out_prefix, ".rds"))
write.table(Assoc_peaks, file = paste0(argv$out_prefix, ".txt"))
}
	

get_peaks(assoc_comb, p_min = argv$p_min, top = argv$top,window = argv$window)




