#Make it so peaks that overlap are combined (before plotting)

library(argparser)
library(magrittr)
library(ggplot2)


argp <- arg_parser("Get peak ranges") %>%
  add_argument("--assoc_file", help = "name of file with association results for all chromosomes") %>%
  add_argument("--peaks_file", help = "name of file with peaks from association results in same format as association results") %>%
  add_argument("--out_prefix", help = "file to save SNP ranges", default = "peak_ranges") %>%
	add_argument("--plot", help = "if TRUE, will plot each peak. Note: this automatically creates a seperate folder for plots to save to, so must run from within workign folder and don't pass a file path to out_prefix", flag = TRUE)

argv <- parse_args(argp)


assoc <- readRDS(argv$assoc_file)

#Peaks should have columns for variant.id, chr, pos, Score.pval (with names that match assoc
peaks <- read.table(argv$peaks_file)

#This function saves the plot as obejct p
plot_the_peak <- function(Chr, Pos, color = NULL){
	p <<- ggplot(assoc[assoc$chr == Chr & assoc$pos > Pos - 1250000 & assoc$pos < Pos + 1250000,], aes(-log10(Score.pval),pos, color = color)) +
	geom_point() +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20), labels  = scales::label_number_auto()) +
	xlab("-log10(P-value)") +
	ylab(paste(Chr, "position"))
}


#Get ranges surrounding each peak

#initialize ranges table
ranges <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(ranges) <- c("variant.id", "chr", "lead_SNP_pos", "start_pos", "end_pos")

i <- 1
for(i in 1:(nrow(peaks))){
	ranges[i,1:3] <- peaks[i, c("variant.id", "chr", "pos"), drop = TRUE]
	#Get SNPs with close-ish magnitude
	p_thresh <- -log10(peaks[i, "Score.pval"])/(-log10(peaks[i, "Score.pval"])/10+2)
	assoc_range <- assoc[assoc$chr == peaks[i, "chr"] & assoc$pos > peaks[i, "pos"] - 1000000 & assoc$pos < peaks[i, "pos"] + 1000000 & assoc$Score.pval < 10^(-p_thresh),]
	ranges[i, "start_pos"] <- min(assoc_range$pos)
	ranges[i, "end_pos"] <- max(assoc_range$pos)
	i <- i + 1 
}



#Combine ranges that overlap
#not done yet



write.table(ranges, paste0(argv$out_prefix, ".txt"))


if(argv$plot == TRUE){	
	if(!(dir.exists("peaks_zoom"))){
		dir.create("peaks_zoom")
	}

	#plot each peak and color range to double check that it was correctly captured
	range_color <- c()
	for(i in 1:nrow(ranges)){
		range_color2 <- assoc[assoc$chr == ranges[i, "chr"] & assoc$pos >= 	ranges[i, "start_pos"]  & assoc$pos <= ranges[i, "end_pos"], "variant.id"]
		range_color <- append(range_color, range_color2)
		assoc$color <- ifelse(assoc$variant.id %in% range_color, "red", "black")
		i <- i + 1
	}


	for(i in 1:nrow(ranges)){
	plot_the_peak(peaks[i, "chr"], peaks[i, "pos"], color = color)
	ggsave(paste0("peaks_zoom/", argv$out_prefix, i, "_peak_zoom.png"), plot=p)
	i <- i+1
	}

}



