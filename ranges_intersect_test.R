#This will look for overlap in the peaks and/or ranges of 2 datasets (ex. 2 association tests or an association test and a ranges table and will return the lead SNP from the first dataset. If there is no lead SNP, it will return the range from that table
#Need a QC step to make sure start pos comes before end_pos (otherwise everything will overlap!

library(argparser)
library(magrittr)
library(ggplot2)


argp <- arg_parser("intersect of ranges table") %>%
  add_argument("ranges_table_a", help = "table with ranges, usually based on assoc results") %>%
  add_argument("ranges_table_b", help = "table with ranges, usually based on assoc results") %>%
  add_argument("--chr_a", help = "name of chromosome varaible in table a", default = "chr") %>%
  add_argument("--start_a", help = "name of start position varaible in table a", default = "start_pos") %>%
  add_argument("--end_a", help = "name of end position variable in table a", default = "end_pos") %>%
  add_argument("--lead_a", help = "name of lead SNP variable in table a", default = "variant.id") %>% #In future, will set up so it doesn't need this and can, instead, save the range
  add_argument("--chr_b", help = "name of chromosome varaible in table a", default = "chr") %>%
  add_argument("--start_b", help = "name of start position variable in table b", default = "start_pos") %>%
  add_argument("--end_b", help = "name of start position varaible in table b", default = "end_pos") %>%
  add_argument("--out_prefix", help = "file to save overlap to", default = "")

argv <- parse_args(argp)

chr_a <- argv$chr_a
start_a <- argv$start_a
end_a <- argv$end_a
lead_a <- argv$lead_a
start_b <- argv$start_b
end_b <- argv$end_b
chr_b <- argv$chr_b


table_a <- read.table(argv$ranges_table_a, header = TRUE)
table_b <- read.table(argv$ranges_table_b, header = TRUE)


shared_SNP <- data.frame(matrix(NA, ncol = ncol(table_a), nrow = nrow(table_a)))
colnames(shared_SNP) <-  colnames(table_a)

#for(i in 1:nrow(table_a)){
#	if(table_a[i, start_a] > table_a[i, end_a]) {
#		print("start position in table a is after end position")
#		break
#}

i <- 1
j <- 1

for(i in 1:nrow(table_a)){
	for(j in 1:nrow(table_b)){
		if(table_a[i, chr_a] == table_b[j, chr_b]) { #If chrosomomes match
			#Check every possible way a and b could overlap
			if (table_a[i, end_a] >= table_b[j, start_b] & table_a[i, end_a] <= table_b[j, end_b]) {
				shared_SNP[i,] <- table_a[i,]
			} else if (table_a[i, start_a] >= table_b[j, start_b] & table_a[i, start_a] <= table_b[j, end_b]) {
				shared_SNP[i,] <- table_a[i,]
			} else if (table_a[i, start_a] < table_b[j, start_b]  & table_a[i, end_a] > table_b[j, end_b]) {
				shared_SNP[i,] <- table_a[i,]
			} else if (table_a[i, start_a] > table_b[j, start_b] & table_a[i, end_a] < table_b[j, end_b]) {
				shared_SNP[i,] <- table_a[i,]
			}
		}
		j <- j + 1
	}
	j <- 1
	i <- i + 1
}


write.table(shared_SNP, file = paste0(argv$out_prefix, "_intersect_ranges.txt"))