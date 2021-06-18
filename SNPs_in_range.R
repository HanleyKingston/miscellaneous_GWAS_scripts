SNPs_in_range <- function(SNP_table, ranges_table){
	SNPs_in_range.df <- cbind(SNP_table[0,c("variant.id", "chr", "pos", "Score.pval")], ranges_table[0,])
	for(i in 1:nrow(ranges_table)){
		SNP_range <- SNP_table[SNP_table$chr == ranges_table$chromosome[i] & SNP_table$pos >= ranges_table$start_pos[i] & SNP_table$pos <= ranges_table$end_pos[i], ]
		min_SNP_in_range <- SNP_range[SNP_range$Score.pval == min(SNP_range$Score.pval),]
		min_SNP_in_range2 <- cbind(min_SNP_in_range[,c("variant.id", "chr", "pos", "Score.pval")], ranges_table[i,])
		SNPs_in_range.df <- rbind(SNPs_in_range.df, min_SNP_in_range2)
	}
	return(SNPs_in_range.df)
}
