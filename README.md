# GWAS_scripts
note: -- indicates an optional argument

## SNP_peaks.R
Finds top peaks from association test results from GENESIS saved by chromsomes in format <prefix>_chr1.RData 
Takes arguemnts: 
--prefix", ex. "<prefix>_chr1.RData" 
--p_min", maximum p-value below which to search for peaks, default = 5e-4 
--top", number of SNPs peaks to save, default = NULL 
--window", bp distance within which SNPs should be grouped into a peak", default = 15000 
--out_file", file to save SNP peaks, default = "SNP_peaks.rds" 
