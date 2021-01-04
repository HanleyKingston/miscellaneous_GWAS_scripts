rand <- SNPs[sample(1:nrow(SNPs),1), "variant.id"]
SNPs[SNPs$variant.id == rand, "pos"] == variants_gds[variants_gds$variant.id == rand, "pos"] & SNPs[SNPs$variant.id == rand, "chrom"] == variants_gds[variants_gds$variant.id == rand, "chr"] 
