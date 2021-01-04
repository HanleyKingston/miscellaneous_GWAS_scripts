library(argparser)
library(tidyverse)
library(glue)
library(SeqVarTools)
library(SNPRelate)
library(dplyr)
library(tibble)

argp <- arg_parser("Correlation of variants with PCs")
argp <- add_argument(argp, "--pca_file", help="path to PCA file on disk")
argp <- add_argument(argp, "--gds_file", help="path to GDS file on disk")
argp <- add_argument(argp, "--outfile", help = "output file")
argp <- add_argument(argp, "--npcs", default = 10, help = "number of PCs to compute correlations for")
argp <- add_argument(argp, "--block-size", default = 10000, help = "genotype block size for computing correlations")
argv <- parse_args(argp)
print(argv)

gds_file <- argv$gds_file
pca_file <- argv$pca_file
n_pcs = argv$npcs
block_size = argv$block_size
outfile = argv$outfile

# Read in PCair and subset to unrelated only.
pcair <- readRDS(pca_file)
pcs_all <- pcair$vectors
pcs_unrel <- pcs_all[rownames(pcs_all) %in% pcair$unrels, ]
dim(pcs_unrel)

gds <- seqOpen(gds_file)

seqSetFilter(gds, sample.id = rownames(pcs_unrel))


# Set up an iterator to iterate over all variants.
iterator <- SeqVarBlockIterator(gds, variantBlock=block_size)
n.iter <- length(variantFilter(iterator))
i <- 1
iterate <- TRUE
corr_list <- list()
while(iterate) {

  # Calculate correlations between genotype and PCs for this block.
  geno <- refDosage(iterator)
  variants <- data.frame(
    variant.id = seqGetData(iterator, "variant.id"),
    position = seqGetData(iterator, "position"),
    chrom = seqGetData(iterator, "chromosome"),
    stringsAsFactors = FALSE
  )

  # Quick check, add code to fix if it fails.
  stopifnot(all.equal(rownames(geno), rownames(pcs_unrel)))

  corr <- cor(geno, pcs_unrel[, 1:n_pcs], use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    setNames(sprintf("PC", 1:n_pcs)) %>%
    rownames_to_column(var = "variant.id") %>%
    mutate(variant.id = as.integer(variant.id))
  corr_list[[i]] <- variants %>% left_join(corr, by = "variant.id")

  message(sprintf("Finished iteration %s/%s", i, n.iter))

  # Iterate.
  i <- i + 1
  iterate <- iterateFilter(iterator)

}
seqClose(gds)

corr <- bind_rows(corr_list) %>%
  select(variant.id, chrom, everything())

corr <- corr[!is.na(corr$PC),]
#Need to fix rownames!

saveRDS(corr, file = outfile)



