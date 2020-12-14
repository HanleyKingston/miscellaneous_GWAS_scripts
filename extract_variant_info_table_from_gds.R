library(SeqArray)
library(GENESIS)
library(SeqVarTools)

gds <- seqOpen("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/LOAD_families_010413_genotyped.gds")

gds_vars <- cbind(variantInfo(gds), variant.name = seqGetData(gds, "annotation/id"))
