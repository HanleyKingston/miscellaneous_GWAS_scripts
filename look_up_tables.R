#For extracting P-values for hits
library(dplyr)

#Hits:
E2_noPCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/SNP_peaks.txt")
E2_4PCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_4PCs/SNP_peaks.txt")
E4_noPCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_noPCs/SNP_peaks.txt")
E4_4PCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_4PCs/SNP_peaks.txt")
AD <- read.table("/acct/hkings/NIAGADS/analyses/AD/AD_4PCs/SNP_peaks.txt")
ADxAPOE <- read.table("/acct/hkings/NIAGADS/analyses/AD/AD_xE4E2_4PCs/SNP_peaks.txt")
AAO <- read.table("/acct/hkings/NIAGADS/analyses/AAO/AAO_4PCs/SNP_peaks.txt")
AAOxAPOE <- read.table("/acct/hkings/NIAGADS/analyses/AAO/AAO_xE4E2_4PCs/SNP_peaks.txt")


#P_values:
E2_noPCs_P <- readRDS("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/assoc.rds")
E2_4PCs_P <- readRDS("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_4PCs/assoc.rds")
E4_noPCs_P <- readRDS("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_noPCs/assoc.rds")
E4_4PCs_P <- readRDS("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_4PCs/assoc.rds")
AD_P <- readRDS("/acct/hkings/NIAGADS/analyses/AD/AD_4PCs/assoc.rds")
ADxAPOE_P <- readRDS("/acct/hkings/NIAGADS/analyses/AD/AD_xE4E2_4PCs/assoc.rds")
AAO_P <- readRDS("/acct/hkings/NIAGADS/analyses/AAO/AAO_4PCs/assoc.rds")
AAOxAPOE_P <- readRDS("/acct/hkings/NIAGADS/analyses/AAO/AAO_xE4E2_4PCs/assoc.rds")


APOE_region <- E4_noPCs_P[E4_noPCs_P$chr == 19 & E4_noPCs_P$pos > 49500000 & E4_noPCs_P$pos < 50750000,] #This is defined very broadly
plot(-log10(APOE_region$Score.pval)~APOE_region$pos)

get_p <- function(hit_from, p_from, thresh){
    p.list <- c()
    hit_from_thresh <- hit_from[hit_from$Score.pval < thresh, ]
    #exclude APOE region:
    hit_from_thresh <- hit_from_thresh[!(hit_from_thresh$chr == 19 & hit_from_thresh$pos > 49500000 & hit_from_thresh$pos < 50750000),]
    if(nrow(hit_from_thresh)>0){
        for(row in 1:nrow(hit_from_thresh)){
            hit_from_row <- hit_from_thresh[row,]
            p.val <- signif(p_from[p_from$chr == hit_from_row$chr & p_from$pos == hit_from_row$pos, "Score.pval", drop = TRUE], 2)
            p.list <- append(p.list, p.val)
            }
        } else {
        p.list <- NA
        }
    paste(p.list, collapse = "\n")
}
      
get_p(E2_noPCs, E2_noPCs, 1.5e-7)


hits.list <- list("E2_noPCs" = E2_noPCs, "E2_4PCs" = E2_4PCs, "E4_noPCs" =  E4_noPCs, "E4_4PCs" = E4_4PCs, "AD" = AD, "ADxAPOE" = ADxAPOE, "AAO" = AAO, "AAOxAPOE" = AAOxAPOE)
p_lookup.list <- list("E2_noPCs" = E2_noPCs_P, "E2_4PCs" = E2_4PCs_P, "E4_noPCs" =  E4_noPCs_P, "E4_4PCs" = E4_4PCs_P, "AD" = AD_P, "ADxAPOE" = ADxAPOE_P, "AAO" = AAO_P, "AAOxAPOE" = AAOxAPOE_P)

i <- 1
lookup.table <- data.frame(matrix(ncol = length(hits.list), nrow = length(p_lookup.list)))
rownames(lookup.table) <- names(p_lookup.list)
colnames(lookup.table) <- names(hits.list)

for(hits.table in hits.list){
     p_values.list <- c()
    j <- 1
    for(p.table in p_lookup.list){
        lookup.table[j,i] <- get_p(hits.table, p.table, 3e-6)
        j <- j + 1
    }
    i <- i + 1
 }

#Note: because there are multiple lines per cell, this isn't easy to read in a text editor. Open in excel
write.csv(lookup.table, "/acct/hkings/NIAGADS/analyses/lookup_table_NIAGADS_APOE_AD_AAO_sig.csv")

write.csv(lookup.table, "/acct/hkings/NIAGADS/analyses/lookup_table_NIAGADS_APOE_AD_AAO_sug_sig.csv")

