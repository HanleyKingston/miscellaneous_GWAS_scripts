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


get_p <- function(hit_from, p_from, thresh){
    p.list <- c()
    hit_from_thresh <- hit_from[hit_from$Score.pval < thresh, ]
    if(nrow(hit_from_thresh)>0){
        for(row in 1:nrow(hit_from_thresh)){
            hit_from_row <- hit_from_thresh[row,]
            p.val <- p_from[p_from$chr == hit_from_row$chr & p_from$pos == hit_from_row$pos, "Score.pval", drop = TRUE]
            p.list <- append(p.list, p.val)
            }
        } else {
        p.list <- NA
        }
    paste(p.list, collapse = "\n")
}
      
get_p(E2_noPCs, E4_noPCs_P, 5e-6)



hits.list <- list("E2_noPCs" = E2_noPCs, "E2_4PCs" = E2_4PCs, "E4_noPCs" =  E4_noPCs, "E4_4PCs" = E4_4PCs, "AD" = AD, "ADxAPOE" = ADxAPOE, "AAO" = AAO, "AAOxAPOE" = AAOxAPOE)

p_lookup.list <- list("E2_noPCs" = E2_noPCs_P, "E2_4PCs" = E2_4PCs_P, "E4_noPCs" =  E4_noPCs_P, "E4_4PCs" = E4_4PCs_P, "AD" = AD_P, "ADxAPOE" = ADxAPOE_P, "AAO" = AAO_P, "AAOxAPOE" = AAOxAPOE_P)

i <- 1
for(hits.table in hits.list){
    j <- 1
    for(p.table in p_lookup.list){
        print(data.frame(hit_from = names(hits.list)[i],
                   lookup = names(p_lookup.list)[j],
                   p_in_lookup = get_p(hits.table, p.table, 5e-6)))
        j <- j + 1
    }
    i <- i + 1
 }
