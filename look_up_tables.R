#For extracting P-values for hits
#Expects files to have columns Score.pval, pos, and chr
library(dplyr)

#Load files to extract hits from (this should already be fitlered to lead SNPs, not a full association test):
E2_noPCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/SNP_peaks.txt")
E2_4PCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_4PCs/SNP_peaks.txt")
E4_noPCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_noPCs/SNP_peaks.txt")
E4_4PCs <- read.table("/acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_4PCs/SNP_peaks.txt")
AD <- read.table("/acct/hkings/NIAGADS/analyses/AD/AD_4PCs/SNP_peaks.txt")
ADxAPOE <- read.table("/acct/hkings/NIAGADS/analyses/AD/AD_xE4E2_4PCs/SNP_peaks.txt")
AAO <- read.table("/acct/hkings/NIAGADS/analyses/AAO/AAO_4PCs/SNP_peaks.txt")
AAOxAPOE <- read.table("/acct/hkings/NIAGADS/analyses/AAO/AAO_xE4E2_4PCs/SNP_peaks.txt")


#Load files we want to look up P values from:
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
    #In itialize empty lists to host chr and pos from hit_from and p-valeus from p_from
    chrpos.list <- c()
    p.list <- c()
    hit_from_thresh <- hit_from[hit_from$Score.pval < thresh, ]
    #exclude APOE region:
    hit_from_thresh <- hit_from_thresh[!(hit_from_thresh$chr == 19 & hit_from_thresh$pos > 49500000 & hit_from_thresh$pos < 50750000),]
    if(nrow(hit_from_thresh)>0){
        for(row in 1:nrow(hit_from_thresh)){
            hit_from_row <- hit_from_thresh[row,]
            #Record the chr and pos of each hit in hit_from:
            hit_chrpos <- paste0(hit_from_row$chr, ": ", hit_from_row$pos)
            chrpos.list <- append(chrpos.list, hit_chrpos)
            #Record p-values associated with this his in p_from:
            p.val <- signif(p_from[p_from$chr == hit_from_row$chr & p_from$pos == hit_from_row$pos, "Score.pval", drop = TRUE], 2)
            p.list <- append(p.list, p.val)
            }
        } else {
        p.list <- NA
        }
    #Collapse both vectors into single, line seperated strings:
    chrpos.list <- paste(chrpos.list, collapse = "\n")
    chrpos.list_save <<- chrpos.list
    paste(p.list, collapse = "\n")
}
      
get_p(E2_noPCs, E2_noPCs, 1.5e-7)

#hit_from_tables.list is a named list of tables to extract chromosomal positions from based on a p-value threshold
#p_from_tables.list is a named list of tables (usually association results) to extract p values from based on matches to hits in hits tables
#thresh is the p-value threshold below which to lookup up hits (gets passed to get_p fuction)
make_hits_table <- function(hit_from_tables.list, p_from_tables.list, thresh){
    print(paste0("p-value-threshold: ", thresh))
    i <- 1
    lookup.table <- data.frame(matrix(ncol = length(hit_from_tables.list), nrow = length(p_from_tables.list)+1))
    rownames(lookup.table) <- c("hits", names(p_from_tables.list))
    colnames(lookup.table) <- names(hit_from_tables.list)

    for(hits.table in hit_from_tables.list){
        p_values.list <- c()
        j <- 2
        for(p.table in p_from_tables.list){
            lookup.table[j,i] <- get_p(hits.table, p.table, 3e-6)
            j <- j + 1
        }
        lookup.table[1,i] <- chrpos.list_save
        i <- i + 1
    }
    return(lookup.table)
}
         
         

hits.list <- list("E2_noPCs" = E2_noPCs, "E2_4PCs" = E2_4PCs, "E4_noPCs" =  E4_noPCs, "E4_4PCs" = E4_4PCs, "AD" = AD, "ADxAPOE" = ADxAPOE, "AAO" = AAO, "AAOxAPOE" = AAOxAPOE)
p_lookup.list <- list("E2_noPCs" = E2_noPCs_P, "E2_4PCs" = E2_4PCs_P, "E4_noPCs" =  E4_noPCs_P, "E4_4PCs" = E4_4PCs_P, "AD" = AD_P, "ADxAPOE" = ADxAPOE_P, "AAO" = AAO_P, "AAOxAPOE" = AAOxAPOE_P)


lookup.table_sug_sig <- make_hits_table(hits.list, p_lookup.list, 3e-6)
write.csv(lookup.table_sug_sig, "/acct/hkings/NIAGADS/analyses/lookup_table_NIAGADS_APOE_AD_AAO_sug_sig.csv")

lookup.table_sig <- make_hits_table(hits.list, p_lookup.list, 1.5e-7)
write.csv(lookup.table_sig, "/acct/hkings/NIAGADS/analyses/lookup_table_NIAGADS_APOE_AD_AAO_sig.csv")
#Note: because there are multiple lines per cell, this isn't easy to read in a text editor. Open in excel


