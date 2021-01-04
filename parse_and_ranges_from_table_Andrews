# Run with: python parse_table.py Andrews_table1.txt
## This is specific to parsing data from table 1 in this paper: References

Andrews, S. J., Fulton-Howard, B., & Goate, A. (2020). Interpretation of risk loci from genome-wide association studies of alzheimer's disease. The Lancet. Neurology, 19(4), 326-335. doi:10.1016/S1474-4422(19)30435-1
## but can be adapted to pull GWAS assoc results from any table and put into a tsv or csv

#This code will extract the loci information and SNP position for the SNPs found in the 4 papers in teh Andrews review paper
#For each SNP, 50K is added to the highest possition and subtracted from the smallest position found in each of the 4 papers
#If the SNP was only found in 1 paper, or the same position was identified in each paper, the final range will be length 100K
#This is done to be able to identify SNPs from my analysis that are in strong LD with those in this paper, given that they (probably?) did not use the same chip

import sys
import re

file = open(sys.argv[1], "r").readlines()[2:]
myFile = open("Andrews_table1_parsed.txt", "w")

position_list = []
loci = open(sys.argv[1], "r").readlines()[1].strip() #initializes the first loci label from the second line

for line in file:
    if re.match("[XS](\d)+", line): #Finds all lines with specific SNPs at that locus from each paper (up to 4 unique positions)
        SNP_info = line.strip().split("\t")[2] #splits line on the tab into a list of elements and takes the first (chr:pos:ref:alt)
        SNP_pos = re.match("((\d)+)(\:)((\d)+)(\S)", SNP_info) #Split the string into a set of regular expressions
        if SNP_pos != None: #This is to prevent an error for lines with no matching syntax, although there shouldn't be any, so need to check why this is happening
            position_list.append(int(SNP_pos.group(4))) #Extracts just the position info (group 3) from the regular expression
    else: #The line defines the SNP locus
            min_SNP = min(position_list)-50000 #find the range to use from position info extracted above (min and max position plus 5000 bp)
            max_SNP = max(position_list)+50000
            myFile.write(f"{line.strip()}:\t{SNP_pos.group(1)}\t{min_SNP}\t{max_SNP}\n")  #assign this range of positions to the previously saved loci and chromosome and print to file
            position_list = [] #re-initialize the position list for next loci
            loci = line.strip() #assign the current line to the next loci (this order is necissary because the file has the loci before the position)
