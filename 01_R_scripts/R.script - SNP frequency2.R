#*******************************************************************************
# An R.script for generating SNP frequencies from sequence alignments
# It takes an alignment from CLC exported in Phylip format
# change the input file at the load data step to match your input file at line 23
#*******************************************************************************
#clear environment
rm(list = ls())

#-------------------------------------------------------------------------------
#load applicable functions and packages
suppressMessages(library(data.table)) # reading in data using fread
suppressMessages(library(tidyverse))  # data aggregation

#-------------------------------------------------------------------------------
# function for checking whether all bases are identical per locus
#-------------------------------------------------------------------------------
# load functions
source("02_R_functions/function_R - identical_bases_per_locus.R")   #checking for polymorphic loci
source("02_R_functions/function_R - snpFreq_per_locus.R")           #snp freq per column/locus

#*******************************************************************************
# load ligase data 
#*******************************************************************************
# skip line 1 as well as empty lines and select colums 1 and 2
fastaIN <- fread("03_data_input/example_alignment.phy", skip =1, header = F, blank.lines.skip = T) 

# split sequence column into multiple columns
fastaIN_seqs <- fastaIN[,tstrsplit(V2, "")]  #require data.table

# select sequence names
fastaIN_names <- fastaIN[,.(V1)]

# bind the two datasets and remame column 1 as it matches column 2
fastaIN <- cbind(fastaIN_names, fastaIN_seqs); names(fastaIN)[1] <- c("V0")

# convert to data.frame
fastaIN <- as.data.frame(fastaIN)

#-------------------------------------------------------------------------------
# replace all gaps with NA to exclude them from the analysis 
fastaIN[fastaIN == "-"] <- NA

# create new dataframe with just the sequence IDs, to be used in the subsequent loop
fastaIN_v1 <- as.data.frame(fastaIN[1])

# keep only columns that represent segregating sites
for(colN in 2:dim(fastaIN)[2]) {
  if(dim(table(fastaIN[colN])) != 1 & dim(table(fastaIN[colN])) != 0) { # !=1 means identical bases, != 0 means loci with NA values
    fastaIN_v1 <- as.data.frame(cbind(fastaIN_v1, fastaIN[colN]))
  }
}

# create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(
                           n = as.numeric(),
                           freq = as.numeric(),
                           base = as.character(),
                           locus = as.character()
                           )

# compute sample size and Func_SnpFrequency for each nucleotide per locus
#-------------------------------------------------------------------------------
for(i in 2:dim(fastaIN_v1)[2]) {
  freqStats_v2 <- Func_SnpFreq(fastaIN_v1[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>%
    mutate(
           base = row.names(freqStats_v2),
           locus = colnames(fastaIN_v1[i]))
  
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
}

# remove V's from locus and filter all locus with 100% freq
#-------------------------------------------------------------------------------
snpFreq <- freqStats_v1 %>%
  mutate(
         locus = gsub("V", "", locus),   #replace "V" with blanks
         locus = as.numeric(locus)) %>%  #convert locust to numeric
  mutate(freq = round(freq, 1)) #round off frequency

# remove temp files
#-------------------------------------------------------------------------------
rm(fastaIN, fastaIN_names, fastaIN_seqs, fastaIN_v1, freqStats_v1, freqStats_v2,
   colN, i, all_identical, Func_SnpFreq)

#*******************************************************************************
# summarise polymorphic data per locus, showing 3D7 ref and non-ref alleles, as well as frequencies
#*******************************************************************************
snpFreq <- snpFreq %>% mutate(locus = as.numeric(locus)) %>%         #convert locus from chr to numeric
  mutate(perc_n = paste(freq, " [",n,"]", sep = "")) %>%             #merge freq and n
  group_by(locus) %>% arrange(locus, desc(freq)) %>%                 #group by locus the sort by locus and freq (descending)
  unite(base_freq, c(base, perc_n), sep = " , ", remove = F) %>%     #merge base and perc_n        
  spread(key = base, value = base_freq) %>%                          #convert from long to wide
  group_by(locus) %>% summarise_all(list(~first(.[!is.na(.)]))) %>%  #group by locus then collapse rows where there are NA's
  gather(NonRef, Freqs, 6:dim(.)[2]) %>%                             #convert from wide to long
  filter(!is.na(Freqs)) %>%                                          #remove rows that have NA's for Freqs variable
  separate("A", c("Wildtype", "Wildtype_Freq"), sep = " , ") %>%     #seperate A into two variables
  separate(Freqs, c("Mutant", "Mutant_Freq"), sep = " , ") %>%       #seperate Freqs into two variables
  select(locus, Wildtype, Wildtype_Freq, Mutant, Mutant_Freq) %>%    #select relevant columns
  arrange(locus)  

# write data
fwrite(snpFreq, "04_data_output/example_snpFreq.csv")
