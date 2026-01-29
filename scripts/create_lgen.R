
## Originally written by Ryan Schubert for the Wheeler Lab, github@RyanSchu Lab github@wheelerlab script to convert illumina final report 
## gentoype data to lgen data. This script will generate the .lgen file and .map file used for plink format. For instructions on how to 
## generate the .fam file please see the .md file
# Adapted by Nathan Ribeiro, UMCG Immunogenetics
library(dplyr)
library(tidyr)
library(argparse)

# Get the input and output settings ----
args <- commandArgs(trailingOnly = TRUE)
illumina_file  <- args[1]
sample_info <- args[2]
lgen_out <- args[3]
map_out <- args[4]
fam_out <- args[5]

# Conversion -----
Finalreport<-as.data.frame(read.table(file=illumina_file, sep='\t', skip = 9, header = T))

Finalreport <- filter(Finalreport, Sample.ID %in% samples) %>%
        filter(Allele1...Top != 'I')
Finalreport["empty"]="0"

# Create mapping file with SNP positions
map<-select(Finalreport, Chr, SNP.Name, empty, Position)
map<-map[!duplicated(map),]
map<-map[complete.cases(map),]

# Create lgen file
lgen<-select(Finalreport, Sample.ID, SNP.Name, Allele1...Top, Allele2...Top)
lgen<-lgen[!duplicated(lgen),]
lgen<-lgen[complete.cases(lgen),]
lgen<-filter(lgen, Allele1...Top != "-" & Allele2...Top != "-")

# Add a dummy FID variable unique to every sample.id
lgen$FID <- as.integer(factor(lgen$Sample.ID))

# Just reordering columns
lgen <- select(lgen, FID, Sample.ID, SNP.Name, Allele1...Top, Allele2...Top)

# Create fam file
fam <- data.frame(FID = lgen$FID, IID = lgen$Sample.ID, Paternal_ID = 0,
                  Maternal_ID = 0)
fam <- fam %>% distinct(IID, .keep_all = TRUE)

# Add sex and phenotype information
fam_info <- read.csv(sample_info, header = TRUE)
fam <- left_join(fam, fam_info)

# Write outputs ----
write.table(map, file = map_out, sep = "\t", col.names = F, row.names = F, quote = F)
write.table(lgen, file = lgen_out, sep = "\t", col.names = F, row.names = F, quote = F)
write.table(fam, file = fam_out, sep = "\t", col.names = F, row.names = F, quote = F)