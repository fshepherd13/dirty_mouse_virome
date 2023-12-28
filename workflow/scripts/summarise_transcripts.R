library('tximport')

################Read in relevant dataframes################
#Read in dataframe that contains the trinity ID's linked to lineage information
df <- read.csv(snakemake@input[["lineages"]], header=TRUE)
#Read in files of from salmon quantification
files <- snakemake@input[["salmon_files"]]
#Link sample ID info to the name of each salmon quant file
names(files) <- sapply(strsplit(dirname(files), "/"), `[`, 4)

#Perform data analysis for accession level lineage assignments, first extracting tximport estimated counts
#Link salmon quant files to the transcript lineage info
txi_sp <- tximport(files, type = "salmon", tx2gen = df[c("trinity_id", "species")], ignoreTxVersion = FALSE)
#Create dataframe out of species counts
sp_counts <- data.frame("species"=row.names(txi_sp$counts), "counts"=txi_sp$counts)
#Edit row names to be sequential numbering
rownames(sp_counts) <- seq_len(nrow(sp_counts))
#Extract the unique combinations of taxonomy lineages
tax <- unique(df[,c("superkingdom", "phylum", "class", "order", "family", "genus", "species")])
#Merge the full lineage information to the species designation
final_sp <- merge(sp_counts, tax, by="species")
write.csv(final_sp, file=file.path(snakemake@output[["species_counts"]]), row.names = FALSE)

#Do the same but output abundances
#Create dataframe out of species abundances (TPM)
sp_tpm <- data.frame("species"=row.names(txi_sp$abundance), "tpm"=txi_sp$abundance)
#Edit row names to be sequential numbering
rownames(sp_tpm) <- seq_len(nrow(sp_tpm))
#Extract the unique combinations of taxonomy lineages
tax <- unique(df[,c("superkingdom", "phylum", "class", "order", "family", "genus", "species")])
#Merge the full lineage information to the species designation
final_sp_tpm <- merge(sp_tpm, tax, by="species")
write.csv(final_sp_tpm, file=file.path(snakemake@output[["species_tpm"]]), row.names = FALSE)


#Perform the same analyses for family level lineage assignments
#First, estimated counts
txi_fam <- tximport(files, type = "salmon", tx2gen = df[c("trinity_id", "family")], ignoreTxVersion = FALSE)
fam_counts <- data.frame("family"=row.names(txi_fam$counts), "counts"=txi_fam$counts)
rownames(fam_counts) <- seq_len(nrow(fam_counts))
tax <- unique(tax[,c("superkingdom", "phylum", "class", "order", "family")])
final_fam <- merge(fam_counts, tax, by="family")
write.csv(final_fam, file=file.path(snakemake@output[["family_counts"]]), row.names=FALSE)

#Next, abundances (TPM)
fam_tpm <- data.frame("family"=row.names(txi_fam$abundance), "counts"=txi_fam$abundance)
rownames(fam_tpm) <- seq_len(nrow(fam_tpm))
tax <- unique(tax[,c("superkingdom", "phylum", "class", "order", "family")])
final_fam_tpm <- merge(fam_tpm, tax, by="family")
write.csv(final_fam_tpm, file=file.path(snakemake@output[["family_tpm"]]), row.names=FALSE)
