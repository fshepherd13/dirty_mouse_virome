## Analysis of metatranscriptomic RNAseq data for pathogen discovery

Snakemake implementation of a pipeline for quantifying viral species and families in dirty mouse RNAseq data. 

### Input
Raw reads from metagenomic sequencing. In our system, the experimental set up has one "dirty" pet store mouse cohoused with one to several clean laboratory mice in a single cage. To identify pathogens circulating in a cage, tissues from animals are harvested, RNA is extracted, and RNA sequencing is performed after polyA selection or rRNA depletion.

### Output
The two main outputs are csv files that show estimated counts for the taxonomic (1) families and (2) species present in each mouse.

### Pipeline overview
The pipeline works by first mapping RNAseq reads to a reference genome (shown here for mouse). Unmapped reads are then combined within a cage to perform a de novo assembly. Contigs < 500 bp are filtered from further analysis to improve accuracy of downstream counts. Taxonomic lineages of the resulting contigs are identified using BLASTn. The contigs are also used to create a reference index that the unmapped reads are then pseudomapped to the reference index with Salmon in order to identify the types/relative quantities of pathogens in the cage. Transcripts are quantified at the species and family level using Salmon and some data wrangling with help of tximport. The final output is the counts per family or species as estimated by tximport. These counts can be used for specific hypothesis testing (such as differential expression) with the help of DESeq or another package, but that is not covered here. 

Because transcripts should be de novo assembled by group, the pipeline is built to be run separately for each transmission group (i.e. one analysis for a pet store mouse and the SPF mice to which the virome is transmitted via cohousing or other method) and any paired SPF animals. I refer to these groups as "cages". To run the pipeline, create a config file that points to a specific sample csv file for the correct group of mice and point snakemake to it via the `--configfile` option.

To run this pipeline, you will need:
1. Snakemake installed on your system
2. An indexed reference genome (used to filter out host reads)
3. A formatted BLAST database


#### Dependencies/versions used
* `STAR` RNAseq read mapper, version 2.7.1a- details [here](https://github.com/alexdobin/STAR)
* `trinity` RNAseq de novo read assembler, version 2.12.0- details [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* `Trimmomatic` read trimmer, version 0.39 - details [here](http://www.usadellab.org/cms/?page=trimmomatic)
* `salmon` RNAseq transcript quantification, version 1.4.0 details [here](https://github.com/COMBINE-lab/salmon) 
* `blast` version 2.13.0, details [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) 
* `snakemake`, latest version used is 6.12.13, details [here](https://snakemake.readthedocs.io/en/stable/)
* `dib-lab/2018-ncbi-lineages` scripts on [github](https://github.com/dib-lab/2018-ncbi-lineages)
* `pyfasta` version 0.5.2, script for splitting fasta files into subsets for a parallelized BLAST, details [here](https://github.com/brentp/pyfasta)
* `tximport`-  details for the tximport command found [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta)

#### Loading dependencies with conda
Most of the dependencies are loaded within the appropriate snakemake rule using an isolated conda environment. Conda environments are found in the `workflow/envs` directory.

#### Remaining dependencies
The `make-lineages-csv.py` script from[ The Lab for Data Intensive Biology](https://github.com/dib-lab)'s `2018-ncbi-lineages` repo has already been downloaded from their github and placed in the `workflow/scripts/` directory. The snakemake rule for assigning lineages based on taxid's found by BLASTn will access the script there. 

#### Indexing reference genome
A STAR index for the reference genome is created automatically via parameters specified in the `config.yaml` file. The fasta and gtf file for the reference is also automatically downloaded via a snakemake wrapper function. 

### Running the pipeline
#### Step 1. Index the reference genome
A STAR index for the reference genome must be created by the user prior to running the pipeline. There is an example script for how to do this at `compile_genome.sh`. It is recommended to save the reference genome in a place where others in your computing group can access it as well (by setting an appropriate `DIR` variable in the script). 

#### Step 2. Adjust the config file to reflect the experimental design.
Adjust the `config.yaml` file to reflect your experimental design. The things you need before running this pipeline are:
* Directory where the raw read files are saved
* Directory where the indexed genome is saved (wherever you saved it in step 1)
* A metadata file that lists sample information (see below)
* An experiment ID
* Locations of the BLAST database you want to use


In the metadata file (the path to which should be listed under the `sample_file` parameter in the `config` file), the only column you definitely need is the 'Sample' column, which should point to the sample IDs given to you by your sequencing core. This is the text up until the lane identification number. For example, if your returned paired raw read files are sample_1_id_L001_R1.fastq.gz sample_1_id_L001_R2.fastq.gz, you should list sample_1_id in your sample.csv file. The pipeline automatically looks for the R1 and R2 files at the directory defined in the `config.yaml` file.

#### Step 3. Install snakemake
Use mamba to install snakemake. At UMN MSI, this is done with:
```
module load mamba
mamba env create -n snakemake -c bioconda snakemake
conda activate snakemake
```

#### Step 4. Run the pipeline
eg:
```
snakemake --cores 50 --use-conda
```

An example slurm script for submitting this to MSI is found at `workflow/snakemake.srun`.

*Notes on resources* This pipeline requires a lot of computing resources to run correctly, particularly Trinity de novo assembly. I run this on a dedicated node at MSI and ask for at least 400 GB and 50 cores. It often takes several days to run as well depending on the complexity of the sequencing reads. 