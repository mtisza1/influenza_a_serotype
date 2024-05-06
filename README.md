# influenza_a_serotype
 Assign reads from a short-read sequencing library to influenza A serotypes with `iav_serotype` command line tool.


 The basic premise is to competitively align reads to an up-to-date database of influenza A sequences tagged with information about segment and serotype. Then, determined if a read (pair) aligns better to a particular serotype.

`Input`: paired-end short reads - OR - long reads

`Output`: 1. files of serotype-specific reads. 2. Read summary table. 3. Plot of serotype distribution in sample
 

# Installation

1) Clone this github repo

2) Create `conda` environment using .yaml file, e.g.

`mamba env create -f influenza_a_serotype/environment/iav_serotype.yaml`

3) activate environment

`conda activate iav_serotype`

4) Download and unpack database files from Zenodo (~115 MB total). `cd` to the directory you want them to live.

`wget https://zenodo.org/records/11123391/files/Influenza_A_segment_sequences.tar.gz`

`tar -xvf Influenza_A_segment_sequences.tar.gz`


You should now have these 2 files:

`DBs/v1.1/Influenza_A_segment_info1.tsv`

`DBs/v1.1/Influenza_A_segment_sequences.fna`

5) (optional) set database as conda environmental variable

`conda env config vars set IAVS_DB=/path/to/DBs/v1.1`

# Run `iav_serotype`

Right now, requirement is either 1 set of paired-end short reads per run, or 1 or more long read files. Any and all reads must be decompressed with the `.fastq` extension.

example:

short paired-end reads:

`iav_serotype -r my_reads/virome.R1.fastq my_reads/virome.R2.fastq -s my_virome_iav -o iav_project --db /path/to/DBs/v1.1`

short paired-end reads:

`iav_serotype -r long_reads/virome1.fastq long_reads/virome2.fastq long_reads/virome3.fastq -s my_lr_virome_iav -o iav_project --db /path/to/DBs/v1.1`

# Database notes

## v1.1

Added sequences and metadata rows to `v1.0`, searches on April 26, 2024:

1) FluDB query Influenza A "complete", collection date(06-01-2022 - 05-31-2023)

2) NCBI virus taxid=2955291, length filter( 3000 > 700 ), collection date(06-01-2023 - 04-26-2024)

Then, sequences and metadata were parsed to remove any duplicate accessions.

## v1.0
using NCBI datasets command line tool, (April 23, 2024), this is how I downloaded Influenza A assmeblies:

`datasets download genome taxon 2955291 --assembly-level complete --exclude-atypical --filename Influenza_A_dataset1 --include gbff`

This retrieved about about 8,000 genomes/64,000 segments. There are about 1,000,000 segments on GenBank, so I am not sure why they weren't all retrieved.

After unzipping the downloaded directory I used the following script to process these data into a .fasta and table (.tsv)

`python parse_iav_genbank/influenza_a_gbf_to_fna_and_table.py ncbi_dataset`

Note: I'm working on getting a more complete database, as I think there is relevant diversity not included here.

