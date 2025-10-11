# influenza_a_serotype
 Assign reads from a short- or long-read sequencing library to influenza A serotypes with `iav_serotype` command line tool.


 The basic premise is to competitively align reads to an up-to-date database of influenza A sequences tagged with information about segment and serotype. Then, determined if a read (pair) aligns better to a particular serotype.

`Input`: paired-end short reads - OR - long reads

`Output`: 1. files of serotype-specific reads. 2. Read summary table. 3. Plot of serotype distribution in sample
 

# Installation

1) Clone this github repo

2) Create `conda` environment using .yaml file, e.g.

`mamba env create -f influenza_a_serotype/environment/iav_serotype.yaml`

3) activate environment

`conda activate iav_serotype`

4) Use `pip` to install command line tool.

`cd influenza_a_serotype`

`pip install .`

(you should now be able to call `iav_serotype` from the command line to bring up the help menu)

5) Download and unpack database files (`v1.25`) from Zenodo (~1.7 GB total). `cd` to the directory you want them to live.

`wget https://zenodo.org/records/11509609/files/Influenza_A_segment_sequences.tar.gz`

`tar -xvf Influenza_A_segment_sequences.tar.gz`


You should now have these 2 files:

`DBs/v1.25/Influenza_A_segment_info1.tsv`

`DBs/v1.25/Influenza_A_segment_sequences.fna`

6) (optional) set database as conda environmental variable

`conda env config vars set IAVS_DB=/path/to/DBs/v1.25`


# Update (if necessary)

Latest version is `v0.2.0`


1) `pip uninstall iav_serotype`

2) `cd influenza_a_serotype`

3) `git pull`

4) `pip install .`

5) `iav_serotype --version`


# Run `iav_serotype`

Right now, requirement is either 1 set of paired-end short reads per run, or 1 or more long read files. Any and all reads must be decompressed with the `.fastq` extension.

## examples:

short paired-end reads:

`iav_serotype -r my_reads/virome.R1.fastq my_reads/virome.R2.fastq -s my_virome_iav -o iav_project --db /path/to/DBs/v1.25`

long reads:

`iav_serotype -r long_reads/virome1.fastq long_reads/virome2.fastq long_reads/virome3.fastq -s my_lr_virome_iav -o iav_project --db /path/to/DBs/v1.25`

# Database notes

## v1.25

NOTE: Use database `v1.25` with `iav_serotype v0.1.3` or later. A small number of added reference sequences have short sequences that were likely assembled into the genome by mistake. These will cause assignment of non-specific reads as "ambiguous" IAV in `iav_serotype v0.1.1`, but these alignments are filtered out starting in `iav_serotype v0.1.2`. Further, in `iav_serotype v0.1.3`, the `minimap2` flag `-f 100000` is added to account for very high prevalence minimizers in reference. Thank you.

Description:

Fresh download of all complete influenza A sequences + metadata tables from NCBI virus. Consists of 981,537 complete segements.

1) Search and download date=2024-06-05, NCBI virus taxid=2955291, length filter( 3000 > 700 )

2) Then, sequences and metadata were parsed to remove any duplicate accessions and any sequences with uninformative serotype labels, e.g. "H3", "mixed", "HxNx".

3) Finally, low complexity regions were masked with `bbmask.sh` with default parameters for low-complexity filtering (available with `bbtools`)

[Zenodo](https://zenodo.org/records/11509609)

## v1.1

NOTE: it was pointed out to me that several sequences with uninformative serotype labels were included in both version of the database, e.g. "H3", "mixed", "HxNx". Since these poorly-labeled sequences only represent ~1,600/~210,000 sequences in `v1.1`, performance should not be substantially effected, but these will be removed in future releases.

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


# Future updates

1) Add Influenza B, C, and D
