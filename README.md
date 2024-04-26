# influenza_a_serotype
 Assign reads from a short-read sequencing library to influenza A serotypes

# Installation

1) Clone this github repo

2) Create `conda` environment using .yaml file, e.g.

`mamba env create -f influenza_a_serotype/environment/iav_serotype.yaml`

3) activate environment

`conda activate iav_serotype`

4) Get databse files from MJT (~115 MB total). Put them in the same directory.

`Influenza_A_segment_info1.tsv`

`Influenza_A_segment_sequences.fna`

5) (optional) set database as conda environmental variable

`conda env config vars set IAVS_DB=/path/to/iav_DB`

# Run `iav_serotype`

Right now, requirement is 1 set of paired-end reads per run. These can be compressed in .bz2 format.

example:

`iav_serotype -r my_reads/virome.R1.fastq my_reads/virome.R2.fastq -s my_virome_iav -o iav_project --db iav_DB`

