#!/bin/bash

# script for generating .slurm scripts aligning metagenome reads to genbank parvovirus b19 ON023021.1
# read stats for preprocessing
#

source /cmmr/prod/envParams/condanewenv.init && conda activate EsViritu

SAMP=$1

READ1=$2

READ2=$3

OUT_DIR=$4

SAMP_DIR="${OUT_DIR}/${SAMP}"

## remember to update this
SCRIPT_DIR="/gpfs2/projects/mjt_repos/influenza_a_serotype"


CPUS=40


if [ ! -d ${OUT_DIR} ] ; then
	mkdir ${OUT_DIR}
fi

if [ ! -d ${OUT_DIR}/${SAMP} ] ; then
	mkdir -p ${OUT_DIR}/${SAMP}
fi

echo "uncompressing bz2 read files"
date

lbzcat -n 40 -c ${READ1} > ${TMPDIR}/${SAMP}.R1.fastq
lbzcat -n 40 -c ${READ2} > ${TMPDIR}/${SAMP}.R2.fastq


echo "runnning fastp for quality"
date

fastp -w 16 -i ${TMPDIR}/${SAMP}.R1.fastq -I ${TMPDIR}/${SAMP}.R2.fastq\
  -o ${TMPDIR}/${SAMP}.R1.fastp.fastq -O ${TMPDIR}/${SAMP}.R2.fastp.fastq


echo "running seqkit stats for read info"
date

seqkit stats -T ${TMPDIR}/${SAMP}.R1.fastp.fastq ${TMPDIR}/${SAMP}.R2.fastp.fastq > ${OUT_DIR}/${SAMP}/${SAMP}_seqstats.tsv


echo "aligning reads to influenza A seq DB with minimap2"
date

minimap2 -t $CPUS -cx sr --secondary=yes\
  /gpfs1/projects/Pools/RD_projects/mjt_TWIST_reference_genomes/influenza_seqs/Influenza_A_segment_sequences.fna\
  ${TMPDIR}/${SAMP}.R1.fastp.fastq ${TMPDIR}/${SAMP}.R2.fastp.fastq > ${OUT_DIR}/${SAMP}/${SAMP}_influenza_A.cigar.paf

if [ -s ${OUT_DIR}/${SAMP}/${SAMP}_influenza_A.cigar.paf ] ; then

	echo "running read-to-serotype assignment R script"
	date

	Rscript ${SCRIPT_DIR}/parse_pafs_influenza_A.R ${SCRIPT_DIR}/Influenza_A_segment_info1.tsv\
	  ${OUT_DIR}/${SAMP}/${SAMP}_influenza_A.cigar.paf ${SAMP} "${OUT_DIR}/${SAMP}"
else
	echo "${OUT_DIR}/${SAMP}/${SAMP}_influenza_A.cigar.paf not found or empty"
fi

READ_NAMES_FILES=$( find ${OUT_DIR}/${SAMP}/ -type f -name "*.txt" )

if [ -n "$READ_NAMES_FILES" ] ; then
	echo "getting reads for each serotype assignment"

	echo "$READ_NAMES_FILES" | while read RN_FILE ; do

		seqkit grep -f $RN_FILE ${TMPDIR}/${SAMP}.R1.fastp.fastq > ${RN_FILE%.txt}.R1.fastq
		seqkit grep -f $RN_FILE ${TMPDIR}/${SAMP}.R2.fastp.fastq > ${RN_FILE%.txt}.R2.fastq
	done
else
	echo "no reads assigned to influenza A"
fi

echo "Influenza A serotype read assignment pipeline ending"
date









