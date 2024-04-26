#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
# uncompress .bz2 fastqs
def unbz2_paired(bzread1:str, bzread2: str, deread1: str, deread2: str, cpu: str):

    return Popen(['lbzcat', '-n', cpu, '-c', bzread1, '>', 
                    deread1],
                    stdout=PIPE, stderr=STDOUT), \
            Popen(['lbzcat', '-n', cpu, '-c', bzread2, '>', 
                    deread2],
                    stdout=PIPE, stderr=STDOUT)
    
# fastp paired-end
def fastp_paired(read1:str, read2: str, fastp_read1: str, 
                 fastp_read2: str, cpus: str,
                 htmlo: str, jsono: str):

    icpus = int(cpus)
    if icpus > 16:
        cpus = str(16)
    
    return Popen(['fastp', '-w', cpus, 
                    '-i', read1, '-I', read2, 
                    '-o', fastp_read1, '-O', fastp_read2,
                    '-h', htmlo, '-j', jsono],
                    stdout=PIPE, stderr=STDOUT)
    

# seqkit stats
def seqkit_stats_paired(read1:str, read2: str, cpus: str, stat_file: str):

    return Popen(['seqkit', 'stats', '-j', cpus, '-T',
                    '-o', stat_file,
                    read1, read2],
                    stdout=PIPE, stderr=STDOUT)
    

