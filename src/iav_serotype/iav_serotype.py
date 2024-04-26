#!/usr/bin/env python

import argparse
import sys, os
import subprocess
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
import time
from datetime import timedelta
import random
import string
import re
import logging
from distutils.spawn import find_executable

####### ####### this tool's functions ####### ####### 
####### ####### ##################### ####### ####### 

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

# run minimap with correct settings
def minimap2(reference: str, read1: str, read2: str, paf_file: str, cpus: str):

    return Popen(['minimap2', '-t', cpus, 
                    '-cx', 'sr', '--secondary=yes', 
                    '-o', paf_file,
                    reference, read1, read2],
                    stdout=PIPE, stderr=STDOUT)

# grep reads of each serotype
def grep_reads(rname_file: str, reads: str, sero_reads: str, cpus: str):
    
    return Popen(['seqkit', 'grep', '-j', cpus,
                    '-o', sero_reads,
                    '-f', rname_file, reads],
                    stdout=PIPE, stderr=STDOUT)

def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None

####### ####### ##################### ####### ####### 
####### ####### ##################### ####### ####### 


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


### entry point function for iav_serotype
def iav_serotype():   
    iavs_starttime = time.perf_counter()
    pathname = os.path.dirname(__file__)  
    iavs_script_path = os.path.abspath(pathname)      
    print(f"this script dir: {str(iavs_script_path)}")

    parentpath = Path(pathname).parents[1]

    __version__ = "0.1.0"

    Def_CPUs = os.cpu_count()

    Def_workdir = os.getcwd()

    parser = argparse.ArgumentParser(description='assign IAV serotype to short reads. \
                                    Version ' + str(__version__))


    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for iav_serotype ')

    required_args.add_argument("-r", "--reads", nargs="+",
                            dest="READS", required=True, 
                            help='read file(s) in .fastq format. May be compressed in .bz2 format. \
                                Specify 2 files of paired end reads.')
    required_args.add_argument("-s", "--sample", 
                            dest="SAMPLE", type=str, required=True, 
                            help='Sample name. No space characters, please.')
    required_args.add_argument("-o", "--output_dir", 
                            dest="OUTPUT_DIR", type=str, required=True, 
                            help='Output directory name. Will be created if it does not exist. \
                                Can be shared with other samples. No space characters, please. ')


    optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for iav_serotype.')

    optional_args.add_argument('--version', action='version', version=str(__version__))

    optional_args.add_argument("-t", "--cpu", dest="CPU", type=int, default=Def_CPUs, 
                            help=f"Default: {Def_CPUs} -- Example: 32 -- Number of CPUs available for iav_serotype. ")
    optional_args.add_argument("-wd", "--working_directory", dest="c_workdir", type=str, default=Def_workdir, 
                            help=f"Default: {Def_workdir} -- Set working directory with absolute or relative path. \
                                run directory will be created within.")
    optional_args.add_argument('-q', "--qual", dest="QUAL", type=str2bool, default='True',
                            help='True or False. Remove low-quality reads with fastp?')
    optional_args.add_argument("--temp", 
                            dest="TEMP_DIR", type=str, default='default',
                            help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/')
    optional_args.add_argument("--keep", 
                            dest="KEEP", type=str2bool, default='False',
                            help='True of False. Keep the intermediate files, located in the temporary directory? \
                                These can add up, so it is not recommended if space is a concern.')
    optional_args.add_argument("-p", "--read_format", 
                            dest="READ_FMT", type=str, default='paired',
                            help='ONLY SUPPORTING PAIRED RIGHT NOW. Must be 2 files.')

    optional_args.add_argument("--db", 
                            dest="DB", type=str, default='default',
                            help='path to sequence database. If not set, iav_serotype looks for environmental \
                                variable IAVS_DB. Then, if this variable is unset, it this is unset, \
                                DB path is assumed to be ' + iavs_script_path.replace("src/iav_serotype", "DBs/v1.0"))

    optional_args.add_argument("--thresh", 
                            dest="THRESH", type=int, default=0.9,
                            help=f'Default = 0.9 -- minimum read score for serotype assignment (ANI*AF). \
                                Maximum possible score = 1')
    args = parser.parse_args()

    ## make directories
    out_directory = os.path.join(str(args.c_workdir), str(args.OUTPUT_DIR))

    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)

    samp_out_dir = os.path.join(out_directory, str(args.SAMPLE))

    if not os.path.isdir(samp_out_dir):
        os.makedirs(samp_out_dir)

    if str(args.TEMP_DIR) == 'default':
        iavs_temp = os.path.join(samp_out_dir, f'{str(args.SAMPLE)}_temp')
    else:
        iavs_temp = str(args.TEMP_DIR)

    if not os.path.isdir(iavs_temp):
        os.makedirs(iavs_temp)

    #### define logger #####
    logger = logging.getLogger("iavs_logger")
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(os.path.join(out_directory, f"{str(args.SAMPLE)}_iavs.log"))
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    #########################

    READS = ' '.join(map(str,args.READS))


    if len(READS.split()) == 2 and str(args.READ_FMT).lower() == "paired":
        logger.info(f'{len(READS.split())} {str(args.READ_FMT).lower()} read files')
    elif len(READS.split()) != 2 and str(args.READ_FMT).lower() == "paired":
        logger.info("if stating --read_format paired, must provide exactly 2 read files. exiting.")
        sys.exit()

    if args.DB == "default" and os.getenv('IAVS_DB') != None:
        args.DB = os.getenv('IAVS_DB')
    elif args.DB == "default":
        args.DB = iavs_script_path.replace("src/iav_serotype", "DBs/v1.0")
    
    DB_file = os.path.join(str(args.DB), 'Influenza_A_segment_sequences.fna')

    if not os.path.isfile(DB_file):
        logger.error(f'database fasta file not found at {DB_file}. exiting.')
        sys.exit()


    print("DB: ", str(args.DB))

    print("version ", str(__version__))

    completedProc = subprocess.run(['Rscript', str(iavs_script_path) + '/check_R_libraries1.R'])

    print(completedProc.returncode)
    if completedProc.returncode != 0 :
        print ("some required R packages are not found. Required:")
        print ("dplyr, data.table, stringr, ggplot2")
        print ("Did you activate the conda environment?")
        print ("see yml. Exiting")
        quit()



    tool_dep_list = ['minimap2', 'fastp', 'seqkit'] #lbzcat
    
    for tool in tool_dep_list:
        if not is_tool(tool):
            logger.warning(f"{tool} is not found. Exiting.")
            sys.exit()   


    #### define logging of subprocess (for parse_pafs_influenza_A.R) ####
    def log_subprocess_output(pipe):
        for line in iter(pipe.readline, b''): # b'\n'-separated lines
            logger.info(line.decode("utf-8").rstrip('\n'))

    ### Actually do stuff

    if not os.path.isfile(str(args.READS[0])) or os.path.getsize(str(args.READS[0])) == 0:
        logger.error(f'before: fastq processing, read 1 file empty or does not exist')
        sys.exit()

    if not os.path.isfile(str(args.READS[1])) or os.path.getsize(str(args.READS[1])) == 0:
        logger.error(f'before: fastq processing, read 2 file empty or does not exist')
        sys.exit()
    
    bz2 = False
    for read in args.READS:
        if read.endswith('bz2'):
            bz2 = True
    
    if bz2 == True:
        lbzcat_starttime = time.perf_counter()

        logger.info(f'running lbzcat to uncompress reads')

        decomp_r1 = f'{iavs_temp}/{str(args.SAMPLE)}_R1.fastq'
        decomp_r2 = f'{iavs_temp}/{str(args.SAMPLE)}_R2.fastq'

        try:

            unbz2run = unbz2_paired(str(args.READS[0]), str(args.READS[1]), 
                                   decomp_r1, decomp_r2, 
                                   str(args.CPU))

            with unbz2run.stdout:
                log_subprocess_output(unbz2run.stdout)
            exitcode = unbz2run.wait()
        except Exception as e:
            logger.error(e)
            sys.exit



        lbzcat_endtime = time.perf_counter()

        time_taken = lbzcat_endtime - lbzcat_starttime

        logger.info(f"> lbzcat took {timedelta(seconds=time_taken)}")
    else:
        decomp_r1 = str(args.READS[0])
        decomp_r2 = str(args.READS[1])

    if not os.path.isfile(decomp_r1) or os.path.getsize(decomp_r1) == 0:
        logger.error(f'before: fastq quality, read 1 file empty or does not exist')
        sys.exit()

    if not os.path.isfile(decomp_r2) or os.path.getsize(decomp_r2) == 0:
        logger.error(f'before: fastq quality, read 2 file empty or does not exist')
        sys.exit()

    if args.QUAL == True:

        fastp_starttime = time.perf_counter()


        ready_r1 = f'{iavs_temp}/{str(args.SAMPLE)}_R1.fastp.fastq'
        ready_r2 = f'{iavs_temp}/{str(args.SAMPLE)}_R2.fastp.fastq'

        htmlo = os.path.join(iavs_temp, f'{str(args.SAMPLE)}_fastp.html')
        jsono = os.path.join(iavs_temp, f'{str(args.SAMPLE)}_fastp.json')

        logger.info(f'running fastp to QC reads')
        try:
            fastprun = fastp_paired(decomp_r1, decomp_r2, 
                                    ready_r1, ready_r2, 
                                    str(args.CPU),
                                    htmlo, jsono)

            with fastprun.stdout:
                log_subprocess_output(fastprun.stdout)
            exitcode = fastprun.wait()

            
        except Exception as e:
            logger.error(e)
            sys.exit

        fastp_endtime = time.perf_counter()

        time_taken = fastp_endtime - fastp_starttime

        logger.info(f"> fastp took {timedelta(seconds=time_taken)}")
    else:
        ready_r1 = decomp_r1
        ready_r2 = decomp_r2
    
    if not os.path.isfile(ready_r1) or os.path.getsize(ready_r1) == 0:
        logger.error(f'before: minimap2, read 1 file empty or does not exist')
        sys.exit()

    if not os.path.isfile(ready_r2) or os.path.getsize(ready_r2) == 0:
        logger.error(f'before: minimap2, read 2 file empty or does not exist')
        sys.exit()
    
    # seqkit stats
    #f'{outdir}/{samp}/{samp}_read_stats.tsv'

    seqkit_starttime = time.perf_counter()

    stat_file = f'{samp_out_dir}/{str(args.SAMPLE)}_read_stats.tsv'

    logger.info(f'running seqkit stats on reads')
    try:

        seqkitrun = seqkit_stats_paired(ready_r1, ready_r2, str(args.CPU), stat_file)

        with seqkitrun.stdout:
            log_subprocess_output(seqkitrun.stdout)
        exitcode = seqkitrun.wait()
        
    except Exception as e:
        logger.error(e)
        sys.exit

    seqkit_endtime = time.perf_counter()

    time_taken = seqkit_endtime - seqkit_starttime

    logger.info(f"> seqkit stats took {timedelta(seconds=time_taken)}")

    # minimap

    mm_starttime = time.perf_counter()

    paf_file = f'{samp_out_dir}/{str(args.SAMPLE)}_influenza_A.cigar.paf'

    logger.info(f'running minimap2 to align reads to influenza A reference sequences')
    try:

        minimaprun = minimap2(DB_file, ready_r1, ready_r2, paf_file, str(args.CPU))

        with minimaprun.stdout:
            log_subprocess_output(minimaprun.stdout)
        exitcode = minimaprun.wait()
    except Exception as e:
        logger.error(e)
        sys.exit

    mm_endtime = time.perf_counter()

    time_taken = mm_endtime - mm_starttime

    logger.info(f"> minimap2 took {timedelta(seconds=time_taken)}")


    ### run the R script for parsing and assignment
    Rprocess = Popen(['Rscript', 
                     str(f'{iavs_script_path}/parse_pafs_influenza_A.R'), 
                     str(f'{args.DB}/Influenza_A_segment_info1.tsv'), 
                     str(paf_file), 
                     str(args.SAMPLE), 
                     samp_out_dir,
                     str(args.THRESH)],
                    stdout=PIPE, stderr=STDOUT)

    with Rprocess.stdout:
        log_subprocess_output(Rprocess.stdout)
    exitcode = Rprocess.wait()


    # retrieve reads by serotype

    assigned_read_list = []
    for readlist in os.listdir(samp_out_dir):
        if readlist.endswith('.txt'):
            f = os.path.join(samp_out_dir, readlist)

            if os.path.isfile(f) and os.path.getsize(f) > 0:
                assigned_read_list.append(f)
    
    if assigned_read_list:
        grep_starttime = time.perf_counter()

        logger.info('getting reads for each serotype assignment')

        for serotype_r in assigned_read_list:
            sero_reads1 = serotype_r.replace('.txt', '.R1.fastq')
            sero_reads2 = serotype_r.replace('.txt', '.R2.fastq')

            try:
                #read R1
                grep1run = grep_reads(serotype_r, ready_r1, sero_reads1, str(args.CPU))

                with grep1run.stdout:
                    log_subprocess_output(grep1run.stdout)
                exitcode = grep1run.wait()

                 #read R2
                grep2run = grep_reads(serotype_r, ready_r2, sero_reads2, str(args.CPU))

                with grep2run.stdout:
                    log_subprocess_output(grep2run.stdout)
                exitcode = grep2run.wait()


            except Exception as e:
                logger.error(e)
            
        grep_endtime = time.perf_counter()

        time_taken = grep_endtime - grep_starttime

        logger.info(f"> seqkit grep of serotyped reads took {timedelta(seconds=time_taken)}")
    else:
        logger.info(f'no reads assigned to influenza A')

    if not args.KEEP == True:
        if len(iavs_temp) > 0:
            for tempfile in os.listdir(iavs_temp):
                if tempfile.startswith(f'{str(args.SAMPLE)}'):
                    subprocess.run(['rm', os.path.join(iavs_temp, tempfile)])


    iavs_endtime = time.perf_counter()

    time_taken = iavs_endtime - iavs_starttime

    time_taken = round(time_taken, 2) 

    logger.info(f"This iav_serotype run finished in {timedelta(seconds=time_taken)}")

    logger.info(f'>>> if successful, run files can be found in: {samp_out_dir}')

if __name__ == "__main__":
    iav_serotype()