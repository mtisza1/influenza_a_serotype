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
import pysam
import polars as pl
try:
    from .iav_funcs import (
        unbz2_file,
        fastp_paired,
        seqkit_stats_paired,
        seqkit_stats_long,
        minimap2_sr,
        minimap2_long,
        grep_reads,
        is_tool,
        bam_list_fastq
    )
    from .assign import assign_serotypes
except:
    from iav_funcs import (
        unbz2_file,
        fastp_paired,
        seqkit_stats_paired,
        seqkit_stats_long,
        minimap2_sr,
        minimap2_long,
        grep_reads,
        is_tool,
        bam_list_fastq
    )
    from assign import assign_serotypes

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

    __version__ = "0.2.0"

    Def_CPUs = os.cpu_count()

    Def_workdir = os.getcwd()

    parser = argparse.ArgumentParser(description='assign IAV serotype to short reads. \
                                    Version ' + str(__version__))


    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for iav_serotype ')

    required_args.add_argument("-r", "--reads", nargs="+",
                            dest="READS", required=True, 
                            help='read file(s) in .fastq format, uncompressed. \
                                Specify 2 files of paired end short reads OR one or more long read files')
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
                            dest="READ_FMT", type=str, choices=['short_paired', 'long'], default='short_paired',
                            help='default = short_paired. short_paired must provide 2 .fastq files.')

    optional_args.add_argument("--db", 
                            dest="DB", type=str, default='default',
                            help='path to sequence database. If not set, iav_serotype looks for environmental \
                                variable IAVS_DB. Then, if this variable is unset, it this is unset, \
                                DB path is assumed to be ' + iavs_script_path.replace("src/iav_serotype", "DBs/v1.0"))

    optional_args.add_argument("--thresh", 
                            dest="THRESH", type=float, default=90,
                            help=f'Default = 0.9 -- minimum read score for serotype assignment (ANI*AF). \
                                Maximum possible score = 1')
    #MM_SET
    optional_args.add_argument("--mm_preset", 
                            dest="MM_SET", type=str, choices=['map-pb', 'map-ont', 'map-hifi', 'lr:hq'], default='lr:hq',
                            help=f'Default = lr:hq -- minimap2 alignemnt preset (for long read data only)')
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


    if len(READS.split()) == 2 and str(args.READ_FMT).lower() == "short_paired":
        logger.info(f'{len(READS.split())} {str(args.READ_FMT).lower()} short read files')
    elif len(READS.split()) != 2 and str(args.READ_FMT).lower() == "short_paired":
        logger.info("if stating --read_format short_paired, must provide exactly 2 read files. exiting.")
        sys.exit()
    elif str(args.READ_FMT).lower() == "long":
        logger.info(f'{len(READS.split())} {str(args.READ_FMT).lower()} long read files')

    if args.DB == "default" and os.getenv('IAVS_DB') != None:
        args.DB = os.getenv('IAVS_DB')
    elif args.DB == "default":
        args.DB = iavs_script_path.replace("src/iav_serotype", "DBs/v1.0")
    
    DB_file = os.path.join(str(args.DB), 'Influenza_A_segment_sequences.fna')

    if not os.path.isfile(DB_file):
        logger.error(f'database fasta file not found at {DB_file}. exiting.')
        sys.exit()


    logger.info(f"DB: {str(args.DB)}")

    logger.info(f"version:  {str(__version__)}")

    logger.info(f"read format : {str(args.READ_FMT)}")

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

    if str(args.READ_FMT).lower() == "short_paired":
        if not os.path.isfile(str(args.READS[0])) or os.path.getsize(str(args.READS[0])) == 0:
            logger.error(f'before: fastq processing, read 1 file empty or does not exist')
            sys.exit()

        if not os.path.isfile(str(args.READS[1])) or os.path.getsize(str(args.READS[1])) == 0:
            logger.error(f'before: fastq processing, read 2 file empty or does not exist')
            sys.exit()
    elif str(args.READ_FMT).lower() == "long":
        for readf in args.READS:
            if not os.path.isfile(str(readf)) or os.path.getsize(str(readf)) == 0:
                logger.error(f'before: fastq processing, read file {readf} is empty or does not exist')
    
    bz2 = False
    ### not developing this option right now.
    #for read in args.READS:
    #    if read.endswith('bz2'):
    #        bz2 = True
    
    if bz2 == True:
        lbzcat_starttime = time.perf_counter()

        logger.info(f'running lbzcat to uncompress reads')

        decomp_r1 = f'{iavs_temp}/{str(args.SAMPLE)}_R1.fastq'
        decomp_r2 = f'{iavs_temp}/{str(args.SAMPLE)}_R2.fastq'

        try:

            unbz2run1 = unbz2_file(str(args.READS[0]), 
                                   decomp_r1, str(args.CPU))

            with unbz2run1.stdout:
                log_subprocess_output(unbz2run1.stdout)
            exitcode = unbz2run1.wait()
        except Exception as e:
            logger.error(e)
            sys.exit

        try:

            unbz2run2 = unbz2_file(str(args.READS[1]), 
                                   decomp_r2, str(args.CPU))

            with unbz2run2.stdout:
                log_subprocess_output(unbz2run2.stdout)
            exitcode = unbz2run2.wait()
        except Exception as e:
            logger.error(e)
            sys.exit


        lbzcat_endtime = time.perf_counter()

        time_taken = lbzcat_endtime - lbzcat_starttime

        logger.info(f"> lbzcat took {timedelta(seconds=time_taken)}")

    stat_file = f'{samp_out_dir}/{str(args.SAMPLE)}_read_stats.tsv'
    aln_stem = f'{samp_out_dir}/{str(args.SAMPLE)}_influenza_A'

    if str(args.READ_FMT).lower() == "short_paired": #######
        if args.QUAL == True:

            fastp_starttime = time.perf_counter()


            ready_r1 = f'{iavs_temp}/{str(args.SAMPLE)}_R1.fastp.fastq'
            ready_r2 = f'{iavs_temp}/{str(args.SAMPLE)}_R2.fastp.fastq'

            htmlo = os.path.join(iavs_temp, f'{str(args.SAMPLE)}_fastp.html')
            jsono = os.path.join(iavs_temp, f'{str(args.SAMPLE)}_fastp.json')

            logger.info(f'running fastp to QC reads')
            try:
                fastprun = fastp_paired(str(args.READS[0]), str(args.READS[1]), 
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
            ready_r1 = str(args.READS[0])
            ready_r2 = str(args.READS[1])
    
        if not os.path.isfile(ready_r1) or os.path.getsize(ready_r1) == 0:
            logger.error(f'before: minimap2, read 1 file empty or does not exist')
            sys.exit()

        if not os.path.isfile(ready_r2) or os.path.getsize(ready_r2) == 0:
            logger.error(f'before: minimap2, read 2 file empty or does not exist')
            sys.exit()
    
        # seqkit stats

        seqkit_starttime = time.perf_counter()

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


        logger.info(f'running minimap2 to align reads to influenza A reference sequences')
        try:
            sorted_bam = minimap2_sr(DB_file, ready_r1, ready_r2, aln_stem, str(args.CPU))
        except Exception as e:
            logger.error(e)
            sys.exit

        mm_endtime = time.perf_counter()

        time_taken = mm_endtime - mm_starttime

        logger.info(f"> minimap2 took {timedelta(seconds=time_taken)}")

        # Python assignment pipeline (replaces R script)
        try:
            assign_outputs = assign_serotypes(
                sorted_bam,
                os.path.join(str(args.DB), 'Influenza_A_segment_info1.tsv'),
                str(args.SAMPLE),
                samp_out_dir,
                float(args.THRESH),
            )
            logger.info(f"Assignment outputs: {assign_outputs}")
        except Exception as e:
            logger.error(f"Assignment step failed: {e}")
            sys.exit()


    ### long reads option ###

    elif str(args.READ_FMT).lower() == "long":
        # seqkit stats

        seqkit_starttime = time.perf_counter()

        logger.info(f'running seqkit stats on reads')
        try:

            seqkitrun = seqkit_stats_long(str(READS), str(args.CPU), stat_file)

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


        logger.info(f'running minimap2 to align reads to influenza A reference sequences')
        try:
            sorted_bam = minimap2_long(DB_file, str(READS), str(args.MM_SET), aln_stem, str(args.CPU))
        except Exception as e:
            logger.error(e)
            sys.exit

        mm_endtime = time.perf_counter()

        time_taken = mm_endtime - mm_starttime

        logger.info(f"> minimap2 took {timedelta(seconds=time_taken)}")

        # Python assignment pipeline (replaces R script)
        try:
            assign_outputs = assign_serotypes(
                sorted_bam,
                os.path.join(str(args.DB), 'Influenza_A_segment_info1.tsv'),
                str(args.SAMPLE),
                samp_out_dir,
                float(args.THRESH),
            )
            logger.info(f"Assignment outputs: {assign_outputs}")
        except Exception as e:
            logger.error(f"Assignment step failed: {e}")
            sys.exit()

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

            ### grep short reads
            if str(args.READ_FMT).lower() == "short_paired":
            
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

            ### grep long reads
            elif str(args.READ_FMT).lower() == "long":
                sero_reads = serotype_r.replace('.txt', '.fastq')
                try:
                    greplongrun = grep_reads(serotype_r, str(READS), sero_reads, str(args.CPU))

                    with greplongrun.stdout:
                        log_subprocess_output(greplongrun.stdout)
                    exitcode = greplongrun.wait()
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