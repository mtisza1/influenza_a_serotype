#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT

def grep_reads(rname_file: str, reads: str, sero_reads: str, cpus: int):
    
    return Popen(['seqkit', 'grep', '-j', cpus,
                    '-o', sero_reads,
                    '-f', rname_file, reads],
                    stdout=PIPE, stderr=STDOUT)
