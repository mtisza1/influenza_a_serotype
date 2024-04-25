#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT

def minimap2(reference: str, read1: str, read2: str, paf_file: str, cpus: int):

    return Popen(['minimap2', '-t', cpus, 
                    '-cx', 'sr', '--secondary=yes', 
                    '-o', paf_file,
                    reference, read1, read2],
                    stdout=PIPE, stderr=STDOUT)
