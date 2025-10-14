import subprocess
from subprocess import Popen, PIPE, STDOUT
import os
import pysam
import shutil

# uncompress .bz2 fastqs
def unbz2_file(bzread1: str, deread1: str, cpu: str):

    return Popen(['lbzcat', '-n', cpu, '-c', bzread1],
                    stdout=open(deread1, "w"), stderr=STDOUT)


# fastp paired-end
def fastp_paired(read1: str, read2: str, fastp_read1: str,
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


# seqkit stats paired short reads
def seqkit_stats_paired(read1: str, read2: str, cpus: str, stat_file: str):

    return Popen(['seqkit', 'stats', '-j', cpus, '-T',
                    '-o', stat_file,
                    read1, read2],
                    stdout=PIPE, stderr=STDOUT)


# seqkit stats all
def seqkit_stats_long(reads: str, cpus: str, stat_file: str):

    return Popen(['seqkit', 'stats', '-j', cpus, '-T',
                    '-o', stat_file,
                    reads],
                    stdout=PIPE, stderr=STDOUT)


# run minimap with correct settings
## should be for paired short reads
def minimap2_sr(reference: str, read1: str, read2: str, file_stem: str, cpus: str):
    """
    Align paired-end short reads with minimap2 to SAM, filter by alignment fraction >= 0.8,
    write BAM, sort it, and return the sorted BAM path. Output paths are derived from `file_stem`.
    """
    sam_file = f"{file_stem}.sam"
    filtered_bam = f"{file_stem}.filtered.bam"
    sorted_bam = f"{file_stem}.sorted.bam"

    # Run minimap2 to produce SAM (-ax sr)
    mini2_command = [
        'minimap2', '-t', str(cpus),
        '-ax', 'sr', '--secondary=yes', '--MD',
        '-f', '100000', '--sam-hit-only',
        '-N', '1000',
        '-p', '0.90',
        reference, read1, read2
    ]

    with open(sam_file, 'w') as sam_out:
        mini2_process = subprocess.Popen(mini2_command, stdout=sam_out, stderr=subprocess.PIPE)
        _, stderr = mini2_process.communicate()
        if mini2_process.returncode != 0:
            raise RuntimeError(f"minimap2 failed with code {mini2_process.returncode}: {stderr.decode().strip()}")

    # write BAM
    with pysam.AlignmentFile(sam_file, "r") as in_sam, \
         pysam.AlignmentFile(filtered_bam, "wb", template=in_sam) as out_bam:
        for aln in in_sam:
            out_bam.write(aln)

    # Sort the BAM
    pysam.sort('-@', str(cpus), '-o', sorted_bam, filtered_bam)

    if os.path.isfile(sam_file):
        os.remove(sam_file)
    if os.path.isfile(filtered_bam):
        os.remove(filtered_bam)

    return sorted_bam


# run minimap with long read correct settings
## should be for long reads
def minimap2_long(reference: str, reads: str, map_set, file_stem: str, cpus: str):
    """
    Align long reads with minimap2 to SAM, filter by alignment fraction >= 0.8,
    write BAM, sort it, and return the sorted BAM path. Output paths are derived from `file_stem`.
    """
    sam_file = f"{file_stem}.sam"
    filtered_bam = f"{file_stem}.filtered.bam"
    sorted_bam = f"{file_stem}.sorted.bam"

    # Support multiple read files when `reads` is a space-delimited string
    read_args = reads.split() if isinstance(reads, str) else list(reads)
    mini2_command = [
        'minimap2', '-t', str(cpus),
        '-ax', str(map_set), '--secondary=yes', '--MD',
        '-f', '100000', '--sam-hit-only',
        '-N', '1000',
        '-p', '0.90',
        reference,
    ] + read_args

    with open(sam_file, 'w') as sam_out:
        mini2_process = subprocess.Popen(mini2_command, stdout=sam_out, stderr=subprocess.PIPE)
        _, stderr = mini2_process.communicate()
        if mini2_process.returncode != 0:
            raise RuntimeError(f"minimap2 failed with code {mini2_process.returncode}: {stderr.decode().strip()}")

    # write BAM
    with pysam.AlignmentFile(sam_file, "r") as in_sam, \
         pysam.AlignmentFile(filtered_bam, "wb", template=in_sam) as out_bam:
        for aln in in_sam:
            out_bam.write(aln)

    # Sort the BAM
    pysam.sort('-@', str(cpus), '-o', sorted_bam, filtered_bam)

    return sorted_bam


# grep reads of each serotype
def grep_reads(rname_file: str, reads: str, sero_reads: str, cpus: str):
    read_args = reads.split() if isinstance(reads, str) else list(reads)
    return Popen(['seqkit', 'grep', '-j', cpus,
                    '-o', sero_reads,
                    '-r',
                    '-f', rname_file, *read_args],
                    stdout=PIPE, stderr=STDOUT)

def bam_list_fastq(rname_file: str, sorted_bam: str, out_fastq_r1: str, out_fastq_r2: str):
    """
    Extract reads listed in `rname_file` from a sorted BAM and write them in paired FASTQ format.
    - `out_fastq_r1` receives read1 mates (flag is_read1)
    - `out_fastq_r2` receives read2 mates (flag is_read2)
    Notes:
    - Only alignments present in the BAM are emitted (e.g., mates filtered out earlier won't be written).
    - The `rname_file` may contain one or more whitespace-separated columns; the read ID is taken
      as the last token of each non-empty line.
    Returns a tuple of output FASTQ paths: (out_fastq_r1, out_fastq_r2).
    """
    # Load read IDs
    wanted = set()
    with open(rname_file, 'r') as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            rid = s.split()[-1]
            wanted.add(rid)

    if not wanted:
        # Nothing to do; create empty files to be consistent with downstream steps
        open(out_fastq_r1, 'w').close()
        open(out_fastq_r2, 'w').close()
        return out_fastq_r1, out_fastq_r2

    written_r1 = set()
    written_r2 = set()
    same_out = out_fastq_r1 == out_fastq_r2
    with pysam.AlignmentFile(sorted_bam, 'rb') as bam:
        if same_out:
            with open(out_fastq_r1, 'w') as fq:
                for aln in bam.fetch(until_eof=True):
                    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                        continue
                    qn = aln.query_name
                    if qn not in wanted:
                        continue
                    seq = aln.query_sequence or ''
                    qual = aln.qual
                    if qual is None and aln.query_qualities is not None:
                        qual = ''.join(chr(q + 33) for q in aln.query_qualities)
                    if qual is None:
                        qual = 'I' * len(seq)
                    # Write all matching primary alignments to the single FASTQ, dedup by qname
                    if qn in written_r1:
                        continue
                    fq.write(f"@{qn}\n{seq}\n+\n{qual}\n")
                    written_r1.add(qn)
        else:
            with open(out_fastq_r1, 'w') as fq1, open(out_fastq_r2, 'w') as fq2:
                for aln in bam.fetch(until_eof=True):
                    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                        continue
                    qn = aln.query_name
                    if qn not in wanted:
                        continue
                    seq = aln.query_sequence or ''
                    qual = aln.qual
                    if qual is None and aln.query_qualities is not None:
                        qual = ''.join(chr(q + 33) for q in aln.query_qualities)
                    if qual is None:
                        qual = 'I' * len(seq)
                    if aln.is_paired:
                        if aln.is_read1 and qn not in written_r1:
                            fq1.write(f"@{qn}\n{seq}\n+\n{qual}\n")
                            written_r1.add(qn)
                        elif aln.is_read2 and qn not in written_r2:
                            fq2.write(f"@{qn}\n{seq}\n+\n{qual}\n")
                            written_r2.add(qn)
                    else:
                        if qn not in written_r1:
                            fq1.write(f"@{qn}\n{seq}\n+\n{qual}\n")
                            written_r1.add(qn)

    return out_fastq_r1, out_fastq_r2

def is_tool(name):
    """Check whether `name` is on PATH."""
    
    return shutil.which(name) is not None