#!/usr/bin/env python

import glob
import multiprocessing
import os
import re
import subprocess

import dxpy

@dxpy.entry_point('main')
def main(fastq, genomeindex_targz):

    fastq_dxfile = dxpy.DXFile(fastq)
    dxpy.download_dxfile(fastq_dxfile.get_id(), "input.fastq")

    genome_dxfile = dxpy.DXFile(genomeindex_targz)
    dxpy.download_dxfile(genome_dxfile.get_id(), "genome.tar.gz")
    os.makedirs("genome")
    tar_cmd = "tar xzvf genome.tar.gz -C genome"
    subprocess.check_call(tar_cmd, shell=True)
    genome_file = glob.glob("genome/*.bwt")[0]
    genome_file = re.sub("\.bwt$", "", genome_file)

    bwa_cmd = ("bwa mem -t {nproc} {genome} {fastq} | "
               "samtools view -u -S - | "
               "samtools sort -m 256M -@ {nproc} - output"
               .format(nproc=multiprocessing.cpu_count(),
                       genome=genome_file,
                       fastq="input.fastq"))
    subprocess.check_call(bwa_cmd, shell=True)

    bam = dxpy.upload_local_file("output.bam")

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    output = {}
    output["bam"] = dxpy.dxlink(bam)

    return output
