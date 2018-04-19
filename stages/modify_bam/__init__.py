#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import subprocess
import cellranger.utils as cr_utils
import os

__MRO__ = '''
stage MODIFY_BAM(
    in bam[] filtered_bams,
    in path reference_path,
    out bam[] modified_filtered_bams,
    src py "stages/modify_bam",
) split using(
    in bam filtered_bam,
    out bam modified_filtered_bam,
)
'''

def split(args):
    return {'chunks': [{'filtered_bam': x, '__mem_gb': 8} for x in args.filtered_bams]}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)


    # Correct the STAR mapping from 255 to 60 and take care of split reads
    star_args = ['gatk-launch', 'SplitNCigarReads',
                 '-R', genome_fasta_path,
                 '-I', args.filtered_bam,
                 '-O', outs.modified_filtered_bam,
                 '--skip-mapping-quality-transform', 'false',
                 '--create-output-bam-index', 'false',
                 '--TMP_DIR', os.getcwd()]

    subprocess.check_call(star_args)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.modified_filtered_bams = [chunk.modified_filtered_bam for chunk in chunk_outs]