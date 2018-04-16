#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import tenkit.bam as tk_bam
import cellranger.utils as cr_utils
import os

__MRO__ = '''
stage MODIFY_BAM(
    in bam[] filtered_bams,
    in path reference_path,
    in int[] clusters,
    out bam[] sorted_bams,
    src py "modify_bam",
) split using(
    in bam filtered_bam,
    in int cluster_id,
    out bam sorted_bam,
)
'''

def split(args):
    return {'chunks': [{'filtered_bam': bam, 'cluster_id': i, '__mem_gb': 8}
                       for bam, i in zip(args.filtered_bams, args.clusters)]}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    tmp_bam = martian.make_path(str(args.cluster_id) + '.unsorted.bam')
    # Correct the STAR mapping from 255 to 60 and take care of split reads
    star_args = ['gatk-launch', 'SplitNCigarReads',
                 '-R', genome_fasta_path,
                 '-I', args.filtered_bam,
                 '-O', tmp_bam,
                 '--skip-mapping-quality-transform', 'false',
                 '--create-output-bam-index', 'false',
                 '--TMP_DIR', os.getcwd()]
    subprocess.check_call(star_args)

    tmp_bam2 = martian.make_path(str(args.cluster_id) + '.fixed.bam')
    in_bam = tk_bam.create_bam_infile(tmp_bam)
    out_bam, _ = tk_bam.create_bam_outfile(tmp_bam2, None, None, template=in_bam)
    for rec in in_bam:
        rec.is_duplicate = False
        out_bam.write(rec)
    out_bam.close()


    outs.modified_filtered_bam = martian.make_path('{}.bam'.format(args.cluster_id))
    subprocess.check_call(['sambamba', 'sort', '-t', str(args.__threads), '-o', outs.sorted_bam, tmp_bam2])
    os.remove(tmp_bam)
    os.remove(tmp_bam2)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.sorted_bams = [chunk.sorted_bam for chunk in chunk_outs]
