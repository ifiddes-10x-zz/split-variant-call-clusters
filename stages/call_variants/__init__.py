#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import os
import martian
import subprocess
import tenkit.tabix as tk_tabix
import cellranger.utils as cr_utils

__MRO__ = '''
stage CALL_VARIANTS(
    in bam[] sorted_bams,
    in int[] clusters,
    in path reference_path,
    in bed target_regions,
    out vcf.gz variants,
    src py "call_variants",
) split using (
    in bam sorted_bam,
    in int cluster_id,
    out vcf.gz variant_subset,
)
'''


def split(args):
    return {'chunks': [{'sorted_bam': bam, 'cluster_id': i, '__threads': 2, '__mem_gb': 2 * 6}
                       for i, bam in zip(args.clusters, args.sorted_bams)],
            'join': {'__mem_gb': 16}}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
    tmp = martian.make_path('{}.vcf'.format(args.cluster_id))

    # Run GATK4
    gatk_args = ['gatk-launch', 'HaplotypeCaller',
                 '-R', genome_fasta_path,
                 '-I', args.sorted_bam,
                 '-O', tmp,
                 '-L', args.target_regions,
                 '--native-pair-hmm-threads', str(args.__threads)]

    subprocess.check_call(gatk_args)

    # fix the header
    recs = [x for x in open(tmp)]
    with open(tmp, 'w') as outf:
        for l in recs:
            if l.startswith('#CHROM'):
                l = l.split()
                l[-1] = str(args.cluster_id)
                l = '\t'.join(l) + '\n'
            outf.write(l)

    tk_tabix.sort_vcf(tmp, outs.variant_subset.replace('.gz', ''))
    tk_tabix.index_vcf(outs.variant_subset.replace('.gz', ''))
    os.remove(tmp)
    os.remove(tmp + '.idx')

def join(args, outs, chunk_defs, chunk_outs):
    # final merge to make one combined VCF
    tmp = martian.make_path('tmp.vcf')
    cmd = ['vcf-merge'] + [chunk.variant_subset for chunk in chunk_outs]
    with open(tmp, 'w') as outf:
        subprocess.check_call(cmd, stdout=outf)
    # Sort and index the files
    tk_tabix.sort_vcf(tmp, outs.variants.replace('.gz', ''))
    tk_tabix.index_vcf(outs.variants.replace('.gz', ''))
    os.remove(tmp)