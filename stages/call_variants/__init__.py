#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import os
import martian
import collections
import subprocess
import pysam
import tenkit.bio_io as tk_io
import tenkit.tabix as tk_tabix
import cellranger.utils as cr_utils

__MRO__ = '''
stage CALL_VARIANTS(
    in bam[] merged_bams,
    in int[] merged_clusters,
    in path reference_path,
    out vcf variants,
    src py "call_variants",
) split using (
    in string locus,
    in bam merged_bam,
    out vcf.gz variant_subset,
)
'''


def split(args):
    # hacky fix to get chrom sizes
    fasta_index = cr_utils.get_reference_genome_fasta(args.reference_path) + '.fai'
    assert os.path.exists(fasta_index)

    loci = []
    chunk_size = 5 * 10 ** 7
    for l in open(fasta_index):
        l = l.split()
        chrom = l[0]
        chrom_length = int(l[1])
        region_start = 0
        while region_start < chrom_length:
            start = region_start
            end = region_start + chunk_size
            if end > chrom_length:
                end = chrom_length
            loci.append('{}:{}-{}'.format(chrom, start, end))
            region_start = end

    chunks = []
    for bam, cluster_id in zip(args.merged_bams, args.merged_clusters):
        for locus in loci:
            chunks.append({'locus': locus, 'merged_bam': bam, 'cluster_id': cluster_id, '__mem_gb': 6})
    return {'chunks': chunks, 'join': {'__mem_gb': 16}}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    #with open(outs.variant_subset, 'w') as outf:
        #subprocess.check_call(['freebayes', '-f', genome_fasta_path, '-r', args.locus, args.merged_bam],
        #                      stdout=outf)

    # get the sample. Making the assumption that this is a single-sample BAM
    s = pysam.Samfile(args.merged_bam)
    sample = s.header['RG'][0]['SM']

    subprocess.check_call(['gatk-launch', 'Mutect2', '-R', genome_fasta_path, '--intervals',
                           args.locus, '-I', args.merged_bam, '-tumor', sample,
                           '-O', outs.variant_subset, '--TMP_DIR', os.getcwd(),
                           '--native-pair-hmm-threads', str(args.__threads)])

    # fix the header
    recs = [x for x in open(outs.variant_subset)]
    with open(outs.variant_subset, 'w') as outf:
        for l in recs:
            if l.startswith('#CHROM'):
                l = l.split()
                l[-1] = str(args.cluster_id)
                l = '\t'.join(l) + '\n'
            outf.write(l)

def join(args, outs, chunk_defs, chunk_outs):
    # mapping of cluster ID -> VCFs
    to_merge = collections.defaultdict(list)
    for o, d in zip(chunk_outs, chunk_defs):
        to_merge[d.cluster_id].append(o.variant_subset)

    # merge each VCF subset for a cluster
    merged_vcfs = []
    for cluster_id, vcf_list in to_merge.iteritems():
        merged_vcf = martian.make_path('{}.vcf'.format(cluster_id))
        tk_io.combine_vcfs(merged_vcf, vcf_list)
        merged_vcfs.append(merged_vcf + '.gz')

    # final merge to make one combined VCF
    tmp = martian.make_path('tmp.vcf')
    cmd = ['vcf-merge'] + merged_vcfs
    with open(tmp, 'w') as outf:
        subprocess.check_call(cmd, stdout=outf)
    # Sort and index the files
    tk_tabix.sort_vcf(tmp, outs.variants.replace('.gz', ''))
    tk_tabix.index_vcf(outs.variants.replace('.gz', ''))
    os.remove(tmp)