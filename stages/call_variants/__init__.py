#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import os
import martian
import collections
import subprocess
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
    out vcf variant_subset,
)
'''


def split(args):
    # hacky fix to get chrom sizes
    ref_dict = cr_utils.get_reference_genome_fasta(args.reference_path).replace('.fa', '.dict')
    assert os.path.exists(ref_dict)
    lines = [x.split() for x in open(ref_dict)][1:]
    seqs = [x[1].split(':')[1] for x in lines]
    sizes = map(int, [x[2].split(':')[1] for x in lines])

    # another hack because of the reference we are using
    def is_primary(chrom):
        chrom = chrom.replace('chr', '')
        return chrom.isdigit() or chrom == 'Y' or chrom == 'X'

    loci = [[chrom, 0, size] for chrom, size in zip(seqs, sizes) if is_primary(chrom)]
    chunks = []
    for bam, cluster_id in zip(args.merged_bams, args.merged_clusters):
        for locus in loci:
            chunks.append({'locus': locus, 'merged_bam': bam, 'cluster_id': cluster_id,
                           '__mem_gb': 40, '__threads': 8})
    return {'chunks': chunks, 'join': {'__mem_gb': 16}}


def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write('\t'.join(map(str, args.locus)))

    # Run GATK4
    gatk_args = ['gatk-launch', 'HaplotypeCaller',
                 '-R', genome_fasta_path,
                 '-I', args.merged_bam,
                 '-O', outs.variant_subset,
                 '-L', bed_path,
                 '--minimum-mapping-quality', '30',
                 '--min-base-quality-score', '20',
                 '--dont-use-soft-clipped-bases', 'true',
                 '--add-output-vcf-command-line', 'false',
                 '--native-pair-hmm-threads', str(args.__threads)]

    subprocess.check_call(gatk_args)

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
        to_merge[d.cluster].append(o.variant_subset)

    # merge each VCF subset for a cluster
    merged_vcfs = []
    for cluster_id, vcf_list in to_merge.iteritems():
        merged_vcf = martian.make_path('{}.vcf.gz'.format(cluster_id))
        tk_io.combine_vcfs(vcf_list, merged_vcf)
        merged_vcfs.append(merged_vcf)

    # final merge to make one combined VCF
    tmp = martian.make_path('tmp.vcf')
    cmd = ['vcf-merge'] + merged_vcfs
    with open(tmp, 'w') as outf:
        subprocess.check_call(cmd, stdout=outf)
    # Sort and index the files
    tk_tabix.sort_vcf(tmp, outs.variants)
    tk_tabix.index_vcf(outs.variants)