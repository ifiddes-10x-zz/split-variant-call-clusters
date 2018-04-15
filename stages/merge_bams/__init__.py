"""
Merges a set of BAMs that are split by locus and barcode into a per-barcode BAM

"""
import martian
import itertools
import os
import subprocess
import tenkit.bam as tk_bam

__MRO__ = """
stage MERGE_BAMS(
    in bam[] modified_filtered_bams,
    in int[] clusters,
    out bam[] merged_bams,
    out int[] merged_clusters,
    src py "merge_bams",
) split using(
    in bam filtered_bam,
    out bam merged_bam,
)
"""

def split(args):
    chunks = []
    for cluster_id, bams in itertools.groupby(zip(args.clusters, args.modified_filtered_bams), key=lambda x: x[0]):
        bams = zip(*list(bams))[1]
        chunks.append({'cluster_bams': bams, 'cluster_id': cluster_id, '__mem_gb': 8*6, '__threads': 8})
    return {'chunks': chunks}


def main(args, outs):
    args.coerce_strings()
    tmp_bam = martian.make_path(str(args.cluster_id) + '.unsorted.bam')
    tk_bam.concatenate(tmp_bam, args.cluster_bams)
    outs.merged_bam = martian.make_path('{}.bam'.format(args.cluster_id))
    subprocess.check_call(['sambamba', 'sort', '-t', str(args.__threads), '-o', outs.merged_bam, tmp_bam])
    os.remove(tmp_bam)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.merged_bams = [str(chunk.merged_bam) for chunk in chunk_outs]
    outs.merged_clusters = [int(chunk.cluster_id) for chunk in chunk_defs]