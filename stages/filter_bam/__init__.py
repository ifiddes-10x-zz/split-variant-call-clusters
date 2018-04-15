"""
Filter BAM. Filters for both read quality as well as barcodes, producing a set of BAMs.

This adapts the cellranger filter_reads stage to also split based on barcode clusters


"""
import itertools
import pandas as pd
import tenkit.bam as tk_bam
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils


__MRO__ = """
stage FILTER_BAM(
    in bam possorted_bam,
    in csv barcode_clusters,
    out bam[] filtered_bams,
    out int[] clusters,
    src py "filter_bam",
) split using(
    in string chunk_start,
    in string chunk_end,
    in int cluster_id,
    in string[] barcodes,
    out bam filtered_bam,
)
"""


def split(args):
    df = pd.read_csv(args.barcode_clusters)
    # construct BAM chunks
    with tk_bam.create_bam_infile(args.possorted_bam) as in_bam:
        chunks = tk_bam.chunk_bam_records(in_bam, chunk_bound_key=cr_utils.pos_sort_key,
                                          chunk_size_gb=cr_constants.BAM_CHUNK_SIZE_GB,
                                          max_chunks=cr_constants.MAX_BAM_CHUNKS)
    # nest BAM chunks with clusters
    bc_chunks = []
    for cluster_id, d in df.groupby('Cluster'):
        for c in chunks:
            bc_chunks.append({'chunk_start': c['chunk_start'], 'chunk_end': c['chunk_end'],
                              'cluster_bcs': d.Barcode.tolist(), 'cluster_id': cluster_id,
                              '__mem_gb': 8})
    return {'chunks': bc_chunks}

def main(args, outs):
    outs.coerce_strings()

    in_bam = tk_bam.create_bam_infile(args.possorted_bam)
    in_bam_chunk = tk_bam.read_bam_chunk(in_bam, (args.chunk_start, args.chunk_end))
    out_bam, _ = tk_bam.create_bam_outfile(outs.filtered_bam, None, None, template=in_bam)
    cluster_bcs = set(args.cluster_bcs)

    for (tid, pos), reads_iter in itertools.groupby(in_bam_chunk, key=cr_utils.pos_sort_key):
        dupe_keys = set()
        for read in reads_iter:
            if cr_utils.get_read_barcode(read) not in cluster_bcs:
                continue

            if cr_utils.is_read_dupe_candidate(read, cr_utils.get_high_conf_mapq({"high_conf_mapq":60})):
                dupe_key = (cr_utils.si_pcr_dupe_func(read), cr_utils.get_read_umi(read))
                if dupe_key in dupe_keys:
                    continue

                dupe_keys.add(dupe_key)
                out_bam.write(read)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    outs.clusters = []
    outs.filtered_bams = []
    for d, o in zip(chunk_defs, chunk_outs):
        outs.clusters.append(d.cluster_id)
        outs.filtered_bams.append(o.filtered_bam)
