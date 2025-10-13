import os
import re
import logging
from typing import Dict, List, Tuple

import pysam
import polars as pl
import matplotlib.pyplot as plt

logger = logging.getLogger("iavs_logger")


def bam_to_score(in_bam: str) -> pl.DataFrame:
    filter_records = []

    with pysam.AlignmentFile(in_bam, 'r') as ministream:
        for record in ministream:
            readLength: int=record.infer_read_length()
            alignLength: int=record.query_alignment_length
            if record.is_paired:
                if record.is_read1:
                    orient: str = "R1"
                else:
                    orient: str = "R2"
            else:
                orient: str = "N"

            # fixes cases when the reads have .1, /1, .2, or /2 at the end
            endtup = (".1", "/1", ".2", "/2")
            if record.query_name.endswith(endtup):
                pairn: str = record.query_name[:-2]
            else:
                pairn: str = record.query_name

            # NM flag reports "edit distance" between read and ref
            try:
                alignEdit: int=dict(record.get_tags(with_value_type=False)).get('NM')
            except:
                alignEdit: int=0
            try:
                alignProp: float=alignLength/readLength
            except:
                alignProp: float=0
            ## assume the edit distance is between aligned part of read and ref
            try:
                alignAcc: float=(alignLength-alignEdit)/alignLength
            except:
                alignAcc: float=0

            if alignLength >= 100 and alignProp >= 0.9 and alignAcc >= 0.8:
                filter_records.append(
                    [
                        record.query_name,
                        record.reference_name,
                        readLength,
                        alignProp,
                        alignLength,
                        alignAcc,
                        orient,
                        pairn
                    ]
                        )
                
    filter_df = pl.DataFrame(filter_records,
                            schema=['qname', 'rname', 
                                    'read_length', 'aln_prop', 
                                    'aln_len', 'aln_acc', 
                                    'orientation', 'pair_name'],
                            orient = "row")

    # return a polars dataframe
    return filter_df

def load_flu_info(db_info_path: str) -> pl.DataFrame:
    df = pl.read_csv(
        db_info_path, 
        separator='\t', 
        has_header=True,
        schema_overrides={"accession": pl.Utf8, "serotype": pl.Utf8, "segment": pl.Utf8}
        )
    # Ensure required columns
    needed = {"accession", "serotype", "segment"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Flu info file missing columns: {missing}")
    return df.select(["accession", "serotype", "segment"]).with_columns(
        pl.col("accession").cast(pl.Utf8),
        pl.col("serotype").cast(pl.Utf8),
        pl.col("segment").cast(pl.Utf8),
    )

def pair_alns(alndf: pl.DataFrame) -> pl.DataFrame:
    mini2_p_df = alndf.with_columns(
        (pl.col("aln_acc") * pl.col("aln_prop")).alias("aln_score")
    ).sort(
        ['qname', 'aln_score'], 
        descending = True
    ).group_by(
        ['qname', 'pair_name', 'rname']
    ).agg(
        pl.col('read_length').first().alias('read_length'),
        pl.col('aln_prop').first().alias('aln_prop'),
        pl.col('aln_len').first().alias('aln_len'),
        pl.col('aln_acc').first().alias('aln_acc'),
        pl.col('aln_score').first().alias('aln_score')
    ).group_by(
        ['pair_name', 'rname']
    ).agg(
        pl.col('read_length').sum().alias('pair_length'),
        pl.col('read_length').count().alias('pair_count'),
        (pl.col("aln_acc").mean() * pl.col('aln_prop').mean() * pl.col('read_length').sum()).alias('pair_score'),  
    )
    return mini2_p_df


def compute_assignment(df_align: pl.DataFrame, flu_info: pl.DataFrame, score_thresh: float = 90) -> pl.DataFrame:
    if df_align.is_empty():
        return df_align

    merged = df_align.join(flu_info, left_on="rname", right_on="accession", how="inner")

    try:
        listed_clf_df = merged.sort(
            ['pair_name', 'pair_score', 'pair_length', 'serotype'], 
            descending = True
        ).group_by(['pair_name', 'serotype']).agg(
            pl.col('pair_score').first().alias('sero_score'),
            pl.col('pair_length').first().alias('sero_al_length'),
            pl.col('pair_count').first().alias('sero_read_count')
        ).sort(
            ['pair_name', 'sero_score'], 
            descending = True
        ).group_by(['pair_name']).agg(
            pl.col('sero_score').first().alias('first_score'),
            pl.col('sero_score').slice(1,1).first().alias('second_score'),
            pl.col('serotype').first().alias('first_sero'),
            pl.col('serotype').slice(1,1).first().alias('second_sero'),
            pl.col('sero_read_count').first().alias('first_nreads'),
            pl.col('sero_read_count').slice(1,1).first().alias('second_nreads'),
            pl.col('sero_al_length').first().alias('first_alength'),
            pl.col('sero_al_length').slice(1,1).first().alias('second_alength'),
        ).with_columns([
            # decision logic: best serotype has to be better than any other serotype
            pl.when(
                (
                    (pl.col("first_score") - 1 >= pl.col("second_score")) | \
                    pl.col("second_score").is_null() | \
                    (pl.col('first_sero') == pl.col('second_sero'))
                )
            ).then(
                pl.col("first_sero")
            ).otherwise(
                pl.lit("ambiguous")
            ).alias("read_assignment")
        ]).filter(
            pl.col("first_score") >= float(score_thresh)
        )
        return listed_clf_df

    except Exception as e:
        logger.warning(f"could not parse taxonomy from alignments")
        logger.warning(e)


def write_outputs(sum_df: pl.DataFrame, sample: str, out_dir: str) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    if sum_df.is_empty():
        return outputs

    summary_path = os.path.join(out_dir, f"{sample}_per_read_summary.tsv")
    sum_df.write_csv(summary_path, separator='\t', include_header=True)
    outputs["summary_tsv"] = summary_path

    # Bar plot of read_assignment (no pandas/pyarrow)
    counts_df = sum_df.group_by("read_assignment").agg(
        pl.col("first_nreads").sum().alias("read_count")
    ).sort(['read_assignment'])
    #counts_df = sum_df.group_by("read_assignment").len().sort("read_assignment")
    labels = counts_df["read_assignment"].to_list()
    values = counts_df["read_count"].to_list()
    plt.figure(figsize=(8, 4))
    plt.bar(labels, values, color="#4C78A8")
    plt.xticks(rotation=90)
    plt.ylabel("read count")
    plt.tight_layout()
    plot_path = os.path.join(out_dir, f"{sample}_read_serotype_assignment.pdf")
    plt.savefig(plot_path)
    plt.close()
    outputs["assignment_plot_pdf"] = plot_path

    #counts_df = counts_df.rename({"len": "read_count"})

    sero_path = os.path.join(out_dir, f"{sample}_per_serotype_summary.tsv")
    counts_df.write_csv(sero_path, separator='\t', include_header=True)
    outputs["sero_tsv"] = sero_path

    # Per-serotype read lists: one row per qname already
    best_per_read = sum_df

    # Write one file per assignment label
    parts = best_per_read.partition_by("read_assignment", as_dict=True, maintain_order=True)
    for label, g in parts.items():
        label = str(label[0])
        rname_path = os.path.join(out_dir, f"{sample}_{label}.txt")
        # Two columns: read_assignment, pair_name
        pl.DataFrame({
            "pair_name": g["pair_name"],
        }).write_csv(rname_path, separator='\t', include_header=False)
        outputs[f"reads_{label}"] = rname_path

    return outputs


def assign_serotypes(bam_path: str, db_info_path: str, sample: str, out_dir: str, score_thresh: float = 90) -> Dict[str, str]:
    #running pairs_alns on the output of bam_to_score
    score_df = bam_to_score(bam_path)
    df_align = pair_alns(
        score_df
    )
    if df_align.is_empty():
        logger.info("No alignments found; skipping assignment outputs.")
        return {}

    score_path = os.path.join(out_dir, f"{sample}_scores_summary.tsv")
    score_df.write_csv(score_path, separator='\t', include_header=True)

    pair_path = os.path.join(out_dir, f"{sample}_pairs_summary.tsv")
    df_align.write_csv(pair_path, separator='\t', include_header=True)
    flu_info = load_flu_info(db_info_path)
    sum_df = compute_assignment(df_align, flu_info, score_thresh)

    outputs = write_outputs(sum_df, sample, out_dir)
    return outputs
