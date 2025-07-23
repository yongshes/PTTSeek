#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyBigWig
import ruptures as rpt
import argparse
import os
from multiprocessing import Pool, cpu_count

# Global: store args for multiprocessing
global_args = {}

def init_worker(args_dict):
    global global_args
    global_args = args_dict

def process_row(row):
    import pyBigWig  # Import inside process
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    gene = row["gene"]
    strand = row["strand"]

    try:
        bw = pyBigWig.open(global_args["bigwig"])
        if end - start < global_args["min_length"]:
            return None

        values = bw.values(chrom, start, end, numpy=True)
        bw.close()

        values = np.nan_to_num(values)
        if np.all(values == 0):
            return None
        if strand == "-":
            values = values[::-1]

        bin_size = global_args["bin"]
        num_bins = (end - start) // bin_size
        if num_bins < 5:
            return None

        binned = np.array([
            np.mean(values[i*bin_size:(i+1)*bin_size])
            for i in range(num_bins)
        ])
        if np.all(binned == 0):
            return None

        # Change point detection
        model = rpt.Binseg(model="l2").fit(binned)
        breakpoints = model.predict(n_bkps=1)

        # Find signal drop at changepoint
        if len(breakpoints) < 1:
            return None

        i = breakpoints[0]
        s1 = binned[:i]
        s2 = binned[i:]
        if len(s1) < 1 or len(s2) < 1:
            return None

        drop = np.mean(s2) - np.mean(s1)
        if drop >= 0:
            return None  # not a drop

        # Map index back to genome coordinate
        if strand == "+":
            tts_pos = start + i * bin_size
        else:
            tts_pos = end - i * bin_size

        return [chrom, tts_pos, tts_pos + 1, gene, drop, strand]

    except Exception as e:
        return None

def main():
    parser = argparse.ArgumentParser(description="Parallel TTS detection using ruptures")
    parser.add_argument("-b", "--bigwig", required=True, help="Input bigWig file (strand-specific)")
    parser.add_argument("-i", "--intron_bed", required=True, help="Intron BED file (with gene in col 4)")
    parser.add_argument("-o", "--output", default="tts_parallel.tsv", help="Output TSV file")
    parser.add_argument("--bin", type=int, default=10, help="Bin size in bp")
    parser.add_argument("--min-length", type=int, default=100, help="Minimum intron length")
    args = parser.parse_args()

    # Load BED
    bed = pd.read_csv(args.intron_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "gene", "score", "strand"])
    bed = bed.dropna(subset=["start", "end"])
    bed = bed[bed["end"] > bed["start"]]

    # Prepare global args for workers
    global_args_dict = {
        "bigwig": args.bigwig,
        "bin": args.bin,
        "min_length": args.min_length
    }

    print(f"Processing {len(bed)} introns with {cpu_count()} cores...")
    with Pool(cpu_count(), initializer=init_worker, initargs=(global_args_dict,)) as pool:
        results = pool.map(process_row, bed.to_dict("records"))

    # Collect best per gene
    gene_best = {}
    for r in results:
        if r is None:
            continue
        chrom, start, end, gene, drop, strand = r
        if gene not in gene_best or drop < gene_best[gene][4]:
            gene_best[gene] = r

    df = pd.DataFrame(list(gene_best.values()), columns=["chrom", "start", "end", "gene", "drop_value", "strand"])
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} TTS calls to {args.output}")

if __name__ == "__main__":
    main()

