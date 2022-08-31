#!/usr/bin/env python


import argparse
import sys
import pysam
from pathlib import Path


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reads from genome",
        epilog="Example: genome2reads genome.fa",
    )
    parser.add_argument(
        "genome",
        type=Path,
        help="Genome in FASTA format.",
    )
    parser.add_argument(
        "--output",
        default="reads.fastq",
        help="Output file in FASTQ format.",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        help="The desired coverage",
        type=int,
        default="30",
    )
    parser.add_argument(
        "-l",
        "--length",
        help="The desired read length",
        type=int,
        default="150",
    )
    parser.add_argument(
        "-t",
        "--tiling",
        help="The desired tiling",
        type=int,
        default="2",
    )

    return parser.parse_args(argv)


def read_genome(genome):
    """Read the genome from a FASTA file."""
    seq = ""
    with pysam.FastxFile(genome) as fh:
        for entry in fh:
            seq += entry.sequence
    return seq


def generate_reads(genome_name, genome_seq, coverage, read_length, tiling):
    """Generate reads from the genome.
    Args:
        genome_name (str): The name of the genome.
        genome_seq (str): The sequence of the genome.
        coverage (int): The desired coverage.
        read_length (int): The desired read length.
        tiling (int): The desired tiling.
    Returns:
        list: A list of reads.
    """
    reads = []
    phred = '@'*read_length

    for c in range(coverage):
        cnt = 0
        i = c * tiling
        while i + read_length < len(genome_seq):
            read_seq = genome_seq[i:i+read_length]
            fq_entry = f"@GENOME2READS.{cnt*(c+1)}\n{read_seq}\n+\n{phred}\n"
            reads.append(fq_entry)
            i += read_length + tiling
            cnt += 1
    return reads


def write_fastq(reads, output):
    """Write the reads to a FASTQ file."""
    with open(output, "w") as fh:
        for read in reads:
            fh.write(read)


def main(argv=None):
    args = parse_args(argv)
    print(args)
    genome_name = args.genome.stem
    genome_seq = read_genome(args.genome)
    reads = generate_reads(genome_name, genome_seq, args.coverage, args.length, args.tiling)
    write_fastq(reads, f"{args.genome.stem}.fastq")


if __name__ == "__main__":
    sys.exit(main())
