#!/usr/bin/env python3

import argparse
from assignment4_utils import (
    parse_hgnc_file,
    extract_chromosome,
    analyze_pseudogenes,
    analyze_locus_types,
    analyze_genes_per_chromosome,
    analyze_gene_groups,
    get_longest_ucsc_genes,
    chromosome_sort_key
)
from io_utils import FileHandler


##############################################
# WRITE OUTPUT FILES (USING FILEHANDLER)
##############################################

def write_gene_statistics_file(output_dir, chr_counts, pseudo_stats, locus_stats, group_stats):
    path = f"{output_dir}/gene_statistics.txt"
    with FileHandler(path, "w") as out:
        out.write("=== Gene Statistics Summary ===\n\n")

        total_genes = sum(count for _, count in chr_counts)
        out.write(f"Total genes: {total_genes}\n")
        out.write(f"Chromosomes found: {len(chr_counts)}\n\n")

        out.write("Genes per chromosome:\n")
        for chrom, count in chr_counts:
            out.write(f"  {chrom}: {count}\n")

        out.write("\n=== Pseudogenes ===\n")
        out.write(f"Total pseudogenes: {pseudo_stats['total_pseudogenes']}\n")
        out.write(f"Top pseudogene chromosome: {pseudo_stats['top_chromosome']}\n\n")

        out.write("=== Locus Types ===\n")
        out.write(f"Number of unique locus types: {locus_stats['num_unique_types']}\n")
        out.write(f"Most common locus type: {locus_stats['most_common_type']} ({locus_stats['most_common_count']} genes)\n\n")

        out.write("=== Gene Groups ===\n")
        out.write(f"Number of unique gene groups: {group_stats['num_unique_groups']}\n\n")

        out.write("Top 3 gene groups:\n")
        for name, count in group_stats["top_3_groups"]:
            out.write(f"  - {name}: {count}\n")

    print(f"Wrote gene_statistics.txt → {path}")


def write_top_gene_groups_file(output_dir, group_stats):
    path = f"{output_dir}/top_gene_groups.txt"
    with FileHandler(path, "w") as out:
        out.write("=== Top Gene Groups ===\n\n")
        out.write("Rank\tGroup Name\tCount\n")

        rank = 1
        for name, count in group_stats["sorted_groups"]:
            out.write(f"{rank}\t{name}\t{count}\n")
            rank += 1

    print(f"Wrote top_gene_groups.txt → {path}")


def write_pseudogenes_by_chromosome(output_dir, pseudo_stats):
    path = f"{output_dir}/pseudogenes_by_chromosome.txt"
    with FileHandler(path, "w") as out:
        out.write("Chromosome\tPseudogene_Count\n")
        for chrom, count in pseudo_stats["by_chromosome"]:
            out.write(f"{chrom}\t{count}\n")

    print(f"Wrote pseudogenes_by_chromosome.txt → {path}")


def write_locus_type_distribution(output_dir, locus_stats):
    path = f"{output_dir}/locus_type_distribution.txt"
    with FileHandler(path, "w") as out:
        out.write("Locus_Type\tCount\n")
        for locus_type, count in locus_stats["counts"].items():
            out.write(f"{locus_type}\t{count}\n")

    print(f"Wrote locus_type_distribution.txt → {path}")


def write_genes_per_chromosome(output_dir, chr_counts):
    path = f"{output_dir}/genes_per_chromosome.txt"
    with FileHandler(path, "w") as out:
        out.write("Chromosome\tGene_Count\n")
        for chrom, count in chr_counts:
            out.write(f"{chrom}\t{count}\n")

    print(f"Wrote genes_per_chromosome.txt → {path}")


def write_longest_ucsc_genes_file(output_dir, longest_genes):
    path = f"{output_dir}/longest_ucsc_genes.txt"
    with FileHandler(path, "w") as out:
        out.write("=== Longest UCSC Genes ===\n\n")
        out.write("Symbol\tUCSC_ID\tChrom\tStart\tEnd\tLength\n")

        for symbol, ucsc, chrom, start, end, length in longest_genes:
            out.write(f"{symbol}\t{ucsc}\t{chrom}\t{start}\t{end}\t{length}\n")

    print(f"Wrote longest_ucsc_genes.txt → {path}")


##############################################
# MAIN
##############################################

def main():
    parser = argparse.ArgumentParser(description="Generate HGNC gene statistics.")
    parser.add_argument("-i", "--input", required=True, help="Input HGNC TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output

    # Load data
    gene_data = parse_hgnc_file(input_file)
    print(f"Loaded {len(gene_data)} genes from {input_file}")

    # Analyses
    pseudo_stats = analyze_pseudogenes(gene_data)
    locus_stats = analyze_locus_types(gene_data)
    counts_raw = analyze_genes_per_chromosome(gene_data)
    chr_counts = sorted(counts_raw.items(), key=chromosome_sort_key)
    group_stats = analyze_gene_groups(gene_data)

    print(f"Chromosomes found: {len(chr_counts)}")
    print(f"Total pseudogenes: {pseudo_stats['total_pseudogenes']}")
    print(f"Top pseudogene chromosome: {pseudo_stats['top_chromosome']}")
    print(f"Most common locus type: {locus_stats['most_common_type']} ({locus_stats['most_common_count']} genes)")
    print(f"Number of unique gene groups: {group_stats['num_unique_groups']}")
    print("Top 3 gene groups:")
    for name, count in group_stats["top_3_groups"]:
        print(f"  - {name}: {count}")

    # Write all required files
    write_pseudogenes_by_chromosome(output_dir, pseudo_stats)
    write_locus_type_distribution(output_dir, locus_stats)
    write_genes_per_chromosome(output_dir, chr_counts)
    write_gene_statistics_file(output_dir, chr_counts, pseudo_stats, locus_stats, group_stats)
    write_top_gene_groups_file(output_dir, group_stats)

    # Longest genes
    longest_genes = get_longest_ucsc_genes(gene_data, top_n=100)
    write_longest_ucsc_genes_file(output_dir, longest_genes)


if __name__ == "__main__":
    main()

