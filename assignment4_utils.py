#!/usr/bin/env python3
import re

############################################################
# 1. parse_hgnc_file
############################################################

def parse_hgnc_file(input_file):
    """
    Parse an HGNC TSV file into a dictionary of gene data.
    Keys = gene symbol
    Values = dict of columns
    """
    gene_data = {}

    with open(input_file, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))

            symbol = row.get("symbol", "").strip()
            if symbol == "":
                continue

            gene_data[symbol] = row

    return gene_data


############################################################
# 2. extract_chromosome
############################################################

def extract_chromosome(location):
    """
    Extract chromosome ID from HGNC 'location' string.
    Examples:
      '7q34' → '7'
      '19q13.43' → '19'
      'Xp22.2' → 'X'
      'Yq11.23' → 'Y'
    """
    if not location:
        return "Unknown"

    match = re.match(r"^([0-9]{1,2}|X|Y)", location)
    if match:
        return match.group(1)

    return "Unknown"


############################################################
# 3. chromosome_sort_key
############################################################

def chromosome_sort_key(item):
    """
    Sorting key for chromosome ordering.
    Works for items like ('1', count) or ('X', count).
    """
    if isinstance(item, tuple):
        chrom = item[0]
    else:
        chrom = item

    if chrom.isdigit():
        return (0, int(chrom))
    if chrom == "X":
        return (1, 23)
    if chrom == "Y":
        return (1, 24)

    return (2, chrom)


############################################################
# 4. analyze_pseudogenes
############################################################

def analyze_pseudogenes(gene_data):
    """
    Count pseudogenes per chromosome.
    """
    chrom_counts = {}

    for symbol, row in gene_data.items():
        locus_type = row.get("locus_type", "")
        if "pseudogene" not in locus_type:
            continue

        chrom = extract_chromosome(row.get("location", ""))
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

    sorted_counts = sorted(chrom_counts.items(), key=chromosome_sort_key)

    top_chrom = None
    if sorted_counts:
        top_chrom = max(sorted_counts, key=lambda x: x[1])

    return {
        "total_pseudogenes": sum(chrom_counts.values()),
        "by_chromosome": sorted_counts,
        "top_chromosome": top_chrom
    }


############################################################
# 5. analyze_locus_types
############################################################

def analyze_locus_types(gene_data):
    """
    Count how many genes appear in each locus_type.
    """
    counts = {}

    for symbol, row in gene_data.items():
        lt = row.get("locus_type", "unknown")
        counts[lt] = counts.get(lt, 0) + 1

    most_common_type = max(counts, key=counts.get)
    most_common_count = counts[most_common_type]

    return {
        "counts": counts,
        "num_unique_types": len(counts),
        "most_common_type": most_common_type,
        "most_common_count": most_common_count
    }


############################################################
# 6. analyze_genes_per_chromosome
############################################################

def analyze_genes_per_chromosome(gene_data):
    """Returns a dict: {chromosome: count}"""
    chrom_counts = {}

    for symbol, row in gene_data.items():
        chrom = extract_chromosome(row.get("location", ""))
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

    return chrom_counts


############################################################
# 7. analyze_gene_groups
############################################################

def analyze_gene_groups(gene_data):
    """
    Analyze gene_group counts.
    """
    counts = {}

    for symbol, row in gene_data.items():
        gg = row.get("gene_group", "")
        groups = gg.split(",") if gg else []

        for g in groups:
            g = g.strip()
            if g:
                counts[g] = counts.get(g, 0) + 1

    sorted_groups = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    top_3 = sorted_groups[:3]

    return {
        "counts": counts,
        "sorted_groups": sorted_groups,
        "top_3_groups": top_3,
        "num_unique_groups": len(counts)
    }


############################################################
# 8. get_longest_ucsc_genes  (FULLY INDENTED)
############################################################

def get_longest_ucsc_genes(gene_data, top_n=100):
    """
    Compute longest genes using location_sortable field.
    Returns list of:
    (symbol, ucsc_id, chromosome, start, end, length)
    """
    results = []

    for symbol, info in gene_data.items():
        loc = info.get("location_sortable", "")
        ucsc = info.get("ucsc_id", "")

        # Require UCSC ID
        if not ucsc:
            continue

        # Require usable coordinates
        if not loc or ":" not in loc or "-" not in loc:
            continue

        try:
            chrom, coords = loc.split(":")
            start_str, end_str = coords.split("-")
            start = int(start_str)
            end = int(end_str)
            length = end - start
        except Exception:
            continue

        results.append((symbol, ucsc, chrom, start, end, length))

    results.sort(key=lambda x: x[5], reverse=True)
    return results[:top_n]

