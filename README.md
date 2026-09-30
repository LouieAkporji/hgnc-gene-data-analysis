# HGNC Gene Data Analysis

A Python-based bioinformatics project for analyzing gene data from the HUGO Gene Nomenclature Committee (HGNC). This project processes HGNC gene records and generates summary statistics describing chromosome distributions, pseudogenes, locus types, gene groups, and genomic lengths.

## Project Overview

The program reads an HGNC TSV dataset and performs several analyses, including:

- Counting genes across chromosomes
- Identifying and counting pseudogenes by chromosome
- Analyzing the distribution of gene locus types
- Identifying common gene groups
- Determining the longest genes using genomic coordinates
- Generating organized text output files containing analysis results

## Technologies & Skills

- Python 3
- Bioinformatics data processing
- TSV file parsing
- Dictionaries and data structures
- Regular expressions
- Object-oriented programming
- File I/O and context managers
- Command-line arguments with `argparse`
- Modular Python programming

## Project Structure

- `gene_stats.py` — Main command-line program that performs the analyses and generates output files.
- `assignment4_utils.py` — Functions for parsing HGNC data and calculating gene statistics.
- `io_utils.py` — Custom context manager for file handling.

## Running the Program

The program accepts an HGNC TSV input file and an output directory:

```bash
python3 gene_stats.py -i <input_file.tsv> -o <output_directory>
