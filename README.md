# ORF Analyzer

A Python-based bioinformatics tool for identifying Open Reading Frames (ORFs) in DNA sequences using Biopython.

## Overview

Open Reading Frames (ORFs) are DNA regions that begin with a start codon and end with an in-frame stop codon. Identifying ORFs is a useful first step in sequence analysis and gene prediction.

This project analyzes DNA sequences from a FASTA file, searches for ORF candidates on both the forward and reverse-complement strands, translates detected ORFs into protein sequences, and exports the results to a CSV file.

## Features

- Reads DNA sequences from FASTA files.
- Searches all three reading frames on the forward strand.
- Searches all three reading frames on the reverse-complement strand.
- Identifies start codons (`ATG`) and stop codons (`TAA`, `TAG`, `TGA`).
- Reports multiple ORF candidates when multiple start codons lead to the same in-frame stop codon.
- Translates detected ORFs into protein sequences.
- Reports strand, reading frame, coordinates, sequence length, DNA sequence, and protein translation.
- Exports detected ORFs to a CSV file.
- Handles sequences in which no ORFs are detected.

## Technologies

- Python
- Biopython
- CSV

## Project Structure

```text
orf-analyzer/
├── data/
│   └── example_sequences.fasta
├── results/
│   └── orfs.csv
├── src/
│   └── orf_analyzer.py
└── README.md
```

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/draganaivanova/orf-analyzer.git
cd orf-analyzer
```

### 2. Install the dependency

Make sure Python is installed, then run:

```bash
python -m pip install biopython
```

## Usage

From the project root directory, run:

```bash
python src/orf_analyzer.py
```

The program reads sequences from `data/example_sequences.fasta`, displays detected ORFs in the terminal, and saves the results to:

```text
results/orfs.csv
```

## Output

The CSV file contains the following columns:

| Column | Description |
|---|---|
| `sequence_id` | Identifier of the input FASTA sequence |
| `strand` | DNA strand (`+` or `-`) |
| `frame` | Reading frame, numbered 1–3 |
| `start` | Start coordinate |
| `end` | End coordinate |
| `length` | ORF length in nucleotides |
| `orf_sequence` | DNA sequence of the ORF |
| `protein` | Translated protein sequence |

Coordinates are zero-based and end-exclusive, following Python's slicing convention. Reverse-strand coordinates are mapped to the original input sequence.

## Current Limitations

- The current implementation uses `ATG` as the start codon and `TAA`, `TAG`, and `TGA` as stop codons.
- ORF candidates shorter than 9 nucleotides are excluded.
- Multiple overlapping ORF candidates may be reported when different start codons share the same downstream in-frame stop codon.
- This is an educational sequence-analysis project, not a complete gene-prediction pipeline.
- The current implementation expects suitable nucleotide sequences and does not yet provide comprehensive input validation.

## Future Improvements

- Add automated unit tests.
- Improve input validation and error handling.
- Add configurable minimum ORF length.
- Improve CSV reporting and command-line options.
- Compare detected ORFs with annotations from reference datasets.

## Author

Dragana Ivanova

Developed as a bioinformatics portfolio project to practice Python programming, sequence analysis, and the use of Biopython.