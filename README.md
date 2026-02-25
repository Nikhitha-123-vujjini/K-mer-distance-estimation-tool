# k-mer Distance Estimation Tool

## Overview

This project implements a Python-based, alignment-free tool for comparing DNA sequences using k-mer frequency profiles. The program extracts overlapping k-mers of a user-defined length from FASTA files, converts raw counts into normalized frequency vectors, and computes pairwise distances between all sequences using Euclidean, Manhattan, or Cosine distance metrics. Results are saved as a CSV distance matrix and visualized as a heatmap image.

## Objectives

- Implement an alignment-free approach to DNA sequence comparison
- Support multiple FASTA file inputs
- Extract k-mers and build normalized frequency vectors using Python only
- Compute pairwise distances using three different metrics
- Output results as both a CSV distance matrix and a heatmap visualization

## Dependencies

- Python 3
- Biopython — FASTA file parsing (Bio.SeqIO)
- Matplotlib — heatmap visualization

## Usage
python kmer_tool.py <fasta_file(s)> -k <kmer_size> -m <metric> -o <output_dir>

### Arguments
| Argument | Flag | Description | Default |
|----------|------|------------|----------|
| FASTA file(s) | positional | One or more input FASTA files | required |
| K-mer size | `-k / --kmer_size` | Length of k-mers to extract | required |
| Distance metric | `-m / --metric` | `euclidean`, `manhattan`, or `cosine` | euclidean |
| Output directory | `-o / --output_dir` | Directory to save results | results |

## Example
python kmer_tool.py seq1.fasta seq2.fasta -k 4 -m cosine -o my_results

## Code Architecture
The pipeline flows as follows:

FASTA input → load_fasta → extract_kmers → count_kmers → unique_kmers → count_to_vectors → pairwise_comparisons → heatmap visualization

- load_fasta() — reads FASTA files using Biopython's SeqIO parser, generates unique labels (basename_header), and stores sequences as strings in a dictionary
- extract_kmers() — implements a sliding-window algorithm to extract all overlapping substrings of size k; handles sequences of any length with no gaps
- count_kmers() — converts the list of k-mers into a frequency dictionary {kmer: count}; forms the basis for vector normalization
- unique_kmers() — uses a Python set() to collect all unique k-mers across all sequences, then converts to a sorted list to ensure all vectors share identical dimensions
- count_to_vectors() — converts k-mer counts into normalized frequency vectors using: frequency = count / total k-mers
- euclidean() / manhattan() / cosine() — manual implementations of each distance metric
- pairwise_comparision() — computes an N × N distance matrix by comparing every vector against every other
- heatmap_plot() — converts the numeric matrix into a color-coded heatmap and saves it as a PNG
- main() — orchestrates the full pipeline from input parsing to output generation


## Distance Metrics
Three distance metrics are available with the -m flag:
### Euclidean
Distance = sqrt( sum( (v1_i - v2_i)^2 ) )

Measures straight-line distance between two vectors. Sensitive to magnitude differences — large differences are amplified by squaring. Useful for detecting overall compositional variation between sequences.
### Manhattan
Distance = sum( abs(v1_i - v2_i) )

Measures total absolute difference across all vector positions. Less sensitive to large individual differences than Euclidean and efficient to compute in high-dimensional k-mer space.
### Cosine
cosine_similarity = dot(v1, v2) / ( sqrt(sum(v1_i^2)) * sqrt(sum(v2_i^2)) )

cosine_distance = 1 - cosine_similarity

Measures the angle between two vectors rather than their magnitude. Focuses on compositional pattern similarity, making it ideal for length-independent comparisons. Sequences with similar k-mer composition yield a small cosine distance.

## Output
All results are saved to the specified output directory:

- distance_matrix.csv — pairwise distance matrix for all input sequences. The diagonal is 0 (each sequence compared to itself) and off-diagonal values reflect k-mer dissimilarity.
- heatmap.png — color-coded heatmap of the distance matrix using a coolwarm color scale. Blue indicates similarity and red indicates dissimilarity.

## Results Summary
The tool was tested on artificial sequences,orthologs and real genomic data. Across all three metrics:

- Euclidean: highlights compositional differences well; smaller absolute values due to the square root operation
- Manhattan: clearly separates each sequence from others; diagonal is distinctly blue (self-similarity) with red off-diagonal blocks
- Cosine: best for cross-species comparison — closely related sequences cluster in blue while sequences from different organisms appear in red; particularly effective for comparing orthologs across species

### Cosine sample results
![Cosine samples results](cos_samples/heatmap.png)

### Cosine orthologs results
![Cosine orthologs results](cos_ortho/heatmap.png)

### Cosine genome results
![Cosine genomes results](cos_genome/heatmap.png)

## Limitations & Future Improvements

- Performance on large sequences — storing all k-mer vectors in memory can be slow or intensive for very large genomes
- Additional distance metrics — future versions could support Bray-Curtis dissimilarity, Jensen-Shannon divergence, or Jaccard similarity
- Phylogenetic tree output — the distance matrix could be used as input to build a neighbour-joining or UPGMA phylogenetic tree
- Ambiguous base handling — sequences containing characters outside A/T/G/C (e.g. N) are not currently filtered


## Challenges

- Handling large FASTA files efficiently
- Implementing all three distance metrics manually without external libraries
- Choosing appropriate data structures for k-mer storage and vector alignment
- Managing heatmap label formatting for readability


### Author

**Nikhitha Vujjini**
Ms. Bioinformatics.
