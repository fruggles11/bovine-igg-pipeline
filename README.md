# Bovine IgG Repertoire Analysis Pipeline

A Nextflow pipeline for analyzing bovine immunoglobulin sequences from Oxford Nanopore amplicon data.

## Overview

This pipeline processes heavy and light chain IgG amplicons from cow PBMCs, clusters reads into consensus sequences, and optionally annotates them with V/D/J gene assignments and CDR3 extraction.

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04+)
- [Docker](https://www.docker.com/) or [Singularity/Apptainer](https://apptainer.org/)

## Quick Start

```bash
nextflow run fruggles11/bovine-igg-pipeline --fastq_dir /path/to/fastq_pass
```

That's it! The pipeline will automatically detect all barcode directories and classify reads as heavy or light chain based on primer sequences.

To use a specific version:
```bash
nextflow run fruggles11/bovine-igg-pipeline -r main --fastq_dir /path/to/fastq_pass
```

## Input Data

The pipeline expects Oxford Nanopore basecalled FASTQ files organized by barcode:

```
fastq_pass/
├── barcode01/
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
├── barcode02/
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
└── barcode03/    # Any number of barcodes supported
    └── ...
```

**Key features:**
- All directories matching `barcode*/` are automatically detected
- Reads are classified as heavy or light chain based on primer sequence matching
- No need to specify which barcode contains which chain type
- Unmatched reads (no primer detected) are saved separately for review

## Setting Up IgBLAST Annotation (Optional)

For V/D/J gene annotation and CDR3 extraction, you need to download bovine germline genes from IMGT:

1. Go to [IMGT/GENE-DB](https://www.imgt.org/genedb/)
2. For each gene type, select:
   - Species: "Bos taurus"
   - Group: IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, or IGLJ
   - Functionality: "functional" (recommended)
3. Click Search, select all genes, and Export as FASTA
4. Save files in `resources/germlines/`:
   - `bovine_IGHV.fasta`
   - `bovine_IGHD.fasta`
   - `bovine_IGHJ.fasta`
   - `bovine_IGKV.fasta`
   - `bovine_IGKJ.fasta`
   - `bovine_IGLV.fasta`
   - `bovine_IGLJ.fasta`

If germline files are not available, run with `--skip_annotation true` to skip the annotation step.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq_dir` | (required) | Path to directory containing barcode subdirectories |
| `--primer_table` | `resources/bovine_primers.csv` | CSV file with primer sequences for chain classification |
| `--primer_mismatch` | `2` | Allowed mismatches when matching primers (max 3) |
| `--min_len` | `400` | Minimum amplicon length |
| `--max_len` | `800` | Maximum amplicon length |
| `--min_qual` | `10` | Minimum quality score |
| `--min_reads` | `100` | Minimum reads per chain per barcode |
| `--similar_consensus` | `98` | Clustering threshold for merging groups |
| `--skip_annotation` | `false` | Skip IgBLAST annotation |
| `--results` | `./results` | Output directory |

## Output

Results are organized by barcode:

```
results/
├── 1_merged_reads/
│   └── barcode01/
│       └── barcode01_merged.fastq.gz
├── 2_classified_reads/
│   └── barcode01/
│       ├── barcode01_heavy.fastq.gz      # Heavy chain reads
│       ├── barcode01_light.fastq.gz      # Light chain reads
│       └── barcode01_unmatched.fastq.gz  # Reads without primer match
├── 3_filtered_reads/
│   └── barcode01/
│       ├── barcode01_heavy_filtered.fastq.gz
│       └── barcode01_light_filtered.fastq.gz
├── 4_consensus_sequences/
│   └── barcode01/
│       ├── barcode01_heavy_consensus/
│       └── barcode01_light_consensus/
├── 5_annotations/
│   └── barcode01/
│       ├── barcode01_heavy_annotations.tsv
│       └── barcode01_light_annotations.tsv
└── 6_reports/
    └── summary_stats.tsv
```

## Pipeline Steps

1. **Merge Reads** - Concatenate reads from each barcode
2. **Classify by Primer** - Sort reads into heavy/light chain based on primer sequences
3. **Quality Filter** - Filter by length and quality
4. **Trim Primers** - Remove primer and adapter sequences
5. **Cluster Reads** - Cluster similar sequences using amplicon_sorter
6. **Annotate** (optional) - V/D/J gene assignment with IgBLAST
7. **Parse CDR3** (optional) - Extract CDR3 sequences
8. **Report** - Generate summary statistics

## Bovine IgG Considerations

This pipeline is specifically tuned for bovine immunoglobulins:

- **Ultralong CDR3H**: Cattle can have CDR3H up to 70+ amino acids (vs ~15 for humans)
- **Limited V gene usage**: Cattle primarily use IGHV1-7 for heavy chains
- **Clustering parameters**: Tuned to distinguish true variants from ONT sequencing errors

## License

MIT License
