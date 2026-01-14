# Bovine IgG Repertoire Analysis Pipeline

A Nextflow pipeline for analyzing bovine immunoglobulin sequences from Oxford Nanopore amplicon data.

## Overview

This pipeline processes heavy and light chain IgG amplicons from cow PBMCs, clusters reads into consensus sequences, and optionally annotates them with V/D/J gene assignments and CDR3 extraction.

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04+)
- [Docker](https://www.docker.com/) or [Singularity/Apptainer](https://apptainer.org/)

## Quick Start

```bash
nextflow run . \
  --fastq_dir /path/to/fastq_pass \
  --heavy_barcode barcode01 \
  --light_barcode barcode02
```

## Input Data

The pipeline expects Oxford Nanopore basecalled FASTQ files organized by barcode:

```
fastq_pass/
├── barcode01/    # Heavy chain reads
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
└── barcode02/    # Light chain reads
    ├── file1.fastq.gz
    └── file2.fastq.gz
```

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
| `--heavy_barcode` | `barcode01` | Name of heavy chain barcode directory |
| `--light_barcode` | `barcode02` | Name of light chain barcode directory |
| `--min_len` | `400` | Minimum amplicon length |
| `--max_len` | `800` | Maximum amplicon length |
| `--min_qual` | `10` | Minimum quality score |
| `--min_reads` | `100` | Minimum reads per sample |
| `--similar_consensus` | `98` | Clustering threshold for merging groups |
| `--skip_annotation` | `false` | Skip IgBLAST annotation |
| `--results` | `./results` | Output directory |

## Output

```
results/
├── 1_merged_reads/           # Merged reads per chain
├── 2_filtered_reads/         # Quality-filtered and trimmed reads
├── 3_consensus_sequences/    # Clustered consensus sequences
├── 4_annotations/            # V/D/J annotations and CDR3 sequences
└── 5_reports/                # Summary statistics
```

## Pipeline Steps

1. **Merge Reads** - Concatenate reads from each barcode
2. **Quality Filter** - Filter by length and quality
3. **Trim Primers** - Remove primer and adapter sequences
4. **Cluster Reads** - Cluster similar sequences using amplicon_sorter
5. **Annotate** (optional) - V/D/J gene assignment with IgBLAST
6. **Parse CDR3** (optional) - Extract CDR3 sequences
7. **Report** - Generate summary statistics

## Bovine IgG Considerations

This pipeline is specifically tuned for bovine immunoglobulins:

- **Ultralong CDR3H**: Cattle can have CDR3H up to 70+ amino acids (vs ~15 for humans)
- **Limited V gene usage**: Cattle primarily use IGHV1-7 for heavy chains
- **Clustering parameters**: Tuned to distinguish true variants from ONT sequencing errors

## License

MIT License
