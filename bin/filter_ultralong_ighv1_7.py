#!/usr/bin/env python3
"""
Filter AIRR-format sequences for IGHV1-7 gene family with ultra-long CDR3H3.

Takes an AIRR-format TSV produced by Change-O MakeDb (from IgBLAST annotation)
and filters for:
  - V gene assignment in the IGHV1-7 family
  - Junction length >= min_cdr3_aa amino acids (CDR3 + 2 anchor residues)

Outputs a filtered TSV and FASTA with V gene and CDR3 length in the headers.

Usage:
    filter_ultralong_ighv1_7.py --db makedb_output.tsv --min_cdr3_aa 50
    filter_ultralong_ighv1_7.py --db makedb_output.tsv --fasta consensus.fasta --min_cdr3_aa 50
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_fasta(fasta_file):
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    return sequences


def cdr3_aa_length(row):
    # Prefer junction_aa (direct amino acid junction sequence)
    junction_aa = row.get('junction_aa', '')
    if junction_aa and junction_aa not in ('NA', '', 'None'):
        return len(junction_aa)
    # Fall back to junction_length in nucleotides (includes 2 anchor codons = 6 nt)
    jl = row.get('junction_length', '') or ''
    try:
        return int(jl) // 3
    except ValueError:
        return 0


def main():
    parser = argparse.ArgumentParser(
        description='Filter AIRR sequences for IGHV1-7 + ultra-long CDR3H3'
    )
    parser.add_argument('--db', required=True,
                        help='AIRR-format TSV from Change-O MakeDb')
    parser.add_argument('--fasta',
                        help='Optional FASTA with input sequences (if not embedded in TSV)')
    parser.add_argument('--min_cdr3_aa', type=int, default=50,
                        help='Minimum CDR3 length in amino acids, default 50')
    parser.add_argument('--out_tsv',
                        help='Output filtered TSV (default: <db>_ighv1_7_ultralong.tsv)')
    parser.add_argument('--out_fasta',
                        help='Output filtered FASTA (default: <db>_ighv1_7_ultralong.fasta)')
    args = parser.parse_args()

    db_path = Path(args.db)
    out_tsv = args.out_tsv or str(db_path.with_suffix('')) + '_ighv1_7_ultralong.tsv'
    out_fasta = args.out_fasta or str(db_path.with_suffix('')) + '_ighv1_7_ultralong.fasta'

    # Load external FASTA if provided
    fasta_seqs = {}
    if args.fasta:
        fasta_seqs = parse_fasta(args.fasta)

    # Filter AIRR TSV
    kept = []
    total = 0
    fieldnames = []

    with open(args.db) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames or []
        for row in reader:
            total += 1
            v_call = row.get('v_call', '') or ''

            is_ighv1_7 = 'IGHV1-7' in v_call
            is_ultralong = cdr3_aa_length(row) >= args.min_cdr3_aa

            if is_ighv1_7 and is_ultralong:
                kept.append(row)

    print(f"Total sequences:                        {total}")
    print(f"IGHV1-7 + CDR3 >= {args.min_cdr3_aa} aa:              {len(kept)}")

    if not kept:
        print("No sequences passed the filter.")
        sys.exit(0)

    # Write filtered TSV
    with open(out_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(kept)
    print(f"Filtered TSV  -> {out_tsv}")

    # Write filtered FASTA
    # Use external FASTA if provided, otherwise fall back to 'sequence' column in AIRR TSV
    written = 0
    with open(out_fasta, 'w') as f:
        for row in kept:
            seq_id = row.get('sequence_id', '')
            seq = fasta_seqs.get(seq_id, '') or row.get('sequence', '') or ''
            if not seq:
                continue
            v_call = row.get('v_call', 'NA')
            cdr3_len = cdr3_aa_length(row)
            junction_aa = row.get('junction_aa', '')
            header = f">{seq_id} v_call={v_call} cdr3_aa={cdr3_len}"
            if junction_aa and junction_aa not in ('NA', '', 'None'):
                header += f" junction_aa={junction_aa}"
            f.write(header + '\n')
            f.write(seq + '\n')
            written += 1

    if written:
        print(f"Filtered FASTA -> {out_fasta} ({written} sequences)")
    else:
        print("Warning: no sequences written to FASTA — sequence column not found in TSV "
              "and no --fasta provided.")


if __name__ == '__main__':
    main()
