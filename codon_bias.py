# Compute codon usage frequencies, RSCU, and CAI from CDS FASTA files.
import argparse
import math
from collections import defaultdict

BASES = ['T', 'C', 'A', 'G']
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def read_fasta(path):
    # Yield (header, sequence) tuples from a FASTA file.
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.upper().replace('U', 'T'))
        if header is not None:
            yield header, ''.join(seq)


def count_codons(seq):
    # Count codons in an in-frame sequence, ignoring incomplete trailing bases.
    counts = defaultdict(int)
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i + 3]
        if codon in CODON_TABLE:
            counts[codon] += 1
    return counts


def synonymous_groups():
    # Map amino acid -> list of codons, excluding stop codons.
    groups = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        if aa != '*':
            groups[aa].append(codon)
    return groups


def compute_rscu(counts):
    # RSCU per codon: observed count scaled by synonymous family mean.
    groups = synonymous_groups()
    rscu = {}
    for aa, codons in groups.items():
        total = sum(counts.get(c, 0) for c in codons)
        n = len(codons)
        for c in codons:
            if total == 0:
                rscu[c] = 0.0
            else:
                rscu[c] = counts.get(c, 0) * n / total
    return rscu


def reference_weights(ref_counts):
    # Relative adaptiveness w = RSCU / max(RSCU) within each synonymous family.
    groups = synonymous_groups()
    rscu = compute_rscu(ref_counts)
    weights = {}
    for aa, codons in groups.items():
        max_rscu = max(rscu.get(c, 0.0) for c in codons)
        for c in codons:
            if max_rscu == 0:
                weights[c] = 0.01
            else:
                weights[c] = max(rscu.get(c, 0.0) / max_rscu, 0.01)
    return weights


def compute_cai(counts, weights):
    # CAI = geometric mean of codon weights, excluding Met, Trp and stops.
    groups = synonymous_groups()
    multi = set()
    for aa, codons in groups.items():
        if len(codons) > 1:
            multi.update(codons)
    log_sum = 0.0
    n = 0
    for codon, c in counts.items():
        if codon in multi and codon in weights:
            log_sum += c * math.log(weights[codon])
            n += c
    if n == 0:
        return 0.0
    return math.exp(log_sum / n)


def main():
    ap = argparse.ArgumentParser(description='Codon usage bias and CAI from CDS FASTA.')
    ap.add_argument('--fasta', required=True, help='CDS FASTA to score')
    ap.add_argument('--reference', required=True, help='Highly expressed reference FASTA')
    ap.add_argument('--out', default='cai_scores.csv', help='Output CSV path')
    args = ap.parse_args()

    ref_counts = defaultdict(int)
    for _, seq in read_fasta(args.reference):
        for codon, c in count_codons(seq).items():
            ref_counts[codon] += c
    weights = reference_weights(ref_counts)

    with open(args.out, 'w') as out:
        out.write('gene,n_codons,cai\n')
        for header, seq in read_fasta(args.fasta):
            counts = count_codons(seq)
            cai = compute_cai(counts, weights)
            out.write('{},{},{:.4f}\n'.format(header, sum(counts.values()), cai))
    print('Wrote CAI scores to', args.out)


if __name__ == '__main__':
    main()
