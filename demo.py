# Self-contained demo: simulate biased CDS, compute RSCU and CAI, save plots.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from codon_bias import (
    count_codons, compute_rscu, reference_weights, compute_cai,
    synonymous_groups, CODON_TABLE,
)


def preferred_codons(rng):
    # Pick one optimal codon per amino acid to simulate translational preference.
    groups = synonymous_groups()
    pref = {}
    for aa, codons in groups.items():
        pref[aa] = codons[rng.integers(len(codons))]
    return pref


def simulate_cds(rng, pref, bias_strength, n_codons):
    # Build a coding sequence drawn from per-amino-acid codon distributions.
    groups = synonymous_groups()
    aas = [aa for aa in groups if len(groups[aa]) > 1]
    seq = ['ATG']
    for _ in range(n_codons):
        aa = aas[rng.integers(len(aas))]
        codons = groups[aa]
        probs = np.ones(len(codons))
        for i, c in enumerate(codons):
            if c == pref[aa]:
                probs[i] += bias_strength
        probs = probs / probs.sum()
        seq.append(codons[rng.choice(len(codons), p=probs)])
    return ''.join(seq)


def main():
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    rng = np.random.default_rng(42)

    pref = preferred_codons(rng)

    # Highly expressed genes: strong bias toward preferred codons.
    high_seqs = [simulate_cds(rng, pref, bias_strength=6.0, n_codons=200) for _ in range(40)]
    # Background genes: weak bias.
    bg_seqs = [simulate_cds(rng, pref, bias_strength=0.3, n_codons=200) for _ in range(40)]

    # Reference weights from the highly expressed set.
    ref_counts = {}
    for seq in high_seqs:
        for codon, c in count_codons(seq).items():
            ref_counts[codon] = ref_counts.get(codon, 0) + c
    weights = reference_weights(ref_counts)

    high_cai = [compute_cai(count_codons(s), weights) for s in high_seqs]
    bg_cai = [compute_cai(count_codons(s), weights) for s in bg_seqs]

    # Aggregate RSCU for each group for the heatmap.
    def group_rscu(seqs):
        agg = {}
        for seq in seqs:
            for codon, c in count_codons(seq).items():
                agg[codon] = agg.get(codon, 0) + c
        return compute_rscu(agg)

    high_rscu = group_rscu(high_seqs)
    bg_rscu = group_rscu(bg_seqs)

    groups = synonymous_groups()
    codon_order = []
    for aa in sorted(groups):
        if len(groups[aa]) > 1:
            codon_order.extend(sorted(groups[aa]))

    heat = np.array([
        [high_rscu.get(c, 0.0) for c in codon_order],
        [bg_rscu.get(c, 0.0) for c in codon_order],
    ])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    im = axes[0].imshow(heat, aspect='auto', cmap='viridis', vmin=0, vmax=heat.max())
    axes[0].set_yticks([0, 1])
    axes[0].set_yticklabels(['highly expressed', 'background'])
    axes[0].set_xticks(range(len(codon_order)))
    axes[0].set_xticklabels(codon_order, rotation=90, fontsize=5)
    axes[0].set_title('RSCU by codon')
    fig.colorbar(im, ax=axes[0], fraction=0.025, label='RSCU')

    axes[1].hist(bg_cai, bins=15, alpha=0.6, label='background', color='grey')
    axes[1].hist(high_cai, bins=15, alpha=0.6, label='highly expressed', color='darkorange')
    axes[1].set_xlabel('CAI')
    axes[1].set_ylabel('gene count')
    axes[1].set_title('CAI distribution')
    axes[1].legend()

    fig.tight_layout()
    fig.savefig('figures/demo.png', dpi=120)

    summary = pd.DataFrame({
        'group': ['highly_expressed', 'background'],
        'n_genes': [len(high_cai), len(bg_cai)],
        'mean_cai': [np.mean(high_cai), np.mean(bg_cai)],
        'median_cai': [np.median(high_cai), np.median(bg_cai)],
        'sd_cai': [np.std(high_cai), np.std(bg_cai)],
    })
    summary.to_csv('results/summary.csv', index=False)
    print(summary)
    print('Saved figures/demo.png and results/summary.csv')


if __name__ == '__main__':
    main()
