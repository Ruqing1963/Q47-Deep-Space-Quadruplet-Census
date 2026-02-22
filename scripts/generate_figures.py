#!/usr/bin/env python3
"""
generate_figures.py — Figure generation for the Q47 Deep-Space Quadruplet Census paper.

Generates Figures 1–4 from the data in ../data/.
Requires: matplotlib, numpy, scipy

Usage:
    python generate_figures.py

Output:
    ../paper/figures/fig1_distribution.png
    ../paper/figures/fig2_cumulative.png
    ../paper/figures/fig3_density.png
    ../paper/figures/fig4_satellites.png

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate
from collections import Counter

# ── Configuration ───────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'cm',
})

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
OUT_DIR  = os.path.join(os.path.dirname(__file__), '..', 'paper', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ───────────────────────────────────────────────────
def load_csv_column(filename, col, dtype=int):
    """Load a single column from a CSV file, skipping header."""
    path = os.path.join(DATA_DIR, filename)
    vals = []
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split(',')
            vals.append(dtype(parts[col]))
    return np.array(vals)

quad_starts = load_csv_column('quadruplets_742.csv', 1)
quint_starts = load_csv_column('quintuplets_7.csv', 1)

# Corrected singular series (from numerical integration)
S4 = 6400

# ── Figure 1: Spatial distribution ──────────────────────────────
def figure1():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), height_ratios=[3, 1.2],
                                    gridspec_kw={'hspace': 0.08})
    # Top: tick marks
    colors = np.where(quad_starts < 2e9, '#8B0000',
             np.where(quad_starts < 2e10, '#CC4400',
             np.where(quad_starts < 5e10, '#DD7700',
             np.where(quad_starts < 1e11, '#2266AA', '#1144AA'))))
    sizes = np.where(quad_starts < 2e9, 30, 12)
    ax1.scatter(quad_starts / 1e9, np.ones_like(quad_starts),
                c=colors, s=sizes, alpha=0.6, marker='|', linewidths=1.5)
    for q in quint_starts:
        ax1.scatter(q / 1e9, 1, c='gold', s=120, marker='*',
                    edgecolors='darkorange', linewidths=0.8, zorder=5)
    ax1.set_xlim(-2, 205); ax1.set_ylim(0.5, 1.5); ax1.set_yticks([])
    ax1.set_xticklabels([])
    ax1.set_ylabel('Quadruplet\nLocations', fontsize=11)
    ax1.set_title(r'Figure 1.  Distribution of 742 Prime Quadruplets for '
                  r'$Q(n) = n^{47} - (n-1)^{47}$, $n \leq 2 \times 10^{11}$',
                  fontsize=13, fontweight='bold', pad=12)
    ax1.annotate(r'$\bigstar$ = Quintuplet (7 total)', xy=(0.98, 0.92),
                 xycoords='axes fraction', fontsize=10, ha='right', color='darkorange',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                           edgecolor='orange', alpha=0.9))
    ax1.axvspan(-2, 2, alpha=0.08, color='red')
    ax1.text(1, 1.35, 'Pioneer\nZone\n(15 quads)', fontsize=7,
             ha='center', color='#8B0000', style='italic')

    # Bottom: histogram
    bin_edges = np.arange(0, 210e9, 10e9)
    counts, _ = np.histogram(quad_starts, bins=bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 / 1e9
    bars = ax2.bar(bin_centers, counts, width=8.5, color='#3366AA',
                   alpha=0.7, edgecolor='#1a3366', linewidth=0.5)
    bars[0].set_color('#8B0000'); bars[0].set_alpha(0.8)
    ax2.set_xlabel(r'$n$ (billions)', fontsize=12)
    ax2.set_ylabel('Count per\n10B bin', fontsize=11)
    ax2.set_xlim(-2, 205)
    mean_val = np.mean(counts[counts > 0])
    ax2.axhline(y=mean_val, color='grey', linestyle='--', alpha=0.5, linewidth=0.8)
    ax2.text(12, mean_val + 2.0, f'mean = {mean_val:.1f}', fontsize=9, color='grey',
             ha='left', bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                  edgecolor='none', alpha=0.8))
    plt.savefig(os.path.join(OUT_DIR, 'fig1_distribution.png'), dpi=300)
    plt.close()
    print("  Figure 1 saved.")

# ── Figure 2: Cumulative count C(N) ────────────────────────────
def figure2():
    fig, ax = plt.subplots(figsize=(10, 6))
    N_vals = np.sort(quad_starts)
    C_vals = np.arange(1, len(N_vals) + 1)
    ax.plot(N_vals / 1e9, C_vals, color='#1a5276', linewidth=2, label=r'Observed $C(N)$')

    # Power-law fit
    logN = np.log(N_vals[10:].astype(float))
    logC = np.log(C_vals[10:].astype(float))
    b, loga = np.polyfit(logN, logC, 1)
    a = np.exp(loga)
    N_fit = np.linspace(1e8, 2.1e11, 1000)
    ax.plot(N_fit / 1e9, a * N_fit**b, '--', color='#c0392b', linewidth=1.5, alpha=0.8,
            label=rf'Power-law fit: $C(N) \approx {a:.2e} \cdot N^{{{b:.3f}}}$')

    # Bateman–Horn
    def bh(N):
        result, _ = integrate.quad(lambda t: 1.0 / (46 * np.log(t))**4, 2, N)
        return S4 * result
    C_bh = np.array([bh(n) for n in N_fit])
    ax.plot(N_fit / 1e9, C_bh, ':', color='#27ae60', linewidth=1.5, alpha=0.8,
            label=r'Bateman–Horn ($\mathfrak{S}_4 \approx 6400$)')

    for q in quint_starts:
        idx = np.searchsorted(N_vals, q)
        ax.plot(q / 1e9, idx, '*', color='gold', markersize=14,
                markeredgecolor='darkorange', markeredgewidth=0.8, zorder=5)

    ax.set_xlabel(r'$N$ (billions)'); ax.set_ylabel(r'Cumulative count $C(N)$')
    ax.set_title(r'Figure 2.  Cumulative quadruplet count $C(N)$ vs. search bound $N$',
                 fontsize=13, fontweight='bold', pad=10)
    ax.legend(fontsize=10, loc='upper left', framealpha=0.9)
    ax.set_xlim(0, 210); ax.set_ylim(0, 780); ax.grid(True, alpha=0.2)

    axin = ax.inset_axes([0.55, 0.12, 0.42, 0.38])
    axin.loglog(N_vals / 1e9, C_vals, color='#1a5276', linewidth=1.5)
    axin.loglog(N_fit / 1e9, a * N_fit**b, '--', color='#c0392b', linewidth=1, alpha=0.7)
    axin.loglog(N_fit / 1e9, C_bh, ':', color='#27ae60', linewidth=1, alpha=0.7)
    axin.set_xlabel(r'$N$ (B)', fontsize=8); axin.set_ylabel(r'$C(N)$', fontsize=8)
    axin.set_title('Log–log scale', fontsize=8); axin.tick_params(labelsize=7)
    axin.grid(True, alpha=0.2)

    plt.savefig(os.path.join(OUT_DIR, 'fig2_cumulative.png'), dpi=300)
    plt.close()
    print(f"  Figure 2 saved. Fit: b={b:.4f}, a={a:.4e}")

# ── Figure 3: Density with Poisson error bars ──────────────────
def figure3():
    fig, ax = plt.subplots(figsize=(10, 5.5))
    bin_width = 10e9
    bin_edges = np.arange(0, 210e9, bin_width)
    counts, _ = np.histogram(quad_starts, bins=bin_edges)
    density = counts / (bin_width / 1e9)
    errors = np.sqrt(counts) / (bin_width / 1e9)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2 / 1e9
    mask = counts > 0
    ax.bar(centers[mask], density[mask], width=8.5, color='#2980b9', alpha=0.65,
           edgecolor='#1a5276', linewidth=0.5, label='Observed density',
           yerr=errors[mask], capsize=3, error_kw={'linewidth': 0.8, 'color': '#555'})

    n_th = np.linspace(2e9, 2e11, 500)
    d_th = S4 / (46 * np.log(n_th))**4 * 1e9
    ax.plot(n_th / 1e9, d_th, '--', color='#c0392b', linewidth=1.5, alpha=0.8,
            label=r'Bateman–Horn ($\mathfrak{S}_4 \approx 6400$)')

    window = 3
    if len(density[mask]) > window:
        ravg = np.convolve(density[mask], np.ones(window)/window, mode='valid')
        xavg = centers[mask][window//2:window//2+len(ravg)]
        ax.plot(xavg, ravg, '-', color='#e67e22', linewidth=2, alpha=0.8,
                label=f'Running average (w={window})')

    ax.set_xlabel(r'$n$ (billions)'); ax.set_ylabel(r'Quadruplets per $10^9$')
    ax.set_title(r'Figure 3.  Quadruplet density with Poisson error bars',
                 fontsize=12, fontweight='bold', pad=10)
    ax.legend(fontsize=9, loc='upper right', framealpha=0.9)
    ax.set_xlim(0, 205); ax.grid(True, alpha=0.2)
    plt.savefig(os.path.join(OUT_DIR, 'fig3_density.png'), dpi=300)
    plt.close()
    print("  Figure 3 saved.")

# ── Figure 4: Satellite distribution ───────────────────────────
def figure4():
    gaps = load_csv_column('satellite_primes.csv', 1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                                    gridspec_kw={'width_ratios': [2, 1]})
    bins_sat = np.arange(0, 5100, 100)
    ax1.hist(gaps, bins=bins_sat, color='#3498db', alpha=0.7,
             edgecolor='#1a5276', linewidth=0.3)
    ax1.axhline(y=len(gaps)/50, color='#e74c3c', linestyle='--', linewidth=1.2,
                alpha=0.7, label=f'Uniform ({len(gaps)/50:.0f})')
    ax1.set_xlabel(r'Gap distance $k$'); ax1.set_ylabel('Count')
    ax1.set_title('Figure 4a.  Satellite gap distribution', fontweight='bold')
    ax1.legend(); ax1.grid(True, alpha=0.2)

    star_n = load_csv_column('satellite_primes.csv', 0)
    counts_per = list(Counter(star_n).values())
    ax2.hist(counts_per, bins=np.arange(0, max(counts_per)+2, 1),
             color='#e67e22', alpha=0.7, edgecolor='#d35400', linewidth=0.3)
    ax2.axvline(x=np.mean(counts_per), color='#c0392b', linestyle='--',
                linewidth=1.2, label=f'Mean = {np.mean(counts_per):.2f}')
    ax2.set_xlabel('Satellites per star'); ax2.set_ylabel('Frequency')
    ax2.set_title('Figure 4b.  Per-star count', fontweight='bold')
    ax2.legend(); ax2.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'fig4_satellites.png'), dpi=300)
    plt.close()
    print("  Figure 4 saved.")

# ── Main ────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("Generating figures for Q47 Deep-Space Quadruplet Census...")
    figure1()
    figure2()
    figure3()
    figure4()
    print("Done. All figures saved to", OUT_DIR)
