# Q47 Deep-Space Quadruplet Census

**A Complete Census of Prime Quadruplets for Q(n) = n<sup>47</sup> − (n−1)<sup>47</sup> in the Range 1 ≤ n ≤ 2×10<sup>11</sup>**

*Part II of the Titan Project*

Ruqing Chen — GUT Geoservice Inc., Montréal, Canada

---

## Summary

This repository contains the data, algorithms, and paper for a large-scale computational census of **prime quadruplets** — runs of four consecutive integers *n, n+1, n+2, n+3* for which all four values of Q(n) = n<sup>47</sup> − (n−1)<sup>47</sup> are probable primes of 340–520 digits (magnitudes up to ~10<sup>520</sup>).

### Key Results

| Quantity | Value |
|----------|-------|
| Search range | 1 ≤ n ≤ 2×10<sup>11</sup> |
| **Prime quadruplets** | **742** |
| **Prime quintuplets** | **7** |
| Sextuplets | 0 |
| Satellite primes (radius 5000) | 12,983 |
| Near-twin pairs (P, P−2) | 7 |
| Inferred singular series 𝔖₄ | **≈ 6,400** |
| Computational cost | ~240 CPU-core-hours |

### Highlights

- **Bateman–Horn validation**: The cumulative count *C(N)* matches the Bateman–Horn prediction via proper numerical integration (not merely the leading asymptotic term), yielding 𝔖₄ ≈ 6,385.
- **Moduli space clustering**: The tightest quadruplet pair (Δn = 146,312 at n ≈ 8.25×10<sup>10</sup>) demonstrates local alignment in the mod-q residue space.
- **Mod-3 symmetric phase shift**: Q(n) ≡ 1 (mod 3) universally, eliminating 1/3 of even offsets in *both* left and right satellite spectra — a symmetric constraint, not a directional chirality.
- **Conductor rigidity**: No anomalous clusters (k ≥ 6) exceed the Hardy–Littlewood statistical envelope across the full 200-billion range.

---

## Repository Structure

```
Q47-Deep-Space-Quadruplet-Census/
├── README.md
├── LICENSE
├── paper/
│   ├── Q47_Complete_Census.tex        # LaTeX source (13 pages)
│   ├── Q47_Complete_Census.pdf        # Compiled paper
│   └── figures/
│       ├── fig1_distribution.png      # Spatial distribution of 742 quadruplets
│       ├── fig2_cumulative.png        # Cumulative count C(N) with BH fit
│       ├── fig3_density.png           # Density with Poisson error bars
│       └── fig4_satellites.png        # Satellite prime gap distribution
├── data/
│   ├── quadruplets_742.csv            # Complete catalog: all 742 starting n
│   ├── quintuplets_7.csv              # All 7 quintuplet starting n
│   ├── pioneer_zone_15.csv            # The 15 pioneer-zone quadruplets
│   ├── satellite_primes.csv           # Left-spectrum satellite primes
│   ├── morphology_pioneer.csv         # Morphological census (n ≤ 2×10⁹)
│   ├── distribution_by_range.csv      # Quadruplet counts by range bin
│   └── closest_pairs_top10.csv        # 10 tightest quadruplet pairs
├── scripts/
│   ├── titan_sweeper.py               # Quadruplet search algorithm
│   ├── titan_radar_ultimate_5000.py   # Satellite prime radar (radius 5000)
│   └── generate_figures.py            # Reproduce all paper figures
└── .gitignore
```

---

## Data Files

### `quadruplets_742.csv`
Complete catalog of all 742 prime quadruplets. Columns: `index, starting_n, approx_digits, n_over_1e9`.

### `quintuplets_7.csv`
All 7 prime quintuplets (5 consecutive prime-generating integers). Each is a pair of overlapping quadruplets.

### `satellite_primes.csv`
Left-spectrum satellite primes within radius 5000 of each quadruplet member. Columns: `main_star_n, gap_k` (meaning P − k is prime).

### `morphology_pioneer.csv`
Complete morphological hierarchy for the pioneer zone (n ≤ 2×10⁹) from Part I: solitary, pair, triplet, and quadruplet counts.

---

## Reproducing the Figures

```bash
cd scripts
pip install matplotlib numpy scipy
python generate_figures.py
```

This regenerates all four figures from the CSV data files.

---

## Running the Search

### Quadruplet search
```bash
pip install gmpy2
python scripts/titan_sweeper.py 1000000 2000000000 --output results.csv
```

### Satellite radar
```bash
python scripts/titan_radar_ultimate_5000.py
```

---

## Compiling the Paper

```bash
cd paper
pdflatex Q47_Complete_Census.tex
pdflatex Q47_Complete_Census.tex   # second pass for cross-references
```

---

## Relation to Part I

This work (Part II) extends the foundational morphological census established in:

> R. Chen, "Statistical Morphology and Geodesic Rigidity of Prime Constellations in Q(n) = n⁴⁷ − (n−1)⁴⁷: A Complete Census of the First 2 Billion Cases," GUT Geoservice Inc., February 2026.
> [GitHub: Ruqing1963/Q47-Prime-Constellation-Census](https://github.com/Ruqing1963/Q47-Prime-Constellation-Census)

Part I established the ground-state statistics (18.47M primes, morphological hierarchy, geodesic rigidity, mod-3 exclusion principle) for n ≤ 2×10⁹. Part II extends the quadruplet search 100× further to n = 2×10¹¹.

---

## Citation

If you use this data or code, please cite:

```bibtex
@misc{Chen2026b,
  author = {Ruqing Chen},
  title  = {A Complete Census of Prime Quadruplets for
            {$Q(n) = n^{47} - (n-1)^{47}$}
            in the Range {$1 \leq n \leq 2 \times 10^{11}$}},
  year   = {2026},
  note   = {GUT Geoservice Inc., Montr\'eal},
  url    = {https://github.com/Ruqing1963/Q47-Deep-Space-Quadruplet-Census}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
