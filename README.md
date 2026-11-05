# Filled-in peeling explorations exhaust every critical Boltzmann planar map

Yitzchak Shmalo

This repository contains the paper and the computations for the following result. Let **q** be an admissible critical weight sequence for bipartite Boltzmann planar maps. For every finite boundary length and every peeling algorithm that never stops, the filled-in peeling exploration of the infinite **q**-Boltzmann map with a boundary discovers the whole map almost surely.

Curien proved this for weight sequences of type a in [3/2, 5/2] (Corollary 27 of his Saint-Flour notes, *Peeling Random Planar Maps*) and asked whether it holds for every critical weight sequence. The answer is yes. The proof replaces the tail asymptotics and the scaling limit of the perimeter process by an exact formula for the negative jumps in terms of the positive ones,

    nu(-k) = d(k)/(2k-1) * sum_{m>=1} m(2m+1) d(m) nu(m) / (m+k),   d(k) = 4^{-k} binom(2k,k),

together with a martingale argument for log h(P_n) along the perimeter process.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/kernel_check.py`: checks the kernel formula on an explicit critical sequence with f''(Z) finite and radius of convergence exactly Z: the critical equations, the formula for k up to 1000, the total mass of nu, and the harmonicity of h. Output in `code/kernel_check_out.txt`. Needs mpmath, numpy and scipy; runs in a few minutes.
- `code/universal_checks.py`: builds transition rows of the perimeter process from the kernel formula for seven positive laws (single atoms and two mixtures) at perimeters up to 16384, and records the ratios of the one-step quantities in the proof to the rate R(p). Output in `code/universal_checks_out.txt` and `code/universal_checks.json`.
- `code/experiment.py`: quadrature check of the closed form of the disk integrals for the same sequence, their j^{-2} asymptotics, and truncated row sums. Output in `code/experiment_out.txt` and `code/experiment_results.json`.
- `lean/PeelingAlgebra.lean`: six elementary identities and inequalities used in the paper, checked with Lean 4 and Mathlib.

## Running

    cd code
    python kernel_check.py > kernel_check_out.txt
    python universal_checks.py > universal_checks_out.txt
    python experiment.py > experiment_out.txt

The Lean file compiles inside any Lean 4 project with Mathlib (`lake env lean PeelingAlgebra.lean`).

## License

Paper: CC BY 4.0. Code: MIT.
