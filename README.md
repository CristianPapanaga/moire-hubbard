# moire-hubbard

Density-matrix renormalization group (DMRG) study of the **extended moiré Hubbard model**
for twisted homobilayers of the transition-metal dichalcogenide WSe₂.

This repository is the companion code to the BSc thesis

> *"Exotic physics in 2d materials: computational study of twisted homobilayers of the
> transition-metal dichalcogenide WSe₂"* — Cristian Papanaga & Ciarán Hickey,
> University College Dublin (2024).

The code uses the [ITensor](https://itensor.org/) library (Julia) to build the Hamiltonian as a
matrix product operator (MPO), find its ground state with DMRG, and then extract and plot
observables such as the electron density, real-space occupation, spin correlations, and
charge/magnetic structure factors.

---

## Physics background

Twisted homobilayers of WSe₂ are well-described by a Hubbard model defined on a moiré
superlattice. The key ingredient that distinguishes the *moiré* Hubbard model from the standard
Hubbard model is a **spin-dependent hopping phase** originating from the interlayer displacement
field. It breaks inversion symmetry while preserving time-reversal symmetry.

### The model

The code solves the extended moiré Hubbard Hamiltonian

$$H = -t\sum_{\langle i,j\rangle}\left[c^\dagger_{i\uparrow}c_{j\uparrow}e^{i\phi} + c^\dagger_{i\downarrow}c_{j\downarrow}e^{-i\phi} + \mathrm{h.c.}\right] + U_0\sum_i n_{i\uparrow}n_{i\downarrow} + U_1\sum_{\langle i,j\rangle} n_i n_j - \mu\sum_i n_i$$

on a **quasi-1D triangular ladder** with open boundary conditions, where

| Symbol | Meaning |
|--------|---------|
| `t`    | nearest-neighbour hopping constant |
| `phi` (`φ`) | spin-dependent hopping phase (tunable via the displacement field) |
| `U_0`  | on-site (`n↑ n↓`) density–density interaction strength |
| `U_1`  | nearest-neighbour density–density interaction strength |
| `mu` (`μ`) | chemical potential (controls the electron density) |
| `L`    | number of unit cells along the ladder (`N = 2L` sites total) |

The triangular ladder is mapped to a 1D chain by translating one sublattice along the transverse
direction, giving on-site, nearest-neighbour and **next-nearest-neighbour** terms — precisely the
`j` and `j+2` terms seen in `src/moire_hubbard.jl`. This 1D geometry is optimal for DMRG.

### Key results

The thesis reports, from DMRG simulations of the interacting model:

- A pronounced **Mott insulating state at half-filling** (ρ = 1), developing for interaction
  strengths `U ≈ 2`–`3` and driven by strong electron correlations.
- Indications of **generalised Wigner crystal (GWC)** states at `1/8` and `3/8` fillings for a
  stronger nearest-neighbour scaling `α = U₁/U₀ = 0.4`.
- An anomalous plateau at `13/24` filling attributed to a finite-size effect of the open
  boundary conditions (two electrons localised in the ladder corners).

---

## Requirements

- [Julia](https://julialang.org/) ≥ 1.10 (developed and tested with Julia 1.12)
- The Julia packages listed in `Project.toml`:
  [`ITensors`](https://github.com/ITensor/ITensors.jl),
  [`ITensorMPS`](https://github.com/ITensor/ITensorMPS.jl),
  [`HDF5`](https://github.com/JuliaIO/HDF5.jl),
  [`Plots`](https://github.com/JuliaPlots/Plots.jl),
  [`PGFPlotsX`](https://github.com/JuliaPlots/PGFPlotsX.jl),
  [`LaTeXStrings`](https://github.com/stevengj/LaTeXStrings.jl),
  [`FFTW`](https://github.com/JuliaMath/FFTW.jl),
  [`NaturalSort`](https://github.com/JuliaString/NaturalSort.jl), and `LinearAlgebra`.

## Installation

```bash
git clone https://github.com/CristianPapanaga/moire-hubbard.git
cd moire-hubbard
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## Repository structure

```
moire-hubbard/
├── Project.toml            # Package environment / dependencies
├── Manifest.toml           # Resolved dependency versions
├── src/
│   ├── moire_hubbard.jl    # Builds the MPO Hamiltonian and runs DMRG
│   ├── utils.jl            # Data loading, lattice construction & plotting helpers
│   ├── testing.jl          # Driver: run DMRG and save ground states to HDF5
│   └── data_anaylsis.jl    # Driver: load HDF5 data and produce plots
├── output/                 # HDF5 files containing the MPS ground states (git-ignored)
├── plots/                  # Generated PDF figures (git-ignored)
└── Thesis.pdf              # The thesis this code accompanies
```

### `src/moire_hubbard.jl`

Defines the core function

```julia
moire_hubbard(L, nsweeps, maxdim, cutoff, t, U_0, U_1, mu, phi)
```

which

1. constructs the `Electron` site indices (`N = 2L` sites),
2. assembles the `μ`-independent part of the Hamiltonian `H₀` as an `OpSum` (nearest- and
   next-nearest-neighbour hopping with the spin-dependent phase `φ`, plus the `U₀` and `U₁`
   interactions),
3. builds the chemical-potential term `H_μ`, and
4. runs DMRG **in parallel over the chemical-potential values** `mu` using `Threads.@threads`
   (each thread writes to its own pre-allocated slot, so the loop is thread-safe).

It returns the ground-state energies and the truncated ground-state MPS for every `μ`.

> Run with multiple threads for a speed-up: `julia --threads=N`.

### `src/utils.jl`

Helper functions for post-processing:

| Function | Purpose |
|----------|---------|
| `lattice_coords_2d(a, b, site_unit_cell_coords, N_x, N_y)` | Generate 2D lattice coordinates from lattice vectors |
| `load_data(file_name)` | Load MPS ground states + parameters from an HDF5 file (parses the file name) |
| `plot_density_vs_mu(...)` | Electron density ρ vs chemical potential μ |
| `plot_density_vs_mu_all_U(...)` | ρ vs μ overlaid for all `U` values in a folder |
| `plot_lattice_occupation(...)` | Real-space occupation expectation value on the lattice |
| `plot_spin_correlations(...)` | Sᶻ–Sᶻ spin correlation map |
| `plot_charge_structure_factor(...)` | Charge structure factor (1D Fourier transform of `Ntot`) |
| `plot_magnetic_structure_factor(...)` | Magnetic structure factor (1D Fourier transform of the spin correlations) |

### `src/testing.jl`

End-to-end driver. It sets the DMRG parameters (number of sweeps, maximum bond dimension,
truncation cutoff), the system parameters (`L`, `t`, `μ`, `U`, `φ`), runs `moire_hubbard`, and
saves the resulting ground states to an HDF5 file in `output/`.

### `src/data_anaylsis.jl`

Loads a previously saved HDF5 file with `load_data` and generates all the figures (density vs μ,
real-space occupation, spin correlations, and the charge/magnetic structure factors).

---

## Usage

### 1. Run DMRG and save the ground states

Edit the `DMRG PARAMETERS` and `SYSTEM PARAMETERS` blocks at the top of `src/testing.jl`, then:

```bash
julia --project=. --threads=4 src/testing.jl
```

This writes a file such as

```
output/(0.1) L = [12]  t = [1]  U_0 = [0.1]  U_1 = [0.05]  phi = [0.5235987755982988].h5
```

### 2. Analyse the data and make plots

Point `src/data_anaylsis.jl` at the HDF5 file you just created, then:

```bash
julia --project=. src/data_anaylsis.jl
```

Figures are saved as PDFs under `plots/<parameters>/`.

---

## File-naming convention

Data and plot file names encode the simulation parameters in a machine-readable way:

```
(<sort-key>) L = [<L>]  t = [<t>]  U_0 = [<U_0>]  U_1 = [<U_1>]  phi = [<phi>].h5
```

- The `(<sort-key>)` prefix (e.g. the interaction strength `U`) is used only for natural-sort
  ordering, because HDF5 returns dataset names alphabetically.
- Values inside `[...]` are parsed back by `load_data` via a regex, so parameters are recovered
  directly from the file name without any sidecar metadata.
- Each `μ` value is stored as a separate dataset named `(i) mu = [<mu_i>]`.

---

## References

- C. Papanaga & C. Hickey, *Exotic physics in 2d materials: computational study of twisted
  homobilayers of the transition-metal dichalcogenide WSe₂*, BSc thesis, University College
  Dublin (2024) — see `Thesis.pdf`.
- ITensor: M. Fishman, S. R. White, E. M. Stoudenmire, *The ITensor Software Library for Tensor
  Network Calculations*, SciPost Phys. Codebases **4** (2022).

