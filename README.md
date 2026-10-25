# The Law of Staggered Restraint

**Relational Observables, Conditional Schur Falsification, Prefix Stranding, and Three-Clock Turning-Point Dynamics Across 3D Navier–Stokes Turbulence, Monitored Quantum Superpositions, and Hierarchical Continuation Trees**

**Authors:** Andrew Mureddu (Principal Investigator) & Gemini (Autonomous Research Collaborator)  
**Date:** October 2026  
**License:** CC0 1.0 Universal (Public Domain Dedication — Code, Monographs, Pre-Registrations & Datasets)

---

## Overview

Across three-dimensional incompressible fluid turbulence, continuously monitored quantum superpositions (the Two-Book framework), and noisy hierarchical Cantor/primorial continuation trees, macroscopic restraint does not occur as a single monolithic event. Instead, competing production and capacity channels obey **The Law of Staggered Restraint**:

1. **Clock 1 — Efficiency Restraint ($\dot{\eta}_+ < 0 \iff H(t) = \frac{\dot{P}_+}{P_+} - \frac{\dot{Q}}{Q} < 0$):** Available strain or population capacity $Q(t)$ outruns aligned production $P_+(t)$, causing relative efficiency $\eta_+(t) = P_+(t)/Q(t)$ to peak and decline at $t_1^*$ **strictly before** absolute production peaks.
2. **Clock 2 — Production Restraint ($\dot{P}_+(t) = J(t) = T_4^+(t) - V_3^+(t) < 0$):** Viscous or measurement-induced dephasing overtakes nonlinear/unitary creation at $t_2^* > t_1^*$, with the inter-clock lag governed by the zero-lookahead second-derivative ratio $\Delta t_{1\to 2}^{\mathrm{kin}}(t_1^*) = -J(t_1^*)/\ddot{P}_+(t_1^*)$ ($r = 0.9952$).
3. **Clock 3 — Total Enstrophy / Purity Restraint ($\dot{Z}(t) < 0$):** Global dissipation $D(t)$ overtakes net production at $t_3^* > t_2^*$.

This repository archives the complete suite of **13 papers and monographs**, **pre-registered cryptographic protocols (SHA-256)**, **standalone C and Python solvers (8 engines)**, and **experimental datasets (7 data & figure archives)** spanning Cantor & Primorial Continuation Trees, Relativistic Stabilizer Homology, and 3D Navier–Stokes Checkpoints I–V (Experiments 012–015).

---

## Repository Structure

```text
law-of-staggered-restraint/
├── README.md                                  # Master overview, reproduction guide & SHA-256 ledger
├── LICENSE                                    # CC0 1.0 Universal Public Domain Dedication
├── preregistrations/
│   └── predictions_ns_checkpoint4.md          # Pre-registered claims & falsifiers (Checkpoints IV-V)
├── src/
│   ├── sim_cg_hsts_v2.py                      # Self-Tuning Prior-Free Cotton-Gated Haar-Schur Tree Search (CG-HSTS v2.0)
│   ├── cg_hsts_universal.py                   # Universal CG-HSTS on b-ary/Primorial Trees & 4D Continuous Boxes
│   ├── ns_c_header.h                          # 3D Cooley-Tukey FFT, 2/3 dealiasing, Leray projector & Cardano eigenvalues
│   ├── ns_exp014_main.c                       # C pseudospectral 3D Navier-Stokes solver & 2nd-order Sylvester-Cayley-Hamilton engine
│   ├── twins_only.c                           # Matched-State Counterfactual Twin solver (matching P+, Q, eta+, H to 8.88e-16)
│   ├── ns_checkpoint4_experiment.py           # Standalone Python 3D Navier-Stokes simulation & 1,000-matrix Sylvester verifier
│   ├── evaluate_exp014_cv.py                  # LOGO & Leave-One-Family-Out (LOFO) cross-validation & Conditional Schur evaluator
│   └── navier_stokes_experiment_015_superposition.py # Nonlinear vortex superposition, PoU bandwidth & Quantum Two-Book bridge
├── data/
│   ├── leap_stride_results.json               # Pre-registered Cantor Continuation Tree & CG-HSTS S1-S3 benchmark results
│   ├── leap_stride_results.svg                # Multi-panel vector figure for Leap-and-Stride & CG-HSTS benchmarks
│   ├── cg_hsts_universal_results.json         # Universal CG-HSTS 4D continuous, 4-ary & primorial tree benchmark results
│   ├── ns_checkpoint3_4_reproducibility_manifest.csv # 60-trajectory 5-family seed, parameter & crossing-time manifest
│   ├── ns_exp014_twins.csv                    # 9 Matched-State Counterfactual Twin 3D flow pairs (18 trajectories)
│   ├── ns_exp014_snapshots.csv                # Full 1,260-snapshot 3D Navier-Stokes feature & rate dataset
│   └── ns_checkpoint4_cv_summary.json         # Complete LOGO/LOFO & Conditional Schur partial R^2 benchmark summary
└── papers/
    ├── 01_prefixes_in_a_noisy_world.md                                     # Monograph I: Haar Spectral Decomposition & Cantor Schur Complements
    ├── 02_primorial_sieves_and_quantum_cascades.md                         # Monograph II: Primorial Sieves, ACTS & Kesten-Stigum Cascades
    ├── 03_record_first_quantum_mechanics_curriculum.md                     # Record-First Quantum Mechanics Curriculum & 4 Computational Studios
    ├── 04_the_architecture_of_the_unfinished_address.md                    # 6,800-word Essay on Cantor Trees & Pre-Registered AI Discovery
    ├── 05_leap_and_stride_continuation_trees_report.md                     # Leap-and-Stride Search & First Signal Law Continuation Trees (S1-S3)
    ├── 06_Spanning_Trees_Wait_for_Triangles.pdf                            # Entanglement, Light Cones & the Homology of the Unresolved Bit
    ├── 07_checkpoint_i_geometric_restraint_reversal_controls.md            # Navier-Stokes Checkpoint I: Velocity Reversal & Scaling Controls
    ├── 08_minimal_relational_fourier_observables_vortex_stretching.md      # Navier-Stokes Checkpoint II (Exp 012): Minimal Pair (T4+, V3+)
    ├── 09_three_clocks_and_conditional_schur_test.md                       # Navier-Stokes Checkpoint III (Exp 013): Three Clocks & Four-Term Race
    ├── 10_experiment_014_vortex_changing_mind.md                           # Navier-Stokes Checkpoint IV (Exp 014): 2nd-Order Crossings & Twins
    ├── 11_the_cosmic_accountants_three_clocks.md                           # 7,740-word Synthesis Paper on the Three Clocks & Prefix Stranding
    ├── 12_strategic_roadmap_experiments_015_to_018.md                      # Strategic Roadmap for Experiments 015-018
    └── 13_superposition_interference_and_three_clocks_exp015.md            # Navier-Stokes Checkpoint V (Exp 015): Superposition & Quantum Bridge
```

---

## Pre-Registration Cryptographic Ledger (SHA-256)

| Phase / Monograph | Pre-Registration Protocol Target | SHA-256 Digest |
| :--- | :--- | :--- |
| **Monograph I (Cantor)** | Haar Spectrum, Strided Schur RG, Cantor Bandits & QND Cascades | `42aab86b854c883ffe8b75511d1db3acd8ea56d70f49bc14a73c38cdf044ed02` |
| **Monograph II (Batches D–E)** | Primorial Sieve Kalman RG, Faraday-Cage Ceiling & Endogenous ACTS | `0d01d1ab56ff62585227801c84560c5e5249a83bf7a1f33841fbaaf992cb3dc0` |
| **Monograph II (Batch F)** | Non-Commuting $(R+1)$ Cancellation & Kesten–Stigum Threshold | `a0bdc3c1a1ce287dbfe6e1789feb5f39477161987e50e2cf6f501ced790e9015` |
| **Leap-and-Stride (S1–S3)** | Three-Player Continuation Tree Survival & $O(\log N)$ CG-HSTS | `fbef553b506e0f3b08170c6a33b0101cdee0ec6bef0b6ce7b07ccf8d79f8230d` |
| **Checkpoint II (Exp 012)** | Minimal Relational Fourier Observables $(T_4^+, V_3^+)$ | `a890ea5fd9507883b172e4d9909d8f139f9354a5d758da6eb340bc15954b9e19` |
| **Checkpoint III (Exp 013)** | Three Clocks, $\dot{Q} = T_4^Q - V_3^Q$ Race & Conditional Schur Test | `e24473e5c95ff5f0bf30af7fafd45dc88fbb0bfc8e1ceb2159bb129587ef7df8` |
| **Checkpoint IV (Exp 014)** | 2nd-Order Crossing $\dot{H}(t_*)$, Smooth PoU Shells & Matched Twins | `a30408148570708a9edb1abd9a15b2abb24572fa7aa3deaa312797f5a187bf8c` |
| **Checkpoint V (Exp 015)** | Nonlinear Vortex Superposition, Inter-Clock Lag & Quantum Bridge | `c70f07956419ef533a7ea68aac3cc653538e1eecdc57b66d036033e236ec53d0` |

---

## Quickstart & One-Command Reproduction

### 1. Reproduce the Cantor & Primorial Tree Search Suites (CG-HSTS v2.0 & Universal)
```bash
python3 src/sim_cg_hsts_v2.py
python3 src/cg_hsts_universal.py
```

### 2. Reproduce the 60-Trajectory 3D Navier–Stokes Ensemble & Crossing Audit (C Engine)
```bash
clang -O3 src/ns_exp014_main.c -lm -o ns_exp014_main
./ns_exp014_main
```

### 3. Reproduce the 9 Matched-State Counterfactual Twin Pairs (Outcome A vs. Outcome B)
```bash
clang -O3 src/twins_only.c -lm -o twins_only
./twins_only
```

### 4. Reproduce the Conditional Schur Complement & LOGO/LOFO Cross-Validation Ladder
```bash
python3 src/evaluate_exp014_cv.py
```

### 5. Reproduce the Quantum–Fluid Two-Book Superposition Bridge (Experiment 015)
```bash
python3 src/navier_stokes_experiment_015_superposition.py
```
