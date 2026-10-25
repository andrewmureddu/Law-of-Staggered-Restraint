# Briefing Document: Geometric Restraint and Reversal Controls in Unforced Navier–Stokes Flow (Checkpoint I)

**Authors:** Andrew Mureddu & Gemini  
**Date:** October 2026  
**License:** CC0 1.0 Universal (Public Domain Dedication)  
**Google Doc Edition:** [Briefing Document: Geometric Restraint and Reversal Controls in Unforced Navier–Stokes Flow](https://docs.google.com/document/d/1dzrEpqAfECvIqFIC6urxvflBuLQXM_e6h29e19TA57k/edit?usp=drivesdk&ouid=112427431233984459040)

---

## Executive Summary

This briefing document summarizes the initial research checkpoint evaluating geometric restraint, vortex stretching, and directional coherence in smooth, incompressible, unforced Navier–Stokes flow on a periodic domain $[0, 2\pi)^3$. The primary objective of this research suite is to establish exact control experiments, implement a validated Fourier solver, and disentangle the diagnostic measurements required to test proposed feedback mechanisms from spatial concentration to suppressed vortex stretching.

### Critical Takeaways

- **Obstruction to Simple Geometric Explanations:** Spatial vorticity concentration and pairwise directional coherence alone cannot determine the sign or magnitude of instantaneous net vortex stretching ($P$). A velocity-reversed initial field ($u \mapsto -u$) shares identical spatial vorticity magnitudes, concentration statistics, velocity/vorticity spectra, viscous dissipation, and pairwise directional coherence, yet exhibits an exact sign reversal in vortex stretching production ($P \mapsto -P$).
- **Amplitude Scaling Constraints:** Proposed analytical estimates that bound stretching using only amplitude-invariant geometric descriptors cannot hold across arbitrary scales. Scaling initial data by a factor $A$ ($u_A = A u_0$) scales stretching as $A^3$ while enstrophy and dissipation scale as $A^2$. Any viable geometric estimate must incorporate amplitude-sensitive information, restricted data classes, or dynamically controlled quantities along the evolution.
- **Competing Geometry in Viscous Damping:** Decomposing vorticity into magnitude $q = |\omega|$ and direction $\xi = \omega/|\omega|$ demonstrates that directional gradients directly participate in viscous damping ($|\nabla \omega|^2 = |\nabla q|^2 + q^2 |\nabla \xi|^2$). Reducing directional gradients (increasing directional coherence) reduces directional viscous dissipation, proving that coherence does not act as an automatically positive damping term.
- **Diagnostic Refinement and Numerical Convergence:** Numerical tests confirm that velocity-field convergence does not guarantee the convergence of sharply thresholded percentile diagnostics. A continuous, intensity-weighted, threshold-free coherence defect ($C_{\text{dir}}$) displays superior stability under grid refinement (relative error $\sim 10^{-8}$ to $10^{-6}$) compared to percentile-masked diagnostics ($\sim 10^{-2}$).

---

## Analytical Foundations and Reversal Controls

### Mathematical Bookkeeping

The unforced, incompressible Navier–Stokes system is analyzed on the periodic domain $[0, 2\pi)^3$. Defining $\langle \cdot \rangle$ as the normalized spatial average over $(2\pi)^3$:

- **Vorticity:** $\omega = \nabla \times u$
- **Strain Tensor:** $S = \frac{1}{2}(\nabla u + \nabla u^{\mathrm{T}})$
- **Kinetic Energy:** $E = \frac{1}{2}\langle |u|^2 \rangle$
- **Enstrophy:** $W = \langle |\omega|^2 \rangle$, with normalized enstrophy $Z = \frac{1}{2}W$
- **Vortex Stretching Production Density:** $p(x) = \omega \cdot S \omega$
- **Net Vortex Stretching Production:** $P = \langle \omega \cdot S \omega \rangle = P_+ - P_-$, where $P_+ = \langle \max(p,0) \rangle$ and $P_- = \langle \max(-p,0) \rangle$
- **Viscous Dissipation:** $D = \nu \langle |\nabla \omega|^2 \rangle$

For smooth periodic solutions, the exact global conservation and evolution equations are:
$$\frac{dE}{dt} = -\nu W, \qquad \frac{dZ}{dt} = P - D.$$

### The Velocity Reversal Control ($u \mapsto -u$)

For any smooth, periodic, divergence-free velocity field $u$, mapping $u \mapsto -u$ transforms the fields as $\omega \mapsto -\omega$, $S \mapsto -S$, $p(x) \mapsto -p(x)$. Under this transformation:
- $|\omega(x)|$, vorticity spectra, energy, enstrophy, spatial concentration statistics ($C$), viscous dissipation ($D$), and pairwise directional angles remain **exact invariants**.
- Integrated stretching production flips sign ($P \mapsto -P$), and positive and negative production integrals interchange ($P_+ \leftrightarrow P_-$).

### Explicit Shear Field Witness

An explicit divergence-free velocity field $u_0$ is constructed:
$$u_0(x,y,z) = \big(\sin y, \, 0, \, \sin x - \sin(x+y)\big)$$
Direct calculation yields:
$$E_0 = \frac{3}{4}, \quad W_0 = 2, \quad Z_0 = 1, \quad P_0 = \frac{1}{4}, \quad \langle |\nabla \omega_0|^2 \rangle = 3.$$
At viscosity $\nu = 0.05$:
$$Z'(0; u_0) = 0.25 - 0.05(3) = +0.10 \quad (\text{Enstrophy Growth}), \qquad Z'(0; -u_0) = -0.25 - 0.05(3) = -0.40 \quad (\text{Enstrophy Decay}).$$

---

## Quantitative Results and Empirical Validation

| Flow Family | $N$ | Viscosity ($\nu$) | Step ($\Delta t$) | Sign | Initial $Z$ | Final $Z$ ($t=1$) | Final $P/D$ ($t=1$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **2D Taylor-Green** | 16 | 0.05 | 0.01 | $+1$ | $0.500000$ | $0.409365$ | $0.000000$ |
| **ABC Beltrami** | 16 | 0.05 | 0.01 | $+1$ | $1.500000$ | $1.357256$ | $0.000000$ |
| **Witness ($u_0$)** | 24 | 0.05 | 0.01 | $+1$ | $1.000000$ | $1.209431$ | $1.727887$ |
| **Witness ($-u_0$)** | 24 | 0.05 | 0.01 | $-1$ | $1.000000$ | $0.843076$ | $1.141326$ |
| **Witness Check** | 32 | 0.05 | 0.005 | $+1$ | $1.000000$ | $1.209431$ | $1.727887$ |
| **3D Taylor-Green** | 24 | 0.05 | 0.01 | $+1$ | $0.375000$ | $0.295480$ | $0.410728$ |
| **3D Taylor-Green** | 32 | 0.05 | 0.01 | $+1$ | $0.375000$ | $0.295480$ | $0.410731$ |
| **3D Taylor-Green** | 48 | 0.05 | 0.01 | $+1$ | $0.375000$ | $0.295480$ | $0.410731$ |
| **3D Taylor-Green** | 32 | 0.01 | 0.005 | $+1$ | $0.375000$ | $0.388428$ | $2.560624$ |
| **3D Taylor-Green** | 48 | 0.01 | 0.005 | $+1$ | $0.375000$ | $0.388428$ | $2.560626$ |
