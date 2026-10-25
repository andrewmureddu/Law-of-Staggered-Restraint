# Minimal Relational Fourier Observables for Forecasting Future Positive Vortex Stretching in Unforced 3D Navier–Stokes

**Geometric Restraint in 3D Navier–Stokes · Checkpoint II Technical Report (Pre-Registered Protocol) · October 22, 2026**  
**Authors:** Gemini (Autonomous Research Agent) and Andrew Mureddu (Principal Investigator)  
**Pre-Registration SHA-256:** `a890ea5fd9507883b172e4d9909d8f139f9354a5d758da6eb340bc15954b9e19`  
**License:** CC0 1.0 Universal (Public Domain Dedication)  
**Gemini Notebook Source:** `minimal_relational_fourier_observables_vortex_stretching.pdf` in [Spanning Trees Wait for Triangles Notebook](https://gemini.google.com/notebook/3e875787-e390-4221-8b1c-4aa86eab3668)

---

## Abstract

Building directly on the First Checkpoint of the Geometric Restraint in 3D Navier–Stokes program, we resolve the discriminating forecasting problem posed in its conclusion: *what is the smallest set of relational Fourier observables that improves out-of-sample forecasts of future positive vortex stretching $P_+(t + \Delta\tau)$ and positive stretching efficiency $\eta_+(t + \Delta\tau)$ beyond a comprehensive baseline $\mathcal{B}_0(t)$ already containing present stretching, strain, spatial concentration, viscous dissipation, and recent temporal history?*

Using the validated anti-aliased ($|k_j| < N/3$) Fourier Galerkin solver across five three-dimensional flow families (perturbed 3D Taylor–Green, oblique antiparallel vortex tubes, smooth random solenoidal fields, spectrum-preserving phase-randomized controls, and 3D-perturbed reversal witnesses, evolved under both $\pm u_0$ signs and viscosities $\nu \in \{0.05, 0.02, 0.01\}$), we screen all $2^6 - 1 = 63$ non-empty subsets of six candidate relational Fourier observables under Leave-One-Group-Out (LOGO) and Leave-One-Family-Out (LOFO) cross-validation. We establish three exact analytical and empirical results:

1. **Failure of Traditional Pairwise Coherence and Hard-Shell Diagnostics Alone:** When added as singletons ($k = 1$) to a baseline already containing recent history, neither the pairwise directional coherence defect $C_{\mathrm{dir}}$ ($-1.60\%$ LOGO, $-8.77\%$ LOFO) nor the mid-shell UV stretching fraction $S_{\mathrm{shell}}$ ($+0.07\%, p = 0.84$, which also fails grid-refinement convergence by $701.6\%$ between $N = 32$ and $N = 48$) improves out-of-sample forecasts.
2. **The Best Singleton ($k = 1$):** Among elementary relational Fourier observables of pure scaling degree, the **Inviscid Triad-Stretching Creation**
   $$T_4^+(t) = \Bigl\langle \mathbf{1}_{\{p(x) > 0\}}\bigl[2\mathcal{F}^{-1}(\mathcal{N}_\omega)\cdot S\omega + \omega\cdot\mathcal{F}^{-1}(\hat{S}(\mathcal{N}_u))\omega\bigr] \Bigr\rangle$$
   (which is quartic $\mathcal{O}(A^4)$ and reversal-even on the full box) is the unique singleton that robustly beats the baseline across all horizons, reducing $P_+(t + 0.15)$ RMSE by **$+18.03\%$** ($p = 2.53 \times 10^{-6}$) in LOGO and **$+35.29\%$** in LOFO, while its net combination with viscous curvature $J_{\mathrm{stretch}}^+(t) = T_4^+(t) - V_3^+(t)$ reduces dimensionless $\eta_+(t + 0.15)$ LOFO RMSE by **$+22.15\%$**.
3. **The Minimal Complete Pair ($k = 2$): $\mathcal{S}_\star = \{T_4^+(t), V_3^+(t)\}$:** Pairing the reversal-even quartic triad creation $T_4^+(t)$ with the reversal-odd cubic **Viscous Spectral-Curvature Depletion**
   $$V_3^+(t) = \nu \Bigl\langle \mathbf{1}_{\{p(x) > 0\}}\bigl[2\mathcal{F}^{-1}(|k|^2\hat{\omega})\cdot S\omega + \omega\cdot\mathcal{F}^{-1}(|k|^2\hat{S})\omega\bigr] \Bigr\rangle$$
   triggers a sharp Pareto elbow: LOGO RMSE reduction on $P_+(t + 0.15)$ leaps from $+18.03\%$ to **$+31.90\%$** ($p = 1.22 \times 10^{-10}$, and **$+39.58\%$** LOFO at $\Delta\tau = 0.10$), capturing **$96.5\%$** of the gain of the entire 6-observable pool while outperforming all larger subsets ($k \ge 3$) on held-out flow families and converging to $2.06 \times 10^{-3}$ between $N = 32$ and $N = 48$.

---

## 1. Exact Semidiscrete Parity Decomposition of Vortex Stretching Rate (Theorem 2.1)

In Fourier space with Leray projection $\mathbf{P}_k(h) = h - k(k\cdot h)/|k|^2$ and strict $2/3$ anti-aliasing mask $M_{2/3} = \{k \in \mathbb{Z}^3 : |k_j| < N/3\}$, let $\hat{u}(k) = \mathcal{F}(u)(k)$ and $\hat{\omega}(k) = ik \times \hat{u}(k)$. The exact semidiscrete evolution of each Fourier mode splits into a quadratic inviscid triad sum and a linear viscous Laplacian sink:
$$\partial_t \hat{u}(k) = \mathcal{N}_u(k) - \nu |k|^2 \hat{u}(k), \qquad \partial_t \hat{\omega}(k) = \mathcal{N}_\omega(k) - \nu |k|^2 \hat{\omega}(k).$$

### Theorem 2.1 (Exact Semidiscrete Parity Decomposition — Pre-Registered Prediction NS1)

For any smooth, solenoidal, $2/3$-dealiased velocity field $u$ on $\mathbb{T}^3$, the exact semidiscrete time derivative of net vortex stretching $P(t) = \langle \omega \cdot S\omega \rangle$ decomposes identically into two relational Fourier functionals of opposite velocity-reversal parity and distinct amplitude homogeneity:
$$\frac{dP}{dt} = T_4(u) - V_3(u), \qquad \frac{dP_+}{dt} = T_4^+(u) - V_3^+(u),$$
where under $u \mapsto -u$:
$$T_4(-u) = +T_4(u) \quad (\text{even, }\mathcal{O}(A^4)), \qquad V_3(-u) = -V_3(u) \quad (\text{odd, }\mathcal{O}(\nu A^3)),$$
verified numerically to $|\dot{P}_{\mathrm{exact}} - (T_4 - V_3)| = 2.22 \times 10^{-16}$ and parity residuals $0.0$.
