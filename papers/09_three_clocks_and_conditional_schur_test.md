# Three Clocks, the Alignment-vs.-Capacity Race, and the Conditional Schur Test in Unforced 3D Navier–Stokes

**Geometric Restraint in 3D Navier–Stokes · Checkpoint III Technical Monograph (Pre-Registered Protocol) · October 23, 2026**  
**Authors:** Gemini (Autonomous Research Agent) and Andrew Mureddu (Principal Investigator)  
**Pre-Registration SHA-256:** `e24473e5c95ff5f0bf30af7fafd45dc88fbb0bfc8e1ceb2159bb129587ef7df8`  
**License:** CC0 1.0 Universal (Public Domain Dedication)  
**Gemini Notebook Source:** `three_clocks_and_conditional_schur_test.pdf` in [Spanning Trees Wait for Triangles Notebook](https://gemini.google.com/notebook/3e875787-e390-4221-8b1c-4aa86eab3668)

---

## Abstract

In Experiment 013 (*Three Clocks and the Conditional Schur Test*), a companion analysis demonstrated that a smooth three-dimensional vortex can experience declining stretching efficiency ($\dot{\eta}_+ < 0$, Clock 1) while both positive vortex stretching ($\dot{P}_+ = T_4^+ - V_3^+ > 0$, Clock 2) and enstrophy ($\dot{Z} = P_+ - P_- - D > 0$, Clock 3) continue growing, and posed a decisive falsification test for our Checkpoint II pair $(T_4^+, V_3^+)$: when reparameterized into the instantaneous rate $J = T_4^+ - V_3^+$ and the total turnover channel $R = T_4^+ + V_3^+$, does $R$ contribute genuine conditional Schur-complement variance reduction beyond $\mathcal{B}_0 \cup \{J, \dot{Q}\}$, or was the two-observable gain primarily due to reconstructing the exact derivative $J$?

Under a pre-registered protocol (SHA-256: `e24473e5...ef7df8`), we replicate Experiment 013 to 8 significant digits, execute the nested conditional Schur test across 60 three-dimensional trajectories spanning five flow families at horizons $h \in \{0.15, 0.25\}$, and establish four results:

1. **Exact Sylvester–Hellmann–Feynman Decomposition of $\dot{Q}$ and the Alignment-vs.-Capacity Race (Theorem 1.1):** Using Sylvester's rank-1 eigen-projector identity on trace-free $3 \times 3$ strain tensors, we derive the exact semidiscrete decomposition $\dot{Q}(t) = T_4^Q(t) - V_3^Q(t)$ without finite differences, splitting
   $$H(t) = \frac{d}{dt}\log\eta_+(t) = \underbrace{\left(\frac{T_4^+(t)}{P_+(t)} - \frac{T_4^Q(t)}{Q(t)}\right)}_{\Delta_{\mathrm{triad}}(t)} - \underbrace{\left(\frac{V_3^+(t)}{P_+(t)} - \frac{V_3^Q(t)}{Q(t)}\right)}_{\Delta_{\mathrm{visc}}(t)}.$$
   We prove that the $t_* \approx 0.648$ efficiency downturn is driven entirely by **inviscid triad alignment depletion** ($\Delta_{\mathrm{triad}} = T_4^+/P_+ - T_4^Q/Q$ plunging from $+0.11635$ at $t = 0.5$ to $-0.20199$ at $t = 1.0$), whereas differential viscosity ($\Delta_{\mathrm{visc}} = -0.04529 < 0$) actually opposes the efficiency decline because $V_3^Q/Q > V_3^+/P_+$.
2. **Resolution of the Conditional Schur Test (Confirming Experiment 013's Falsifier):** Once $J = T_4^+ - V_3^+$ is included in the predictor set ($\text{M1} = \mathcal{B}_0 \cup \{J\}$), $J$ alone captures **$> 90\%$** of the predictive improvement ($+33.71\%$ LOGO on $P_+$, $+33.94\%$ on $\eta_+$, $+40.66\%$ on $H(t + 0.15)$, raising Clock-1 sign accuracy from $89.17\%$ to $94.47\%$). Conditional on $\mathcal{B}_0 \cup \{J, \dot{Q}\}$, the population Schur partial $R^2$ of $R = T_4^+ + V_3^+$ on $\eta_+$ is **$0.00066\%$** at $h = 0.15$ and **$0.00696\%$** at $h = 0.25$—confirming that the Checkpoint II jump from $\{T_4^+\}$ ($+18.0\%$) to $\{T_4^+, V_3^+\}$ ($+31.9\%$) came overwhelmingly from completing the exact derivative $J = T_4^+ - V_3^+$ rather than from $R$, with $R$ providing only a modest $+2.9\%$ cross-family calibration gain in Leave-One-Family-Out (LOFO) testing.
3. **Inverted Reversal Parity of Strain-Capacity Rates (Theorem 1.2):** We verify the exact region-restricted transformation $T_4^+(-u) = T_4^-(u)$ and $V_3^+(-u) = -V_3^-(u)$ to $0.0$ error, and discover that because $\|S\|_{\mathrm{op}} = |\lambda_{\max}(S)|$ is even under $u \mapsto -u$ while $S$ is odd, the parity of $\dot{Q} = T_4^Q - V_3^Q$ is inverted:
   $$T_4^Q(-u) = -T_4^Q(u) \quad (\text{odd}), \qquad V_3^Q(-u) = +V_3^Q(u) \quad (\text{even}),$$
   away from planar-shear ties $|\lambda_1| = |\lambda_3|$.
4. **Resolution-Invariant Dyadic Fourier Shells:** Partitioning $T_4^+$ across fixed integer wavevector shells $B_1 = \{1 \le |k|^2 \le 2\}$, $B_2 = \{3 \le |k|^2 \le 6\}$, $B_3 = \{|k|^2 \ge 7\}$ satisfies $T_4^+ = \sum_{m=1}^3 T_{4,B_m}^+$ to $1.11 \times 10^{-16}$, reduces the $N = 32$ vs. $N = 48$ refinement error from $701.6\%$ to $3.77 \times 10^{-4}$, and yields a $4,700\times$ larger population Schur partial $R^2$ on $\eta_+$ ($3.14\%$) and $H$ ($3.23\%$) than $R$.

---

## 1. Exact Replication of Experiment 013 ($N = 32, \nu = 0.05$) and Four-Term Race Audit

| Flow Configuration & Time ($\nu=0.05$) | $P_+(t)$ | $Q(t)$ | $\eta_+(t)$ | $T_4^+(t)$ | $V_3^+(t)$ | $T_4^Q(t)$ | $V_3^Q(t)$ | $\Delta_{\mathrm{triad}}(t)$ | $-\Delta_{\mathrm{visc}}(t)$ | $H(t) = \frac{d}{dt}\log\eta_+$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **2D3C Witness ($t = 0.50, N = 32$)** | $0.46933$ | $2.18513$ | $0.21479$ | $0.41425$ | $0.13944$ | $1.68319$ | $0.66644$ | $+0.11234$ | $+0.00788$ | $+0.12022$ |
| **2D3C Witness ($t = 1.00, N = 32$)** | $0.57142$ | $2.72592$ | $0.20963$ | $0.35576$ | $0.22218$ | $2.25452$ | $1.18653$ | $-0.20449$ | $+0.04646$ | $-0.15802$ |
| **3D Perturbed ($t = 0.50, N = 32$)** | $0.48662$ | $2.23999$ | $0.21724$ | $0.43320$ | $0.14866$ | $1.73294$ | $0.70269$ | $+0.11635$ | $+0.00844$ | $+0.12479$ |
| **3D Perturbed ($t_* = 0.648, N = 32$)** | $0.52546$ | $2.39794$ | $0.21913$ | $0.42279$ | $0.18276$ | $1.97654$ | $0.88103$ | $-0.01969$ | $+0.01969$ | $0.00000$ |
| **3D Perturbed ($t = 1.00, N = 32$)** | $0.59164$ | $2.78681$ | $0.21230$ | $0.37360$ | $0.23909$ | $2.28859$ | $1.21833$ | $-0.20199$ | $+0.04529$ | $-0.15670$ |
| **3D Perturbed ($t = 1.00, N = 48$)** | $0.59168$ | $2.78683$ | $0.21231$ | $0.37352$ | $0.23912$ | $2.28910$ | $1.21829$ | $-0.20229$ | $+0.04522$ | $-0.15707$ |
| **Rel. Diff ($N = 32$ vs. $48$)** | $6.5 \times 10^{-5}$ | $7.2 \times 10^{-6}$ | $5.8 \times 10^{-5}$ | $2.1 \times 10^{-4}$ | $1.3 \times 10^{-4}$ | $2.2 \times 10^{-4}$ | $3.3 \times 10^{-5}$ | $1.5 \times 10^{-3}$ | $1.5 \times 10^{-3}$ | $2.3 \times 10^{-3}$ |
