# Pre-Registration Protocol: Navier–Stokes Checkpoint IV (Experiment 014 — The Vortex That Knew It Was About to Change Its Mind)

**Date:** October 23, 2026  
**Authors:** Gemini (Autonomous Research Agent) & Andrew Mureddu (Principal Investigator)  
**Target:** Second-Order Crossing Dynamics $\dot{H}(t_*)$, Smooth Scale-Normalized Partition-of-Unity Fourier Shells (Prefix-Stranding Cure), Pre-Crossing Event Forecasting $\mathcal{E}_h(t)$, and Matched-State Counterfactual Twins (Outcome A vs. Outcome B).

## 1. Pre-Registered Claims, Equations, and Falsifiers

### Claim C4.1 (Exact Second-Order Sylvester–Cayley–Hamilton Operator for $\ddot{Q}$, $\ddot{P}_+$, and the Crossing Acceleration $\dot{H}(t_*)$)
1. For any $3\times 3$ symmetric trace-free strain tensor $S(x,t)$ with simple dominant eigenvalue $\lambda_{\max}(x,t)$ ($|\lambda_{\max}| = \|S\|_{\mathrm{op}}$, $3\lambda_{\max}^2 \neq J_2$), differentiating $\lambda_{\max}^3 - J_2 \lambda_{\max} - J_3 = 0$ twice yields the exact eigenvector-free second derivative of the operator norm:
   $$\partial_t^2 \|S\|_{\mathrm{op}} = \mathrm{sign}(\lambda_{\max})\,\ddot{\lambda}_{\max}, \qquad \ddot{\lambda}_{\max} = \frac{\lambda_{\max}\ddot{J}_2 + \ddot{J}_3 + 2\dot{J}_2\dot{\lambda}_{\max} - 6\lambda_{\max}(\dot{\lambda}_{\max})^2}{3\lambda_{\max}^2 - J_2},$$
   where $\dot{J}_2 = \mathrm{Tr}(S\dot{S})$, $\ddot{J}_2 = \mathrm{Tr}(\dot{S}^2 + S\ddot{S})$, $\dot{J}_3 = \mathrm{Tr}(S^2\dot{S})$, and $\ddot{J}_3 = \mathrm{Tr}(2S\dot{S}^2 + S^2\ddot{S})$.
2. Substituting the exact single-snapshot Navier–Stokes field accelerations decomposes $\ddot{P}_+(t)$ and $\ddot{Q}(t)$ into a Quintic-Quartic-Cubic hierarchy plus two positive distributional corner-flux integrals ($\mathcal{B}_{\partial\Omega_+}$ on the zero-stretching surface $p=0$ and $\mathcal{B}_{\Sigma_0}^Q$ on the planar-shear surface $\lambda_2(S)=0$):
   $$\ddot{P}_+(t) = T_5^+(t) - M_4^+(t) + V_{3,2}^+(t) + \mathcal{B}_{\partial\Omega_+}(t), \qquad \ddot{Q}(t) = T_5^Q(t) - M_4^Q(t) + V_{3,2}^Q(t) + \mathcal{B}_{\Sigma_0}^Q(t).$$
3. At the Clock-1 efficiency turning point $H(t_*) = 0$, the squared first-derivative terms cancel identically:
   $$\dot{H}(t_*) = \frac{\ddot{P}_+(t_*)}{P_+(t_*)} - \frac{\ddot{Q}(t_*)}{Q(t_*)} = \Delta_{T_5}(t_*) - \Delta_{M_4}(t_*) + \Delta_{V_{3,2}}(t_*) + \frac{\mathcal{B}_{\partial\Omega_+}(t_*)}{P_+(t_*)} - \frac{\mathcal{B}_{\Sigma_0}^Q(t_*)}{Q(t_*)}.$$

### Claim C4.2 (Curing "Prefix Stranding" via Smooth Scale-Normalized Partition-of-Unity Fourier Shells)
Replacing hard integer shells with smooth, overlapping, scale-normalized Partition-of-Unity (PoU) filters $w_m(\xi)$ in $\xi(k; u) = \log_2(|k| / k_{\mathrm{char}}(u))$ where $k_{\mathrm{char}}(u) = \sqrt{D(u)/(\nu W(u))}$ preserves $\sum_{m=1}^3 \widetilde{T}_{4,m}^+(u) \equiv T_4^+(u)$ to machine precision ($< 10^{-15}$) while curing the Leave-One-Family-Out (LOFO) generalization collapse ($+15.05\% \to > +33\%$).

### Claim C4.3 (Experiment 014 Pre-Crossing Event Detection $\mathcal{E}_h(t)$ and Remaining Time-to-Turn $\tau_{\mathrm{turn}}(t)$)
Restricting strictly to pre-crossing states where $H(t) > 0$, adding the exact first-derivative rates $(J, \dot{Q})$, smooth scale-normalized PoU shells, and the exact single-snapshot crossing acceleration $\dot{H}(t)$ improves out-of-sample $\mathcal{E}_h(t)$ crossing detection and $\tau_{\mathrm{turn}}(t) = t_* - t$ RMSE under both LOGO and LOFO.

### Claim C4.4 (Outcome A vs. Outcome B — Matched-State Counterfactual Twin Flows)
Pairs of 3D velocity fields $(u^{(A)}(t_0), u^{(B)}(t_0))$ matched to $< 10^{-11}$ on $(P_+, Q, \eta_+, J, \dot{Q}, H)$ at $t_0$ diverge in their subsequent Clock-1 turning time $t_*^{(A)} \neq t_*^{(B)}$ due to differences in their high-wavenumber triad shell distribution and initial second derivative $\dot{H}(t_0)$.
