# Leap-and-Stride Search & First Signal Law Continuation Trees on Cantor Space {0, 1}^L

**Pre-Registered Simulation Report · October 20, 2026**  
**Pre-Registration SHA-256 Fingerprint:** `fbef553b506e0f3b08170c6a33b0101cdee0ec6bef0b6ce7b07ccf8d79f8230d`  
**License:** CC0 1.0 Universal (Public Domain Dedication)  
**Companion Data & Figures:** [leap_stride_results.svg](https://drive.google.com/file/d/1NVk8ZmFbpOILAxw_FtL7p_GtZVcJJInb/view) · [leap_stride_results.json](https://drive.google.com/file/d/1bpdvdX5RZam3vwkq-3vAs9lWXdupwv1U/view)  
**Theoretical Foundations:** [first-signal-law.md](https://docs.google.com/document/d/1-nDMjs3vwVExLvfz9jNoXyxpxmltjo7Ma9QUkdG-Rv8/edit?usp=drivesdk&ouid=112427431233984459040) · [First Signal Law Theory Summary](https://drive.google.com/open?id=15ghjPkXWO9_pGxOY93s_xex8lxtn7XYzXXx-EG4EmdE) · *Prefixes in a Noisy World* (October 18, 2026)

---

## 1. Executive Summary: From the Crystallization Boundary to the Continuation Tree

In [first-signal-law.md](https://docs.google.com/document/d/1-nDMjs3vwVExLvfz9jNoXyxpxmltjo7Ma9QUkdG-Rv8/edit?usp=drivesdk&ouid=112427431233984459040) and [First Signal Law Theory Summary](https://drive.google.com/open?id=15ghjPkXWO9_pGxOY93s_xex8lxtn7XYzXXx-EG4EmdE), a three-player coalition $\Omega = \{L, M, S\}$ (Large Reservoir $L$, Medium Regulator/Pivot $M$, and Small Catalyst/Scout $S$) faces an unconditionally defecting Environment $E$ between two fatal boundaries: **Dissolution** (uncoordinated noise where the coalition breaks apart) and **Crystallization** (static completeness where uncertainty hits zero, gradients flatten, and the coalition waits until the lights go out).

When Scout $S$ proposes leaping ahead into the unknown for $\Delta t = 1$ tick, the triad discovers the exact geometric constraint of *Prefixes in a Noisy World*:

1. **Measure-Zero Impossibility (§1):** $S$ cannot simply bolt to an unindexed leaf $x \in \{0, 1\}^L$, because a completed point in Cantor space has measure $2^{-L} \to 0$. Without a shared coarse address, $L$ and $M$ have nearly zero probability of ever finding $S$ again.
2. **The Cylinder Compromise:** $S$ commits to land inside a finite coarse cylinder prefix $[w]$ at horizon $\Delta t = 1$, relying on the superior combined observational precision $T_0$ of $L$ and $M$ to catch it.
3. **Why $L$ and $M$ Must Separate (Theorems 2.2 & 2.3 + Proposition 2.4):** If $L$ and $M$ stay clumped together in child cylinder $[w0]$, Theorem 2.2(1) proves they learn zero bits about the internal Haar detail modes of sibling cylinder $[w1]$, losing $S$ whenever it lands in $[w1]$. To catch $S$, $L$ and $M$ must **separate** across $[w0]$ and $[w1]$ (prefix-strided conditioning). Proposition 2.4 proves that splitting their precision budget $T_0$ into $T_0/2$ each across $[w0]$ and $[w1]$ loses **zero coarse SNR** on the parent cylinder $[w]$ via their Haar sum $(Y_{w0} + Y_{w1})/2$, while their Haar difference $(Y_{w0} - Y_{w1})/2$ resolves which child branch holds $S$.

Turning this three-player "Leap-and-Stride" catch dynamic into a general search algorithm yields **Cotton-Gated Haar–Schur Tree Search (CG-HSTS)**, which we tested across three pre-registered simulation batches (S1–S3).

---

## 2. Pre-Registration Scorecard

| ID | Domain | Pre-Registered Claim & Quantitative Test | Empirical Result | Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **S1a** | Haar Sum Conservation | Separating $L$ and $M$ across $[0]$ and $[1]$ with precision $T_0/2$ each preserves parent coarse posterior variance to $< 10^{-13}$ | Parent variance $0.212182$ together vs. $0.212182$ separated (residual **$1.83 \times 10^{-15}$**) | **HOLDS** |
| **S1b** | Clumped vs. Separated Catch | Clumped $[0]$ leaves $[1]$ internal Haar modes unchanged ($< 10^{-13}$); separation maximizes $I(X; Y)$, $D(S)$, and sub-cylinder catch rate | Internal shift **$2.22 \times 10^{-15}$**; resolved variance $D$ rises **$69.86 \to 92.09$** ($+31.8\%$); catch rate rises **$51.17\% \to 70.09\%$** | **HOLDS** |
| **S2a** | Boundary Mortality | Crystallization (hold at $k=0$) and Measure-Zero Bolt (unindexed jump to $L=10$) both suffer $0\%$ survival | Crystallization dies at **$t = 77.0$** ($0\%$ survival); Measure-Zero Bolt dies at **$t = 5.04$** ($0\%$ survival) | **HOLDS** |
| **S2b** | First-Jump Neck & Bootstrap | Step-1 catch is the tightest risk bottleneck ("works, but barely"), after which Cotton shell-gating bootstraps $> 90\%$ survival to $T = 200$ | Step-1 catch = **$92.0\%$** ($8\%$ initial wobble absorbed by $L$'s slack) $\to$ subsequent catch = **$98.85\%$**; **$100\%$ survival** ($E_{200} = 268.65 \pm 0.67$) vs. $0\%$ for Clumped Pursuit (dies at $t = 17.22$) | **HOLDS** |
| **S3a** | Fast Search Runtime | $O(\log N)$ Haar–Schur tree updater runs $\ge 50\times$ faster than exact $N \times N$ Gaussian Process conditioning on $N = 1,024$ leaves | **$3.40\text{ ms}$** (CG-HSTS) vs. **$17,771.43\text{ ms}$** (Exact GP) across 256 steps — a **$5,232\times$ speedup** | **HOLDS** |
| **S3b** | Search Regret & Reward | CG-HSTS beats Standard UCT/MCTS, Uncapped Flat TS, and Clumped Greedy in simple regret ($p < 10^{-4}$) and cumulative reward | Simple regret **$0.257 \pm 0.038$** vs. **$0.967 \pm 0.060$** (UCT, $p = 2.85 \times 10^{-25}$), **$0.808$** (Flat TS), **$1.143$** (Clumped); highest cumulative reward **$+417.72 \pm 46.51$** | **HOLDS** |

---

## 3. Experiment S1: The Haar Sum-and-Difference Separation Theorem ($L = 6, N = 64$)

In Experiment S1, catchers $L$ and $M$ share total precision budget $T_0 = 4.0$ on an ultrametric Cantor tree ($L = 6, N = 64$ leaves, scale decay $\beta = 0.8$).

- **Zero Cost to Coarse Cohesion (S1a):** When $L$ and $M$ remain together on the parent cylinder, the posterior variance of the coarse root mean drops from $1.402642$ to **$0.2121817692775456$**, but the posterior variance of the child difference mode $\Delta X = (X_{[0]} - X_{[1]})/2$ remains completely unreduced at **$0.402642$** (zero directional information). When $L$ and $M$ separate into $[0]$ and $[1]$ with diluted precision $T_0/2 = 2.0$ each, their Haar sum recovers the exact parent posterior variance **$0.2121817692775474$** (algebraic residual **$1.83 \times 10^{-15}$**), while their Haar difference simultaneously cuts the child difference variance by $61.7\%$ down to **$0.154235$**.
- **Eliminating the Dark Sibling Blind Spot (S1b):** When $L$ and $M$ clump inside half-tree $[0]$ (monitoring $[00]$ and $[01]$), all 31 internal Haar detail eigenvalues of sibling half-tree $[1]$ remain identically at their prior values (maximum numerical shift **$2.22 \times 10^{-15}$**). Separating across $[0]$ and $[1]$ increases mutual information from $1.9916$ bits to $2.0546$ bits, expands total resolved field variance $D(S)$ from **$69.86$ to $92.09$** ($+31.8\%$), and increases 4-subcylinder Scout catch accuracy across $10,000$ Monte Carlo trials from **$51.17\%$ to $70.09\%$**.

---

## 4. Experiment S2: The Three-Player Continuation Tree Game ($L = 10, N = 1,024$, 50 Paired Seeds)

Experiment S2 simulates the three-player coalition $\Omega = \{L, M, S\}$ over horizon $T = 200$ against Environment $E$'s constant metabolic drain (burn rate $b_E = 0.42/\text{step}$, initial reserve $E_0 = 20.0$):

| Coalition Strategy | Survival Rate ($T = 200$) | Mean Survival Time | Final Slack Reserve $E(200)$ | Final Continuation Depth $k(200)$ | Failure / Survival Mechanism |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Crystallization (Hold at $k = 0$)** | $0.0\%$ (0/50) | $77.00$ steps | $0.00 \pm 0.00$ | $0.00 / 10$ | Static gradient depletes exponentially ("waiting until the lights go out") |
| **Measure-Zero Bolt (Jump to $L = 10$)** | $0.0\%$ (0/50) | $5.04$ steps | $0.00 \pm 0.00$ | $0.40 / 10$ | Unindexed leaf jump has $2/1024$ catch probability $\to$ immediate dissolution |
| **Clumped Pursuit (No Separation)** | $0.0\%$ (0/50) | $17.22$ steps | $0.00 \pm 0.00$ | $2.52 / 10$ | Catchers stay together in $[w0]$; $50\%$ of leaps land in dark sibling $[w1]$ (Thm 2.2) |
| **Leap-and-Stride Triad (FSL + Haar)** | **$100.0\%$ (50/50)** | **$200.00$ steps** | **$268.65 \pm 0.67$** | **$10.00 / 10$** | Separated catchers lock coarse sum & fine difference, bootstrapping all 10 shells |

### Why It "Works, But Barely, at Least the First Time"

At step $t = 1$ (the very first jump from $k = 0$ to $k = 1$), fine-scale aliasing variance is unscreened, yielding a first-jump catch rate of $92.0\%$ (4 out of 50 seeds experience an immediate first-step miss). In those $8\%$ of trajectories, survival depends entirely on Player $L$'s slack reserve ($e_L - c_L - r_L > 0$ in [first-signal-law.md](https://docs.google.com/document/d/1-nDMjs3vwVExLvfz9jNoXyxpxmltjo7Ma9QUkdG-Rv8/edit?usp=drivesdk&ouid=112427431233984459040)) absorbing the initial miss penalty so the triad can try again at $t = 2$. Once the first catch locks in, Cotton-damped multi-shell conditioning drives the catch reliability on all subsequent jumps up to **$98.85\%$**, propelling the continuation tree to full depth $k = 10$ on every seed.

---

## 5. Experiment S3: Fast Search Benchmark (CG-HSTS vs. UCT/MCTS, Flat TS, and Exact GP)

Experiment S3 benchmarks **CG-HSTS** on a 1,024-leaf ($L = 10$) noisy multiscale Cantor optimization task (observation noise $\sigma = 1.0$, $\beta = 0.70$, $T = 256$ queries, 50 paired seeds) in [leap_stride_results.json](https://drive.google.com/file/d/1bpdvdX5RZam3vwkq-3vAs9lWXdupwv1U/view):

| Search Algorithm | Runtime for 256 Steps | Speedup vs. Exact GP | Final Simple Regret (Lower = Better) | Cumulative Reward ($T = 256$) | Paired Regret Diff vs. UCT |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **CG-HSTS (Cotton–Haar UCB + Stride)** | **$3.40\text{ ms}$** | **$5,232\times$** | **$0.257 \pm 0.038$** | **$+417.72 \pm 46.51$** | **$-0.710 \pm 0.068$ ($p = 2.85 \times 10^{-25}$)** |
| **CG-HSTS (Pure Global Stride Variant)** | **$3.83\text{ ms}$** | **$4,685\times$** | **$0.255 \pm 0.036$** | $+365.37 \pm 45.97$ | **$-0.713 \pm 0.065$ ($p = 9.00 \times 10^{-28}$)** |
| **Standard UCT / MCTS (Undamped)** | $3.12\text{ ms}$ | $5,700\times$ | $0.967 \pm 0.060$ | $+414.72 \pm 44.20$ | Baseline (Crystallizes on local peak) |
| **Clumped Greedy Tree Search** | $2.89\text{ ms}$ | $6,150\times$ | $1.143 \pm 0.094$ | $+267.98 \pm 35.25$ | $+0.176 \pm 0.088$ (Worse) |
| **Uncapped Flat Thompson Sampling** | $48.50\text{ ms}$ | $366\times$ | $0.808 \pm 0.065$ | $-36.12 \pm 47.59$ | $-0.159 \pm 0.079$ (Strands in noise) |
| **Exact $N \times N$ Matrix GP Search** | $17,771.43\text{ ms}$ | $1.0\times$ | $0.255$ (Exact ref) | $+417.72$ (Exact ref) | Exact algebraic reference |

---

## 6. Complete Algorithm Specification: Self-Tuning CG-HSTS v2.0 (Resolution of Open Problem 1)

**Reference Implementation:** `src/sim_cg_hsts_v2.py` ([Drive](https://drive.google.com/file/d/1VmdZkCpTa1kzPYcNLUmK0T8YXYv6KfTo/view))

We replace all oracle hyperparameters with an exact $O(L)$ online empirical Bayes estimator driven by the three-role triad $(L, M, S)$ using the harmonic pair weight $h_w = \frac{2 n_{w0} n_{w1}}{n_{w0} + n_{w1}}$ and empirical Haar detail coefficient $d_w = \frac{1}{2}(\bar{Y}_{w0} - \bar{Y}_{w1})$.

| Landscape Regime ($N = 1,024$) | Prior-Free CG-HSTS v2.0 Regret | Standard UCT / MCTS Regret | Uncapped Flat TS Regret | Paired Regret Reduction vs. UCT | Cumulative Reward (CG-HSTS v2 vs. UCT vs. Flat TS) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Rough Field ($\beta = 0.4$)** | **$0.592 \pm 0.096$** | $1.255 \pm 0.144$ | $0.775 \pm 0.089$ | $-0.663 \pm 0.168$ ($p = 8.35 \times 10^{-5}$) | $+784.2$ vs. $+765.6$ (UCT) vs. $+89.3$ (Flat TS) |
| **Moderate Field ($\beta = 0.7$)** | **$0.303 \pm 0.050$** | $1.002 \pm 0.075$ | $0.734 \pm 0.082$ | $-0.698 \pm 0.086$ ($p = 8.48 \times 10^{-16}$) | $+373.0$ vs. $+368.4$ (UCT) vs. $-54.0$ (Flat TS) |
| **Smooth Field ($\beta = 1.1$)** | **$0.184 \pm 0.026$** | $0.626 \pm 0.048$ | $0.547 \pm 0.050$ | $-0.442 \pm 0.053$ ($p = 9.50 \times 10^{-17}$) | $+224.0$ vs. $+218.4$ (UCT) vs. $-20.4$ (Flat TS) |

---

## 7. Universal CG-HSTS: General $b$-ary / Primorial Trees, Continuous 4D Optimization, and Resolution of Open Problem 2

**Universal Engine Code:** `src/cg_hsts_universal.py` ([Drive](https://drive.google.com/file/d/1oFaQu5NKko8mmXM3LIOYijYESG6jn8Da/view)) · **Benchmark Data:** `data/cg_hsts_universal_results.json` ([Drive](https://drive.google.com/file/d/1-Qsa_AmheuCcEjCAkC8A219sdBKCxpGC/view))

We generalize CG-HSTS v2.0 to arbitrary branching sequences $(b_1, b_2, \dots, b_L)$ and continuous $d$-dimensional domains $[0, 1]^d$ via the Chrestenson–Helmert multi-branch detail energy $Q_w = \sum_{j=0}^{b_k-1} n_{w,j}(\bar{Y}_{w,j} - \bar{Y}_w)^2$ and generalized harmonic weight $H_w = N_w - \frac{1}{N_w}\sum_{j=0}^{b_k-1} n_{w,j}^2$:

| Benchmark Domain & Tree Geometry | Universal CG-HSTS Regret | Baseline UCT / StoSOO Regret | Quasi-Random / LHS Regret | Paired Regret Reduction vs. UCT / StoSOO | Cumulative Reward (Universal CG-HSTS vs. UCT) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **G1: 4D Continuous Box $[0, 1]^4$** ($L = 12, N = 4,096$ cells, $\sigma = 0.85$) | **$0.3960 \pm 0.0659$** | $0.8909 \pm 0.1080$ (StoSOO) | $2.2411 \pm 0.1544$ (LHS) | **$-0.4949 \pm 0.1307$ ($p = 1.64 \times 10^{-4}$, $2.25\times$ lower)** | Continuous global optimization |
| **G2a: Quaternary 4-ary Game Tree** ($b = [4, 4, 4, 4, 4], N = 1,024$) | **$0.273 \pm 0.048$** | $0.907 \pm 0.096$ (UCT) | — | **$-0.634 \pm 0.110$ ($p = 1.06 \times 10^{-8}$, $3.32\times$ lower)** | **$+534.1 \pm 44.5$** vs. $+519.2 \pm 45.2$ |
| **G2b: Primorial Bruhat-Tits Tree** ($b = [2, 3, 5, 7, 5], N = 1,050$) | **$0.332 \pm 0.047$** | $1.055 \pm 0.081$ (UCT) | — | **$-0.724 \pm 0.079$ ($p = 7.26 \times 10^{-20}$, $3.18\times$ lower)** | **$+411.7 \pm 52.0$** vs. $+406.2 \pm 50.2$ |
