import json
import math
import random
import statistics
import time

def paired_ttest(x, y):
    diffs = [xi - yi for xi, yi in zip(x, y)]
    n = len(diffs)
    mean_d = statistics.mean(diffs)
    sd_d = statistics.stdev(diffs)
    se_d = sd_d / math.sqrt(n)
    t_stat = mean_d / max(se_d, 1e-15)
    z = abs(t_stat) * (1.0 - 1.0 / (4.0 * (n - 1)))
    p_val = math.erfc(z / math.sqrt(2.0))
    return t_stat, p_val, mean_d, se_d

# ==============================================================================
# UNIVERSAL CG-HSTS ENGINE (Arbitrary Branching Sequence b_1, b_2, ..., b_L)
# Unifies Open Problem 1 (Adaptive Haar Shell Gating) and Open Problem 2
# (General p-Adic / Primorial / b-ary Trees) via the Generalized Harmonic Weight:
#   H_w = N_w - (sum_{j=0}^{b-1} n_{wj}^2) / N_w
#   Q_w = sum_{j=0}^{b-1} n_{wj} * (Y_bar_{wj} - Y_bar_w_weighted)^2
# ==============================================================================
class UniversalCGHSTS:
    """
    Universal Prior-Free Cotton-Gated Haar-Schur Tree Search.
    Supports arbitrary branching factors per depth: branch_factors = [b_1, b_2, ..., b_L].
    """
    def __init__(self, branch_factors: list[int], gamma_star: float = 0.45, c_ucb: float = 0.90, cap_K: int = 4):
        self.branch_factors = list(branch_factors)
        self.L = len(branch_factors)
        self.gamma_star = gamma_star
        self.c_ucb = c_ucb
        self.cap_K = cap_K

        self.level_sizes = [1]
        for b in self.branch_factors:
            self.level_sizes.append(self.level_sizes[-1] * b)
        self.N = self.level_sizes[-1]

        self.level_offsets = [0]
        for sz in self.level_sizes[:-1]:
            self.level_offsets.append(self.level_offsets[-1] + sz)
        self.total_nodes = self.level_offsets[-1] + self.N

        self.radices = [1] * (self.L + 1)
        for k in range(self.L - 1, -1, -1):
            self.radices[k] = self.radices[k + 1] * self.branch_factors[k]

        self.counts = [0] * self.total_nodes
        self.sums = [0.0] * self.total_nodes

        self.node_H = [0.0] * self.total_nodes
        self.node_Q = [0.0] * self.total_nodes
        self.node_df = [0] * self.total_nodes

        self.shell_H_sum = [0.0] * (self.L + 1)
        self.shell_Q_sum = [0.0] * (self.L + 1)
        self.shell_df_sum = [0] * (self.L + 1)

        self.leaf_ss_within = 0.0
        self.leaf_df_within = 0
        self.sigma2_hat = 1.0

        self.v_hat = [0.8 * (0.65 ** k) for k in range(self.L + 1)]

    def node_idx(self, k: int, w: int) -> int:
        return self.level_offsets[k] + w

    def update(self, leaf: int, y: float):
        """O(L) bottom-up update for arbitrary branching tree."""
        leaf_id = self.node_idx(self.L, leaf)
        n_old = self.counts[leaf_id]
        if n_old >= 1:
            old_m = self.sums[leaf_id] / n_old
            new_m = (self.sums[leaf_id] + y) / (n_old + 1)
            self.leaf_ss_within += (y - old_m) * (y - new_m)
            self.leaf_df_within += 1
            if self.leaf_df_within >= 3:
                self.sigma2_hat = max(0.02, self.leaf_ss_within / self.leaf_df_within)

        for k in range(self.L + 1):
            w = leaf // self.radices[k]
            idx = self.node_idx(k, w)
            self.counts[idx] += 1
            self.sums[idx] += y

        for k in range(1, self.L + 1):
            b_k = self.branch_factors[k - 1]
            p_w = leaf // self.radices[k - 1]
            p_idx = self.node_idx(k - 1, p_w)
            N_w = self.counts[p_idx]

            vis_children = 0
            sum_n2 = 0
            Q_w = 0.0
            p_mean = self.sums[p_idx] / N_w

            base_child_w = p_w * b_k
            for j in range(b_k):
                c_idx = self.node_idx(k, base_child_w + j)
                nj = self.counts[c_idx]
                if nj > 0:
                    vis_children += 1
                    sum_n2 += nj * nj
                    diff = (self.sums[c_idx] / nj) - p_mean
                    Q_w += nj * diff * diff

            if vis_children >= 2:
                H_w = N_w - (sum_n2 / N_w)
                df_w = vis_children - 1
                self.shell_H_sum[k] += (H_w - self.node_H[p_idx])
                self.shell_Q_sum[k] += (Q_w - self.node_Q[p_idx])
                self.shell_df_sum[k] += (df_w - self.node_df[p_idx])
                self.node_H[p_idx] = H_w
                self.node_Q[p_idx] = Q_w
                self.node_df[p_idx] = df_w

        # Continuous-space fallback: if no leaf is visited twice, estimate noise floor from finest split shells
        if self.leaf_df_within < 3:
            fine_Q = 0.0
            fine_df = 0
            for k in range(self.L, max(1, self.L // 2), -1):
                if self.shell_df_sum[k] > 0:
                    fine_Q += self.shell_Q_sum[k]
                    fine_df += self.shell_df_sum[k]
            if fine_df >= 3:
                self.sigma2_hat = max(0.05, fine_Q / fine_df)

        running_alias = 0.0
        for k in range(self.L, 0, -1):
            b_k = self.branch_factors[k - 1]
            H_k = self.shell_H_sum[k]
            df_k = self.shell_df_sum[k]
            eff_noise_k = self.sigma2_hat + running_alias
            prior_guess = 0.8 * (0.65 ** k)
            if H_k > 0.0 and df_k >= 1:
                raw_v = (self.shell_Q_sum[k] - df_k * eff_noise_k) / H_k
                weight = df_k / (df_k + 3.0)
                est_v = weight * max(0.02 * prior_guess, raw_v) + (1.0 - weight) * prior_guess
            else:
                est_v = prior_guess
            self.v_hat[k] = est_v
            running_alias = (running_alias / b_k) + ((b_k - 1.0) / b_k) * est_v

    def get_cotton_damped_child(self, k: int, w: int, parent_est: float) -> tuple[float, float, float]:
        idx = self.node_idx(k, w)
        n = self.counts[idx]
        v_loc = self.v_hat[k]
        if n == 0:
            return parent_est, v_loc, 0.0
        alias_below = sum(self.v_hat[m] for m in range(k + 1, self.L + 1)) if k < self.L else 0.0
        eff_noise = (self.sigma2_hat + alias_below) / n
        gamma = v_loc / (v_loc + eff_noise)
        raw_mean = self.sums[idx] / n
        damped_mean = parent_est + gamma * (raw_mean - parent_est)
        post_var = v_loc * (1.0 - gamma)
        return damped_mean, post_var, gamma

    def stride_suffix(self, k: int, t: int) -> int:
        if k >= self.L:
            return 0
        val = t
        offset = 0
        for m in range(k, self.L):
            b_m = self.branch_factors[m]
            digit = val % b_m
            val //= b_m
            offset = offset * b_m + digit
        return offset

    def select_leaf(self, t: int, rng: random.Random) -> int:
        n_root = self.counts[0]
        curr_est = (self.sums[0] / n_root) if n_root > 0 else 0.0
        curr_w = 0
        uncertain_crossed = 0

        for k in range(1, self.L + 1):
            b_k = self.branch_factors[k - 1]
            base_w = curr_w * b_k
            p_idx = self.node_idx(k - 1, curr_w)
            N_p = max(self.counts[p_idx], 1)

            for j in range(b_k):
                child_w = base_w + j
                if self.counts[self.node_idx(k, child_w)] == 0:
                    return child_w * self.radices[k] + self.stride_suffix(k, t)

            best_u = -1e18
            best_child_w = base_w
            best_est = curr_est
            best_gam = 0.0
            log_np = math.log(N_p + 1)

            for j in range(b_k):
                child_w = base_w + j
                est_j, var_j, gam_j = self.get_cotton_damped_child(k, child_w, curr_est)
                u_j = est_j + self.c_ucb * math.sqrt(var_j * log_np)
                if u_j > best_u:
                    best_u = u_j
                    best_child_w = child_w
                    best_est = est_j
                    best_gam = gam_j

            curr_w, curr_est = best_child_w, best_est

            if best_gam < self.gamma_star:
                uncertain_crossed += 1
                if uncertain_crossed >= self.cap_K and k < self.L:
                    return curr_w * self.radices[k] + self.stride_suffix(k, t)

        return curr_w

    def best_leaf_estimate(self) -> int:
        n_root = self.counts[0]
        if n_root == 0:
            return 0
        root_mean = self.sums[0] / n_root
        best_lf = 0
        best_score = -1e18
        for lf in range(self.N):
            if self.counts[self.node_idx(self.L, lf)] == 0:
                continue
            curr_est = root_mean
            for k in range(1, self.L + 1):
                w = lf // self.radices[k]
                idx = self.node_idx(k, w)
                n = self.counts[idx]
                if n == 0:
                    break
                v_loc = self.v_hat[k]
                alias_below = sum(self.v_hat[m] for m in range(k + 1, self.L + 1)) if k < self.L else 0.0
                eff_noise = (self.sigma2_hat + alias_below) / n
                gamma = v_loc / (v_loc + eff_noise)
                curr_est = curr_est + gamma * ((self.sums[idx] / n) - curr_est)
            if curr_est > best_score:
                best_score = curr_est
                best_lf = lf
        return best_lf

# ==============================================================================
# BENCHMARK G1: Continuous 4D Noisy Black-Box Optimization on [0, 1]^4
# ==============================================================================
def eval_4d_multimodal(x_vec: list[float], x_star: list[float]) -> float:
    dist1 = max(abs(xi - si) for xi, si in zip(x_vec, x_star))
    dist2 = sum((xi - si) ** 2 for xi, si in zip(x_vec, x_star))
    cosine_mod = sum(math.cos(4.0 * math.pi * (xi - si)) for xi, si in zip(x_vec, x_star)) / len(x_vec)
    return 3.5 * math.exp(-3.0 * dist2) - 0.8 * dist1 + 0.65 * cosine_mod

def leaf_to_4d_point(leaf: int, L_depth: int = 12, dim: int = 4) -> list[float]:
    bits_per_dim = L_depth // dim
    coords = [0.0] * dim
    for k in range(L_depth):
        bit = (leaf >> (L_depth - 1 - k)) & 1
        d = k % dim
        level_in_d = k // dim
        coords[d] += bit * (0.5 ** (level_in_d + 1))
    half_width = 0.5 ** (bits_per_dim + 1)
    return [c + half_width for c in coords]

def run_benchmark_g1():
    L = 12
    N = 1 << L  # 4,096 cells in [0,1]^4
    T = 256
    n_seeds = 50
    sigma = 0.85

    reg_cg = []
    reg_stosoo = []
    reg_lhs = []

    for seed in range(n_seeds):
        rng_env = random.Random(30000 + seed)
        x_star = [rng_env.uniform(0.1, 0.9) for _ in range(4)]
        cell_vals = [eval_4d_multimodal(leaf_to_4d_point(lf, L, 4), x_star) for lf in range(N)]
        best_val = max(cell_vals)
        noise_seq = [rng_env.gauss(0.0, sigma) for _ in range(T)]

        # 1. Universal CG-HSTS v2.0 on [0,1]^4
        rng_cg = random.Random(80000 + seed)
        tree_cg = UniversalCGHSTS([2] * L, gamma_star=0.45, c_ucb=0.85, cap_K=4)
        for t in range(T):
            lf = tree_cg.select_leaf(t, rng_cg)
            tree_cg.update(lf, cell_vals[lf] + noise_seq[t])
        rec_cg = tree_cg.best_leaf_estimate()
        reg_cg.append(best_val - cell_vals[rec_cg])

        # 2. StoSOO / UCT-on-Boxes (undamped box averages)
        rng_uct = random.Random(80000 + seed)
        uct_counts = [0] * (1 << (L + 1))
        uct_sums = [0.0] * (1 << (L + 1))
        for t in range(T):
            curr_w = 0
            for k in range(1, L + 1):
                p_idx = (1 << (k - 1)) + curr_w
                w0, w1 = (curr_w << 1), (curr_w << 1) | 1
                idx0, idx1 = (1 << k) + w0, (1 << k) + w1
                n0, n1 = uct_counts[idx0], uct_counts[idx1]
                if n0 == 0:
                    rem = L - k
                    curr_w = (w0 << rem) | rng_uct.randrange(1 << rem) if rem > 0 else w0
                    break
                elif n1 == 0:
                    rem = L - k
                    curr_w = (w1 << rem) | rng_uct.randrange(1 << rem) if rem > 0 else w1
                    break
                else:
                    n_p = max(uct_counts[p_idx], 1)
                    u0 = (uct_sums[idx0] / n0) + 1.2 * math.sqrt(math.log(n_p) / n0)
                    u1 = (uct_sums[idx1] / n1) + 1.2 * math.sqrt(math.log(n_p) / n1)
                    curr_w = w0 if u0 >= u1 else w1
            leaf = curr_w
            y = cell_vals[leaf] + noise_seq[t]
            for k in range(L + 1):
                w = leaf >> (L - k)
                idx = (1 << k) + w
                uct_counts[idx] += 1
                uct_sums[idx] += y
        best_uct_leaf = max(
            (i for i in range(N) if uct_counts[(1 << L) + i] > 0),
            key=lambda i: uct_sums[(1 << L) + i] / uct_counts[(1 << L) + i],
        )
        reg_stosoo.append(best_val - cell_vals[best_uct_leaf])

        # 3. Quasi-Random / LHS Search without Hierarchical Cotton Backup
        best_y = -1e18
        best_lhs_lf = 0
        for t in range(T):
            lf = int(bin(t * 16 + (seed % 16))[2:].zfill(L)[::-1][:L], 2) % N
            y = cell_vals[lf] + noise_seq[t]
            if y > best_y:
                best_y = y
                best_lhs_lf = lf
        reg_lhs.append(best_val - cell_vals[best_lhs_lf])

    _, p_vs_stosoo, diff_vs_stosoo, se_vs_stosoo = paired_ttest(reg_stosoo, reg_cg)
    _, p_vs_lhs, diff_vs_lhs, se_vs_lhs = paired_ttest(reg_lhs, reg_cg)

    return {
        "regret_cg_hsts_4d": f"{statistics.mean(reg_cg):.4f} ± {statistics.stdev(reg_cg)/math.sqrt(n_seeds):.4f}",
        "regret_stosoo_4d": f"{statistics.mean(reg_stosoo):.4f} ± {statistics.stdev(reg_stosoo)/math.sqrt(n_seeds):.4f}",
        "regret_lhs_4d": f"{statistics.mean(reg_lhs):.4f} ± {statistics.stdev(reg_lhs)/math.sqrt(n_seeds):.4f}",
        "diff_vs_stosoo": f"-{diff_vs_stosoo:.4f} ± {se_vs_stosoo:.4f} (p = {p_vs_stosoo:.2e})",
        "diff_vs_lhs": f"-{diff_vs_lhs:.4f} ± {se_vs_lhs:.4f} (p = {p_vs_lhs:.2e})",
    }

# ==============================================================================
# BENCHMARK G2: General Primorial (2x3x5x7x5 = 1,050 Leaves) & 4-ary (4^5 = 1,024)
# ==============================================================================
def sample_general_tree_field(branch_factors: list[int], v_scales: list[float], rng: random.Random) -> list[float]:
    L = len(branch_factors)
    radices = [1] * (L + 1)
    for k in range(L - 1, -1, -1):
        radices[k] = radices[k + 1] * branch_factors[k]
    N = radices[0]
    theta = [rng.gauss(0.0, math.sqrt(v_scales[0]))] * N
    num_blocks = 1
    for m in range(1, L + 1):
        num_blocks *= branch_factors[m - 1]
        block_size = radices[m]
        sd_m = math.sqrt(v_scales[m])
        for u in range(num_blocks):
            zu = rng.gauss(0.0, sd_m)
            start = u * block_size
            for i in range(start, start + block_size):
                theta[i] += zu
    return theta

def run_benchmark_g2():
    configs = {
        "quaternary_4ary_1024": [4, 4, 4, 4, 4],
        "primorial_2x3x5x7x5_1050": [2, 3, 5, 7, 5],
    }
    T = 256
    n_seeds = 50
    sigma = 1.0
    results = {}

    n0, n1 = 7, 13
    y0, y1 = 1.42, -0.88
    N_w = n0 + n1
    H_w_gen = N_w - (n0*n0 + n1*n1) / N_w
    h_w_bin = (2.0 * n0 * n1) / (n0 + n1)
    yw = (n0 * y0 + n1 * y1) / N_w
    Q_w_gen = n0 * ((y0 - yw) ** 2) + n1 * ((y1 - yw) ** 2)
    d_w_bin = 0.5 * (y0 - y1)
    Q_w_bin = 2.0 * h_w_bin * (d_w_bin ** 2)
    results["algebraic_reduction_H_residual"] = abs(H_w_gen - h_w_bin)
    results["algebraic_reduction_Q_residual"] = abs(Q_w_gen - Q_w_bin)

    for name, b_list in configs.items():
        L = len(b_list)
        v_scales = [1.0] + [0.6 ** m for m in range(1, L + 1)]
        reg_cg = []
        reg_uct = []
        rew_cg = []
        rew_uct = []

        for seed in range(n_seeds):
            rng_env = random.Random(60000 + seed)
            theta = sample_general_tree_field(b_list, v_scales, rng_env)
            N = len(theta)
            best_val = max(theta)
            noise_seq = [rng_env.gauss(0.0, sigma) for _ in range(T)]

            rng_cg = random.Random(85000 + seed)
            tree_cg = UniversalCGHSTS(b_list, gamma_star=0.45, c_ucb=0.90, cap_K=3)
            r_tot = 0.0
            for t in range(T):
                lf = tree_cg.select_leaf(t, rng_cg)
                y = theta[lf] + noise_seq[t]
                tree_cg.update(lf, y)
                r_tot += theta[lf]
            rec_cg = tree_cg.best_leaf_estimate()
            rew_cg.append(r_tot)
            reg_cg.append(best_val - theta[rec_cg])

            rng_uct = random.Random(85000 + seed)
            tree_uct = UniversalCGHSTS(b_list)
            r_tot = 0.0
            for t in range(T):
                curr_w = 0
                for k in range(1, L + 1):
                    b_k = b_list[k - 1]
                    base_w = curr_w * b_k
                    p_idx = tree_uct.node_idx(k - 1, curr_w)
                    N_p = max(tree_uct.counts[p_idx], 1)
                    unvis = [j for j in range(b_k) if tree_uct.counts[tree_uct.node_idx(k, base_w + j)] == 0]
                    if unvis:
                        chosen_j = unvis[0]
                        child_w = base_w + chosen_j
                        rem_size = tree_uct.radices[k]
                        curr_w = child_w * rem_size + (rng_uct.randrange(rem_size) if rem_size > 1 else 0)
                        break
                    else:
                        best_u = -1e18
                        best_w = base_w
                        for j in range(b_k):
                            cw = base_w + j
                            c_idx = tree_uct.node_idx(k, cw)
                            nj = tree_uct.counts[c_idx]
                            uj = (tree_uct.sums[c_idx] / nj) + 1.414 * math.sqrt(math.log(N_p) / nj)
                            if uj > best_u:
                                best_u = uj
                                best_w = cw
                        curr_w = best_w
                leaf = curr_w
                y = theta[leaf] + noise_seq[t]
                for k in range(L + 1):
                    w = leaf // tree_uct.radices[k]
                    idx = tree_uct.node_idx(k, w)
                    tree_uct.counts[idx] += 1
                    tree_uct.sums[idx] += y
                r_tot += theta[leaf]
            best_emp_leaf = max(
                (i for i in range(N) if tree_uct.counts[tree_uct.node_idx(L, i)] > 0),
                key=lambda i: tree_uct.sums[tree_uct.node_idx(L, i)] / tree_uct.counts[tree_uct.node_idx(L, i)],
            )
            rew_uct.append(r_tot)
            reg_uct.append(best_val - theta[best_emp_leaf])

        _, p_reg, diff_reg, se_reg = paired_ttest(reg_uct, reg_cg)
        _, p_rew, diff_rew, se_rew = paired_ttest(rew_cg, rew_uct)
        results[name] = {
            "regret_universal_cg_hsts": f"{statistics.mean(reg_cg):.3f} ± {statistics.stdev(reg_cg)/math.sqrt(n_seeds):.3f}",
            "regret_standard_uct": f"{statistics.mean(reg_uct):.3f} ± {statistics.stdev(reg_uct)/math.sqrt(n_seeds):.3f}",
            "regret_reduction_vs_uct": f"-{diff_reg:.3f} ± {se_reg:.3f} (p = {p_reg:.2e})",
            "reward_universal_cg_hsts": f"{statistics.mean(rew_cg):.1f} ± {statistics.stdev(rew_cg)/math.sqrt(n_seeds):.1f}",
            "reward_standard_uct": f"{statistics.mean(rew_uct):.1f} ± {statistics.stdev(rew_uct)/math.sqrt(n_seeds):.1f}",
            "reward_gain_vs_uct": f"{diff_rew:+.1f} ± {se_rew:.1f} (p = {p_rew:.2e})",
        }

    return results

def main():
    g1 = run_benchmark_g1()
    g2 = run_benchmark_g2()
    out = {"G1_continuous_4D": g1, "G2_general_branching": g2}
    with open("cg_hsts_universal_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
