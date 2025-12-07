# American_Option_Pricing_LSMC
American option pricing using Longstaff–Schwartz Least Squares Monte Carlo with binomial benchmark and volatility scenario analysis.
# American Option Pricing using Least Squares Monte Carlo (LSMC)

This repository implements **American option pricing** using the
**Longstaff–Schwartz Least Squares Monte Carlo (LSMC)** algorithm,
and benchmarks it against a high-resolution **Cox–Ross–Rubinstein (CRR) binomial tree**.

The goals of this project are:

- To show a full, working implementation of the LSMC method for American options.
- To validate the Monte Carlo estimator by:
  - Comparing prices to a binomial benchmark.
  - Studying convergence as the number of simulated paths increases.
  - Measuring relative pricing error across **20+ volatility scenarios**.
- To demonstrate a clear understanding of the **theory + numerical methods** behind American option pricing.

---

## 1. Background

### 1.1 American vs European options

- A **European** option can only be exercised at maturity.
- An **American** option can be exercised at any time up to maturity.

Because of this early-exercise feature, there is **no simple Black–Scholes closed form**
for general American options. Numerical methods are required.

Two commonly used approaches:

- **Binomial tree (CRR model)** – a discrete-time lattice where the underlying price
  moves up/down at each step; option value is computed by backward induction. :contentReference[oaicite:0]{index=0}  
- **Least Squares Monte Carlo (LSMC)** – simulate many price paths, then use **regression**
  to approximate the continuation value and decide whether to exercise or continue
  at each time step. :contentReference[oaicite:1]{index=1}  

This project implements **both**:

- CRR binomial → used as a **reference “ground truth”**
- LSMC → fast, flexible Monte Carlo method that scales to path-dependent / multi-factor models

---

## 2. Model Setup

We assume the underlying asset \( S_t \) follows **risk-neutral GBM**:

\[
dS_t = (r - q) S_t\, dt + \sigma S_t\, dW_t
\]

where:

- \( r \) – risk-free rate  
- \( q \) – continuous dividend yield (0 in this simple version)  
- \( \sigma \) – volatility  
- \( W_t \) – standard Brownian motion  

We work under the risk-neutral measure so the discounted option price is a martingale.

In code, GBM paths are simulated by the log-Euler scheme:

\[
S_{t+\Delta t} = S_t \exp\big((r - q - \tfrac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}\, Z\big),
\quad Z \sim \mathcal{N}(0,1).
\]

---

## 3. Methods

### 3.1 CRR Binomial Tree (Benchmark)

The **Cox–Ross–Rubinstein** binomial model builds a recombining tree with
up/down movements per step: :contentReference[oaicite:2]{index=2}

- Time step: \( \Delta t = T / N \)
- Up factor: \( u = e^{\sigma\sqrt{\Delta t}} \)
- Down factor: \( d = 1/u \)
- Risk-neutral probability:

\[
p = \frac{e^{r\Delta t} - d}{u - d}.
\]

Algorithm:

1. Compute terminal asset prices \( S_T \) at all nodes.
2. Compute terminal payoffs \( \max(K - S_T, 0) \) for a **put** (or \( \max(S_T - K, 0) \) for a call).
3. Step backwards:
   - **Continuation value** at each node:

     \[
     C = e^{-r\Delta t}\big(p V_{\text{up}} + (1-p) V_{\text{down}}\big)
     \]

   - **Exercise value**:

     \[
     E = \max(K - S, 0)
     \]

   - American option value at that node: \( V = \max(E, C) \).

Using a large number of steps (e.g. 400–500), the binomial price acts as a **high-accuracy reference**.

---

### 3.2 Longstaff–Schwartz LSMC

LSMC prices American options by approximating the **continuation value** via regression. :contentReference[oaicite:3]{index=3}  

Key idea: at any exercise date, the holder compares:

- Immediate exercise payoff \( h(S_t) \)
- Expected discounted value of continuation \( \mathbb{E}[V_{t+1} \mid S_t] \)

We approximate this conditional expectation using least squares.

Algorithm (put option):

1. **Simulate Monte Carlo paths** for \( S_t \), \( t = 0, \dots, N \).  
2. Compute intrinsic payoff at each time:

   \[
   h(S_t) = \max(K - S_t, 0).
   \]

3. Initialize cashflows at maturity: \( CF_i(T) = h(S_i(T)) \).
4. For each time step \( t = N-1, \dots, 1 \) (backward):
   - Discount cashflows one step: \( CF_i(t) = e^{-r\Delta t} CF_i(t+1) \).
   - Restrict to **in-the-money** paths \( \{ i : h(S_i(t)) > 0 \} \).
   - Regress discounted CF on basis functions of \( S_t \):

     \[
     CF_i(t) \approx \beta_0 + \beta_1 S_i(t) + \beta_2 S_i(t)^2 + \dots
     \]

     giving an estimate of continuation value \( \hat{C}(S_t) \).
   - For ITM paths, compare:

     - exercise: \( E_i = h(S_i(t)) \)
     - continuation: \( \hat{C}_i = \hat{C}(S_i(t)) \)

     If \( E_i > \hat{C}_i \), exercise → set \( CF_i(t) = E_i \);  
     otherwise keep continuation \( CF_i(t) \).
5. After stepping back to \( t=0 \), discount once more and average:

   \[
   \text{Price} \approx \frac{1}{M} \sum_{i=1}^M e^{-r\Delta t} CF_i(1).
   \]

In this implementation, the **basis functions** are simple polynomials  
\[ [1, S_t, S_t^2] \]  
which work surprisingly well for a single-asset American put.

---

## 4. Implementation Details

All logic lives inside a single notebook:

- `American_Option_Pricing_LSMC.ipynb`

Core components:

1. `simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, q, seed)`  
   - Vectorized GBM path simulator (no Python loops over paths).  
2. `price_american_binomial(...)`  
   - Clean CRR tree with early-exercise check at each node.  
3. `lsmc_american_option(...)`  
   - Longstaff–Schwartz pricer with:
     - in-the-money filtering,
     - polynomial basis regression,
     - backward induction over time.

Vectorization + NumPy linear algebra make the method efficient
even with **20,000–50,000 Monte Carlo paths**.

---

## 5. Experiments & Results

All results below are generated under a simple base configuration:

- \( S_0 = 100 \)
- \( K = 100 \) (at-the-money)
- \( T = 1 \) year
- \( r = 5\% \)
- \( q = 0 \)
- American **put** option

### 5.1 Single-Case Comparison (Binomial vs LSMC)

- Binomial American put (500 steps): **~6.09**  
- LSMC (20,000 paths, 50 time steps): **close to binomial**,  
  with relative error on the order of **a few tenths of a percent**.

This confirms that the LSMC implementation is correctly capturing early-exercise behavior.

---

### 5.2 Convergence with Number of Paths

We evaluate LSMC for path counts:

\[
M \in \{2000, 5000, 10000, 20000, 50000\}.
\]

The **convergence plot** (`lsmc_convergence_paths.png`) shows:

- As \( M \) increases, the LSMC price **converges towards the binomial reference**.
- For \( M \ge 10{,}000 \), the prices are very close to the binomial benchmark
  (error ~0.3–0.5% in my runs).

This is exactly what Monte Carlo theory predicts: more paths → lower variance of the estimator.

---

### 5.3 Volatility Scenario Analysis (20+ sigmas)

To test robustness across market regimes, we vary volatility:

\[
\sigma \in [0.10, 0.60] \text{ in 20 evenly spaced steps}.
\]

For each \( \sigma \):

1. Compute **binomial price** (400 steps).  
2. Compute **LSMC price** (20,000 paths, 50 time steps).  
3. Record **relative error**:

\[
\text{error}(\sigma) = \frac{|P_{\text{LSMC}} - P_{\text{binomial}}|}{P_{\text{binomial}}} \times 100\%.
\]

The **error plot** (`vol_error_lsmc_vs_binomial.png`) shows:

- Errors are extremely small, typically in the range **0.1% – 0.5%**.
- Errors remain low **across all 20 volatility scenarios**.
- Slightly higher error at very high vol (σ ≈ 0.5–0.6) is expected,
  since the early-exercise boundary becomes more nonlinear and harder to approximate
  with simple polynomial basis functions.

Overall, this demonstrates that the LSMC implementation is:

- **Numerically stable**, and  
- **Accurate** across different volatility regimes.

---

## 6. How to Run

### 6.1 Requirements

- Python 3.x
- `numpy`
- `pandas`
- `matplotlib`
- `tqdm`

You can install via:

```bash
pip install -r requirements.txt
