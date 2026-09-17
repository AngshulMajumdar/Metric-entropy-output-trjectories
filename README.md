# Metric entropy rates of output trajectories — code and certificates

Code, integer certificates and reference outputs for the paper

> A. Majumdar, *Metric Entropy Rates of Output Trajectories for Input-Injective, Incrementally
> Contracting Systems with Finite Memory*.

Every certified number in the paper can be re-derived from this repository without an eigensolver.

## Requirements

Python ≥ 3.10 and the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

All scripts are run from the repository root; they read and write files relative to it.

## Result categories

The paper labels every numerical result as one of

* **(C) certified** – rigorous enclosures. Interval arithmetic (`mpmath.iv`, outward rounding) is used for all
  enclosures, and Collatz–Wielandt inequalities are checked in exact integer/rational arithmetic
  against the stored certificate vectors in `certs/`.
  *Exception:* `run_E2.py` (Channel E, 2-D state) uses vectorized double precision with every enclosure
  widened by 1e-12 and every ratio by a factor 1 ± 1e-9; it is rigorous under the platform's ≤ 1 ulp
  error bound for `tanh`, `sin`, `cos` (see Appendix A of the paper).
* **(N) numerical estimates** – Feynman–Kac particle estimates (not certified).
* **(V) empirical validation** – Monte Carlo checks of identities and constants.

All values are in bits per sample unless stated otherwise.

## Scripts and the results they produce

| Script | Category | Paper | Output |
|---|---|---|---|
| `run_S.py v` | C | Table 1(a) | `S_v.pkl` |
| `run_S.py eps` | C | Table 1(b), Fig. 1 left | `S_eps.pkl` |
| `run_fam.py a` | C | Table 2 (Channel S_α) | `fam_a.pkl` |
| `run_M.py` | C | Table 3, Fig. 1 right | `M.pkl` |
| `run_fam.py beta 256` | C | Table 4 (Channel M_β) | `fam_beta.pkl` |
| `run_fam.py Mfine` | C | Channel M, ρ = 1/2, 512 cells | `fam_Mfine.pkl` |
| `thm_consts.py` | formula | Table 5 | `thm_consts.pkl` |
| `run_fam.py kx2 32`, `run_fam.py kx2 48` | C | Channel K | `fam_kx2.pkl` |
| `run_fam.py quad` | C | Table 6 (Channel Q) | `fam_quad.pkl` |
| `run_E.py` | C | Table 7 (relaxations of Channel E) | `E_rel.pkl` |
| `run_E2.py n m` | C* | Table 7 (2-D enclosures; n m = 32 128, 64 256, 96 384, 128 512) | `E2.pkl` |
| `volcheck.py` | V | Table 8 | `volcheck.pkl` |
| `constcheck.py` | V | check of the constants in Lemmas 3 and 5 on Channel E | stdout |
| `particles.py M`, `particles.py E` | N | Section 8, Table 9 (last row) | `particles.pkl` |
| `pstudy.py N`, `pstudy.py T` | N | Table 9 | `pstudy.pkl` |
| `pk2.py` | N | Channel K particle estimate | `pk2.pkl` |
| `fig.py` | – | Figure 1 | `figs.pdf` |

Shared modules: `core.py` (grids, step kernels, exact Collatz–Wielandt checks, outward-rounded
logarithms) and `ivlib.py` (interval enclosures of `tanh`, `|·|`, `sech²`).

Without the environment variable `VERIFY`, the certified scripts recompute Perron vectors and
**overwrite** the certificate files in `certs/`. With `VERIFY=1` they load the stored certificates
and use no eigensolver.

## Verification

```
python3 verify.py            # all 14 certified jobs (about 20 minutes)
python3 verify.py 0 1 2 3    # selected jobs
```

`verify.py` runs each certified script with `VERIFY=1`. Each run rebuilds every step matrix from
interval enclosures, rechecks every inequality from `certs/*.npy`, and compares the result with the
reference output in `ref/`. The job indices are:

| # | Job | # | Job |
|---|---|---|---|
| 0 | `run_S.py v` | 7 | `run_fam.py kx2 48` |
| 1 | `run_S.py eps` | 8 | `run_fam.py beta 256` |
| 2 | `run_M.py` | 9 | `run_fam.py Mfine` |
| 3 | `run_E.py` | 10 | `run_E2.py 32 128` |
| 4 | `run_fam.py quad` | 11 | `run_E2.py 64 256` |
| 5 | `run_fam.py a` | 12 | `run_E2.py 96 384` |
| 6 | `run_fam.py kx2 32` | 13 | `run_E2.py 128 512` |

`MANIFEST.sha256` lists the SHA-256 hashes of all certificates, reference outputs and scripts:

```
sha256sum -c MANIFEST.sha256
```

## Layout

```
certs/     integer certificate vectors (*.npy): Perron vectors scaled by 2^40, rounded up (upper
           bounds) or down (lower bounds)
ref/       reference outputs compared by verify.py
*.pkl      outputs of the most recent run of each script
*.py       scripts (see table above)
```

## Channels

| Name | System | Input set |
|---|---|---|
| S_α | y⁺ = α y + 0.5(1 + 0.5 cos y) sin 3x (Channel S: α = 0.6) | [−1, 1] |
| M_β | y⁺ = 0.6 y + tanh(2x + β x₋₁) (Channel M: β = 0.3) | [−ρ, ρ] |
| K | y⁺ = 0.6 y + tanh(2x + 0.2 x₋₁ + 0.1 x₋₂) | [−1/2, 1/2] |
| Q | y⁺ = g(y) + x + 0.4 x x₋₁ + 0.3 x² | [−ρ, ρ] |
| E | y⁺ = 0.5 y + 0.2 sin y + (1 + 0.3 cos y) tanh(x + 0.1 x₋₁) | [−1/2, 1/2] |
