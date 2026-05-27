# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

M1 research project (Université Paris Cité, MMAS 2025–2026) by Arthur Conche and Abdel Niang, supervised by Camille Pouchol. Subject: *Estimation in a linear inverse problem with hidden variables — likelihood, EM algorithm and error analysis* — a toy model of cryo-EM reconstruction.

The repo has two deliverables that evolve in parallel:

- **`doc.tex` → `doc.pdf`** — the written report (French, ~1127 lines of LaTeX).
- **`code python/`** — three notebooks for Cas 1 (1D translations), each self-contained with its own data generation:
  - `graddesc_cas1.ipynb` — gradient descent with Armijo backtracking; produces `signaux_1D.pdf` and `simulation_parametres.pdf`.
  - `em_cas1.ipynb` — EM algorithm (loop-based, readable); validates the E/M-step derivations.
  - `analyse_erreur.ipynb` — Monte-Carlo study of ε_n(σ); uses a fully vectorised EM; produces `erreur_vs_n.pdf` and `erreur_vs_sigma.pdf`.

`Projet_M1.pdf` is the assignment brief from the supervisor. `Exemples_Rapports_M1_MMAS/` contains six PDFs of prior years' reports as style/scope references — they are **not** part of this project's source.

## Build & run

LaTeX (full rebuild with bibliography, run from repo root):

```bash
latexmk -pdf doc.tex          # build doc.pdf
latexmk -c                    # clean aux files (.aux .fls .fdb_latexmk .log .out .toc .synctex.gz)
latexmk -C                    # also remove doc.pdf
```

The `.tex` references `\bibliography{references}`. `references.bib` exists at the repo root with entries for: Dempster/Laird/Rubin 1977 (EM), McLachlan & Krishnan, Wu 1983 (convergence EM), Natterer 2001 (Radon), Frank 2006 (cryo-EM), Singer & Shkolnisky 2011 (cryo-EM), Cavalier 2008 (inverse problems), Boyd & Vandenberghe (optimization).

Python / notebooks (venv lives at `.venv/`, Python 3.14):

```bash
source .venv/bin/activate
jupyter lab                   # or: jupyter notebook
```

Installed: `numpy`, `scipy` (via deps), `matplotlib`, `pandas`, `jupyter`/`jupyterlab`, `ipywidgets`. No `requirements.txt` or `pyproject.toml` — packages were added ad-hoc with `pip install`.

Note the directory name `code python/` contains a space — quote it in shell commands (`cd "code python"`).

## Mathematical architecture (needed before editing `doc.tex` or implementing algorithms)

The whole report is built around one model, referenced everywhere as **(M)** with `\eqref{eq:modele}`:

$$Y_i = A_{Z_i}\theta + \sigma E_i, \quad i = 1,\dots,n$$

with $Z_i$ i.i.d. hidden, $E_i \sim \mathcal{N}(0, \mathrm{Id}_p)$, and $A_z \in L(E, \mathbb{R}^p)$. Two concrete cases instantiate $A_z$:

- **Cas 1 — 1D translations on the torus.** $E = \mathrm{Vect}(e^{ikx})$, $A_z$ = shifted-grid evaluation. Diagonal in the Fourier basis (translation by $z$ = phase modulation $e^{-ikz}$). This is the primary testbed — every algorithm is first developed and validated here.
- **Cas 2 — 2D rotation + Radon projection.** Analog of cryo-EM proper; treated as an extension (Section 7), not the main object.

Section flow in `doc.tex` mirrors the methodology: (§1) setup → (§2) reference case with $A_z = A$ constant, closed-form MLE → (§3) marginal likelihood over $Z$, identifiability obstruction → (§4) gradient descent **and** EM as two solvers, EM's M-step is a quadratic linear system → (§5) numerical study of $\varepsilon_n(\sigma) = \mathbb{E}\|\hat\theta_n - \theta\|^2$ — **note the metric is translation-invariant** (`\inf_{a \in [0,2\pi]} \|\tau_a \hat\theta - \vartheta\|`), because the standard norm is unidentifiable in Cas 1 → (§6) regularisation for ill-posed problems → (§7) extensions.

When implementing algorithms in notebooks, the Fourier diagonalization of Cas 1 is the key efficiency lever: do not rotate in the spatial domain.

## Notebook architecture (`code python/`)

All three notebooks share the same global parameters: `n` (observations), `p` (grid points, 32), `N` (Fourier frequencies, 3), `d = 2N+1` (signal dimension, 7), `sigma` (noise), `n_z` (quadrature points). Core primitives are duplicated across notebooks (each is self-contained).

### Core primitives (all three notebooks)

| Function | Role |
|---|---|
| `build_phi(z)` | Returns the $(p \times d)$ evaluation matrix $\Phi(z)$; `build_phi(z) @ u` gives $A_z\theta$ on the grid |
| `build_phiu(z, u)` | Thin wrapper returning `(p,)` signal vector |
| `qi(z, y, u)` | Computes $\frac{1}{2\sigma^2}\|y - \Phi(z)u\|^2$, clipped to 745 to avoid float64 overflow |
| `log_likelihood(u)` | Marginal log-likelihood via Riemann quadrature + `scipy.logsumexp` |
| `translation_error(u_hat, u_true)` | Translation-invariant relative error via FFT cross-correlation |

### `graddesc_cas1.ipynb`

| Function | Role |
|---|---|
| `compute_gradient(u)` | Gradient $\nabla_\theta \ell_n$ using posterior weights $w_{i,z}$ (loop over `n × n_z`) |
| `gradient_ascent_armijo(u0, ...)` | Gradient ascent with Armijo backtracking; returns `(u_hat, hist_ll, hist_grad)` |

### `em_cas1.ipynb`

| Function | Role |
|---|---|
| `e_step(u)` | Returns posterior weight matrix $(n \times n_z)$; uses softmax stabilisation |
| `m_step(weights)` | Solves the weighted least-squares system $Au = b$ via `np.linalg.solve`; returns $u^{k+1}$ |
| `em(u0, ...)` | Full EM loop; returns `(u_hat, hist_ll, hist_u)`; converges when ll gain < `tol` |

**Initialisation note:** start EM from `u0 = 0.5 * np.random.randn(d)` (reduced amplitude). A large `u0` saturates all quadrature weights at the first E-step (all `qi` hit the clip ceiling), making the weights uniform and locking EM on a bad fixed point.

### `analyse_erreur.ipynb`

Fully vectorised EM for Monte-Carlo speed. At import time it pre-computes:
- `Phi_stack` `(n_z, p, d)` — all $\Phi(z_q)$ at once
- `PhiTPhi_stack` `(n_z, d, d)` — all $\Phi(z_q)^T \Phi(z_q)$ via `einsum`

| Function | Role |
|---|---|
| `e_step(u, observations, sigma)` | Vectorised E-step; residuals computed as `(n, n_z, p)` broadcast |
| `m_step(weights, observations)` | Vectorised M-step via `einsum`; no Python loops |
| `em_run(observations, sigma, ...)` | Multi-start EM; returns best `u` by log-likelihood |
| `generate_data(n, sigma)` | Generates `n` observations $Y_i = \Phi(Z_i)u^\star + \sigma E_i$ |
| `estimate_eps(n, sigma, M)` | Monte-Carlo estimate of $\varepsilon_n(\sigma)$: $\sqrt{\frac{1}{M}\sum \text{err}^2}$ |
| `ref_error(n, sigma)` | Section 2 reference: $\frac{\sigma}{\sqrt{n}}\sqrt{\mathrm{Tr}((A^TA)^{-1})}$ |

Generated figures (all saved to `../figures/`, i.e. `figures/` at repo root):
- `signaux_1D.pdf`, `simulation_parametres.pdf` — from `graddesc_cas1.ipynb`
- `erreur_vs_n.pdf`, `erreur_vs_sigma.pdf` — from `analyse_erreur.ipynb`

## LaTeX conventions in `doc.tex`

- French (`babel` + `inputenc utf8` + `fontenc T1`). Math commentary, comments, and TODOs are all in French.
- Custom shortcuts defined at the preamble (lines 48–61): `\R`, `\N`, `\E`, `\Prob`, `\Zcal`, `\Fcal`, `\Ncal`, `\norm{}`, `\abs{}`, `\inner{}{}`, `\Id`, `\Tr`, `\argmin`, `\argmax`. Use these rather than reinventing.
- Theorems share a single counter (`theorem` numbered by section, with `proposition`/`lemma`/`corollary`/`definition`/`example`/`remark` sharing `[theorem]`).
- Algorithms use `algorithm2e` with `[ruled, vlined, french]`.
- Section bodies frequently end with a French `%`-prefixed outline of TODO content — that outline is the spec for what still needs to be written. Don't delete those scaffolds unless the prose that replaces them covers the same points.

## Things easy to get wrong

- Don't commit LaTeX build artifacts (`*.aux`, `*.fdb_latexmk`, `*.fls`, `*.log`, `*.out`, `*.synctex.gz`, `*.toc`) — they're currently sitting in the working tree because there's no `.gitignore` at repo root, but they should not be tracked.
- `calendar.css` is unrelated to the project (a local FullCalendar tweak that ended up here). Leave it alone.
