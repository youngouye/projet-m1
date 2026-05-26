# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

M1 research project (Université Paris Cité, MMAS 2025–2026) by Arthur Conche and Abdel Niang, supervised by Camille Pouchol. Subject: *Estimation in a linear inverse problem with hidden variables — likelihood, EM algorithm and error analysis* — a toy model of cryo-EM reconstruction.

The repo has two deliverables that evolve in parallel:

- **`doc.tex` → `doc.pdf`** — the written report (French, ~1127 lines of LaTeX).
- **`code python/1DD_optimise.ipynb`** — numerical experiments for Cas 1 (1D translations); produces figures and validates the algorithms discussed in the report.

`Projet_M1.pdf` is the assignment brief from the supervisor. `Exemples_Rapports_M1_MMAS/` contains six PDFs of prior years' reports as style/scope references — they are **not** part of this project's source.

## Build & run

LaTeX (full rebuild with bibliography, run from repo root):

```bash
latexmk -pdf doc.tex          # build doc.pdf
latexmk -c                    # clean aux files (.aux .fls .fdb_latexmk .log .out .toc .synctex.gz)
latexmk -C                    # also remove doc.pdf
```

The `.tex` references `\bibliography{references}` but **`references.bib` does not yet exist** — first `latexmk` run will warn about undefined citations until that file is created.

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

## Notebook architecture (`code python/1DD_optimise.ipynb`)

Global parameters at the top of the notebook: `n` (observations), `p` (grid points), `N` (Fourier frequencies), `d = 2N+1` (signal dimension), `sigma` (noise), `n_z` (quadrature points for the $z$-integral).

Key functions:

| Function | Role |
|---|---|
| `build_phi(z)` | Returns the $(p \times d)$ evaluation matrix $\Phi(z)$; `build_phi(z) @ u` gives $A_z\theta$ on the grid |
| `build_phiu(z, u)` | Thin wrapper: `build_phi(z) @ u`; returns the signal vector on the grid as `(p,)` |
| `qi(z, y, u)` | Computes $\frac{1}{2\sigma^2}\|y - \Phi(z)u\|^2$, clipped to 745 to avoid float64 overflow |
| `log_likelihood(u)` | Marginal log-likelihood via Riemann quadrature + `scipy.logsumexp` for numerical stability |
| `compute_gradient(u)` | Gradient $\nabla_\theta \ell_n$ using posterior weights $w_{i,z}$ (eq. 3.5 in the report) |
| `gradient_ascent_armijo(u0, ...)` | Gradient ascent with Armijo backtracking line search; returns `(u_hat, hist_ll, hist_grad)` |
| `translation_error(u_hat, u_true, xgrid)` | Translation-invariant relative error via FFT cross-correlation (aligns signals before comparing) |

**Performance note:** `compute_gradient` contains a Python loop over `n × n_z` iterations, making it the bottleneck for large `n` or `n_z`. Vectorising with a stacked Phi tensor `(n_z, p, d)` is the natural next step.

Generated figures are saved to `../figures/` (i.e. `figures/` at the repo root): `signaux_1D.pdf` and `simulation_parametres.pdf`.

## LaTeX conventions in `doc.tex`

- French (`babel` + `inputenc utf8` + `fontenc T1`). Math commentary, comments, and TODOs are all in French.
- Custom shortcuts defined at the preamble (lines 48–61): `\R`, `\N`, `\E`, `\Prob`, `\Zcal`, `\Fcal`, `\Ncal`, `\norm{}`, `\abs{}`, `\inner{}{}`, `\Id`, `\Tr`, `\argmin`, `\argmax`. Use these rather than reinventing.
- Theorems share a single counter (`theorem` numbered by section, with `proposition`/`lemma`/`corollary`/`definition`/`example`/`remark` sharing `[theorem]`).
- Algorithms use `algorithm2e` with `[ruled, vlined, french]`.
- Section bodies frequently end with a French `%`-prefixed outline of TODO content — that outline is the spec for what still needs to be written. Don't delete those scaffolds unless the prose that replaces them covers the same points.

## Things easy to get wrong

- Don't commit LaTeX build artifacts (`*.aux`, `*.fdb_latexmk`, `*.fls`, `*.log`, `*.out`, `*.synctex.gz`, `*.toc`) — they're currently sitting in the working tree because there's no `.gitignore` at repo root, but they should not be tracked.
- `calendar.css` is unrelated to the project (a local FullCalendar tweak that ended up here). Leave it alone.
