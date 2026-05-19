# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

M1 research project (Université Paris Cité, MMAS 2025–2026) by Arthur Conche and Abdel Niang, supervised by Camille Pouchol. Subject: *Estimation in a linear inverse problem with hidden variables — likelihood, EM algorithm and error analysis* — a toy model of cryo-EM reconstruction.

The repo has two deliverables that evolve in parallel:

- **`doc.tex` → `doc.pdf`** — the written report (French, ~550 lines of LaTeX).
- **`code python/*.ipynb`** — numerical experiments that produce the figures and validate the algorithms discussed in the report. As of the latest snapshot, `générarion1D.ipynb` is empty (0 bytes) — notebooks are still to be written.

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
- **Cas 2 — 2D rotation + Radon projection.** Analog of cryo-EM proper; treated as an extension (Section 6), not the main object.

Section flow in `doc.tex` mirrors the methodology: (§1) setup → (§2) reference case with $A_z = A$ constant, closed-form MLE → (§3) marginal likelihood over $Z$, identifiability obstruction → (§4) gradient descent **and** EM as two solvers, EM's M-step is a quadratic linear system → (§5) numerical study of $\varepsilon_n(\sigma) = \mathbb{E}\|\hat\theta_n - \theta\|^2$ — **note the metric is translation-invariant** (`\inf_{a \in [0,2\pi]} \|\tau_a \hat\theta - \vartheta\|`), because the standard norm is unidentifiable in Cas 1 → (§6) extensions.

When implementing algorithms in notebooks, the Fourier diagonalization of Cas 1 is the key efficiency lever: do not rotate in the spatial domain.

## LaTeX conventions in `doc.tex`

- French (`babel` + `inputenc utf8` + `fontenc T1`). Math commentary, comments, and TODOs are all in French.
- Custom shortcuts defined at the preamble (lines ~46–67): `\R`, `\N`, `\E`, `\Prob`, `\Zcal`, `\Fcal`, `\Ncal`, `\norm{}`, `\abs{}`, `\inner{}{}`, `\Id`, `\Tr`, `\argmin`, `\argmax`. Use these rather than reinventing.
- Theorems share a single counter (`theorem` numbered by section, with `proposition`/`lemma`/`corollary`/`definition`/`example`/`remark` sharing `[theorem]`).
- Algorithms use `algorithm2e` with `[ruled, vlined, french]`.
- Section bodies frequently end with a French `%`-prefixed outline of TODO content — that outline is the spec for what still needs to be written. Don't delete those scaffolds unless the prose that replaces them covers the same points.

## Things easy to get wrong

- Don't commit LaTeX build artifacts (`*.aux`, `*.fdb_latexmk`, `*.fls`, `*.log`, `*.out`, `*.synctex.gz`, `*.toc`) — they're currently sitting in the working tree because there's no `.gitignore` at repo root, but they should not be tracked.
- `calendar.css` is unrelated to the project (a local FullCalendar tweak that ended up here). Leave it alone.
- The repo is **not** a git repository at the moment (`Is a git repository: false`). Don't assume `git log` works.
