# A Convergent Moment Hierarchy for the Variance of Weighted Spatial-Mean Estimators under Missing-at-Random Sampling

Code and data accompanying the paper:
> Seshadri, A. K. and Pal Majumder, A. (2026). A convergent moment hierarchy for the variance of weighted spatial-mean estimators under missing-at-random sampling.

## Overview

The paper studies the weighted spatial-mean ratio estimator `rhat = R/S`, with
`R = sum_i beta_i s_i r_i` and `S = sum_i beta_i s_i`, when sites report at
random (`s_i ~ Bernoulli(alpha)`, missing-at-random). Because the denominator
is random, the estimator has no tractable closed-form variance; the paper
expands `R/S` about `(E R, E S)` with the delta method and builds a hierarchy of
closed-form variance approximations `mu2^{(n1,n2)}`, indexed by a coupling order
`n1` and a denominator order `n2`. A single dimensionless quantity, the relative
denominator fluctuation `rho = sqrt(Var(S))/E[S]`, controls the truncation
accuracy. This repository contains the closed forms, an exact
(enumeration-based) variance check, the Monte-Carlo validation, and the scripts
that reproduce every figure, applied to India Meteorological Department (IMD)
gridded daily rainfall.

## Repository Structure

```
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── data/
│   ├── indiadat.mat          # 1° × 1° gridded daily rainfall (1901–2011, 357 locations)
│   ├── indialatlon.mat       # Latitude/longitude coordinates for the 1° grid
│   └── latlonmat.txt         # Same coordinates as plain text (for reference)
└── src/
    ├── run_all_figures.m         # Master driver: runs the three section scripts
    ├── PFigs_Section2.m          # Figures 1–2: rho governs convergence and accuracy
    ├── PFigs_Section3.m          # Figures 3–6: closed forms and their limits
    ├── PFigs_SI.m                # Figures S1–S2: source decomposition (SI)
    ├── truncations_scenarioI.m   # Uniform-weight (Scenario I) closed forms
    ├── truncations_general.m     # Arbitrary-weight closed forms and rho
    ├── ratio_var_exact_uniform.m # Exact Var(rhat) for uniform weights (no Monte Carlo)
    ├── mc_truevar.m              # Monte-Carlo moments of rhat under MAR
    ├── apply_fig_style.m         # Consistent figure style
    └── panel_label.m             # Bold (a)/(b)/(c) panel letters
```

The figure scripts do not write files; each opens its results as MATLAB figure
windows.

## Requirements

- **MATLAB** R2016b or later (uses `yyaxis` and `histogram`).
- **Statistics and Machine Learning Toolbox** (for `binornd`, `binopdf`).
- **Parallel Computing Toolbox** is optional: `mc_truevar` uses `parfor`, which
  runs serially if no pool is available.

## Usage

1. Open MATLAB and set the working directory to the `src/` folder:

   ```matlab
   cd('/path/to/this-repository/src')
   ```

2. Reproduce all figures at once:

   ```matlab
   run_all_figures
   ```

   or run any section script on its own, for example:

   ```matlab
   PFigs_Section3
   ```

   Each script loads the rainfall data from `../data/` and calls the shared
   functions. The Monte-Carlo ensemble size `Kensemb` defaults to `1e5`; reduce
   it for quicker previews.

### Figure Descriptions

| Script            | Paper Figure | Description                                                                                       |
| ----------------- | ------------ | ------------------------------------------------------------------------------------------------- |
| `PFigs_Section2.m`| Figs. 1–2    | `rho`-collapse, sharp `O(rho^2)`/`O(rho^4)` rates, and `1/N` convergence; good-event/rare-event split |
| `PFigs_Section3.m`| Figs. 3–6    | Monte Carlo vs the four truncations; mixed moments vs theory; source decomposition; stress regimes |
| `PFigs_SI.m`      | Figs. S1–S2  | Term-by-term source decomposition: raw field and reduced-variance case                            |

## Data Sources

- **1° × 1° gridded rainfall** (`indiadat.mat`, `indialatlon.mat`): Rajeevan, M., Bhate, J., Kale, J. D., & Lal, B. (2006). High resolution daily gridded rainfall data for the Indian region: Analysis of break and active monsoon spells. *Current Science*, 91, 296–306.

The rainfall data are provided by the India Meteorological Department (IMD) and
are redistributed here for reproducibility; please observe IMD's terms of use.

## Citation

If you use this code, please cite both the paper and this archive:

```bibtex
@article{SeshadriPalMajumder2026,
  title   = {A convergent moment hierarchy for the variance of weighted spatial-mean estimators under missing-at-random sampling},
  author  = {Seshadri, Ashwin K. and Pal Majumder, Abhishek},
  year    = {2026}
}

@software{SeshadriPalMajumder2026code,
  title     = {A convergent moment hierarchy for the variance of weighted spatial-mean estimators under missing-at-random sampling (code)},
  author    = {Seshadri, Ashwin K. and Pal Majumder, Abhishek},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22097271}
}
```

## License

This project is licensed under the [MIT License](LICENSE). If you use this code,
please reference the associated paper. The IMD rainfall data are subject to
IMD's own terms of use.
