# Weighted Spatial Averages with Missing Observations: Convergence and Accuracy of Delta-Method Variance Approximations

[![DOI](https://zenodo.org/badge/1346177850.svg)](https://doi.org/10.5281/zenodo.22097270)

Code and data accompanying the paper:

> Seshadri, A. K. and Pal Majumder, A. (2026). Weighted spatial averages with missing observations: Convergence and accuracy of delta-method variance approximations. Submitted to *Spatial Statistics*.

This paper follows:

> Seshadri, A. K. (2018). Statistics of spatial averages and optimal averaging in the presence of missing data. *Spatial Statistics*, 25, 1–21. [doi:10.1016/j.spasta.2018.04.002](https://doi.org/10.1016/j.spasta.2018.04.002) — code at [akseshadri/spatial-averages-missing-data](https://github.com/akseshadri/spatial-averages-missing-data)

## Overview

When observations are missing at random, a weighted spatial average is a ratio of two random quantities, the weighted sum of the reported values and the total weight of the reporting sites. Its variance is usually approximated by truncating a Taylor series of the ratio (the delta method). The variance estimator of Seshadri (2018) is the first-order member of a two-index sequence of such approximations. This paper examines when the series converges, how accurate a given truncation is, and whether higher-order terms help. The code computes the approximations, the exact variance for uniform weights (as a single sum over the number of reporting sites), and the exact variance for small networks with arbitrary weights (by enumerating all reporting patterns). The numerical examples use a homogeneous benchmark, small synthetic networks, and the daily gridded rainfall over India from the India Meteorological Department (IMD).

## Repository Structure

```
├── README.md
├── LICENSE
├── .gitignore
├── CITATION.cff
├── data/
│   └── indiadat.mat          # 1° × 1° gridded daily rainfall (1901–2011, 357 locations)
├── figures/                  # EPS output of the figure scripts
└── src/
    ├── PFig1.m               # Figure 1: error scaling and denominator concentration
    ├── PFig2.m               # Figure 2: validation with the IMD rainfall field
    ├── PFig3.m               # Figure 3: common shift of the field
    ├── PFig4.m               # Figure 4: stress tests
    ├── PFigS1.m              # Figure S1: homogeneous two-parameter scaling
    ├── PFigS2.m              # Figure S2: exact check of the mixed-moment formulas
    ├── PFigS3.m              # Figure S3: term-by-term decomposition for the IMD field
    ├── PTables.m             # Values in Table 4 and Table S3
    ├── loadimd.m             # Helper: site means and covariance of the IMD field
    ├── varexact.m            # Helper: exact variance for uniform weights
    ├── truncvar.m            # Helper: variance approximation mu^(n1,n2) for uniform weights
    ├── enumvar.m             # Helper: exact variance and approximations by enumeration (small N)
    └── mixedmom.m            # Helper: closed forms of the mixed moments mu12, mu21, kappa
```

## Requirements

- **MATLAB** (R2016b or later). No toolboxes are required.
- The scripts also run in **GNU Octave**.

## Usage

1. Set the working directory to the `src/` folder:
   ```matlab
   cd('/path/to/repository/src')
   ```

2. Run any script directly. For example, to reproduce Figure 2:
   ```matlab
   PFig2
   ```

   Each script loads the required data from `../data/`, calls the helper functions, and writes an EPS file to `../figures/`. `PTables` prints its values to the command window.

### Script Descriptions

| Script | Paper item | Description |
|--------|-----------|-------------|
| `PFig1.m` | Fig. 1 | Homogeneous benchmark: errors of the (0,0) and (0,1) approximations against ρ and N, and concentration of the reporting denominator |
| `PFig2.m` | Fig. 2 | IMD field, uniform weights: exact variance and four approximations, with errors relative to the total variance and to the contribution of missing data |
| `PFig3.m` | Fig. 3 | IMD field with a common shift added to all sites: only the (0,1) approximation is unaffected, like the exact variance |
| `PFig4.m` | Fig. 4 | Stress tests: strong spatial correlation, small networks, and concentrated weights |
| `PFigS1.m` | Fig. S1 | Homogeneous benchmark with both availability and network size varied |
| `PFigS2.m` | Fig. S2 | Closed-form mixed moments compared with exact enumeration for a small correlated network |
| `PFigS3.m` | Fig. S3 | Term-by-term decomposition of the low-order approximations for the IMD field |
| `PTables.m` | Table 4, Table S3 | Convergence at fixed network size, and errors along the shift-invariant truncations |

All calculations are exact averages over the reporting process for the given field moments; no Monte Carlo simulation is used.

## Data Source

- **1° × 1° gridded rainfall** (`indiadat.mat`): Rajeevan, M., Bhate, J., Kale, J. D., & Lal, B. (2006). High resolution daily gridded rainfall data for the Indian region: Analysis of break and active monsoon spells. *Current Science*, 91, 296–306.

The rainfall data are provided by the India Meteorological Department (IMD) and
are redistributed here for reproducibility; please observe IMD's terms of use.

## Citation

If you use this code, please cite the paper, this archive, and Seshadri (2018):

```bibtex
@article{SeshadriPalMajumder2026,
  title   = {Weighted spatial averages with missing observations: Convergence and accuracy of delta-method variance approximations},
  author  = {Seshadri, Ashwin K. and Pal Majumder, Abhishek},
  year    = {2026}
}

@software{SeshadriPalMajumder2026code,
  title     = {Weighted spatial averages with missing observations: Convergence and accuracy of delta-method variance approximations (code)},
  author    = {Seshadri, Ashwin K. and Pal Majumder, Abhishek},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22097270}
}

@article{Seshadri2018,
  title   = {Statistics of spatial averages and optimal averaging in the presence of missing data},
  author  = {Seshadri, Ashwin K.},
  journal = {Spatial Statistics},
  volume  = {25},
  pages   = {1--21},
  year    = {2018},
  doi     = {10.1016/j.spasta.2018.04.002}
}
```

## License

This project is licensed under the [MIT License](LICENSE). If you use this code,
please reference the associated paper. The IMD rainfall data are subject to
IMD's own terms of use.
