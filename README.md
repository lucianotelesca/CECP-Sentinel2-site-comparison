# CECP-Sentinel2-site-comparison
This repository contains all data and code necessary to fully reproduce the analysis presented in the manuscript. 
---
Repository Structure
```
/raw_data/
    data_input.mat                  → raw spectral time series (input to pre-processing)

/processed_data/
    data_preprocessed.mat           → pre-processed data (input to CECP analysis)

/code/
    modifyStructWithGap7_sx3.m      → Step 1: outlier removal & temporal alignment
                                       (includes divide_myTableFlexible_sx as local function)
    CECP_compare_sites_allbands.m   → Step 2: complexity-entropy (H, C) analysis
    complexity_entropy_universal.m  → computes H and C for a time series
    maximum_complexity_entropy.m    → upper theoretical limit curve of the CECP
    minimum_complexity_entropy.m    → lower theoretical limit curve of the CECP

README.md
```
---
Sites
The analysis compares two sites:

CP — infected site

FOL — healthy site

---
Input Data — `data_input.mat`

The file `data_input.mat` (MATLAB v5 format) contains the raw spectral time series used as input to the pre-processing pipeline. It stores two top-level structures:

`input_CP` — data for the infected site

`input_FOL` — data for the healthy site

Each structure contains the following fields:

Field	Description	Value used in this study

`originalStruct`	MATLAB struct of tables, one per spectral band	6 bands (see below)

`gap`	Target time step for temporal alignment (days)	10

`flexibleAlignment`	Enable flexible date matching	true

`n`	Day tolerance for flexible alignment (days)	5

`sogliaiqr_up`	Upper IQR multiplier for outlier removal	2

`sogliaiqr_do`	Lower IQR multiplier for outlier removal	2

The `originalStruct` field contains one MATLAB `table` per spectral band. Each table has a `datetime` column (acquisition dates) followed by one column per pixel. The six spectral bands included are:

Band	Description

`blue`	Blue spectral band

`green`	Green spectral band

`red`	Red spectral band

`nir`	Near-infrared band

`swir1`	Short-wave infrared band 1

`swir2`	Short-wave infrared band 2


---
Code

Step 1 — Pre-processing: `modifyStructWithGap7_sx3.m`

Performs outlier removal and temporal alignment on the raw spectral time series. 

Pipeline:

Outlier removal using a median ± k·IQR threshold (applied column-wise, i.e. per pixel)

Temporal alignment to a regular date grid using left-side flexible matching: each observation is assigned to the nearest target date on the left within a tolerance of `n` days

Summary statistics computed before and after processing (valid observations, NaN counts, min/max, outliers removed)

[struct_CP,  stats_CP]  = modifyStructWithGap7_sx3(input_CP);

[struct_FOL, stats_FOL] = modifyStructWithGap7_sx3(input_FOL);

```
---
Step 2 — Complexity-Entropy Analysis: `CECP_compare_sites_allbands.m`

Computes permutation entropy (H) and statistical complexity (C) for every pixel and every spectral band, then compares CP vs FOL using the Complexity-Entropy Causality Plane (CECP), ROC analysis, and Cohen's d effect size.


results = CECP_compare_sites_allbands(struct_CP, struct_FOL);

results = CECP_compare_sites_allbands(struct_CP, struct_FOL, dx);

results = CECP_compare_sites_allbands(struct_CP, struct_FOL, dx, options);

```

---

---
Citation
If you use this code or data, please cite the associated manuscript:
> [Authors], [Year]. [Title]. [Journal]. DOI: [DOI]
---
License
[CC BY 4.0]
