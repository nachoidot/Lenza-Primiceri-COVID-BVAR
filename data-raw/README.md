# Data Preparation Instructions

This directory contains the script to prepare the macroeconomic dataset for the covidBVAR package.

## Files

- `prepare_macro_data.R`: R script that loads the Excel file and creates the `macro_data_may2021.rda` file
- Original data: `../inst/extdata/dataMLprojectMay2021.xlsx`

## How to Run

From the package root directory:

```r
source("data-raw/prepare_macro_data.R")
```

This will:
1. Load data from the Excel file
2. Apply transformations (matching Baseline_May2021.m lines 22-25):
   - Unemployment: exp(x/100)
   - GZ spread: exp(x/100)
   - Real PCE: nominal / deflator
   - Real PCE services: nominal / deflator
3. Create the `macro_data_may2021` list object
4. Save to `data/macro_data_may2021.rda`

## Data Structure

The resulting dataset includes:
- **data**: 390 x 16 matrix of monthly observations (Dec 1988 - May 2021)
- **time**: Date vector
- **variable_names**: Short names for 16 variables
- **variable_descriptions**: Full descriptions
- **baseline_indices**: Indices of 7 variables for baseline model
- **baseline_names**: Names of baseline variables
- **key_dates**: Important observation indices (T0, T1estim, Tfeb2020, Tcovid)
- **notes**: Dataset description

## Variable List

1. CPI - Consumer Price Index
2. ppcedg - PCE: Durable Goods (Price)
3. ppceg - PCE: Goods (Price)
4. ppcendg - PCE: Non-Durable Goods (Price)
5. corepce - Core PCE (Price)
6. ppces - PCE: Services (Price)
7. ip - Industrial Production
8. empl - Employment
9. pce - PCE (Real)
10. pcedg - PCE: Durable Goods (Nominal)
11. pcendg - PCE: Non-Durable Goods (Nominal)
12. ppce - PCE Deflator
13. coreppce - Core PCE (Alternative)
14. pces - PCE: Services (Real)
15. unem - Unemployment Rate
16. GZspread - Gilchrist-Zakrajsek Spread

## Baseline Model Variables (indices: 15, 8, 9, 14, 12, 6, 13)

1. Unemployment
2. Employment
3. PCE
4. PCE: Services
5. PCE (price)
6. PCE: Services (price)
7. Core PCE (price)

## Notes

- Data transformations match the original MATLAB code exactly
- The Tcovid parameter (376) represents March 2020 relative to estimation start (Dec 1988)
- All variables will be log-transformed and multiplied by 100 before VAR estimation
