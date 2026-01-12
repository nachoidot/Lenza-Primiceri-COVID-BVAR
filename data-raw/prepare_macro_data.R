# prepare_macro_data.R
# Data preparation script for COVID-19 BVAR package
# Loads and transforms dataMLprojectMay2021.xlsx
# Exact replication of Baseline_May2021.m lines 16-36

library(readxl)

# Load raw data from Excel
raw_data <- read_excel(
  "../inst/extdata/dataMLprojectMay2021.xlsx",
  sheet = "Monthly"
)

# Extract date column (first column contains Excel serial dates)
# Excel serial dates start from 12/31/1899
excel_dates <- raw_data[[1]]
TimeMACRO <- as.Date(excel_dates, origin = "1899-12-30")

# Extract data columns (all except first)
DataMACRO <- as.matrix(raw_data[, 2:ncol(raw_data)])

# Data transformations (matching MATLAB lines 22-25)
# DataMACRO(:,15)=exp(DataMACRO(:,15)/100);  % unemployment
# DataMACRO(:,16)=exp(DataMACRO(:,16)/100);  % GZ spread
# DataMACRO(:,9)=DataMACRO(:,9)./DataMACRO(:,12);  % real pce
# DataMACRO(:,14)=DataMACRO(:,14)./DataMACRO(:,6);  % real pce services

DataMACRO[, 15] <- exp(DataMACRO[, 15] / 100)  # unemployment
DataMACRO[, 16] <- exp(DataMACRO[, 16] / 100)  # GZ spread
DataMACRO[, 9] <- DataMACRO[, 9] / DataMACRO[, 12]  # real pce
DataMACRO[, 14] <- DataMACRO[, 14] / DataMACRO[, 6]  # real pce services

# Variable names (matching MATLAB line 28-29)
# [ 1      2      3       4        5       6     7   8     9    10      11     12      13      14   15     16    ]
# [cpi  ppcedg  ppceg  ppcendg  corepce  ppces  ip  empl  pce  pcedg  pcendg  ppce  coreppce  pces unem  GZspread]

variable_names <- c(
  "cpi",               # 1: CPI
  "ppcedg",            # 2: PCE: durable goods (price)
  "ppceg",             # 3: PCE: goods (price)
  "ppcendg",           # 4: PCE: non-durable goods (price)
  "corepce",           # 5: core PCE (price)
  "ppces",             # 6: PCE: services (price)
  "ip",                # 7: industrial production
  "empl",              # 8: employment
  "pce",               # 9: PCE (real)
  "pcedg",             # 10: PCE: durable goods (nominal)
  "pcendg",            # 11: PCE: non-durable goods (nominal)
  "ppce",              # 12: PCE (price deflator)
  "coreppce",          # 13: core PCE (price)
  "pces",              # 14: PCE: services (real)
  "unem",              # 15: unemployment rate
  "GZspread"           # 16: Gilchrist-Zakrajsek spread
)

variable_descriptions <- c(
  "Consumer Price Index",
  "PCE: Durable Goods (Price Deflator)",
  "PCE: Goods (Price Deflator)",
  "PCE: Non-Durable Goods (Price Deflator)",
  "Core PCE (Price Deflator)",
  "PCE: Services (Price Deflator)",
  "Industrial Production Index",
  "Total Nonfarm Employment",
  "Personal Consumption Expenditures (Real)",
  "PCE: Durable Goods (Nominal)",
  "PCE: Non-Durable Goods (Nominal)",
  "PCE (Price Deflator)",
  "Core PCE (Price Deflator, Alternative)",
  "PCE: Services (Real)",
  "Unemployment Rate (%)",
  "Gilchrist-Zakrajsek Credit Spread (bps)"
)

# Key dates (matching MATLAB lines 40-47)
T0 <- which(format(TimeMACRO, "%Y") == "1988" & format(TimeMACRO, "%m") == "12")  # Dec 1988
T1estim <- which(format(TimeMACRO, "%Y") == "2021" & format(TimeMACRO, "%m") == "05")  # May 2021
Tfeb2020 <- which(format(TimeMACRO, "%Y") == "2020" & format(TimeMACRO, "%m") == "02")  # Feb 2020
Tcovid <- Tfeb2020 - T0 + 2  # Position of March 2020 (Tcovid = 376 in Baseline model)

# Baseline model variable indices (MATLAB line 32)
# indmacro=[15 8 9 14 12 6 13]
# = [unemployment, employment, PCE, PCE services, PCE deflator, PCE services deflator, core PCE]
baseline_indices <- c(15, 8, 9, 14, 12, 6, 13)
baseline_names <- c(
  "unemployment",
  "employment",
  "PCE",
  "PCE: services",
  "PCE (price)",
  "PCE: services (price)",
  "core PCE (price)"
)

# Create dataset object
macro_data_may2021 <- list(
  data = DataMACRO,
  time = TimeMACRO,
  variable_names = variable_names,
  variable_descriptions = variable_descriptions,
  baseline_indices = baseline_indices,
  baseline_names = baseline_names,
  key_dates = list(
    T0 = T0,              # Start of estimation sample (Dec 1988)
    T1estim = T1estim,    # End of estimation sample (May 2021)
    Tfeb2020 = Tfeb2020,  # Position of Feb 2020
    Tcovid = Tcovid       # Position of March 2020 (for volatility break)
  ),
  notes = paste(
    "Macroeconomic data for COVID-19 BVAR estimation.",
    "Monthly U.S. data from December 1988 to May 2021.",
    "Data transformations applied: unemployment and GZ spread exponentiated,",
    "real PCE computed as nominal/deflator.",
    "Source: Federal Reserve Economic Data (FRED).",
    "Used in Lenza & Primiceri (2022) JAE."
  )
)

# Save as .rda file
usethis::use_data(macro_data_may2021, overwrite = TRUE)

# Print summary
cat("\nData preparation complete!\n")
cat("---------------------------\n")
cat("Total observations:", nrow(DataMACRO), "\n")
cat("Number of variables:", ncol(DataMACRO), "\n")
cat("Time range:", format(TimeMACRO[1], "%b %Y"), "to", format(TimeMACRO[length(TimeMACRO)], "%b %Y"), "\n")
cat("Baseline model uses", length(baseline_indices), "variables:\n")
for (i in seq_along(baseline_names)) {
  cat("  -", baseline_names[i], "\n")
}
cat("\nKey dates:\n")
cat("  - Estimation start (T0):", format(TimeMACRO[T0], "%b %Y"), "(obs", T0, ")\n")
cat("  - Estimation end (T1):", format(TimeMACRO[T1estim], "%b %Y"), "(obs", T1estim, ")\n")
cat("  - Feb 2020:", format(TimeMACRO[Tfeb2020], "%b %Y"), "(obs", Tfeb2020, ")\n")
cat("  - March 2020 (COVID break):", "obs", Tcovid, "relative to T0\n")
cat("\nSaved as: data/macro_data_may2021.rda\n")
