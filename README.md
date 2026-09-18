# brazil_downscalr_usocob.R

An R script that downscales Brazil's national/state-level FABLE Calculator
land-use targets (2000–2050) to a fine-grained spatial grid, using
**MapBiomas UsoCobertura** as the base-year land stock (in place of the
older HILDA+ 2015 stock) and the `downscalr` package's Bayesian
multinomial-logit (MNL) allocation model. It produces gridded land-use
transitions, diagnostic plots, and area-conservation checks.

---

## What it does, at a high level

1. Reads a 2020 MapBiomas land-cover stock table and collapses ~28 detailed
   land-cover columns into 5 FABLE land-use classes (Forest, OtherLand,
   Cropland, Pasture, Urban).
2. Loads national/state FABLE Calculator targets (future period-by-period
   transition targets) and historical HILDA+ transitions (used to train the
   allocation model).
3. Builds a covariate matrix per grid cell (land shares, altitude, slope,
   travel time, livestock, population, crop yield) and estimates an MNL
   model of land-use transition propensity for each source class.
4. Builds cell-level restrictions (protected areas + a manually justified
   geographic exclusion for remote interior Amazon states) that block
   certain land conversions in certain cells.
5. Runs `downscalr::downscale()` to allocate the national/state targets down
   to grid cells for each 5-year period from 2000–2050.
6. Saves results, plots a Cropland-2050 choropleth (change land use and year accordingly, in arguments), runs diagnostics on
   specific "offender" cells, verifies area conservation across periods,
   and produces summary tables and covariate maps.

---

## Requirements

**R packages:** `dplyr`, `tidyr`, `tibble`, `readxl`, `downscalr`, `nloptr`
(used internally by the patch), `sf`, `terra`, `RColorBrewer`, `ggplot2`,
`openxlsx`/`writexl` (for `write.xlsx`).

**Working directory:** the script calls
`setwd("C:/Users/User/Desktop/FABLE/downscalr")` — this is a hardcoded,
machine-specific path and must be changed before running elsewhere. Some
shapefile paths later in the script (`br_states.shp`, `br_biomes.shp`) are
also hardcoded to an absolute Windows path and will need to be updated.

**Input files (expected in the working directory):**

| File | Purpose |
|---|---|
| `UsoCobertura__fable_v2.xlsx` (sheet `usoecob2020`/`usoecob2000`) | MapBiomas 2020 land-cover stock, wide format, one row per grid cell |
| `hildaluc_br_2015_2019.xlsx` | HILDA+ 2015→2019 transition flows, used to train the MNL model |
| `brazil_fable_UP50_CT.xlsx` | FABLE Calculator land-use targets by period |
| `altitude_br.xlsx`, `slope_br.xlsx`, `travel_time_br.xlsx` | Per-cell physical covariates |
| `livestock_br.xlsx`, `pop_2020.xlsx`, `crop_yield_br.xlsx` | Per-cell socioeconomic covariates (`pop_2020.xlsx` also supplies cell geometries via a `.geo` column) |
| `Areas_Protegidas_FABLE.xls` (sheet `FABLE`) | Protected-area (ha) overlap per cell |
| `br_states.shp`, `br_biomes.shp` | State and biome boundary shapefiles used for the geographic exclusion rule and for map overlays |


## Section-by-section description

### Header / setup (lines 1–29)
File banner describing purpose, inputs, and the rule that the LU-class
mapping should only be edited in Section A. Loads packages and sets the
working directory (hardcoded path — must be edited per machine).

### Runtime patch: `solve_biascorr.mnl` (lines 30–296)
A monkey-patch for a bug in `downscalr::solve_biascorr.mnl()`. When exactly
one (`lu.from`, `lu.to`) target triggers the bias-correction fallback,
`restr.mat[, error_ind]` silently drops to a 1-D vector, and the subsequent
`error_restrictions[, ccc]` indexing then fails with a dimensions error.
The patch redefines the whole function with `drop = FALSE` added at that
point and installs it into the `downscalr` namespace via
`assignInNamespace()`, so it takes effect for the rest of the session
without modifying or reinstalling the package. Must run after
`library(downscalr)`.

### LU class list (line 302)
Defines `LU_CLASSES <- c("Forest", "OtherLand", "Cropland", "Pasture", "Urban")`,
the five target classes used throughout the script.

### A. Column → LU class mapping (lines 304–351)
`col_to_lu`: a named vector mapping ~28 raw MapBiomas land-cover column
names (e.g. `Savanna`, `Soybean`, `Mining`) to one of the five FABLE
classes. This is called out as the single place to edit if the
classification scheme needs to change. Columns not listed (e.g. `NoData`)
are dropped.

### B. Load `usoecob2020` and build LU stock table (lines 353–390)
Reads the MapBiomas sheet, keeps `id_c` and all `(ha)` columns, strips the
`(ha)` suffix from column names, then sums the mapped columns per cell into
one column per FABLE class (`lu_Forest`, `lu_Pasture`, etc.), producing
`lu_stock`. Prints total hectares per class as a sanity check.

### C. `brazil_FABLE` (lines 392–417)
Reads the FABLE Calculator target workbook, renames columns to the
`downscalr` schema (`lu.from`, `times`), pivots the wide `To<Class>` target
columns into long format, strips the `To` prefix, converts values from
1000 ha (kha, the FABLE Calculator's native unit) to raw ha to match the
rest of the script, and drops the `NewForest` pseudo-class and any period
before 2000.

### D. `brazil_luc` — HILDA+ transitions (lines 419–445)
Reads historical 2015→2019 land-use transition counts from HILDA+, decodes
two-digit `from`/`to` land-cover codes into FABLE classes via a lookup
table, and aggregates to (cell, lu.from, lu.to, year) totals. This is the
observed-transition data used to train the MNL model in Section G; it is
explicitly noted as unchanged from the original HILDA-based version of the
script.

### E. `brazil_xmat` — covariates (lines 447–537)
Builds the covariate matrix used in the MNL model:
- Reads altitude, slope, travel time, livestock, population, and crop-yield
  tables per cell, and joins them to the land-use stock table.
- Imputes missing values with the column median (`impute_median()`).
- Log-transforms livestock, population, and crop totals (`log1p`).
- Drops zero-variance columns and standardizes (z-scores) all covariates
  globally (`xmat_std`), producing both a wide (`xmat_full`) and long
  (`X_long`) version.
- Saves diagnostic RDS snapshots (`xmat_wide_diag.rds`, `X_long_diag.rds`,
  `lu_stock_diag.rds`) for later inspection.

### F. `restrictions_br` — land-conversion restrictions (lines 539–639)
Builds a table of cells where certain conversions (into Pasture, Urban,
Cropland, or OtherLand) are forbidden. Two sources are combined:
1. **Protected areas:** any cell with nonzero area in `Areas_Protegidas_FABLE.xls`
   (Proteção Integral, Terra Indígena, or Uso Sustentável) is restricted.
2. **Geographic exclusion:** a documented, ground-truthed rule excluding
   deep-interior Amazon cells in Acre, Amazonas, and Roraima (identified via
   a spatial join against `br_states.shp`, confirmed to lie within the
   Amazon biome via `br_biomes.shp`). This targets remote cells the MNL
   model was otherwise assigning implausibly large Cropland gains (driven
   by extreme standardized travel-time values), while deliberately leaving
   out Rondônia, Pará, and Mato Grosso, which contain the real,
   well-documented "arc of deforestation" frontier and should not be
   blocked.

The union of both cell sets becomes `restrictions_br`, a long-format table
of (cell, lu.from, lu.to, value = 1) restriction flags.

### G. MNL estimation (lines 642–726)
For each source land-use class:
- Selects "active" cells (those with an observed nonzero transition in
  `brazil_luc`) and pivots their observed 2015 transitions to shares.
- Builds a per-covariate prior-variance vector (`A0`): most covariates get
  a diffuse prior (`A0 = 1e4`), but `travel_time`, `log_livestock`,
  `log_pop`, and `log_crop` are given a much tighter prior (`A0 = 100`).
  This is explained in a long comment: with a flat diffuse prior, a
  handful of extreme-covariate cells produced runaway posterior betas
  (e.g. `travel_time → Urban` beta ≈ 9.19), making those same cells look
  attractive across multiple destination classes and cascading
  implausible land conversion through them across chained periods.
  Tightening these covariates' priors shrinks the posterior toward zero
  unless the data strongly supports the estimate.
- Fits `downscalr::mnlogit()` (Bayesian MCMC MNL, 100 iterations / 50
  burn-in) and averages posterior coefficients across draws.
- Classes with zero active training cells get placeholder zero-value betas.

Result: `betas_all`, a long-format table of transition-propensity
coefficients per (source class, destination class, covariate).

### H. Starting areas (lines 729–798)
- Converts `lu_stock` into the long-format `br_start_areas` table
  `downscale()` expects.
- **Harmonization:** rescales each class's per-cell values so the national
  class total matches manually specified targets in `HARMONIZE_TARGETS`
  (each cell keeps its existing share of the class total; classes omitted
  from the vector are left untouched).
- **Flat-prior blending:** builds `flat_priors`, a uniform prior applied
  with weight `PRIOR_WEIGHT = 0.7`. This blends 70% uniform / 30%
  econometric propensity per cell, redistributing the "pull" that would
  otherwise concentrate in a few extreme-covariate cells more evenly across
  the landscape.

### I. Downscale (lines 801–815)
Calls `downscalr::downscale()` with the targets, starting areas, covariates,
estimated betas, blended priors, and restrictions, chaining allocation
across all periods from 2000–2050.

### J. Save (lines 818–830)
Saves the raw result object, the extracted transition table
(`downscaled_LUC`), the betas, and the starting areas to `.rds`/`.csv`
files (`results_DS_usocob_states.rds`, `downscaled_LUC_usocob_states.rds/csv`,
etc.).

### K. Plot — Forest/Cropland 2050 map (lines 833–918)
Builds a spatial grid raster from cell geometries (`pop_2020.xlsx`'s `.geo`
column), extracts Cropland stock (diagonal `lu.from == lu.to`) for 2050,
converts to kha, and renders a classed choropleth PNG
(`cropland_2050_ct_states_breaks.png`) using fixed reference color breaks,
with optional state/biome boundary overlays if the shapefiles are present.

### L. Diagnostic: oscillation / offender cells (lines 921–942)
Prints per-period Forest/Cropland/Urban/total area for a specific
problem cell (`68690_a`) and a hardcoded list of "offender" cells known
from earlier runs to show large or oscillating 2020→2050 area changes,
reporting the maximum 2050/2020 ratio and count of cells with ratio > 2 —
a quick check on whether the priors/restrictions changes tamed the
runaway-cell problem described in Section G.

### M. Verification: area conservation (lines 945–1132)
Reconstructs, per class and period, `AreaStart`, `AreaEnd`, `TotalGains`,
and `TotalLosses` from `downscaled_LUC` (mirroring the FABLE Calculator
workbook's own verification layout), then checks two invariants across all
class-period combinations, within a tolerance of `1e-6` kha:
- **Continuity:** one period's closing area equals the next period's
  opening area.
- **Balance:** `AreaStart + TotalGains − TotalLosses == AreaEnd`.

Prints pass/fail summaries and any breaking rows, and saves the full
verification table (note: the code writes with `write.xlsx(...,
"area_continuity_check_usocob_states.xlsx")` but the accompanying message
refers to a `.csv` filename — a minor inconsistency to be aware of).

### N. Total area by land use per time step (lines 1134–1164)
Builds a wide summary table (one row per 5-year time step 2000–2050, one
column per class, plus a `Total` column) from the verification table, as a
quick check that total area stays constant over time (since `downscale()`
should conserve total area), and saves it to
`area_by_landuse_timestep_kha.xlsx`.

### O. Covariate maps (lines 1167–1251)
Rebuilds full cell polygons (rather than centroids) and produces one PNG
choropleth per covariate — both the modeled (log-transformed, standardized)
version and its raw/non-log counterpart — plus a combined faceted panel
with covariates z-scored onto a shared color scale, all saved under a
`covariate_maps/` subdirectory.

---

## Key outputs

| File | Content |
|---|---|
| `results_DS_usocob_states.rds` | Full `downscale()` result object (allocations + solver diagnostics) |
| `downscaled_LUC_usocob_states.rds` / `.csv` | Long-format gridded land-use transitions by period |
| `betas_all_usocob_states.rds` | Estimated MNL coefficients |
| `start_areas_usocob_states.rds` | Harmonized starting areas per cell/class |
| `cropland_2050_ct_states_breaks.png` | Choropleth of 2050 Cropland stock |
| `area_continuity_check_usocob_states.xlsx` | Area-conservation verification table |
| `area_by_landuse_timestep_kha.xlsx` | National area by class, per time step |
| `covariate_maps/*.png` | Per-covariate and combined choropleths |

---

## Things to check before running elsewhere

- Update `setwd(...)` and the two hardcoded shapefile paths
  (`br_states.shp`, `br_biomes.shp`) for your own machine.
- Confirm whether the MapBiomas sheet should be `usoecob2000` or
  `usoecob2020` (the header comment, a `message()`, and the actual
  `read_excel()` call don't all agree).
- `HARMONIZE_TARGETS` and `PRIOR_WEIGHT` are manually tuned constants —
  revisit them if the input data changes.
- Section L's offender-cell list and Section K's color breaks are
  hardcoded from a specific prior run and are diagnostic/reference aids
  rather than general-purpose code.
