# =============================================================================
# brazil_downscalr_usocob.R
# Downscaling with MapBiomas UsoCobertura land-use stock
# (UsoCobertura__fable_v2.xlsx, sheet usoecob2000) replacing hilda_br_2015.xlsx.
#
# The sheet is already in wide format (one row per cell, one column per
# land-use class in ha) — no pivot or comma conversion needed.
# Only the (ha) columns and id_c are used. 
#
# TO CHANGE THE COLUMN -> LU CLASS MAPPING: edit only `col_to_lu` in Section A.
#
# Inputs (all in working directory):
#   UsoCobertura__fable_v2.xlsx   <- MapBiomas stock (replaces hilda_br_2015)
#   hildaluc_br_2015_2019.xlsx    <- HILDA+ transition flows (unchanged)
#   brazil_fable_ct.xlsx          <- FABLE targets (unchanged)
#   altitude_br.xlsx, slope_br.xlsx, travel_time_br.xlsx,
#   livestock_br.xlsx, pop_2020.xlsx, crop_yield_br.xlsx
#   pa_br.xlsx
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(downscalr)
})

# EDIT THIS to point at your own project folder before running.
# (All relative inputs/outputs in this script are read/written from here.)
setwd("C:/Users/User/Desktop/FABLE/downscalr")
# ==============================================================================
# downscalr_patch.R
#
# Runtime monkey-patch for a bug in downscalr::solve_biascorr.mnl() (R/solve_biascorr.R).
#
# Bug: 'error_restrictions = restr.mat[,error_ind]' drops to a 1D vector when
# exactly ONE (lu.from,lu.to) target triggers the bias-correction fallback
# (i.e. sum(error_ind) == 1). The subsequent 'error_restrictions[,ccc]' then
# fails with: 'Error in error_restrictions[, ccc] : incorrect number of
# dimensions'. Fix: add drop = FALSE so it stays a (1-column) matrix.
#
# This block redefines solve_biascorr.mnl with that one-line fix and installs
# it into the downscalr namespace via assignInNamespace(). Must run AFTER
# library(downscalr). No package reinstall or separate file needed.
#
# Upstream fix: change line ~207 of R/solve_biascorr.R in tkrisztin/downscalr
# from 'restr.mat[,error_ind]' to 'restr.mat[,error_ind,drop=FALSE]'.
# ==============================================================================

solve_biascorr.mnl = function(targets,areas,xmat,betas,priors = NULL,restrictions=NULL,
                              options = downscale_control()) {
  lu.from <- unique(targets$lu.from)
  lu.to <- unique(targets$lu.to)
  ks = unique(betas$ks)
  
  out.solver <- list()
  curr.lu.from <- lu.from[1]
  full.out.res = NULL
  for(curr.lu.from in lu.from){
    err.txt = paste0(curr.lu.from," ",options$err.txt)
    
    # Extract targets
    curr.targets = dplyr::filter(targets,lu.from == curr.lu.from)$value
    names(curr.targets) <- targets$lu.to[targets$lu.from == curr.lu.from]
    curr.lu.to = names(curr.targets)
    
    # Extract betas
    curr.betas = dplyr::filter(betas,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
      tidyr::pivot_wider(names_from = "lu.to",values_from = "value",id_cols = "ks") %>%
      tibble::column_to_rownames(var = "ks")
    curr.betas = as.matrix(curr.betas)
    
    # Extract xmat
    ## IMPORTANT CHECK ORDER OF VARIABLES FIXED
    curr.xmat = dplyr::filter(xmat,ks %in% row.names(curr.betas)) %>%
      tidyr::pivot_wider(names_from = "ks",values_from = "value",id_cols = "ns")  %>%
      tibble::column_to_rownames(var = "ns") %>%
      dplyr::select(rownames(curr.betas))
    curr.xmat = as.matrix(curr.xmat)
    
    # Extract areas
    curr.areas = dplyr::filter(areas,lu.from == curr.lu.from)$value
    names(curr.areas) <- areas$ns[areas$lu.from == curr.lu.from]
    # BUGFIX: MW, order of curr.areas wrong, need to re-arrange based on xmat
    if (nrow(curr.xmat) > 0) {
      curr.areas = curr.areas[match(rownames(curr.xmat),names(curr.areas))]
    }
    
    # Extract priors
    ## IMPORTANT CHECK ORDER OF VARIABLES SIMILARLY TO XMAT
    if (!is.null(priors) && any(priors$lu.from == curr.lu.from)) {
      curr.priors = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.prior_weights = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "weight",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      name_match = match(names(curr.areas),row.names(curr.priors))
      curr.priors = curr.priors[name_match,,drop = FALSE]
      curr.prior_weights = curr.prior_weights[name_match,,drop = FALSE]
      # check if betas have been provided for priors already
      mixed_priors =c()
      nonmixed_priors = colnames(curr.priors)
      if (any(colnames(curr.betas) %in% colnames(curr.priors))) {
        mixed_priors = colnames(curr.betas)[which(colnames(curr.betas) %in% colnames(curr.priors))]
        nonmixed_priors = nonmixed_priors[!nonmixed_priors %in% mixed_priors]
        #warning(paste0(err.txt,
        #               "Priors provided for lu.from/lu.to combinations for which betas exist.\n These will be overwritten."))
        #curr.betas = curr.betas[,-which(colnames(curr.betas) %in% colnames(curr.priors)),drop = FALSE]
      }
      curr.priors = as.matrix(curr.priors)
    } else {curr.priors = NULL}
    
    # Extract restrictions
    if (!is.null(restrictions) && any(restrictions$lu.from == curr.lu.from)) {
      curr.restrictions = dplyr::filter(restrictions,lu.from == curr.lu.from) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.restrictions = curr.restrictions[match(names(curr.areas),row.names(curr.restrictions)),,drop = FALSE]
      curr.restrictions = as.matrix(curr.restrictions)
    } else {curr.restrictions = NULL}
    
    p = length(curr.targets)
    n = length(curr.areas)
    p1 = ncol(curr.betas)
    k = nrow(curr.betas)
    
    if (p1 > 0 ) {
      if (ncol(curr.xmat)!=k || nrow(curr.xmat)!=n) {
        stop(paste0(err.txt,"Dimensions of xmat, areas and betas do not match."))
      }
    }
    if (!is.null(curr.priors)) {
      #p2 = ncol(curr.priors)
      p2 = length(nonmixed_priors)
      p2_mixed = length(mixed_priors)
      if (any(curr.priors<0)) {stop(paste0(err.txt,"Priors must be strictly non-negative."))}
    } else {p2 = 0;p2_mixed = 0}
    
    # check restrictions for consistency
    if (!is.null(curr.restrictions) & any(colnames(curr.restrictions) %in% names(curr.targets))) {
      restr.mat = matrix(0,n,p); colnames(restr.mat) = names(curr.targets)
      restr.mat[,colnames(curr.restrictions)] = curr.restrictions
    } else {restr.mat = NULL}
    
    # out.res contains downscaled estimates; priors.mu econometric & other priors for estimation
    out.res = priors.mu = matrix(0,n,p)
    colnames(out.res) = colnames(priors.mu) = names(curr.targets)
    # match econometric priors
    priors.mu[,colnames(curr.betas)] = curr.xmat %*% curr.betas
    # make sure the priors are numerically well behaved
    priors.mu[priors.mu > options$MAX_EXP] = options$MAX_EXP
    priors.mu[priors.mu < -options$MAX_EXP] = -options$MAX_EXP
    priors.mu[,colnames(curr.betas)] = exp(priors.mu[,colnames(curr.betas)])
    # match other priors (if they exist)
    if (p2 > 0) {
      priors.mu[,nonmixed_priors] = curr.priors[,nonmixed_priors]
    }
    if (p2_mixed > 0) {
      w1 = curr.prior_weights[,mixed_priors,drop = FALSE]
      #re-scale exogeneous prior to priors.mu
      eco.priors_sum = apply(priors.mu[,mixed_priors,drop = FALSE],c(2),sum)
      exo.priors = curr.priors[,mixed_priors,drop = FALSE]
      exo.priors_sum = apply(exo.priors, 2, sum)
      exo.priors = t(
        (t(exo.priors) / exo.priors_sum) * eco.priors_sum   )
      priors.mu[,mixed_priors] = as.matrix((1-w1)*priors.mu[,mixed_priors] + w1*exo.priors)
    }
    # remove targets that are all zero
    not.zero = (curr.targets != 0)
    if (all(curr.targets == 0)) {
      
      #catch case if all targets are equal zero
      out.solver[[curr.lu.from]] = NULL
    } else {
      
      #cut out zero targets from targets and priors
      if (any(curr.targets == 0) && !all(curr.targets == 0)) {
        curr.targets = curr.targets[not.zero]
        priors.mu = priors.mu[,not.zero,drop = FALSE]
        if (!is.null(curr.restrictions)) {restr.mat = restr.mat[,not.zero,drop = FALSE]}
      }
      
      #proceed with bias correction
      x0 = curr.targets / sum(curr.targets + 1)
      opts <- list(algorithm = options$algorithm,
                   xtol_rel = options$xtol_rel,
                   xtol_abs = options$xtol_abs,
                   maxeval = options$maxeval
      )
      # check if the optimiser uses gradient or not
      if (grepl("_LD_",opts$algorithm)) {
        eval_grad_f = grad_sqr_diff.mnl
      } else {
        eval_grad_f = NULL
      }
      res.x = nloptr::nloptr(x0 = x0,
                             eval_f = sqr_diff.mnl,
                             eval_grad_f = eval_grad_f,
                             lb = rep(exp(-options$MAX_EXP),length(x0)),
                             ub = rep(exp(options$MAX_EXP),length(x0)),
                             opts=opts,
                             mu = priors.mu,areas = curr.areas,targets = curr.targets,
                             restrictions = restr.mat,cutoff = options$cutoff)
      res.x$par = res.x$solution
      
      out.mu = mu.mnl(res.x$solution[1:length(curr.targets)],priors.mu,curr.areas,restr.mat,options$cutoff)
      
      # Check where out.mu deviates more from the target than max_diff
      error_ind = (curr.targets - colSums(out.mu))^2 > options$max_diff
      # If any diff is larger, do individual logit models boosted by grid search
      if (any(error_ind)) {
        error_targets = curr.targets[error_ind]
        error_restrictions = restr.mat[,error_ind,drop=FALSE]
        
        # Calculate residual areas -  res_areas
        res_areas = curr.areas
        if (any(!error_ind)) {
          res_areas = res_areas - rowSums(out.mu[,!error_ind,drop=FALSE])
        }
        
        # Loop over remaining targets
        for (ccc in 1:length(error_targets)) {
          curr_error_target = error_targets[ccc]
          curr_error_restrictions = error_restrictions[,ccc]
          curr_error_mu = priors.mu[,names(curr_error_target),drop=FALSE]
          
          # Do an iterated grid search to find correct scaling coefficient for priors
          curr_scaling = iterated_grid_search(min_param = -options$MAX_EXP,
                                              max_param = options$MAX_EXP,
                                              func = sqr_diff.mnl,
                                              max_iterations = 10,
                                              precision_threshold = 1e-3,
                                              exp_transform = TRUE,
                                              mu = curr_error_mu,areas = res_areas,
                                              targets = curr_error_target,
                                              restrictions = curr_error_restrictions,
                                              cutoff = options$cutoff)
          
          # Re-scale prior
          curr_error_mu = curr_error_mu * curr_scaling$best_param
          
          # Optimize with scaled priors
          res.x = nloptr::nloptr(x0 = 1,
                                 eval_f = sqr_diff.mnl,
                                 eval_grad_f = eval_grad_f,
                                 lb = exp(-options$MAX_EXP),
                                 ub = exp(options$MAX_EXP),
                                 opts=opts,
                                 mu = curr_error_mu,areas = res_areas,targets = curr_error_target,
                                 restrictions = curr_error_restrictions,cutoff = options$cutoff)
          
          # Calculate areas of current target with mu.mnl
          curr_error_out.mu =
            mu.mnl(res.x$solution,
                   curr_error_mu,res_areas,curr_error_restrictions,
                   options$cutoff)
          
          # Add calculated areas to out.mu
          out.mu[,names(curr_error_target)] = curr_error_out.mu
          
          # Substract areas from res.areas
          res_areas = res_areas - curr_error_out.mu
          
          # Add note to res.x
          res.x$message = "INDIVIDUAL LOGIT BOOSTED BY GRID SEARCH: Standard optimization failed to converge"
        }
      }
      
      if (all(not.zero)) {out.res = out.mu
      } else {out.res[,not.zero] = out.mu}
      out.solver[[curr.lu.from]] = res.x
    }
    
    # add residual own flows in output
    out.res2 = data.frame(ns = names(curr.areas),
                          curr.areas - rowSums(out.res),out.res)
    colnames(out.res2)[2] = paste0(curr.lu.from)
    # pivot into long format
    res.agg <- out.res2 %>%
      pivot_longer(cols = -c("ns"),names_to = "lu.to") %>%
      bind_cols(lu.from = curr.lu.from)
    
    # aggregate results over dataframes
    if(curr.lu.from==lu.from[1]){
      full.out.res <- res.agg
    } else {
      full.out.res = bind_rows(full.out.res,res.agg)
    }
  }
  return(list(out.res = full.out.res, out.solver = out.solver))
}

environment(solve_biascorr.mnl) <- asNamespace("downscalr")
assignInNamespace("solve_biascorr.mnl", solve_biascorr.mnl, ns = "downscalr")
message("Patched downscalr::solve_biascorr.mnl (error_restrictions drop=FALSE fix)")


# ── Inline downscalr patch (solve_biascorr.R drop=FALSE bug) ──────────────────
# Paste the full patch body here, or:
# source("downscalr_patch.R")

LU_CLASSES <- c("Forest", "OtherLand", "Cropland", "Pasture", "Urban")

# EDIT THIS to point at the folder containing br_states.shp / br_biomes.shp
# on your own machine (used in Sections F, K, and O).
SHAPEFILE_DIR <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE"

# =============================================================================
# A.  COLUMN -> LU CLASS MAPPING    <-- edit this to adjust classification
# =============================================================================
# Keys = exact column names from the (ha) columns in usoecob2000 (without " (ha)").
# Values = one of: Forest | OtherLand | Cropland | Pasture | Urban
# Columns not listed here are excluded (NoData, Aquaculture, etc.).

col_to_lu <- c(
  # Forest
  "Forest"         = "Forest",
  "Savanna"        = "Forest",
  "Mangrove"       = "Forest",
  "FloodabForest"  = "Forest",
  "WoodSandB"      = "Forest",
  "Wetland"        = "Forest",
  "Grassland"      = "Forest",
  "HerbaceSand"    = "Forest",
  "Other non Forest Formations" = "Forest",

  # OtherLand
  "BeachSand"      = "OtherLand",
  "RockyOutcrop"   = "OtherLand",
  "HypersalineTF"  = "OtherLand",
  "Aquatic"        = "OtherLand",
  "Mining"         = "OtherLand",
  "OtherNonVeg"    = "OtherLand",
  "PhotovoltaicPP" = "OtherLand",
  "Aquaculture"    = "OtherLand",
  "ComForestry"    = "OtherLand",

  # Cropland
  "Sugarcane"      = "Cropland",
  "MosaicRural"    = "Cropland",
  "PalmOil"        = "Cropland",
  "Soybean"        = "Cropland",
  "Rice"           = "Cropland",
  "OthTemCrop"     = "Cropland",
  "Coffee"         = "Cropland",
  "Citrus"         = "Cropland",
  "OthPereCrop"    = "Cropland",
  "Cotton"         = "Cropland",

  # Pasture
  "Pasture"        = "Pasture",

  # Urban
  "Urban"          = "Urban"
)

# =============================================================================
# B.  Load usoecob2000 and build LU stock table
# =============================================================================

USOCOB_SHEET <- "usoecob2000"  # matches the file header comment (base year 2000)
message("Reading UsoCobertura__fable_v2.xlsx | sheet: ", USOCOB_SHEET, " ...")
raw <- read_excel("UsoCobertura__fable_v2.xlsx", sheet = USOCOB_SHEET)

# Select id_c + all (ha) columns
ha_cols <- names(raw)[grepl("\\(ha\\)", names(raw))]
message("Found ", length(ha_cols), " (ha) columns: ",
        paste(head(ha_cols, 5), collapse = ", "), " ...")

usocob <- raw %>%
  select(id_c, all_of(ha_cols)) %>%
  mutate(id_c = as.character(id_c)) %>%
  rename(ns = id_c)

# Strip " (ha)" suffix from column names for matching against col_to_lu keys
names(usocob) <- gsub(" \\(ha\\)$", "", names(usocob))

# Aggregate columns into 5 LU classes
# For each LU class, sum all mapped columns per cell
lu_stock <- tibble(ns = usocob$ns)

for (lu in LU_CLASSES) {
  mapped_cols <- names(col_to_lu)[col_to_lu == lu]
  present     <- intersect(mapped_cols, names(usocob))
  missing     <- setdiff(mapped_cols, names(usocob))
  if (length(missing) > 0)
    warning("Columns not found in sheet (skipped): ", paste(missing, collapse = ", "))
  lu_stock[[paste0("lu_", lu)]] <- rowSums(usocob[, present, drop = FALSE],
                                            na.rm = TRUE)
}

message("LU stock built | totals (ha):")
for (lu in LU_CLASSES)
  message(sprintf("  %-10s: %12.0f ha", lu, sum(lu_stock[[paste0("lu_", lu)]])))


# =============================================================================
# C.  brazil_FABLE
# =============================================================================

read_br <- function(path) {
  df <- read_excel(path)
  df[ , !names(df) %in% c(".geo", "system.index")]
}

brazil_FABLE <- read_br("brazil_fable_UP50_CT.xlsx") %>%
  rename(lu.from = LandCoverInit, times = YearEnd) %>%
  select(-YearStart) %>%
  pivot_longer(starts_with("To"), names_to = "lu.to", values_to = "value") %>%
  mutate(
    lu.to = sub("^To", "", lu.to),
    # FABLE Calculator outputs are in 1000 ha (kha); start.areas, xmat, and
    # betas are all in raw ha throughout this script. Convert here so
    # downscale() reconciles targets against the stock in matching units.
    value = value * 1000
  ) %>%
  select(lu.from, times, lu.to, value) %>%
  filter(lu.from != "NewForest", lu.to != "NewForest", times >= 2000)

message("brazil_FABLE: ", nrow(brazil_FABLE), " rows | periods: ",
        paste(sort(unique(brazil_FABLE$times)), collapse = ", "))


# =============================================================================
# D.  brazil_luc (HILDA+ transitions — unchanged)
# =============================================================================

hilda_labels <- c(
  "11" = "Urban",     "22" = "Cropland",  "33" = "Pasture",
  "44" = "Forest",    "55" = "OtherLand", "66" = "OtherLand",
  "77" = "OtherLand", "88" = "OtherLand", "99" = "OtherLand"
)

brazil_luc <- read_br("hildaluc_br_2015_2019.xlsx") %>%
  rename(ns = id_c) %>%
  mutate(ns = as.character(ns)) %>%
  pivot_longer(-ns, names_to = "code", values_to = "value") %>%
  mutate(
    code    = sub("^X", "", code),
    lu.from = hilda_labels[substr(code, 1, 2)],
    lu.to   = hilda_labels[substr(code, 3, 4)],
    Ts      = 2015L
  ) %>%
  filter(!is.na(lu.from), !is.na(lu.to)) %>%
  group_by(ns, lu.from, lu.to, Ts) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

message("brazil_luc: ", nrow(brazil_luc), " rows | ",
        n_distinct(brazil_luc$ns), " cells")


# =============================================================================
# E.  brazil_xmat  (UsoCobertura LU shares as covariates)
# =============================================================================

impute_median <- function(df, cols) {
  for (cl in cols) {
    x <- df[[cl]]
    if (anyNA(x)) {
      med <- median(x, na.rm = TRUE)
      message(sprintf("  Imputed %d NA in '%s' with median = %.4f",
                       sum(is.na(x)), cl, med))
      x[is.na(x)] <- med; df[[cl]] <- x
    }
  }
  df
}

altitude_df <- read_br("altitude_br.xlsx") %>%
  rename(ns = id_c, altitude_mean = MeanAltitude) %>%
  mutate(ns = as.character(ns)) %>% select(ns, altitude_mean)

slope_df <- read_br("slope_br.xlsx") %>%
  rename(ns = id_c, slope_mean = MeanSlope) %>%
  mutate(ns = as.character(ns)) %>% select(ns, slope_mean)

travel_df <- read_br("travel_time_br.xlsx") %>%
  rename(ns = id_c, travel_time = MeanTravelTime) %>%
  mutate(ns = as.character(ns)) %>% select(ns, travel_time)

livestock_df <- read_br("livestock_br.xlsx") %>%
  rename(ns = id_c, livestock_all = alllivestock) %>%
  mutate(ns = as.character(ns)) %>% select(ns, livestock_all)

pop_df <- read_br("pop_2020.xlsx") %>%
  rename(ns = id_c, pop_2020 = TotalPop) %>%
  mutate(ns = as.character(ns)) %>% select(ns, pop_2020)

crop_df <- read_br("crop_yield_br.xlsx") %>%
  rename(ns = id_c, crop_total = total) %>%
  mutate(ns = as.character(ns)) %>% select(ns, crop_total)

xmat_wide <- lu_stock %>%
  left_join(altitude_df,  by = "ns") %>%
  left_join(slope_df,     by = "ns") %>%
  left_join(travel_df,    by = "ns") %>%
  left_join(livestock_df, by = "ns") %>%
  left_join(pop_df,       by = "ns") %>%
  left_join(crop_df,      by = "ns")

numeric_cols <- setdiff(names(xmat_wide), "ns")
xmat_wide <- impute_median(xmat_wide, numeric_cols)
stopifnot(!anyNA(xmat_wide))

xmat_wide <- xmat_wide %>%
  mutate(
    log_livestock = log1p(livestock_all),
    log_pop       = log1p(pop_2020),
    log_crop      = log1p(crop_total)
  ) %>%
  select(-livestock_all, -pop_2020, -crop_total)

xmat_numeric_cols <- setdiff(names(xmat_wide), "ns")
stopifnot(!anyNA(xmat_wide),
          all(vapply(xmat_wide[xmat_numeric_cols],
                     function(x) all(is.finite(x)), logical(1))))

# Global standardisation
xmat_full <- xmat_wide %>% column_to_rownames("ns")
ok_var    <- vapply(xmat_full, function(col) !anyNA(col) && var(col) > 0, logical(1))
xmat_full <- xmat_full[ , ok_var, drop = FALSE]
if (any(!ok_var))
  message("Dropped zero-variance covariates: ",
          paste(names(ok_var)[!ok_var], collapse = ", "))

col_means <- colMeans(xmat_full)
col_sds   <- apply(xmat_full, 2, sd); col_sds[col_sds == 0] <- 1
xmat_std  <- scale(xmat_full, center = col_means, scale = col_sds)
stopifnot(all(is.finite(xmat_std)))

X_long <- as.data.frame(xmat_std) %>%
  rownames_to_column("ns") %>%
  pivot_longer(-ns, names_to = "ks", values_to = "value")

# DIAGNOSTIC: save intermediate objects for area-conservation investigation
saveRDS(xmat_wide, "xmat_wide_diag.rds")
saveRDS(X_long,    "X_long_diag.rds")
saveRDS(lu_stock,  "lu_stock_diag.rds")

message("brazil_xmat: ", n_distinct(X_long$ks), " covariates | ",
        n_distinct(X_long$ns), " cells")


# =============================================================================
# F.  restrictions_br  (from Areas_Protegidas_FABLE.xls, sheet "FABLE")
# =============================================================================
# Three protected-area categories (all in ha):
#   Proteção Integral  — strict protection (no extractive use)
#   Terra Indígena     — indigenous territory
#   Uso Sustentável    — sustainable use (limited extractive allowed)
# Any cell with total PA area > 0 is restricted: Cropland, Pasture, and Urban
# expansion into it is forbidden (max.change = 0 / value = 1).

message("Reading Areas_Protegidas_FABLE.xls | sheet: FABLE ...")
pa_raw <- readxl::read_xls("Areas_Protegidas_FABLE.xls", sheet = "FABLE")

# Select id_c + all (ha) columns
pa_ha_cols <- names(pa_raw)[grepl("\\(ha\\)", names(pa_raw))]
message("PA (ha) columns: ", paste(pa_ha_cols, collapse = ", "))

pa_df <- pa_raw %>%
  select(id_c...1, all_of(pa_ha_cols)) %>%
  mutate(
    id_c     = as.character(id_c...1),
    pa_total = rowSums(across(all_of(pa_ha_cols)), na.rm = TRUE)
  )

protected_cells <- pa_df %>% filter(pa_total > 0) %>% pull(id_c)
message("Protected cells (PA overlap): ", length(protected_cells), " / ", nrow(pa_df))

# ── Geographic exclusion: deep-interior western Amazon ───────────────────────
# The PA dataset alone under-restricts a specific set of cells: deep-interior
# Amazon cells (Acre, Amazonas, Roraima -- hundreds of km from any real road
# network or agricultural frontier) that have ZERO or negligible official PA
# overlap, yet the MNL model (driven by extreme standardised travel_time
# values at these remote cells) assigns them large Cropland stocks by 2050 --
# 25,000-55,000+ ha in cells that should be near-pristine forest.
#
# This is geographically distinct from the LEGITIMATE Rondonia/Mato Grosso
# "arc of deforestation" frontier, which has real, well-documented agricultural
# expansion and should NOT be blocked. A pure lon/lat bounding box was tried
# first, but it's an arbitrary cutoff with no ground truth. Using br_states.shp
# to check ACTUAL STATE membership is far more precise and defensible: a
# spatial join confirms ALL 7 known deep-interior offenders fall in Acre (AC)
# or Roraima (RR), while ALL 6 known legitimate frontier cells fall in
# Rondonia (RO) -- a perfect, ground-truthed separation with zero manual
# coordinate tuning. Amazonas (AM) is included in the exclusion state list on
# the same logic (its interior is equally remote/pristine), while Rondonia,
# Para, and Mato Grosso are deliberately left OUT since they contain the real,
# actively-monitored deforestation frontier.
#
# Final rule: a cell is excluded if it is BOTH (a) inside Acre, Amazonas, or
# Roraima (br_states.shp spatial join) AND (b) genuinely inside the Amazon
# biome polygon (br_biomes.shp, GID0==1) -- the state check alone already
# implies biome membership here (verified: all 843 AC/AM/RR cells are inside
# the Amazon biome), but the biome check is kept as a defensive redundancy.

message("Reading br_states.shp and br_biomes.shp for state/biome boundaries ...")
sf::sf_use_s2(FALSE)  # avoid spherical-geometry validity errors on these shapefiles

states_shp <- sf::st_read(file.path(SHAPEFILE_DIR, "br_states.shp"), quiet = TRUE) %>% sf::st_make_valid()
biomes     <- sf::st_read(file.path(SHAPEFILE_DIR, "br_biomes.shp"), quiet = TRUE)
amazon_biome <- biomes %>% filter(GID0 == 1)  # Amazonia

DEEP_INTERIOR_STATES <- c("AC", "AM", "RR")  # Acre, Amazonas, Roraima

geo_raw_excl <- read_excel("pop_2020.xlsx") %>% select(id_c, .geo)
cell_geoms <- lapply(seq_len(nrow(geo_raw_excl)), function(i) {
  g <- tryCatch(sf::st_read(geo_raw_excl$.geo[i], quiet = TRUE), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  sf::st_centroid(g$geometry[[1]])
})
valid_geom <- !sapply(cell_geoms, is.null)

cell_pts <- sf::st_sf(id_c = geo_raw_excl$id_c[valid_geom],
                       geometry = sf::st_sfc(cell_geoms[valid_geom]))
sf::st_crs(cell_pts) <- 4326
cell_pts_states <- sf::st_transform(cell_pts, sf::st_crs(states_shp))
cell_pts_biome  <- sf::st_transform(cell_pts, sf::st_crs(amazon_biome))

cell_info <- sf::st_join(cell_pts_states, states_shp["SIGLAUF3"]) %>%
  sf::st_drop_geometry() %>%
  rename(state = SIGLAUF3) %>%
  mutate(in_amazon_biome = sf::st_intersects(cell_pts_biome, amazon_biome, sparse = FALSE)[, 1])

deep_interior_cells <- cell_info %>%
  filter(state %in% DEEP_INTERIOR_STATES, in_amazon_biome) %>%
  pull(id_c)

message("Deep-interior exclusion cells (state+biome, AC/AM/RR): ", length(deep_interior_cells),
        " / ", nrow(cell_info))

# Union of PA-protected cells and the geographic exclusion zone
all_restricted_cells <- union(protected_cells, deep_interior_cells)
message("Total restricted cells (PA union geographic): ", length(all_restricted_cells))

restrictions_br <- expand.grid(lu.from = LU_CLASSES, lu.to = c("Pasture", "Urban", "Cropland", "OtherLand"),
                                stringsAsFactors = FALSE) %>%
  filter(lu.from != lu.to) %>%
  tidyr::expand_grid(ns = all_restricted_cells) %>%
  mutate(value = 1) %>%
  select(ns, lu.from, lu.to, value)

message("restrictions_br: ", nrow(restrictions_br), " rows")


# =============================================================================
# G.  MNL estimation
# =============================================================================

active_cells <- brazil_luc %>%
  filter(lu.from != lu.to, value > 0) %>%
  distinct(ns) %>% pull(ns)

message("active_cells: ", length(active_cells))

# ── Per-covariate prior variance (A0) ────────────────────────────────────────
# mnlogit() builds beta_prior_var <- diag(k) * A0 (see downscalr::mnlogit).
# A0 as a scalar gives every covariate the same diffuse prior (1e4, near-flat).
# Passing a length-k VECTOR gives each covariate its own prior variance on the
# diagonal. Near-separation in a handful of extreme-covariate training cells
# combined with a flat A0=1e4 lets posterior betas explode (observed:
# travel_time->Urban beta=9.19, log_livestock->Urban beta=-7.79,
# travel_time->Pasture beta=8.27 -- all from cells with extreme standardised
# travel_time, ~3+ SD out). These extreme betas make the same handful of
# cells look overwhelmingly attractive across MULTIPLE destination classes,
# and because xmat is static across all 7 chained periods in downscale(),
# the same distortion re-applies every period, cascading land through those
# cells (Forest->Cropland one period, that Cropland->Urban the next, etc.)
# rather than spreading it plausibly. Tightening these covariates' priors
# shrinks the posterior toward zero unless the data strongly supports it.

A0_DEFAULT <- 1e4
A0_PER_COVARIATE <- c(
  travel_time    = 100,
  log_livestock  = 100,
  log_pop        = 100,
  log_crop       = 100
)
build_A0_vector <- function(covariate_names) {
  a0 <- rep(A0_DEFAULT, length(covariate_names))
  names(a0) <- covariate_names
  override <- intersect(names(A0_PER_COVARIATE), covariate_names)
  a0[override] <- A0_PER_COVARIATE[override]
  a0
}

betas_list <- list()

for (cl in LU_CLASSES) {

  Yraw   <- brazil_luc %>% filter(lu.from == cl & Ts == 2015)
  Y_full <- Yraw %>%
    pivot_wider(id_cols = ns, names_from = lu.to, values_from = value,
                values_fill = 0) %>%
    column_to_rownames("ns")
  for (mc in setdiff(LU_CLASSES, colnames(Y_full))) Y_full[[mc]] <- 0
  Y_full <- Y_full[, LU_CLASSES]

  baseline <- which(colnames(Y_full) == cl)
  keep     <- rownames(Y_full) %in% active_cells & rowSums(Y_full) > 0
  Y        <- Y_full[keep, , drop = FALSE]

  if (nrow(Y) == 0) {
    betas_list[[cl]] <- expand.grid(ks = colnames(xmat_full),
                                     lu.to = setdiff(LU_CLASSES, cl),
                                     stringsAsFactors = FALSE) %>%
      mutate(lu.from = cl, value = 0) %>% select(ks, lu.to, value, lu.from)
    message(sprintf("%-10s: 0 active cells -> placeholder betas", cl)); next
  }

  Y <- Y / rowSums(Y)
  X <- xmat_std[rownames(Y), , drop = FALSE]
  stopifnot(all(is.finite(X)), all(is.finite(as.matrix(Y))))
  message(sprintf("%-10s: %d training cells", cl, nrow(Y)))

  set.seed(42)
  A0_vec <- build_A0_vector(colnames(X))
  res_mnl <- mnlogit(as.matrix(X), as.matrix(Y), baseline = baseline,
                      niter = 100, nburn = 50, A0 = A0_vec)

  pred_coeff <- apply(res_mnl$postb[, -baseline, , drop = FALSE], c(1, 2), mean)
  betas_list[[cl]] <- as.data.frame(pred_coeff) %>%
    rownames_to_column("ks") %>%
    pivot_longer(-ks, names_to = "lu.to", values_to = "value") %>%
    mutate(lu.from = cl)
}

betas_all <- bind_rows(betas_list) %>% select(ks, lu.to, value, lu.from)
stopifnot(all(is.finite(betas_all$value)))
message("betas_all: ", nrow(betas_all), " rows")


# =============================================================================
# H.  Starting areas from UsoCobertura 2020 stock
# =============================================================================

br_start_areas <- lu_stock %>%
  pivot_longer(-ns, names_to = "lu.from", values_to = "value") %>%
  mutate(lu.from = sub("^lu_", "", lu.from))

message("br_start_areas: total = ",
        round(sum(br_start_areas$value) / 1e6, 2), " M ha")

# ── Harmonize toward manually-specified, per-class target totals ───────────
# Each class in HARMONIZE_TARGETS is rescaled independently: every cell's
# value for that class is multiplied by (target / current class total), so
# the class hits the target total exactly while every cell keeps the same
# SHARE of that class's total it had before (no reshuffling across cells,
# no change to other classes). Units: raw ha, matching br_start_areas.
# Leave a class out of the vector to leave it untouched.
HARMONIZE_TARGETS <- c(
  "Forest"   = 589855000,
  "Pasture"  = 170839000,
  "Cropland" = 51749000,
  "OtherLand" = 35744000,
  "Urban" = 2465000
)

if (length(HARMONIZE_TARGETS) > 0) {
  current_totals <- br_start_areas %>%
    filter(lu.from %in% names(HARMONIZE_TARGETS)) %>%
    group_by(lu.from) %>%
    summarise(current = sum(value), .groups = "drop")
  
  scale_tbl <- current_totals %>%
    mutate(
      target = HARMONIZE_TARGETS[lu.from],
      factor = target / current
    )
  
  message("\nHarmonizing br_start_areas to manual per-class targets:")
  print(as.data.frame(scale_tbl %>% select(lu.from, current, target, factor)),
        row.names = FALSE)
  
  br_start_areas <- br_start_areas %>%
    left_join(scale_tbl %>% select(lu.from, factor), by = "lu.from") %>%
    mutate(value = ifelse(!is.na(factor), value * factor, value)) %>%
    select(-factor)
  
  message("br_start_areas: total AFTER harmonization = ",
          round(sum(br_start_areas$value) / 1e6, 2), " M ha")
}

# ── OPTION 2: flat/uniform priors, blended with the econometric (beta) mu ────
# PRIOR_WEIGHT is the fraction of each cell's predicted transition propensity
# that comes from a UNIFORM prior instead of the MNL/beta model, blended as:
#   priors.mu = (1 - weight) * econometric_mu + weight * flat_prior
# (flat_prior is rescaled internally by solve_biascorr.mnl to preserve each
# class's total propensity mass, so this REDISTRIBUTES the same implied pull
# more evenly across cells rather than concentrating it in a handful of
# extreme-covariate "winner" cells). Applied as a constant weight across all
# rows -- same intent as a global "prior_weights" scalar, implemented via the
# per-row `weight` column that complete_priors()/solve_biascorr.mnl support.
PRIOR_WEIGHT <- 0.7

flat_priors <- betas_all %>%
  distinct(lu.from, lu.to) %>%
  tidyr::crossing(ns = unique(X_long$ns)) %>%
  mutate(value = 1, weight = PRIOR_WEIGHT) %>%
  select(ns, lu.from, lu.to, value, weight)

message("flat_priors: ", nrow(flat_priors), " rows | weight = ", PRIOR_WEIGHT)


# =============================================================================
# I.  Downscale
# =============================================================================

results_DS <- downscale(
  targets      = brazil_FABLE,
  start.areas  = br_start_areas,
  xmat         = X_long,
  betas        = betas_all,
  priors       = flat_priors,
  restrictions = restrictions_br
)

message("\n=== Downscaling complete ===")
print(results_DS)


# =============================================================================
# J.  Save
# =============================================================================

downscaled_LUC <- results_DS$out.res

saveRDS(results_DS,       "results_DS_usocob_states.rds")
saveRDS(downscaled_LUC,   "downscaled_LUC_usocob_states.rds")
write.csv(downscaled_LUC, "downscaled_LUC_usocob_states.csv", row.names = FALSE)
saveRDS(betas_all,        "betas_all_usocob_states.rds")
saveRDS(br_start_areas,   "start_areas_usocob_states.rds")

message("Saved: results_DS_usocob_states.rds, downscaled_LUC_usocob_states.rds/csv")


# =============================================================================
# K.  Plot -- Forest 2050, terra, reference-style breaks (1000 ha)
# =============================================================================

library(sf); library(terra); library(RColorBrewer)

geo_raw    <- read_excel("pop_2020.xlsx") %>% select(id_c, .geo)
suffix_off <- 1e6
raw_id     <- geo_raw$id_c
is_suf     <- grepl("_a$", raw_id)
numeric_id <- suppressWarnings(
  ifelse(is_suf, as.numeric(sub("_a$","",raw_id))+suffix_off, as.numeric(raw_id)))

geoms     <- lapply(geo_raw$.geo, function(g) sf::st_read(g,quiet=TRUE)$geometry[[1]])
grid_sf   <- sf::st_sf(id_c=numeric_id, geometry=sf::st_sfc(geoms,crs=4326))
grid_vect <- terra::vect(grid_sf)
tmpl      <- terra::rast(terra::ext(grid_vect), resolution=0.05,
                          crs=terra::crs(grid_vect))
ns_raster <- terra::rasterize(grid_vect, tmpl, field="id_c")

# Stock (diagonal, lu.from == lu.to) = area that stayed in Cropland by 2050,
# converted from raw ha (downscale() output units, matching start.areas)
# to 1000 ha (kha) for the reference plotting scale.
cropland_2050 <- downscaled_LUC %>%
  filter(lu.from == lu.to, lu.to == "Cropland", times == "2020") %>%
  group_by(ns) %>%
  summarise(value_kha = sum(value, na.rm = TRUE) / 1000, .groups = "drop") %>%
  mutate(id_c = suppressWarnings(
    ifelse(grepl("_a$", ns), as.numeric(sub("_a$","",ns))+suffix_off, as.numeric(ns))))

total_Mha <- sum(cropland_2050$value_kha, na.rm = TRUE) / 1000
message(sprintf("Cropland 2050 (CT) total = %.2f Mha", total_Mha))

# Reference breaks (1000 ha), scaled proportionally to this data's max
BREAKS_KHA <- c(0, 0.001, 5.9, 16.8, 37.3, 71.1, 130.4, 211.5, 260, 310)
REDS <- c("#FFFFFF","#fee0d2","#fcbba1","#fc9272","#fb6a4a",
          "#ef3b2c","#cb181d","#a50f15","#67000d")

reclass_mat <- as.matrix(cropland_2050[, c("id_c","value_kha")])
r_cont      <- terra::classify(ns_raster, reclass_mat, others = NA)
mx          <- terra::global(r_cont, "max", na.rm = TRUE)[[1]]
breaks      <- BREAKS_KHA / 310 * mx
breaks[length(breaks)] <- mx * 1.0001

rcl_mat <- cbind(c(-Inf, breaks[-c(1,length(breaks))]), breaks[-1], seq_along(breaks[-1]))
r_class <- terra::classify(r_cont, rcl = rcl_mat, right = TRUE)

fmt <- function(x) if (x==0) "0" else if (x<1) formatC(x,format="f",digits=3) else formatC(x,format="f",digits=1)
LEG_LABELS <- paste0(sapply(BREAKS_KHA[-length(BREAKS_KHA)], fmt), " - ",
                     sapply(BREAKS_KHA[-1], fmt))

BRAZIL_EXT <- terra::ext(-74, -33.5, -35, 6)

# Optional shapefiles -- skip cleanly if not present (sandbox has none)
shp_states_path <- file.path(SHAPEFILE_DIR, "br_states.shp")
shp_biomes_path <- file.path(SHAPEFILE_DIR, "br_biomes.shp")
shp_states <- if (file.exists(shp_states_path)) terra::vect(shp_states_path) else NULL
shp_biomes <- if (file.exists(shp_biomes_path)) terra::vect(shp_biomes_path) else NULL

png("cropland_2050_ct_states_breaks.png", width = 760, height = 820, res = 105, bg = "white")

terra::plot(
  terra::crop(r_class, BRAZIL_EXT),
  col    = REDS,
  type   = "classes",
  legend = FALSE,
  axes   = TRUE,
  mar    = c(3, 3, 3, 1),
  main   = sprintf("CT | Cropland | 2050 | Total: %.2f Mha", total_Mha),
  cex.main = 0.95
)

if (!is.null(shp_states)) terra::plot(shp_states, border="lightgray", lwd=0.3, add=TRUE)
if (!is.null(shp_biomes)) terra::plot(shp_biomes, border="black",     lwd=0.8, add=TRUE)

usr <- par("usr")
legend(
  x = usr[2] - (usr[2]-usr[1])*0.02, y = usr[3] + (usr[4]-usr[3])*0.02,
  legend = LEG_LABELS, fill = REDS, border = NA,
  title = "1000 ha", title.adj = 0, bty = "n", cex = 0.72, xjust = 1, yjust = 0
)

dev.off()
message("Saved: cropland_2050_ct_states_breaks.png")

saveRDS(cropland_2050, "cropland_2050_ct_states_kha.rds")


# =============================================================================
# L.  Diagnostic: 68690_a oscillation + offender cells, WITH flat priors blend
# =============================================================================

cat("\n=== 68690_a: total area per period, WITH priors (weight=0.7) ===\n")
tmp <- downscaled_LUC %>% filter(ns=="68690_a", lu.from==lu.to) %>%
  group_by(times) %>% summarise(Forest=sum(value[lu.to=="Forest"]),
                                  Cropland=sum(value[lu.to=="Cropland"]),
                                  Urban=sum(value[lu.to=="Urban"]),
                                  total=sum(value))
print(as.data.frame(tmp), row.names=FALSE)

offenders <- c("70554","70251","70561","70869","70868","68692_a","70560","71180",
               "71798","68691_a","70871","69626","72418","69000","68690_a")
totals <- downscaled_LUC %>% filter(lu.from==lu.to) %>% group_by(ns,times) %>%
  summarise(total=sum(value), .groups="drop")
wide <- totals %>% filter(times %in% c("2020","2050")) %>%
  tidyr::pivot_wider(names_from=times, values_from=total, names_prefix="y") %>%
  mutate(diff=y2050-y2020, ratio=y2050/y2020)
cat("\nMax ratio (WITH priors):", round(max(wide$ratio, na.rm=TRUE),2), "\n")
cat("Cells with ratio > 2:", sum(wide$ratio > 2, na.rm=TRUE), "\n")
print(as.data.frame(wide %>% filter(ns %in% offenders) %>% arrange(desc(diff))), row.names=FALSE)


# =============================================================================
# M.  Verification: AreaEnd(period t) == AreaStart(period t+1), per class
#     (same continuity logic as the FABLE verification workbook:
#      LandCoverInit | YearStart | YearEnd | AreaStart | AreaEnd | Check)
# =============================================================================
# In downscaled_LUC, each row is a (ns, lu.from, lu.to, times, value) transition
# for the period ENDING in `times`. For a given class c and period t:
#   AreaStart(t, c) = sum(value | lu.from == c, times == t)  -- everything c
#                     had at the start of the period, however it ends up
#   AreaEnd(t, c)   = sum(value | lu.to   == c, times == t)  -- everything
#                     that is class c at the end of the period, from any source
# Continuity requires AreaEnd(t, c) == AreaStart(t+1, c) for consecutive
# periods, i.e. the closing stock downscale() produced for a period must be
# exactly the opening stock it fed into the next period.

area_start <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  group_by(lu.from, times) %>%
  summarise(AreaStart = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.from, YearEnd = times)

area_end <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  group_by(lu.to, times) %>%
  summarise(AreaEnd = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.to, YearEnd = times)

# TotalGains(t, c)  = area that flowed INTO c from every OTHER class this period
# TotalLosses(t, c) = area that flowed OUT of c into every OTHER class this period
# (mirrors the FABLE workbook's TotalGains/TotalLosses columns; the diagonal
#  lu.from == lu.to entries are "stayed as c" and are excluded from both, since
#  they're neither a gain nor a loss for c)
total_gains <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  filter(lu.from != lu.to) %>%
  group_by(lu.to, times) %>%
  summarise(TotalGains = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.to, YearEnd = times)

total_losses <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  filter(lu.from != lu.to) %>%
  group_by(lu.from, times) %>%
  summarise(TotalLosses = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.from, YearEnd = times)

verify_tbl <- full_join(area_start, area_end, by = c("LandCoverInit", "YearEnd")) %>%
  full_join(total_gains,  by = c("LandCoverInit", "YearEnd")) %>%
  full_join(total_losses, by = c("LandCoverInit", "YearEnd")) %>%
  mutate(
    TotalGains  = tidyr::replace_na(TotalGains, 0),
    TotalLosses = tidyr::replace_na(TotalLosses, 0)
  ) %>%
  # downscaled_LUC / out.res is in raw ha throughout the script; convert to
  # 1000 ha (kha) here so this table matches the FABLE workbook's units.
  mutate(
    AreaStart   = AreaStart   / 1000,
    AreaEnd     = AreaEnd     / 1000,
    TotalGains  = TotalGains  / 1000,
    TotalLosses = TotalLosses / 1000
  ) %>%
  arrange(LandCoverInit, YearEnd) %>%
  group_by(LandCoverInit) %>%
  mutate(
    YearStart      = YearEnd - 5L,
    NextAreaStart  = lead(AreaStart),                       # AreaStart (kha) of the following period
    CheckContinuity = NextAreaStart - AreaEnd,              # should be 0: this period's close == next period's open
    CheckBalance    = AreaEnd - (AreaStart + TotalGains - TotalLosses)  # should be 0: gains/losses reconcile start->end
  ) %>%
  ungroup() %>%
  select(LandCoverInit, YearStart, YearEnd, AreaStart, AreaEnd,
         TotalGains, TotalLosses, NextAreaStart, CheckContinuity, CheckBalance)

# tolerance is in kha now (1e-6 kha = 1e-3 ha) -- still far tighter than the
# floating-point residuals (~1e-10 ha) seen in the FABLE workbook itself
TOL <- 1e-6
continuity_breaks <- verify_tbl %>% filter(abs(CheckContinuity) > TOL & !is.na(CheckContinuity))
balance_breaks    <- verify_tbl %>% filter(abs(CheckBalance)    > TOL & !is.na(CheckBalance))

cat("\n=== Area verification (per class, per period) ===\n")
print(as.data.frame(verify_tbl), row.names = FALSE)

if (nrow(continuity_breaks) == 0) {
  cat("\nOK: all", nrow(verify_tbl), "class-period joins are continuous (|diff| <=", TOL, ").\n")
} else {
  cat("\nFAILED:", nrow(continuity_breaks), "class-period joins are NOT continuous (|diff| >", TOL, "):\n")
  print(as.data.frame(continuity_breaks), row.names = FALSE)
}

if (nrow(balance_breaks) == 0) {
  cat("OK: AreaStart + TotalGains - TotalLosses == AreaEnd for all", nrow(verify_tbl), "rows (|diff| <=", TOL, ").\n")
} else {
  cat("FAILED:", nrow(balance_breaks), "rows where AreaStart + TotalGains - TotalLosses != AreaEnd:\n")
  print(as.data.frame(balance_breaks), row.names = FALSE)
}

write.xlsx(verify_tbl, "area_continuity_check_usocob_states.xlsx", row.names = FALSE)
message("Saved: area_continuity_check_usocob_states.xlsx")

# =============================================================================
# N.  Total area by land use, per time step (wide format, 1000 ha)
# =============================================================================
# One row per time step (2000, 2005, ..., 2050), one column per land-use class,
# plus a Total column (sum across classes -- should be constant/near-constant
# across time steps, since downscale() conserves total area).
#
# For 2000 (base year) we use AreaStart of the first period (2000-2005), since
# downscaled_LUC itself only carries transitions for periods ending 2005+.
# For 2005-2050 we use AreaEnd of the period ending in that year, which is the
# total stock of each class at that point in time.

area_2000 <- verify_tbl %>%
  filter(YearStart == min(YearStart)) %>%
  select(LandCoverInit, Year = YearStart, Area = AreaStart)

area_rest <- verify_tbl %>%
  select(LandCoverInit, Year = YearEnd, Area = AreaEnd)

area_by_timestep <- bind_rows(area_2000, area_rest) %>%
  distinct(LandCoverInit, Year, .keep_all = TRUE) %>%
  arrange(Year, LandCoverInit) %>%
  tidyr::pivot_wider(names_from = LandCoverInit, values_from = Area) %>%
  mutate(Total = rowSums(across(-Year), na.rm = TRUE)) %>%
  arrange(Year)

cat("\n=== Total area by land use, per time step (1000 ha) ===\n")
print(as.data.frame(area_by_timestep), row.names = FALSE)

write.xlsx(area_by_timestep, "area_by_landuse_timestep_kha.xlsx", row.names = FALSE)
message("Saved: area_by_landuse_timestep_kha.csv (values in 1000 ha)")


# =============================================================================
# O.  Covariate maps (sf)
# =============================================================================
# One choropleth per covariate (raw, unstandardized values from xmat_wide,
# Section E) plus a combined faceted panel. Reuses geo_raw_excl (built in the
# state/biome-restriction section) but keeps the full cell polygon instead of
# collapsing it to a centroid, so this needs to run AFTER that section (and
# after Section E, where xmat_wide is created).

library(ggplot2)

# Build full-polygon geometries per cell (same .geo parsing as the
# restriction-cell section, but keeping the polygon, not st_centroid())
cell_polys_list <- lapply(seq_len(nrow(geo_raw_excl)), function(i) {
  g <- tryCatch(sf::st_read(geo_raw_excl$.geo[i], quiet = TRUE), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  g$geometry[[1]]
})
valid_poly <- !sapply(cell_polys_list, is.null)

cell_polys <- sf::st_sf(
  ns       = as.character(geo_raw_excl$id_c[valid_poly]),
  geometry = sf::st_sfc(cell_polys_list[valid_poly])
)
sf::st_crs(cell_polys) <- 4326

# Join both the modeled (log-transformed) covariates and the raw,
# non-log versions onto the polygons. The raw versions are rebuilt here
# straight from the original per-source data frames (livestock_df, pop_df,
# crop_df -- read in Section E, before the log1p transform), so nothing in
# Section E needs to change.
covariate_cols <- setdiff(names(xmat_wide), "ns")

xmat_raw_cols <- livestock_df %>%
  left_join(pop_df,  by = "ns") %>%
  left_join(crop_df, by = "ns")
xmat_raw_cols <- impute_median(xmat_raw_cols, setdiff(names(xmat_raw_cols), "ns"))

map_data <- cell_polys %>%
  dplyr::left_join(xmat_wide, by = "ns") %>%
  dplyr::left_join(xmat_raw_cols, by = "ns")

# covariate_cols_all = the 6 modeled covariates (3 log-transformed + altitude/
# slope/travel_time) plus the 3 raw, non-log counterparts (livestock_all,
# pop_2020, crop_total) -- 9 maps total
covariate_cols_all <- union(covariate_cols, setdiff(names(xmat_raw_cols), "ns"))

out_dir <- "covariate_maps"
dir.create(out_dir, showWarnings = FALSE)

# One PNG per covariate (log-transformed AND their non-log/raw counterparts)
for (cov in covariate_cols_all) {
  p <- ggplot(map_data) +
    geom_sf(aes(fill = .data[[cov]]), color = NA) +
    scale_fill_viridis_c(name = cov, na.value = "grey80") +
    labs(title = paste("Covariate:", cov)) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
  
  fname <- file.path(out_dir, paste0("covariate_", cov, ".png"))
  ggsave(fname, p, width = 8, height = 7, dpi = 200)
  message("Saved: ", fname)
}

# Combined faceted panel, all covariates in one figure.
# facet_wrap + coord_sf can't use free fill scales, so all covariates share
# one legend here -- z-score each covariate first so that shared scale is
# actually comparable (raw units still used in the individual PNGs above).
map_long <- map_data %>%
  tidyr::pivot_longer(all_of(covariate_cols_all), names_to = "covariate", values_to = "value") %>%
  dplyr::group_by(covariate) %>%
  dplyr::mutate(value_z = as.numeric(scale(value))) %>%
  dplyr::ungroup()

p_facet <- ggplot(map_long) +
  geom_sf(aes(fill = value_z), color = NA) +
  scale_fill_viridis_c(name = "z-score", na.value = "grey80") +
  facet_wrap(~ covariate) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "covariate_maps_all.png"), p_facet, width = 14, height = 10, dpi = 200)
message("Saved: ", file.path(out_dir, "covariate_maps_all.png"))