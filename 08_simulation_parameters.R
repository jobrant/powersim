# Overdispersion --------------------------------------------
#' Draw counts from a Beta-binomial (or Binomial) at a given mean rate
#'
#' For each element, draws `Binomial(size, p)` where `p = mean_rate` when `s` is
#' `NULL`/`Inf`, or `p ~ Beta(mean_rate * s, (1 - mean_rate) * s)` when `s` is
#' finite (overdispersion, precision `s = alpha + beta`). Boundary means (0 or 1)
#' carry no dispersion and use `p = mean_rate` directly.
#'
#' @param mean_rate Numeric vector of target mean methylation rates.
#' @param size Numeric vector of coverages (recycled/aligned to `mean_rate`).
#' @param s Optional Beta precision; `NULL`/`Inf` gives a plain Binomial.
#' @return Integer vector of methylated counts.
#' @keywords internal
.betabinom_counts <- function(mean_rate, size, s = NULL) {
  mean_rate <- pmin(1, pmax(0, mean_rate))
  p <- mean_rate
  if (!is.null(s) && is.finite(s)) {
    mid <- mean_rate > 0 & mean_rate < 1
    p[mid] <- rbeta(sum(mid), mean_rate[mid] * s, (1 - mean_rate[mid]) * s)
  }
  rbinom(length(p), size = pmax(0L, as.integer(round(size))), prob = p)
}


#' Pooled per-site true methylation rate across replicates
#'
#' Coverage-weighted: `p = sum(mc) / sum(cov)` over the replicate list. Used as
#' the ground-truth rate for parametric pseudo-group construction.
#'
#' @param replicates List of aligned replicate data.tables (need `mc`, `cov`).
#' @return Numeric vector of pooled rates (0 where total coverage is 0).
#' @keywords internal
.pooled_rate <- function(replicates) {
  tot_mc  <- Reduce(`+`, lapply(replicates, function(dt) as.numeric(dt$mc)))
  tot_cov <- Reduce(`+`, lapply(replicates, function(dt) as.numeric(dt$cov)))
  ifelse(tot_cov > 0, tot_mc / tot_cov, 0)
}


#' Draw one parametric replicate from a per-site true rate
#'
#' Reuses the template's real coordinates and coverage, but regenerates `mc`
#' from scratch as `BetaBinom(cov, p * efficiency, dispersion_s)`. Seeded
#' independently so replicates and groups share no realized noise.
#'
#' @param template Source data.table supplying coordinates and `cov`.
#' @param p Per-site true rate vector (e.g. from [.pooled_rate()]).
#' @param efficiency Scalar enzyme efficiency in (0, 1].
#' @param seed_i Integer seed for this replicate.
#' @param dispersion_s Optional Beta precision (`NULL`/`Inf` = binomial).
#' @return data.table copy of `template` with regenerated `mc` and `rate`.
#' @keywords internal
.draw_parametric_rep <- function(template, p, efficiency, seed_i,
                                 dispersion_s = NULL) {
  set.seed(seed_i)
  out <- copy(template)
  mc_new <- .betabinom_counts(p * efficiency, out$cov, dispersion_s)
  set(out, j = "mc", value = mc_new)
  set(out, j = "rate", value = ifelse(out$cov > 0, mc_new / out$cov, 0))
  out
}


#' Simulate enzyme efficiency distortion on a single sample
#'
#' For each site, resample methylated counts to reflect a simulated enzyme
#' efficiency. Two noise models are supported:
#' \describe{
#'   \item{Binomial (default)}{`new_mc ~ Binomial(original_mc, efficiency)`.
#'     This is the minimum-variance model and reproduces the original behavior
#'     exactly. Note `Binomial(mc, eff)` has the same marginal as
#'     `Binomial(cov, rate * eff)`, the form the beta-binomial path generalizes.}
#'   \item{Beta-binomial (`dispersion_s` finite)}{Draws a latent per-site rate
#'     `p ~ Beta(mu * s, (1 - mu) * s)` with `mu = rate * efficiency` and
#'     `s = dispersion_s`, then `new_mc ~ Binomial(cov, p)`. This injects
#'     overdispersion calibrated to real replicate variability (smaller `s` =
#'     more dispersion; `s = Inf` reduces to the Binomial path).}
#' }
#' Coverage stays constant in both cases (efficiency affects methylation
#' labeling, not sequencing depth).
#'
#' @param dt data.table with columns: chr, pos, strand, site, mc, cov, rate
#' @param efficiency Numeric scalar in (0, 1]. The simulated enzyme efficiency.
#'   Values < 1 reduce observed methylation proportional to accessibility.
#' @param seed Optional integer seed for reproducibility.
#' @param dispersion_s Optional Beta precision (`alpha + beta`) for beta-binomial
#'   overdispersion. `NULL` or `Inf` (default) uses the original pure-binomial
#'   model. Calibrate with `inst/scripts/diagnose_variance.R` (e.g. ~26 for the
#'   M-series data, ~101 for PrEC).
#' @return data.table with distorted mc and recalculated rate (cov unchanged).
simulate_efficiency <- function(dt, efficiency, seed = NULL, dispersion_s = NULL) {
  stopifnot(
    is.data.table(dt),
    all(c("mc", "cov", "rate") %in% names(dt)),
    is.numeric(efficiency), length(efficiency) == 1,
    efficiency > 0, efficiency <= 1
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- copy(dt)
  n <- nrow(out)
  
  if (is.null(dispersion_s) || !is.finite(dispersion_s)) {
    # Binomial resampling: each methylated read retained with prob = efficiency
    mc_new <- rbinom(n,
                     size = pmax(0L, as.integer(round(out$mc))),
                     prob = efficiency)
  } else {
    stopifnot(is.numeric(dispersion_s), length(dispersion_s) == 1,
              dispersion_s > 0)
    # Beta-binomial resample at mean rate mu = rate * efficiency, precision s.
    mc_new <- .betabinom_counts(out$rate * efficiency, out$cov, dispersion_s)
  }
  
  set(out, j = "mc", value = mc_new)
  set(out, j = "rate", value = ifelse(out$cov > 0, out$mc / out$cov, 0))
  
  return(out)
}

# Read depth(coverage, rate) --------------------------------------------
#' Simulate enzyme efficiency distortion on a single sample
#'
#' For each site, resample methylated counts to reflect a simulated enzyme
#' efficiency. Two noise models are supported:
#' \describe{
#'   \item{Binomial (default)}{`new_mc ~ Binomial(original_mc, efficiency)`.
#'     This is the minimum-variance model and reproduces the original behavior
#'     exactly. Note `Binomial(mc, eff)` has the same marginal as
#'     `Binomial(cov, rate * eff)`, the form the beta-binomial path generalizes.}
#'   \item{Beta-binomial (`dispersion_s` finite)}{Draws a latent per-site rate
#'     `p ~ Beta(mu * s, (1 - mu) * s)` with `mu = rate * efficiency` and
#'     `s = dispersion_s`, then `new_mc ~ Binomial(cov, p)`. This injects
#'     overdispersion calibrated to real replicate variability (smaller `s` =
#'     more dispersion; `s = Inf` reduces to the Binomial path).}
#' }
#' Coverage stays constant in both cases (efficiency affects methylation
#' labeling, not sequencing depth).
#'
#' @param dt data.table with columns: chr, pos, strand, site, mc, cov, rate
#' @param efficiency Numeric scalar in (0, 1]. The simulated enzyme efficiency.
#'   Values < 1 reduce observed methylation proportional to accessibility.
#' @param seed Optional integer seed for reproducibility.
#' @param dispersion_s Optional Beta precision (`alpha + beta`) for beta-binomial
#'   overdispersion. `NULL` or `Inf` (default) uses the original pure-binomial
#'   model. Calibrate with `inst/scripts/diagnose_variance.R` (e.g. ~26 for the
#'   M-series data, ~101 for PrEC).
#' @return data.table with distorted mc and recalculated rate (cov unchanged).
simulate_efficiency <- function(dt, efficiency, seed = NULL, dispersion_s = NULL) {
  stopifnot(
    is.data.table(dt),
    all(c("mc", "cov", "rate") %in% names(dt)),
    is.numeric(efficiency), length(efficiency) == 1,
    efficiency > 0, efficiency <= 1
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- copy(dt)
  n <- nrow(out)
  
  if (is.null(dispersion_s) || !is.finite(dispersion_s)) {
    # Binomial resampling: each methylated read retained with prob = efficiency
    mc_new <- rbinom(n,
                     size = pmax(0L, as.integer(round(out$mc))),
                     prob = efficiency)
  } else {
    stopifnot(is.numeric(dispersion_s), length(dispersion_s) == 1,
              dispersion_s > 0)
    # Beta-binomial resample at mean rate mu = rate * efficiency, precision s.
    mc_new <- .betabinom_counts(out$rate * efficiency, out$cov, dispersion_s)
  }
  
  set(out, j = "mc", value = mc_new)
  set(out, j = "rate", value = ifelse(out$cov > 0, out$mc / out$cov, 0))
  
  return(out)
}

# Dmr length --------------------------------------------

#' Select random genomic regions for spike-in DMRs (legacy)
#'
#' Samples contiguous regions from the data to serve as ground-truth DMRs.
#' For sparse data, prefer select_dense_spikein_regions() instead.
#'
#' @param dt Reference data.table (any sample; needs chr, pos columns).
#' @param n_regions Number of DMR regions to select.
#' @param region_width_bp Width of each region in base pairs.
#' @param min_sites Minimum number of sites a region must contain to be valid.
#' @param seed Random seed.
#' @return data.table with columns: region_id, chr, start, end, n_sites
select_spikein_regions <- function(dt,
                                   n_regions = 300,
                                   region_width_bp = 500,
                                   min_sites = 5,
                                   seed = 123) {
  set.seed(seed)
  
  # Get unique chromosomes and their site ranges
  chr_ranges <- dt[, .(min_pos = min(pos), max_pos = max(pos), n = .N), by = chr]
  chr_ranges <- chr_ranges[max_pos - min_pos > region_width_bp * 2]
  
  regions <- data.table(
    region_id = integer(),
    chr = character(),
    start = integer(),
    end = integer(),
    n_sites = integer()
  )
  
  attempts <- 0
  max_attempts <- n_regions * 20
  
  while (nrow(regions) < n_regions && attempts < max_attempts) {
    attempts <- attempts + 1
    
    # Pick a random chromosome weighted by number of sites
    chr_pick <- sample(chr_ranges$chr, 1, prob = chr_ranges$n)
    chr_info <- chr_ranges[chr == chr_pick]
    
    # Pick a random start position
    start_pos <- sample(chr_info$min_pos:(chr_info$max_pos - region_width_bp), 1)
    end_pos <- start_pos + region_width_bp
    
    # Count sites in this region
    n_sites_in_region <- dt[chr == chr_pick & pos >= start_pos & pos <= end_pos, .N]
    
    if (n_sites_in_region < min_sites) next
    
    # Check for overlap with existing regions
    if (nrow(regions) > 0) {
      overlaps <- regions[chr == chr_pick & start <= end_pos & end >= start_pos]
      if (nrow(overlaps) > 0) next
    }
    
    regions <- rbind(regions, data.table(
      region_id = nrow(regions) + 1L,
      chr = chr_pick,
      start = start_pos,
      end = end_pos,
      n_sites = n_sites_in_region
    ))
  }
  
  if (nrow(regions) < n_regions) {
    warning(sprintf("Only found %d non-overlapping regions (requested %d)",
                    nrow(regions), n_regions))
  }
  
  message(sprintf("Selected %d spike-in regions (median %d sites/region)",
                  nrow(regions), median(regions$n_sites)))
  
  return(regions)
}

# Effect size --------------------------------------------
#' Inject spike-in DMRs into a sample
#'
#' For sites within spike-in regions, shift the methylation rate by a specified
#' effect size, then resample mc from Binomial(cov, new_rate).
#'
#' Direction is chosen based on the mean baseline rate of each region to
#' maximize the achievable effect. Regions with high baseline rates receive
#' decreasing effects; regions with low baseline rates receive increasing
#' effects. This avoids wasting injections on clamped sites that cannot
#' move in the chosen direction.
#'
#' @param dt data.table of methylation data for one sample.
#' @param regions data.table of spike-in regions (from select_spikein_regions
#'   or select_dense_spikein_regions).
#' @param effect_size Numeric scalar: amount to shift rates.
#' @param direction Character: "auto" (default, rate-aware), "increase",
#'   "decrease", or "both" (random per region, legacy behavior).
#' @param low_rate_threshold Numeric: regions with mean rate below this get
#'   "increase" direction. Default 0.4.
#' @param high_rate_threshold Numeric: regions with mean rate above this get
#'   "decrease" direction. Default 0.6.
#' @param seed Random seed.
#' @return List with:
#'   \item{data}{Modified data.table}
#'   \item{truth}{data.table mapping region_id to injected effect per site}
#'   \item{stats}{Summary of achievable effects across regions}
#' @noRd
inject_spikein_dmrs <- function(dt,
                                regions,
                                effect_size = 0.15,
                                direction = "auto",
                                low_rate_threshold = 0.4,
                                high_rate_threshold = 0.6,
                                seed = NULL) {
  stopifnot(
    is.data.table(dt),
    is.data.table(regions),
    direction %in% c("auto", "increase", "decrease", "both")
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- copy(dt)
  
  # Track ground truth
  truth_list <- list()
  n_injected <- 0L
  n_skipped <- 0L
  achievable_effects <- numeric(nrow(regions))
  
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    
    # Find sites in this region
    site_idx <- which(out$chr == r$chr & out$pos >= r$start & out$pos <= r$end)
    
    if (length(site_idx) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }
    
    old_rates <- out$rate[site_idx]
    mean_rate <- mean(old_rates)
    
    # Determine direction
    dir_sign <- if (direction == "auto") {
      if (mean_rate < low_rate_threshold) {
        1L   # room to increase
      } else if (mean_rate > high_rate_threshold) {
        -1L  # room to decrease
      } else {
        # Middle range: pick direction with most room to move
        room_up   <- 1 - mean_rate
        room_down <- mean_rate
        if (room_up >= room_down) 1L else -1L
      }
    } else if (direction == "increase") {
      1L
    } else if (direction == "decrease") {
      -1L
    } else {
      # "both" — legacy random behavior
      sample(c(1L, -1L), 1L)
    }
    
    # Compute new rates
    new_rates <- pmin(1, pmax(0, old_rates + dir_sign * effect_size))
    
    # Track achievable effect (actual mean shift after clamping)
    actual_mean_effect <- mean(new_rates - old_rates)
    achievable_effects[i] <- actual_mean_effect
    
    # Resample mc from Binomial(cov, new_rate)
    covs <- out$cov[site_idx]
    new_mc <- rbinom(length(site_idx), size = as.integer(covs), prob = new_rates)
    set(out, i = site_idx, j = "mc", value = new_mc)
    set(out, i = site_idx, j = "rate", 
        value = ifelse(out$cov[site_idx] > 0, new_mc / out$cov[site_idx], 0))
    
    truth_list[[i]] <- data.table(
      region_id     = r$region_id,
      chr           = r$chr,
      pos           = out$pos[site_idx],
      original_rate = old_rates,
      injected_rate = new_rates,
      actual_effect = new_rates - old_rates,
      direction     = ifelse(dir_sign == 1L, "increase", "decrease"),
      mean_baseline = mean_rate
    )
    
    n_injected <- n_injected + 1L
  }
  
  truth <- rbindlist(truth_list)
  
  # Summary statistics
  achievable_effects <- achievable_effects[achievable_effects != 0]
  stats <- list(
    n_regions_injected = n_injected,
    n_regions_skipped  = n_skipped,
    n_sites_affected   = nrow(truth),
    mean_achievable_effect = mean(abs(achievable_effects)),
    median_achievable_effect = median(abs(achievable_effects)),
    pct_regions_gt_half_effect = mean(abs(achievable_effects) >= effect_size / 2) * 100
  )
  
  message(sprintf(
    "Injected DMRs in %d regions (%d sites, effect=%.2f, median achievable=%.3f, %.0f%% >%.2f)",
    n_injected, nrow(truth), effect_size,
    stats$median_achievable_effect,
    stats$pct_regions_gt_half_effect,
    effect_size / 2
  ))
  
  return(list(data = out, truth = truth, stats = stats))
}


#' Inject spike-in DMRs into the true rate (parametric mode)
#'
#' Unlike `inject_spikein_dmrs()`, which resamples realized counts, this shifts
#' the per-site TRUE rate `p` within spike-in regions, returning a modified rate
#' vector for group B. Replicates are then drawn from this rate by
#' [.draw_parametric_rep()], so the spike-in is a genuine biological difference
#' that precedes efficiency and sampling. Direction is rate-aware ("auto"),
#' matching `inject_spikein_dmrs()`.
#'
#' @param ref_dt Reference data.table aligned to `p` (needs `chr`, `pos`).
#' @param p Per-site true rate vector (e.g. from [.pooled_rate()]).
#' @param regions Spike-in regions (region_id, chr, start, end).
#' @param effect_size Numeric scalar: amount to shift rates.
#' @param direction "auto" (default), "increase", "decrease", or "both".
#' @param low_rate_threshold Below this mean rate, regions get "increase".
#' @param high_rate_threshold Above this mean rate, regions get "decrease".
#' @param seed Optional seed (only affects direction = "both").
#' @return List with `p` (modified group-B rate vector) and `truth` (per-site
#'   injected effects).
#' @noRd
inject_spikein_parametric <- function(ref_dt, p, regions,
                                      effect_size = 0.15, direction = "auto",
                                      low_rate_threshold = 0.4,
                                      high_rate_threshold = 0.6,
                                      seed = NULL) {
  stopifnot(
    is.data.table(ref_dt), is.data.table(regions),
    length(p) == nrow(ref_dt),
    direction %in% c("auto", "increase", "decrease", "both")
  )
  if (!is.null(seed)) set.seed(seed)
  
  chr <- ref_dt$chr
  pos <- ref_dt$pos
  p_B <- p
  truth_list <- list()
  
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    site_idx <- which(chr == r$chr & pos >= r$start & pos <= r$end)
    if (length(site_idx) == 0) next
    
    old_rates <- p[site_idx]
    mean_rate <- mean(old_rates)
    
    dir_sign <- if (direction == "auto") {
      if (mean_rate < low_rate_threshold) {
        1L
      } else if (mean_rate > high_rate_threshold) {
        -1L
      } else {
        if ((1 - mean_rate) >= mean_rate) 1L else -1L
      }
    } else if (direction == "increase") {
      1L
    } else if (direction == "decrease") {
      -1L
    } else {
      sample(c(1L, -1L), 1L)
    }
    
    new_rates <- pmin(1, pmax(0, old_rates + dir_sign * effect_size))
    p_B[site_idx] <- new_rates
    
    truth_list[[i]] <- data.table(
      region_id     = r$region_id,
      chr           = r$chr,
      pos           = pos[site_idx],
      original_rate = old_rates,
      injected_rate = new_rates,
      actual_effect = new_rates - old_rates,
      direction     = ifelse(dir_sign == 1L, "increase", "decrease"),
      mean_baseline = mean_rate
    )
  }
  
  truth <- rbindlist(truth_list)
  message(sprintf(
    "Parametric spike-in: %d regions, %d sites, effect=%.2f (median achievable=%.3f)",
    length(truth_list), nrow(truth), effect_size,
    if (nrow(truth)) median(abs(truth$actual_effect)) else 0
  ))
  
  list(p = p_B, truth = truth)
}