# =============================================================================
# feature_cta.R - Confirmatory Tetrad Analysis for PLS-SEM (CTA-PLS)
# =============================================================================
# Implements CTA-PLS for empirically testing whether a construct's measurement
# model is reflective or formative.
#
# Under a reflective (common factor) model, all model-implied vanishing tetrads
# equal zero. If any tetrad is significantly non-zero, the reflective
# specification is rejected in favour of a formative specification.
#
# References:
# - Gudergan, S. P., Ringle, C. M., Wende, S. & Will, A. (2008). Confirmatory
#   Tetrad Analysis in PLS Path Modeling. Journal of Business Research, 61(12),
#   1238-1249.
# - Cefis, M., Angelelli, M., Carpita, M. & Ciavolino, E. (2025). Confirmatory
#   Tetrad Analysis in PLS-SEM: A Multiple Testing Correction Perspective.
#   Social Indicators Research (advance online).
# - Bollen, K. A. & Ting, K.-F. (2000). A Tetrad Test for Causal Indicators.
#   Psychological Methods, 5(1), 3-22.
# =============================================================================

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Enumerate non-redundant vanishing tetrads for a set of indicators.
#'
#' For each combination of 4 indicators (i,j,k,l), there are 3 possible tetrads
#' but only 2 are algebraically independent (tau1 - tau2 + tau3 = 0 always).
#' We return 2 per 4-tuple.
#'
#' @param indicators Character vector of indicator names (length >= 4).
#' @return A data.frame with columns: i, j, k, l (indicator names) and
#'   tetrad_id (1 or 2) identifying which of the two non-redundant tetrads.
#' @noRd
enumerate_tetrads <- function(indicators) {
  combos <- utils::combn(indicators, 4)
  n_combos <- ncol(combos)

  # Pre-allocate: 2 tetrads per 4-tuple
  result <- data.frame(
    i = character(2 * n_combos),
    j = character(2 * n_combos),
    k = character(2 * n_combos),
    l = character(2 * n_combos),
    tetrad_id = integer(2 * n_combos),
    stringsAsFactors = FALSE
  )

  for (idx in seq_len(n_combos)) {
    vars <- combos[, idx]
    row_base <- (idx - 1) * 2

    # Tetrad 1: sigma_ij * sigma_kl - sigma_ik * sigma_jl
    result[row_base + 1, ] <- list(vars[1], vars[2], vars[3], vars[4], 1L)
    # Tetrad 2: sigma_ij * sigma_kl - sigma_il * sigma_jk
    result[row_base + 2, ] <- list(vars[1], vars[2], vars[3], vars[4], 2L)
  }

  result
}

#' Compute tetrad values from a covariance matrix.
#'
#' @param cov_mat Covariance matrix with named rows/columns.
#' @param tetrad_df Data.frame from enumerate_tetrads().
#' @return Numeric vector of tetrad values.
#' @noRd
compute_tetrads <- function(cov_mat, tetrad_df) {
  n <- nrow(tetrad_df)
  values <- numeric(n)

  for (r in seq_len(n)) {
    i <- tetrad_df$i[r]
    j <- tetrad_df$j[r]
    k <- tetrad_df$k[r]
    l <- tetrad_df$l[r]

    # Both tetrads share the first product: sigma_ij * sigma_kl
    prod_common <- cov_mat[i, j] * cov_mat[k, l]

    if (tetrad_df$tetrad_id[r] == 1L) {
      # tau1 = sigma_ij * sigma_kl - sigma_ik * sigma_jl
      values[r] <- prod_common - cov_mat[i, k] * cov_mat[j, l]
    } else if (tetrad_df$tetrad_id[r] == 2L) {
      # tau2 = sigma_ij * sigma_kl - sigma_il * sigma_jk
      values[r] <- prod_common - cov_mat[i, l] * cov_mat[j, k]
    } else {
      # tau3 (tau_1342) = sigma_ik * sigma_jl - sigma_il * sigma_jk
      values[r] <- cov_mat[i, k] * cov_mat[j, l] - cov_mat[i, l] * cov_mat[j, k]
    }
  }

  values
}

#' Format a tetrad label for display.
#' @noRd
format_tetrad_label <- function(i, j, k, l, tetrad_id) {
  if (tetrad_id == 1L) {
    paste0("s(", i, ",", j, ")s(", k, ",", l, ") - s(",
           i, ",", k, ")s(", j, ",", l, ")")
  } else if (tetrad_id == 2L) {
    paste0("s(", i, ",", j, ")s(", k, ",", l, ") - s(",
           i, ",", l, ")s(", j, ",", k, ")")
  } else {
    # tau_1342: i,j own; k,l borrowed
    paste0("s(", i, ",", k, ")s(", j, ",", l, ") - s(",
           i, ",", l, ")s(", j, ",", k, ")")
  }
}

#' Resolve indicators for a construct, handling HOC chaining.
#'
#' For standard constructs, returns raw indicator names and the data source.
#' For HOC constructs, returns LOC names as "indicators" and uses construct
#' scores as the data source (testing whether the HOC is reflective over LOCs).
#'
#' @param construct Construct name.
#' @param model seminr model object.
#' @return List with: indicators (character vector), data (matrix),
#'   is_hoc (logical), needs_reestimation (logical).
#' @noRd
resolve_indicators <- function(construct, model) {
  items <- seminr:::items_of_construct(construct, model)
  construct_cols <- colnames(model$outer_weights)

  # Detect HOC: items that are themselves construct names in outer_weights
  loc_names <- items[items %in% construct_cols]

  if (length(loc_names) > 0 && length(loc_names) == length(items)) {
    # This is a HOC -- use LOC construct scores as "indicators"
    # LOC scores are in model$data (NOT model$construct_scores)
    list(
      indicators = loc_names,
      data = model$data[, loc_names, drop = FALSE],
      is_hoc = TRUE,
      needs_reestimation = TRUE
    )
  } else {
    # Standard construct -- use raw indicator data
    list(
      indicators = items,
      data = model$data[, items, drop = FALSE],
      is_hoc = FALSE,
      needs_reestimation = FALSE
    )
  }
}


# =============================================================================
# BORROWING HELPERS (Gudergan et al., 2008, Table 1)
# =============================================================================

#' Get measurement mode ("A" or "B") for a construct.
#' @noRd
get_construct_mode <- function(construct, model) {
  mm_sub <- model$mmMatrix[model$mmMatrix[, "construct"] == construct, , drop = FALSE]
  unique(mm_sub[, "type"])[1]
}

#' Get constructs structurally connected (as source or target) to a construct.
#' @noRd
get_structurally_connected <- function(construct, model) {
  sm <- model$smMatrix
  targets <- unname(sm[sm[, "source"] == construct, "target"])
  sources <- unname(sm[sm[, "target"] == construct, "source"])
  unique(c(targets, sources))
}

#' Find the best donor construct for borrowing indicators.
#'
#' Selects the adjacent construct that maximises testable vanishing tetrads.
#' Rules from Gudergan et al. (2008), Table 1:
#' - Formative focal: no vanishing tetrads regardless -> NULL
#' - 3 own (refl.) + 1 borrowed from reflective: ALL tetrads vanish (score 2)
#' - 3 own (refl.) + 1 borrowed from formative: NONE vanish (score 0)
#' - 2 own (refl.) + 2 borrowed: only tau_1342 vanishes (score 1, any donor mode)
#'
#' @return List with donor info, or NULL if no suitable donor found.
#' @noRd
find_donor <- function(focal, focal_info, model) {
  focal_mode <- get_construct_mode(focal, model)
  if (focal_mode != "A") return(NULL)

  n_own <- length(focal_info$indicators)
  if (n_own < 2 || n_own >= 4) return(NULL)

  n_borrow <- 4L - n_own

  neighbors <- get_structurally_connected(focal, model)
  neighbors <- neighbors[!grepl("*", neighbors, fixed = TRUE)]
  if (length(neighbors) == 0) return(NULL)

  best <- NULL
  best_score <- 0L


  for (donor_name in neighbors) {
    donor_info <- resolve_indicators(donor_name, model)

    # Don't mix HOC (LOC scores) with standard (raw indicators)
    if (focal_info$is_hoc != donor_info$is_hoc) next
    if (length(donor_info$indicators) < n_borrow) next

    donor_mode <- get_construct_mode(donor_name, model)

    if (n_own == 3L) {
      if (donor_mode == "A") {
        score <- 2L; pattern <- "all"
      } else {
        score <- 0L; pattern <- "none"
      }
    } else {
      # n_own == 2: tau_1342 vanishes regardless of donor mode
      score <- 1L; pattern <- "tau_1342"
    }

    if (score > best_score ||
        (score == best_score && !is.null(best) &&
         length(donor_info$indicators) > best$n_donor_indicators)) {
      best_score <- score
      best <- list(
        construct = donor_name,
        mode = donor_mode,
        donor_indicators = donor_info$indicators,
        donor_data = donor_info$data,
        is_hoc = donor_info$is_hoc,
        needs_reestimation = donor_info$needs_reestimation,
        borrowed = donor_info$indicators[seq_len(n_borrow)],
        vanishing_pattern = pattern,
        n_vanishing = score,
        n_donor_indicators = length(donor_info$indicators)
      )
    }
  }

  if (is.null(best) || best$n_vanishing == 0L) return(NULL)
  best
}

#' Enumerate vanishing tetrads for the borrowing case.
#'
#' @param own_indicators Character vector of focal construct's own indicators.
#' @param borrowed_indicators Character vector of borrowed indicator names.
#' @param vanishing_pattern "all" (standard enumeration) or "tau_1342".
#' @return Data frame with columns i, j, k, l, tetrad_id.
#' @noRd
enumerate_borrowed_tetrads <- function(own_indicators, borrowed_indicators,
                                       vanishing_pattern) {
  if (vanishing_pattern == "all") {
    return(enumerate_tetrads(c(own_indicators, borrowed_indicators)))
  }

  # tau_1342: sigma_{ik}*sigma_{jl} - sigma_{il}*sigma_{jk}
  # i, j from own; k, l from borrowed
  own_combos <- utils::combn(own_indicators, 2, simplify = FALSE)
  bor_combos <- utils::combn(borrowed_indicators, 2, simplify = FALSE)

  results <- vector("list", length(own_combos) * length(bor_combos))
  idx <- 1L
  for (oc in own_combos) {
    for (bc in bor_combos) {
      results[[idx]] <- data.frame(
        i = oc[1], j = oc[2], k = bc[1], l = bc[2],
        tetrad_id = 3L, stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, results)
}

# =============================================================================
# MAIN EXPORTED FUNCTION
# =============================================================================

#' Confirmatory Tetrad Analysis for PLS-SEM (CTA-PLS)
#'
#' `assess_cta` tests whether each construct's measurement model is consistent
#' with a reflective (common factor) specification. Under a reflective model,
#' all vanishing tetrads equal zero. If any tetrad is significantly non-zero
#' (after multiple testing correction), the reflective specification is rejected.
#'
#' @param seminr_model A PLS-SEM model estimated by [seminr::estimate_pls()].
#' @param constructs Character vector of construct names to test. If `NULL`
#'   (default), all constructs in the model are tested.
#' @param nboot Number of bootstrap subsamples (default 5000).
#' @param seed Seed for reproducibility (default 123).
#' @param alpha Significance level (default 0.05).
#' @param correction Multiple testing correction method: `"BH"`
#'   (Benjamini-Hochberg, default), `"bonferroni"`, or `"none"`.
#' @param borrow Logical. If `TRUE` (default), constructs with 2--3
#'   indicators borrow indicators from structurally connected constructs
#'   to reach the minimum of 4 needed for tetrad testing (Gudergan et al.,
#'   2008, Table 1). If `FALSE`, constructs with fewer than 4 indicators
#'   are skipped.
#'
#' @return An S3 object of class `cta_analysis` containing:
#' \describe{
#'   \item{construct_results}{Data frame summarising each construct: mode,
#'     number of indicators, tetrads tested, significant tetrads, and verdict.}
#'   \item{tetrad_details}{Named list of per-construct data frames with
#'     individual tetrad estimates, bootstrap CIs, and adjusted p-values.}
#'   \item{nboot}{Number of bootstrap subsamples used.}
#'   \item{alpha}{Significance level used.}
#'   \item{correction}{Multiple testing correction method used.}
#'   \item{skipped}{Character vector of constructs that were skipped.}
#'   \item{borrowing}{Named list of borrowing details for constructs that
#'     borrowed indicators (donor name, mode, and vanishing pattern).}
#' }
#'
#' @details
#' **Minimum indicator requirement:** Without borrowing, constructs with fewer
#' than 4 indicators are skipped. With `borrow = TRUE` (default), constructs
#' with 2--3 indicators borrow from structurally adjacent constructs.
#' A minimum of 2 own indicators is always required.
#'
#' **Borrowing rules** (Gudergan et al., 2008, Table 1): For a reflective
#' focal construct with 3 own indicators, 1 indicator is borrowed from an
#' adjacent reflective construct; all tetrads vanish under H0. For 2 own
#' indicators, 2 are borrowed from any adjacent construct; only the
#' tau_1342 tetrad vanishes. Formative focal constructs cannot be tested
#' via borrowing.
#'
#' For HOC constructs, the lower-order constructs
#' (LOCs) serve as "indicators"; at least 4 LOCs are required.
#'
#' **Interaction constructs** (moderation terms with `*` in the name) are
#' automatically excluded because their measurement specification is determined
#' by the interaction method, not the data.
#'
#' **Multiple testing correction:** The number of tetrads grows combinatorially
#' with indicators. Benjamini-Hochberg (BH) correction is recommended as the
#' default (Cefis et al., 2025). Bonferroni is more conservative. With only
#' 4 indicators (2 tetrads), correction has minimal impact.
#'
#' **Sample size:** CTA-PLS requires adequate sample sizes for reliable results.
#' A warning is issued if N < 200.
#'
#' @seealso [congruence_test()] for congruence coefficient testing
#'
#' @references
#' Gudergan, S. P., Ringle, C. M., Wende, S. & Will, A. (2008). Confirmatory
#' Tetrad Analysis in PLS Path Modeling. *Journal of Business Research*, 61(12),
#' 1238-1249.
#'
#' Cefis, M., Angelelli, M., Carpita, M. & Ciavolino, E. (2025). Confirmatory
#' Tetrad Analysis in PLS-SEM: A Multiple Testing Correction Perspective.
#' *Social Indicators Research* (advance online).
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' # Specify measurement model
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Value",        multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#'
#' # Specify structural model
#' mobi_sm <- relationships(
#'   paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
#'   paths(from = "Expectation", to = c("Value", "Satisfaction")),
#'   paths(from = "Value",       to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#'
#' # Estimate the model
#' pls_model <- estimate_pls(data = mobi, measurement_model = mobi_mm,
#'                            structural_model = mobi_sm)
#'
#' # Run CTA-PLS with borrowing (low nboot for example speed)
#' cta <- assess_cta(pls_model, nboot = 50, seed = 123)
#' print(cta)
#' summary(cta)
#'
#' # Without borrowing -- only constructs with >= 4 indicators are tested
#' cta_no_borrow <- assess_cta(pls_model, nboot = 50, borrow = FALSE)
#'
#' @export
assess_cta <- function(seminr_model,
                       constructs = NULL,
                       nboot = 5000,
                       seed = 123,
                       alpha = 0.05,
                       correction = "BH",
                       borrow = TRUE) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_cta")) {
    return(NULL)
  }

  correction <- match.arg(correction, c("BH", "bonferroni", "none"))

  n_obs <- nrow(seminr_model$data)
  if (n_obs < 200) {
    message("Note: CTA-PLS has limited power with N < 200 (current N = ", n_obs,
            "). Interpret results with caution.")
  }

  # ---------------------------------------------------------------------------
  # Step 2: Determine which constructs to test
  # ---------------------------------------------------------------------------
  all_constructs <- unique(seminr_model$mmMatrix[, "construct"])

  if (is.null(constructs)) {
    constructs <- all_constructs
  } else {
    invalid <- setdiff(constructs, all_constructs)
    if (length(invalid) > 0) {
      warning("Constructs not found in model: ",
              paste(invalid, collapse = ", "), call. = FALSE)
      constructs <- intersect(constructs, all_constructs)
    }
    if (length(constructs) == 0) {
      warning("No valid constructs to test.", call. = FALSE)
      return(NULL)
    }
  }

  # Exclude interaction constructs
  is_intxn <- grepl("*", constructs, fixed = TRUE)
  excluded_interactions <- constructs[is_intxn]
  if (length(excluded_interactions) > 0) {
    message("Excluding interaction constructs (measurement determined by method): ",
            paste(excluded_interactions, collapse = ", "))
    constructs <- constructs[!is_intxn]
  }

  # ---------------------------------------------------------------------------
  # Step 3: Resolve indicators per construct and check minimum count
  # ---------------------------------------------------------------------------
  skipped <- character(0)
  test_constructs <- character(0)
  construct_info <- list()
  borrowing_details <- list()

  for (construct in constructs) {
    info <- resolve_indicators(construct, seminr_model)

    if (length(info$indicators) >= 4) {
      # Standard case: enough indicators
      test_constructs <- c(test_constructs, construct)
      info$borrowing <- NULL
      construct_info[[construct]] <- info
    } else if (borrow && length(info$indicators) >= 2) {
      # Attempt borrowing
      donor <- find_donor(construct, info, seminr_model)
      if (!is.null(donor)) {
        test_constructs <- c(test_constructs, construct)
        borrowed_data <- donor$donor_data[, donor$borrowed, drop = FALSE]
        construct_info[[construct]] <- list(
          indicators = c(info$indicators, donor$borrowed),
          data = cbind(info$data, borrowed_data),
          is_hoc = info$is_hoc || donor$is_hoc,
          needs_reestimation = info$needs_reestimation || donor$needs_reestimation,
          own_indicators = info$indicators,
          borrowed_indicators = donor$borrowed,
          borrowing = list(
            donor = donor$construct,
            donor_mode = donor$mode,
            vanishing_pattern = donor$vanishing_pattern,
            n_vanishing = donor$n_vanishing
          )
        )
        borrowing_details[[construct]] <- construct_info[[construct]]$borrowing
        message("Borrowing ", length(donor$borrowed), " indicator(s) from '",
                donor$construct, "' for '", construct, "'.")
      } else {
        message("Skipping '", construct, "': ", length(info$indicators),
                " indicator(s), no suitable donor found for borrowing.")
        skipped <- c(skipped, construct)
      }
    } else {
      n_ind <- length(info$indicators)
      if (borrow) {
        message("Skipping '", construct, "': only ", n_ind,
                " indicator(s); CTA-PLS requires >= 2 for borrowing.")
      } else {
        message("Skipping '", construct, "': only ", n_ind,
                " indicator(s); CTA-PLS requires >= 4. Set borrow = TRUE ",
                "to test constructs with 2-3 indicators.")
      }
      skipped <- c(skipped, construct)
    }
  }

  if (length(test_constructs) == 0) {
    message("No constructs testable. ",
            if (!borrow) "Set borrow = TRUE to test constructs with 2-3 indicators."
            else "All constructs have too few indicators or no suitable donors.")
    result <- list(
      construct_results = data.frame(
        Construct = character(0), Mode = character(0),
        Indicators = integer(0), Tetrads = integer(0),
        Significant = integer(0), Verdict = character(0),
        stringsAsFactors = FALSE
      ),
      tetrad_details = list(),
      nboot = nboot, alpha = alpha, correction = correction,
      skipped = c(skipped, excluded_interactions),
      borrowing = borrowing_details
    )
    class(result) <- c("cta_analysis", class(result))
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Step 4: Set seed and run CTA-PLS for each construct
  # ---------------------------------------------------------------------------
  set.seed(seed)

  # Pre-build summary data frame
  summary_df <- data.frame(
    Construct = character(0), Mode = character(0),
    Indicators = integer(0), Tetrads = integer(0),
    Significant = integer(0), Verdict = character(0),
    stringsAsFactors = FALSE
  )
  details_list <- list()

  # Determine if any construct needs model re-estimation in bootstrap
  any_hoc <- any(vapply(construct_info, function(x) x$needs_reestimation, logical(1)))

  # ---------------------------------------------------------------------------
  # Bootstrap: resample once, compute tetrads for all constructs
  # ---------------------------------------------------------------------------

  # Pre-compute: enumerate tetrads and compute original values per construct
  original_tetrads <- list()
  tetrad_indices <- list()

  for (construct in test_constructs) {
    info <- construct_info[[construct]]
    if (!is.null(info$borrowing)) {
      tetrad_indices[[construct]] <- enumerate_borrowed_tetrads(
        info$own_indicators, info$borrowed_indicators, info$borrowing$vanishing_pattern
      )
    } else {
      tetrad_indices[[construct]] <- enumerate_tetrads(info$indicators)
    }
    cov_mat <- stats::cov(info$data)
    original_tetrads[[construct]] <- compute_tetrads(cov_mat, tetrad_indices[[construct]])
  }

  # Bootstrap array: store tetrad values per bootstrap iteration
  boot_tetrads <- list()
  for (construct in test_constructs) {
    n_tetrads <- nrow(tetrad_indices[[construct]])
    boot_tetrads[[construct]] <- matrix(NA_real_, nrow = nboot, ncol = n_tetrads)
  }

  for (b in seq_len(nboot)) {
    # Resample rows
    idx <- sample(n_obs, n_obs, replace = TRUE)

    # If any HOC constructs need re-estimation, do it once per bootstrap
    if (any_hoc) {
      resampled_data <- seminr_model$rawdata[idx, ]
      boot_model <- tryCatch(
        suppressMessages(seminr::rerun(seminr_model, data = resampled_data)),
        error = function(e) NULL
      )
    }

    for (construct in test_constructs) {
      info <- construct_info[[construct]]

      if (info$needs_reestimation) {
        # HOC: need LOC construct scores from re-estimated model (in $data)
        if (is.null(boot_model)) {
          # Re-estimation failed; leave NAs
          next
        }
        boot_data <- boot_model$data[, info$indicators, drop = FALSE]
      } else {
        # Standard: just resample the indicator data
        boot_data <- info$data[idx, , drop = FALSE]
      }

      boot_cov <- stats::cov(boot_data)
      boot_tetrads[[construct]][b, ] <- compute_tetrads(boot_cov,
                                                         tetrad_indices[[construct]])
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5: Compute p-values and CIs per construct
  # ---------------------------------------------------------------------------
  alpha_half <- alpha / 2

  for (construct in test_constructs) {
    info <- construct_info[[construct]]
    td <- tetrad_indices[[construct]]
    orig <- original_tetrads[[construct]]
    boot_mat <- boot_tetrads[[construct]]
    n_tetrads <- length(orig)

    # Get measurement mode from mmMatrix
    mm_sub <- seminr_model$mmMatrix[seminr_model$mmMatrix[, "construct"] == construct, ]
    if (is.matrix(mm_sub)) {
      mode_type <- unique(mm_sub[, "type"])[1]
    } else {
      mode_type <- mm_sub["type"]
    }
    mode_label <- if (mode_type == "A") "Mode A (reflective)" else "Mode B (formative)"

    # For HOC: mode is determined by the higher_composite specification
    if (info$is_hoc && is.null(info$borrowing)) {
      hoc_weights <- seminr_model$outer_weights[info$indicators, construct]
      mode_label <- paste0(mode_label, " [HOC]")
    } else if (info$is_hoc) {
      mode_label <- paste0(mode_label, " [HOC]")
    }

    # Add borrowing annotation
    if (!is.null(info$borrowing)) {
      mode_label <- paste0(mode_label, " [borrowed from ", info$borrowing$donor, "]")
    }

    # Compute per-tetrad statistics
    tetrad_labels <- character(n_tetrads)
    t_values <- numeric(n_tetrads)
    boot_means <- numeric(n_tetrads)
    boot_sds <- numeric(n_tetrads)
    ci_lower <- numeric(n_tetrads)
    ci_upper <- numeric(n_tetrads)
    p_values <- numeric(n_tetrads)

    for (t_idx in seq_len(n_tetrads)) {
      tetrad_labels[t_idx] <- format_tetrad_label(
        td$i[t_idx], td$j[t_idx], td$k[t_idx], td$l[t_idx], td$tetrad_id[t_idx]
      )

      boot_vals <- boot_mat[, t_idx]
      boot_vals <- boot_vals[!is.na(boot_vals)]

      if (length(boot_vals) < 10) {
        t_values[t_idx] <- NA_real_
        boot_means[t_idx] <- NA_real_
        boot_sds[t_idx] <- NA_real_
        ci_lower[t_idx] <- NA_real_
        ci_upper[t_idx] <- NA_real_
        p_values[t_idx] <- NA_real_
        next
      }

      boot_means[t_idx] <- mean(boot_vals)
      boot_sds[t_idx] <- sd(boot_vals)

      # Bootstrap t-value: tau / se (Gudergan et al., 2008, p. 1242)
      if (boot_sds[t_idx] < .Machine$double.eps) {
        t_values[t_idx] <- NA_real_
      } else {
        t_values[t_idx] <- orig[t_idx] / boot_sds[t_idx]
      }

      # Percentile confidence interval
      ci_lower[t_idx] <- stats::quantile(boot_vals, probs = alpha_half)
      ci_upper[t_idx] <- stats::quantile(boot_vals, probs = 1 - alpha_half)

      # Two-sided p-value: proportion of bootstrap distribution on opposite
      # side of zero from the observed value (or crossing zero)
      # H0: tau = 0. Reject if CI does not include 0.
      # P-value = 2 * min(proportion >= 0, proportion <= 0)
      p_above <- mean(boot_vals >= 0)
      p_below <- mean(boot_vals <= 0)
      p_values[t_idx] <- 2 * min(p_above, p_below)
    }

    # Apply multiple testing correction
    if (correction == "none" || all(is.na(p_values))) {
      adj_p <- p_values
    } else {
      adj_p <- stats::p.adjust(p_values, method = correction)
    }

    # Significance determination
    significant <- !is.na(adj_p) & adj_p <= alpha

    # Build detail data.frame
    ci_lower_name <- paste0(alpha_half * 100, "% CI")
    ci_upper_name <- paste0((1 - alpha_half) * 100, "% CI")

    detail_df <- data.frame(
      Tetrad = tetrad_labels,
      Estimate = orig,
      T_Value = t_values,
      Boot_Mean = boot_means,
      Boot_SD = boot_sds,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      P_Value = p_values,
      Adj_P = adj_p,
      Significant = significant,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(detail_df)[6:7] <- c(ci_lower_name, ci_upper_name)

    details_list[[construct]] <- detail_df

    # Verdict
    n_sig <- sum(significant, na.rm = TRUE)
    if (n_sig == 0) {
      verdict <- "Reflective supported"
    } else {
      verdict <- "Reflective rejected"
    }

    summary_row <- data.frame(
      Construct = construct,
      Mode = mode_label,
      Indicators = length(info$indicators),
      Tetrads = n_tetrads,
      Significant = n_sig,
      Verdict = verdict,
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, summary_row)
  }

  # ---------------------------------------------------------------------------
  # Step 6: Compile result
  # ---------------------------------------------------------------------------
  result <- list(
    construct_results = summary_df,
    tetrad_details = details_list,
    nboot = nboot,
    alpha = alpha,
    correction = correction,
    skipped = c(skipped, excluded_interactions),
    borrowing = borrowing_details
  )
  class(result) <- c("cta_analysis", class(result))

  result
}


# =============================================================================
# S3 METHODS
# =============================================================================

#' @export
print.cta_analysis <- function(x, ...) {
  cat("Confirmatory Tetrad Analysis (CTA-PLS)\n")
  cat("=======================================\n")
  cat("Bootstrap samples:", x$nboot, "| Alpha:", x$alpha,
      "| Correction:", x$correction, "\n\n")

  if (nrow(x$construct_results) == 0) {
    cat("No constructs tested (all had < 4 indicators).\n")
  } else {
    # Print summary table
    df <- x$construct_results
    col_widths <- pmax(nchar(colnames(df)),
                       apply(df, 2, function(col) max(nchar(as.character(col)))))

    # Header
    header <- paste(mapply(function(name, w) formatC(name, width = w, flag = "-"),
                           colnames(df), col_widths), collapse = "  ")
    cat(header, "\n")
    cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

    # Rows
    for (r in seq_len(nrow(df))) {
      row_str <- paste(mapply(function(val, w) formatC(as.character(val), width = w, flag = "-"),
                              df[r, ], col_widths), collapse = "  ")
      cat(row_str, "\n")
    }
  }

  if (length(x$borrowing) > 0) {
    cat("\nBorrowing:\n")
    for (cname in names(x$borrowing)) {
      b <- x$borrowing[[cname]]
      cat("  ", cname, ": borrowed from ", b$donor,
          " (", b$vanishing_pattern, " pattern, ",
          b$n_vanishing, " vanishing tetrad(s))\n", sep = "")
    }
  }

  if (length(x$skipped) > 0) {
    cat("\nSkipped:", paste(x$skipped, collapse = ", "), "\n")
  }

  invisible(x)
}

#' @export
summary.cta_analysis <- function(object, ...) {
  result <- list(
    construct_results = object$construct_results,
    tetrad_details = object$tetrad_details,
    nboot = object$nboot,
    alpha = object$alpha,
    correction = object$correction,
    skipped = object$skipped,
    borrowing = object$borrowing
  )
  class(result) <- c("summary.cta_analysis", class(result))
  result
}

#' @export
print.summary.cta_analysis <- function(x, digits = 4, ...) {
  cat("Confirmatory Tetrad Analysis (CTA-PLS) \u2014 Detailed Results\n")
  cat("==========================================================\n")
  cat("Bootstrap samples:", x$nboot, "| Alpha:", x$alpha,
      "| Correction:", x$correction, "\n\n")

  for (construct in names(x$tetrad_details)) {
    # Construct header with verdict
    row_idx <- which(x$construct_results$Construct == construct)
    verdict <- x$construct_results$Verdict[row_idx]
    mode <- x$construct_results$Mode[row_idx]

    cat("--- ", construct, " (", mode, ") ---\n", sep = "")
    if (construct %in% names(x$borrowing)) {
      b <- x$borrowing[[construct]]
      cat("Borrowed ", b$n_vanishing, " vanishing tetrad(s) from ", b$donor,
          " (", b$vanishing_pattern, ")\n", sep = "")
    }
    cat("Verdict: ", verdict, "\n\n", sep = "")

    detail <- x$tetrad_details[[construct]]

    # Format numeric columns
    numeric_cols <- c("Estimate", "T_Value", "Boot_Mean", "Boot_SD",
                      colnames(detail)[6], colnames(detail)[7],
                      "P_Value", "Adj_P")
    for (col in numeric_cols) {
      if (col %in% colnames(detail)) {
        detail[[col]] <- formatC(as.numeric(detail[[col]]),
                                 format = "f", digits = digits)
      }
    }
    detail$Significant <- ifelse(detail$Significant, "*", "")

    print(detail, row.names = FALSE, right = FALSE)
    cat("\n")
  }

  if (length(x$skipped) > 0) {
    cat("Skipped constructs: ", paste(x$skipped, collapse = ", "), "\n")
  }

  invisible(x)
}

#' @export
plot.cta_analysis <- function(x, ...) {
  if (nrow(x$construct_results) == 0) {
    message("No constructs were tested \u2014 nothing to plot.")
    return(invisible(x))
  }

  # Collect all adjusted p-values per construct
  constructs <- names(x$tetrad_details)
  n_constructs <- length(constructs)

  # Determine layout
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (n_constructs <= 4) {
    n_cols <- min(n_constructs, 2)
    n_rows <- ceiling(n_constructs / n_cols)
  } else {
    n_cols <- 3
    n_rows <- ceiling(n_constructs / n_cols)
  }
  par(mfrow = c(n_rows, n_cols), mar = c(4, 5, 3, 1))

  for (construct in constructs) {
    detail <- x$tetrad_details[[construct]]
    adj_p <- as.numeric(detail$Adj_P)
    n_tetrads <- length(adj_p)

    row_idx <- which(x$construct_results$Construct == construct)
    verdict <- x$construct_results$Verdict[row_idx]

    # Colour: green if supported, red if rejected
    cols <- ifelse(adj_p <= x$alpha, "firebrick", "steelblue")

    # Dot plot of adjusted p-values
    plot(seq_len(n_tetrads), adj_p,
         pch = 19, col = cols, cex = 1.2,
         xlab = "Tetrad index", ylab = "Adjusted p-value",
         main = construct,
         ylim = c(0, max(1, max(adj_p, na.rm = TRUE))),
         xaxt = "n", ...)
    axis(1, at = seq_len(n_tetrads))
    abline(h = x$alpha, lty = 2, col = "grey40")
    legend("topright", legend = verdict, bty = "n", cex = 0.8,
           text.col = if (grepl("supported", verdict)) "steelblue" else "firebrick")
  }

  invisible(x)
}
