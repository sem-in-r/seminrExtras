# =============================================================================
# feature_micom.R - Measurement Invariance of Composites (MICOM) for two groups
# =============================================================================
# Implements steps 2 and 3 of MICOM (Henseler, Ringle & Sarstedt, 2016) as in
# cSEM::testMICOM() and SmartPLS:
#   Step 2: c = cor(X w1, X w2) on the pooled indicator data X, with each group's
#           weights; one-sided permutation test (p = share of permuted c <= c).
#   Step 3: composite scores from the pooled-data model; difference in means
#           and log ratio of variances; two-sided permutation tests and
#           permutation intervals.
#
# Re-estimation runs on correlation matrices with an internal PLS core that
# reproduces seminr's estimate_pls() weights (Mode A, Mode B; path weighting or
# factorial scheme), so the permutation refits stay fast.
#
# The file also holds INTERNAL, unexported research functions (bottom of the
# file). They are not part of the public interface.
# =============================================================================

# ---- PLS core on a correlation matrix -----------------------------------------

# PLS path model weights from a correlation (or covariance) matrix.
#
# @param S Correlation matrix of all indicators.
# @param blocks List of integer indicator positions, one element per construct.
# @param inner J x J 0/1 matrix, inner[k, j] = 1 if construct k -> construct j.
# @param modes Character vector, "A" or "B" per construct.
# @param scheme "path" (path weighting) or "factorial".
# @return J x P weight matrix; each composite has unit variance.
# @noRd
micom_pls_core <- function(S, blocks, inner, modes, scheme, tol = 1e-10, maxit = 300) {
  J <- length(blocks)
  W <- matrix(0, J, ncol(S))
  for (j in seq_len(J)) {
    b <- blocks[[j]]
    w <- rep(1, length(b))
    W[j, b] <- w / sqrt(drop(t(w) %*% S[b, b, drop = FALSE] %*% w))
  }
  adj <- (inner + t(inner)) > 0
  for (it in seq_len(maxit)) {
    W0 <- W
    C <- W %*% S %*% t(W)
    E <- matrix(0, J, J)
    for (j in seq_len(J)) {
      if (scheme == "path") {
        pred <- which(inner[, j] == 1)
        succ <- which(inner[j, ] == 1)
        if (length(pred)) E[j, pred] <- solve(C[pred, pred, drop = FALSE], C[pred, j])
        if (length(succ)) E[j, succ] <- C[j, succ]
      } else {
        E[j, adj[j, ]] <- C[j, adj[j, ]]
      }
    }
    for (j in seq_len(J)) {
      if (!any(adj[j, ])) next
      b <- blocks[[j]]
      w <- drop(S[b, , drop = FALSE] %*% t(W) %*% E[j, ])
      if (modes[j] == "B" && length(b) > 1) w <- solve(S[b, b, drop = FALSE], w)
      W[j, ] <- 0
      W[j, b] <- w / sqrt(drop(t(w) %*% S[b, b, drop = FALSE] %*% w))
    }
    if (max(abs(W - W0)) < tol) break
  }
  W
}

# Specification of a fitted seminr model in the form the core needs.
# Returns NULL (having warned) when MICOM cannot be applied.
# @noRd
micom_spec <- function(seminr_model, fn_name) {
  if (!validate_seminr_model(seminr_model, fn_name)) return(NULL)
  if (has_higher_order(seminr_model)) {
    warning(fn_name, "() does not support higher-order models.", call. = FALSE)
    return(NULL)
  }
  constructs <- seminr_model$constructs
  if (any(grepl("*", constructs, fixed = TRUE))) {
    warning(fn_name, "() does not support interaction terms: their measurement is fixed by the ",
            "product method, so compositional invariance is not defined for them.", call. = FALSE)
    return(NULL)
  }
  scheme <- if (identical(seminr_model$inner_weights, seminr::path_weighting)) "path"
            else if (identical(seminr_model$inner_weights, seminr::path_factorial)) "factorial"
            else NA_character_
  if (is.na(scheme)) {
    warning(fn_name, "() supports path_weighting and path_factorial inner weights only.", call. = FALSE)
    return(NULL)
  }
  mm <- seminr_model$mmMatrix
  items <- lapply(constructs, function(cn) mm[mm[, "construct"] == cn, "measurement"])
  indicators <- unlist(items)
  blocks <- lapply(items, function(it) match(it, indicators))
  modes <- vapply(constructs, function(cn) {
    ty <- unique(mm[mm[, "construct"] == cn, "type"])
    if (any(ty == "B")) "B" else "A"
  }, character(1))
  sm <- seminr_model$smMatrix
  inner <- matrix(0, length(constructs), length(constructs), dimnames = list(constructs, constructs))
  for (r in seq_len(nrow(sm))) inner[sm[r, "source"], sm[r, "target"]] <- 1
  list(constructs = constructs, indicators = indicators, blocks = blocks, modes = modes,
       inner = inner, scheme = scheme, data = as.matrix(seminr_model$data[, indicators]))
}

# Validated two-level grouping vector, its labels, and the logical index of group 1.
# @noRd
micom_groups <- function(seminr_model, group, n) {
  if (is.character(group) && length(group) == 1) {
    if (!group %in% colnames(seminr_model$rawdata)) {
      stop("Column '", group, "' is not in the model's data.", call. = FALSE)
    }
    group <- seminr_model$rawdata[[group]]
  }
  if (length(group) != n) {
    stop("`group` must have one value per row of the model's data (", n, ").", call. = FALSE)
  }
  if (anyNA(group)) stop("`group` has missing values.", call. = FALSE)
  lev <- if (is.factor(group)) levels(droplevels(group)) else sort(unique(as.character(group)))
  if (length(lev) != 2) {
    stop("MICOM here compares exactly two groups; `group` has ", length(lev), " distinct values.", call. = FALSE)
  }
  list(levels = lev, g1 = as.character(group) == lev[1])
}

#' Measurement invariance of composites (MICOM) for two groups
#'
#' `assess_micom()` runs steps 2 and 3 of the MICOM procedure (Henseler, Ringle &
#' Sarstedt, 2016) for an estimated seminr PLS model and a two-level grouping
#' variable, following the implementation in `cSEM::testMICOM()`.
#'
#' @section Step 1:
#' Configural invariance (same indicators, same model, same algorithm settings in
#' both groups) holds by construction: both groups are estimated with the
#' specification of `seminr_model`. It still has to be argued substantively
#' (identical questionnaires, coding and data treatment).
#'
#' @section Step 2 -- compositional invariance:
#' For each construct, `c` is the correlation between composite scores computed
#' on the pooled indicator data with the first group's weights and with the
#' second group's weights (Henseler et al., 2016, Eq. 4). Group labels are
#' permuted `nperm` times; each time both groups are re-estimated and c is
#' recomputed. `p_value` is the share of permuted c at or below the observed c,
#' and `c_quantile` is the `alpha` quantile of the permuted c. Compositional
#' invariance is supported when `p_value >= alpha` (equivalently, c is not below
#' the quantile).
#'
#' @section Step 3 -- equal composite means and variances:
#' Composite scores come from the model estimated on the pooled data. For each
#' construct the difference in means (first minus second group) and the log
#' ratio of variances are compared with their permutation distributions:
#' two-sided p-values, and the `alpha / 2` and `1 - alpha / 2` permutation
#' quantiles. Equality is supported when the p-value is at least `alpha`.
#'
#' @section Verdict:
#' Following Henseler et al. (2016): `"full"` measurement invariance when step 2
#' and both step-3 tests support invariance; `"partial"` when step 2 does but
#' step 3 does not; `"none"` when step 2 does not.
#'
#' @param seminr_model A PLS model estimated with [seminr::estimate_pls()]
#'   (path weighting or factorial inner scheme; no interaction terms or
#'   higher-order constructs).
#' @param group A vector with exactly two distinct values, one per row of the
#'   model's data, or the name of a column in the data passed to
#'   `estimate_pls()`. The first group is the first factor level (or the first
#'   value in sorted order); differences are first minus second.
#' @param nperm Number of permutations.
#' @param alpha Significance level.
#' @param seed Optional random seed, for reproducible permutations.
#'
#' @return An object of class `micom_analysis`: a list with `step2` and `step3`
#'   data frames (one row per construct), the group labels and sizes, and the
#'   settings used.
#'
#' @references Henseler, J., Ringle, C. M., & Sarstedt, M. (2016). Testing
#' measurement invariance of composites using partial least squares.
#' International Marketing Review, 33(3), 405-431.
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#' mobi_sm <- relationships(
#'   paths(from = c("Image", "Expectation"), to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#' mobi_pls <- estimate_pls(mobi, mobi_mm, mobi_sm)
#'
#' # Illustrative grouping: respondents above vs below the median of CUEX1
#' grp <- ifelse(mobi$CUEX1 > stats::median(mobi$CUEX1), "high", "low")
#' assess_micom(mobi_pls, group = grp, nperm = 99, seed = 123)
#'
#' @export
assess_micom <- function(seminr_model, group, nperm = 1000, alpha = 0.05, seed = NULL) {
  spec <- micom_spec(seminr_model, "assess_micom")
  if (is.null(spec)) return(NULL)
  grp <- micom_groups(seminr_model, group, nrow(spec$data))
  if (!is.null(seed)) set.seed(seed)
  X <- spec$data
  g1 <- grp$g1
  n1 <- sum(g1); N <- nrow(X); J <- length(spec$constructs)
  est <- function(rows) micom_pls_core(stats::cor(X[rows, , drop = FALSE]), spec$blocks, spec$inner, spec$modes, spec$scheme)

  # ---- step 2: c on the pooled indicator data ----
  cstat <- function(idx) {
    Wa <- est(idx); Wb <- est(!idx)
    vapply(seq_len(J), function(j) {
      b <- spec$blocks[[j]]
      if (length(b) < 2) return(1)
      stats::cor(drop(X[, b, drop = FALSE] %*% Wa[j, b]), drop(X[, b, drop = FALSE] %*% Wb[j, b]))
    }, numeric(1))
  }
  c_obs <- cstat(g1)
  c_perm <- replicate(nperm, { idx <- rep(FALSE, N); idx[sample.int(N, n1)] <- TRUE; cstat(idx) })
  if (is.null(dim(c_perm))) c_perm <- matrix(c_perm, nrow = 1)
  p2 <- vapply(seq_len(J), function(j) mean(c_perm[j, ] <= c_obs[j]), numeric(1))
  q2 <- apply(c_perm, 1, stats::quantile, probs = alpha, names = FALSE)

  # ---- step 3: pooled-model composite scores ----
  Wp <- est(rep(TRUE, N))
  Z <- scale(X)
  scores <- vapply(seq_len(J), function(j) { b <- spec$blocks[[j]]; drop(Z[, b, drop = FALSE] %*% Wp[j, b]) }, numeric(N))
  m3 <- function(idx) colMeans(scores[idx, , drop = FALSE]) - colMeans(scores[!idx, , drop = FALSE])
  v3 <- function(idx) log(apply(scores[idx, , drop = FALSE], 2, stats::var)) - log(apply(scores[!idx, , drop = FALSE], 2, stats::var))
  m_obs <- m3(g1); v_obs <- v3(g1)
  perm3 <- replicate(nperm, { idx <- rep(FALSE, N); idx[sample.int(N, n1)] <- TRUE; c(m3(idx), v3(idx)) })
  if (is.null(dim(perm3))) perm3 <- matrix(perm3, ncol = nperm)
  p_two <- function(obs, sims) mean(sims > abs(obs)) + mean(sims < -abs(obs))
  q3 <- function(sims) stats::quantile(sims, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  pm <- vapply(seq_len(J), function(j) p_two(m_obs[j], perm3[j, ]), numeric(1))
  pv <- vapply(seq_len(J), function(j) p_two(v_obs[j], perm3[J + j, ]), numeric(1))
  qm <- t(vapply(seq_len(J), function(j) q3(perm3[j, ]), numeric(2)))
  qv <- t(vapply(seq_len(J), function(j) q3(perm3[J + j, ]), numeric(2)))

  step2 <- data.frame(construct = spec$constructs, c = c_obs, c_quantile = q2, p_value = p2,
                      compositional_invariance = ifelse(p2 >= alpha, "yes", "no"),
                      row.names = NULL, stringsAsFactors = FALSE)
  step3 <- data.frame(construct = spec$constructs,
                      mean_diff = m_obs, mean_lower = qm[, 1], mean_upper = qm[, 2], mean_p = pm,
                      equal_mean = ifelse(pm >= alpha, "yes", "no"),
                      log_var_ratio = v_obs, var_lower = qv[, 1], var_upper = qv[, 2], var_p = pv,
                      equal_var = ifelse(pv >= alpha, "yes", "no"),
                      row.names = NULL, stringsAsFactors = FALSE)
  step3$invariance <- ifelse(step2$compositional_invariance == "no", "none",
                             ifelse(step3$equal_mean == "yes" & step3$equal_var == "yes", "full", "partial"))
  out <- list(step2 = step2, step3 = step3, groups = grp$levels, n = c(n1, N - n1),
              settings = list(nperm = nperm, alpha = alpha, scheme = spec$scheme, seed = seed))
  class(out) <- c("micom_analysis", class(out))
  out
}

#' @param x A `micom_analysis` object.
#' @param digits Number of decimal places to print.
#' @param ... Not used.
#' @rdname assess_micom
#' @export
print.micom_analysis <- function(x, digits = 3, ...) {
  rnd <- function(df) { num <- vapply(df, is.numeric, logical(1)); df[num] <- lapply(df[num], round, digits); df }
  cat("MICOM: measurement invariance of composites (Henseler, Ringle & Sarstedt, 2016)\n")
  cat("==============================================================================\n")
  cat(sprintf("Groups: %s (n = %d) vs %s (n = %d); %d permutations; alpha = %s\n\n",
              x$groups[1], x$n[1], x$groups[2], x$n[2], x$settings$nperm, format(x$settings$alpha)))
  cat("Step 2 -- compositional invariance\n")
  print(rnd(x$step2), row.names = FALSE)
  cat("\nStep 3 -- equal composite means and variances (first minus second group)\n")
  print(rnd(x$step3), row.names = FALSE)
  invisible(x)
}

# =============================================================================
# INTERNAL RESEARCH FUNCTIONS -- NOT EXPORTED
# =============================================================================
# Not part of the public interface; reachable as seminrExtras:::micom_research().
# They implement methods under development (a studentised test of compositional
# invariance and Welch's test for step 3a) and are kept here so the authors can
# use them on seminr models before they are published and documented.

# Unit-length, sign-anchored weight directions for every construct.
# @noRd
micom_directions <- function(X, spec) {
  W <- micom_pls_core(stats::cor(X), spec$blocks, spec$inner, spec$modes, spec$scheme)
  lapply(seq_along(spec$blocks), function(j) {
    w <- W[j, spec$blocks[[j]]]
    if (sum(w) < 0) w <- -w
    w / sqrt(sum(w^2))
  })
}

# Rank-r pseudo-inverse of a symmetric matrix.
# @noRd
micom_pinv <- function(V, r) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  E <- e$vectors[, seq_len(r), drop = FALSE]
  E %*% diag(1 / e$values[seq_len(r)], r) %*% t(E)
}

# Research version of MICOM step 2 / 3a (internal).
#
# Studentised permutation test of each Mode A construct's weight direction in the
# full PLS model (grouped-jackknife covariance, rank-(K - 1) Wald statistic,
# Chung & Romano 2013), the angle between the groups' weight directions (degrees),
# and Welch's t-test for equal composite means on the pooled-model scores.
#
# @param seminr_model,group As in assess_micom().
# @param nperm Permutations for the studentised test.
# @param jackknife_groups Jackknife groups per sample.
# @param seed Optional random seed.
# @return data.frame: construct, mode, items, p_studentized, angle, p_mean_welch.
# @noRd
micom_research <- function(seminr_model, group, nperm = 999, jackknife_groups = 20, seed = NULL) {
  spec <- micom_spec(seminr_model, "micom_research")
  if (is.null(spec)) return(NULL)
  grp <- micom_groups(seminr_model, group, nrow(spec$data))
  if (!is.null(seed)) set.seed(seed)
  X <- spec$data; g1 <- grp$g1; n1 <- sum(g1); N <- nrow(X); J <- length(spec$constructs)
  K <- vapply(spec$blocks, length, integer(1))
  tested <- spec$modes == "A" & K >= 2
  XA <- scale(X[g1, , drop = FALSE]); XB <- scale(X[!g1, , drop = FALSE]); G <- jackknife_groups
  jk <- function(Xg) {
    h <- sample(rep_len(seq_len(G), nrow(Xg)))
    U <- lapply(seq_len(G), function(k) micom_directions(Xg[h != k, , drop = FALSE], spec))
    lapply(seq_len(J), function(j) {
      M <- vapply(U, `[[`, numeric(K[j]), j)
      if (is.null(dim(M))) M <- matrix(M, nrow = 1)
      ((G - 1) / G) * tcrossprod(M - rowMeans(M))
    })
  }
  wstat <- function(a, b) {
    ua <- micom_directions(a, spec); ub <- micom_directions(b, spec); va <- jk(a); vb <- jk(b)
    vapply(seq_len(J), function(j) {
      if (!tested[j]) return(NA_real_)
      d <- ua[[j]] - ub[[j]]
      drop(t(d) %*% micom_pinv(va[[j]] + vb[[j]], K[j] - 1) %*% d)
    }, numeric(1))
  }
  ua <- micom_directions(XA, spec); ub <- micom_directions(XB, spec)
  angle <- vapply(seq_len(J), function(j) {
    if (K[j] < 2) return(NA_real_)
    acos(min(1, max(-1, sum(ua[[j]] * ub[[j]])))) * 180 / pi
  }, numeric(1))
  p_stud <- rep(NA_real_, J)
  if (any(tested)) {
    w_obs <- wstat(XA, XB); P <- rbind(XA, XB)
    w_perm <- replicate(nperm, { i <- sample.int(N, n1); wstat(P[i, , drop = FALSE], P[-i, , drop = FALSE]) })
    if (is.null(dim(w_perm))) w_perm <- matrix(w_perm, nrow = 1)
    p_stud[tested] <- vapply(which(tested), function(j) (1 + sum(w_perm[j, ] >= w_obs[j])) / (nperm + 1), numeric(1))
  }
  Wp <- micom_pls_core(stats::cor(X), spec$blocks, spec$inner, spec$modes, spec$scheme); Z <- scale(X)
  p_welch <- vapply(seq_len(J), function(j) { b <- spec$blocks[[j]]; s <- drop(Z[, b, drop = FALSE] %*% Wp[j, b])
    stats::t.test(s[g1], s[!g1])$p.value }, numeric(1))
  data.frame(construct = spec$constructs, mode = spec$modes, items = K, p_studentized = p_stud,
             angle = angle, p_mean_welch = p_welch, row.names = NULL, stringsAsFactors = FALSE)
}
