# =============================================================================
# feature_coa.R - Composite Overfit Analysis (COA)
# =============================================================================
# This file implements the COA framework for PLS-SEM models, as described in:
#
# - Ray, S., Danks, N. P., & Valdez, A. C. (2022). "A framework for
#   Composite Overfit Analysis." Information Systems Research.
#
# COA identifies cases whose predictions are unstable (overfit) and groups
# them to reveal which model parameters are affected. The analysis proceeds
# in three steps:
# 1. Compute Predictive Deviance (PD) via cross-validation
# 2. Build a deviance tree to identify groups of deviant cases
# 3. Assess parameter instability by removing deviant groups
# =============================================================================

# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

#' Composite Overfit Analysis (COA)
#'
#' `assess_coa` conducts a full Composite Overfit Analysis on a PLS-SEM model.
#' It computes predictive deviance for a focal construct, identifies groups of
#' deviant cases using a decision tree, and assesses parameter instability.
#'
#' @param seminr_model An estimated SEMinR model from `estimate_pls()`.
#' @param focal_construct Name of the endogenous construct to analyze (character).
#' @param deviance_bounds Two-element numeric vector of quantile bounds for
#'   identifying deviant cases (default `c(0.025, 0.975)`).
#' @param params Character vector of model parameters to assess for instability.
#'   Valid values include `"path_coef"`, `"outer_weights"`, `"outer_loadings"`,
#'   `"rSquared"` (default `"path_coef"`).
#' @param technique Prediction technique: `predict_DA` (default) or `predict_EA`.
#' @param noFolds Number of folds for k-fold cross-validation (default 10).
#' @param reps Number of CV repetitions (default 1).
#' @param cores Number of cores for parallel CV (default NULL = sequential).
#' @param seed Random seed for reproducibility (default 123).
#' @param predict_model Optional pre-computed `predict_pls_model` object. If
#'   provided, skips the cross-validation step (saves time when predictions
#'   have already been computed).
#'
#' @return An object of class `coa_analysis` containing:
#'   \item{pls_model}{The original estimated model}
#'   \item{focal_construct}{Name of the analyzed construct}
#'   \item{deviance_bounds}{Quantile bounds used}
#'   \item{predictive_deviance}{Predictive deviance results (class `coa_deviance`)}
#'   \item{deviance_tree}{Deviance tree and deviant groups (class `coa_dtree`)}
#'   \item{unstable}{Parameter instability results (class `coa_unstable`)}
#'
#' @seealso [predictive_deviance()], [deviance_tree()], [unstable_params()],
#'   [group_rules()], [competes()]
#'
#' @examples
#' \donttest{
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",  multi_items("IMAG", 1:5)),
#'   composite("Value",  multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty", multi_items("CUSL", 1:3))
#' )
#'
#' mobi_sm <- relationships(
#'   paths(from = c("Image", "Value"), to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#'
#' mobi_pls <- estimate_pls(
#'   data = mobi,
#'   measurement_model = mobi_mm,
#'   structural_model  = mobi_sm
#' )
#'
#' coa_result <- assess_coa(mobi_pls,
#'                           focal_construct = "Loyalty",
#'                           noFolds = 10,
#'                           cores = 1,
#'                           seed = 123)
#' print(coa_result)
#' summary(coa_result)
#' plot(coa_result, type = "pd")
#' }
#'
#' @export
assess_coa <- function(seminr_model,
                        focal_construct,
                        deviance_bounds = c(0.025, 0.975),
                        params = "path_coef",
                        technique = predict_DA,
                        noFolds = 10,
                        reps = 1,
                        cores = NULL,
                        seed = 123,
                        predict_model = NULL) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_for_prediction(seminr_model, "assess_coa")) {
    return(NULL)
  }

  construct_names <- colnames(seminr_model$construct_scores)
  if (!(focal_construct %in% construct_names)) {
    stop("focal_construct '", focal_construct,
         "' not found in model constructs: ",
         paste(construct_names, collapse = ", "),
         call. = FALSE)
  }

  # Validate deviance_bounds: must be two probabilities with lower < upper
  if (length(deviance_bounds) != 2 ||
      any(deviance_bounds < 0) || any(deviance_bounds > 1) ||
      deviance_bounds[1] >= deviance_bounds[2]) {
    stop("deviance_bounds must be two values in [0, 1] with first < second, e.g. c(0.025, 0.975)",
         call. = FALSE)
  }

  # Validate params against known model slots
  valid_params <- c("path_coef", "outer_weights", "outer_loadings", "rSquared")
  invalid_params <- setdiff(params, valid_params)
  if (length(invalid_params) > 0) {
    stop("Invalid params: ", paste(invalid_params, collapse = ", "),
         ". Valid options: ", paste(valid_params, collapse = ", "),
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Compute predictive deviance via cross-validated predictions
  # ---------------------------------------------------------------------------
  pd <- predictive_deviance(seminr_model, focal_construct,
                            technique = technique,
                            noFolds = noFolds, reps = reps,
                            cores = cores, seed = seed,
                            predict_model = predict_model)

  # ---------------------------------------------------------------------------
  # Step 3: Build deviance tree and identify deviant case groups
  # ---------------------------------------------------------------------------
  dt <- deviance_tree(pd, deviance_bounds = deviance_bounds)

  # ---------------------------------------------------------------------------
  # Step 4: Assess parameter instability by removing deviant groups
  # ---------------------------------------------------------------------------
  unstable <- unstable_params(seminr_model,
                              deviant_groups = dt$deviant_groups,
                              params = params)

  analysis <- list(
    pls_model = seminr_model,
    focal_construct = focal_construct,
    deviance_bounds = deviance_bounds,
    predictive_deviance = pd,
    deviance_tree = dt,
    unstable = unstable
  )

  class(analysis) <- c("coa_analysis", class(analysis))
  analysis
}

# =============================================================================
# STEP 1: PREDICTIVE DEVIANCE
# =============================================================================

#' Compute Predictive Deviance
#'
#' Computes per-case predictive deviance (PD) for a focal construct using
#' k-fold cross-validated predictions from `predict_pls()`.
#'
#' PD = in-sample fitted score - out-of-sample predicted score.
#' High PD indicates a case that is well-fitted in-sample but poorly predicted
#' out-of-sample, suggesting overfit.
#'
#' @param seminr_model An estimated SEMinR model.
#' @param focal_construct Name of the construct to analyze.
#' @param technique Prediction technique (default `predict_DA`).
#' @param noFolds Number of CV folds (default 10).
#' @param reps Number of CV repetitions (default 1).
#' @param cores Number of cores for parallel CV (default NULL).
#' @param seed Random seed (default 123).
#' @param predict_model Optional pre-computed `predict_pls_model` object.
#'
#' @return An object of class `coa_deviance` containing PD values, metrics,
#'   and the data frame used for tree construction.
#'
#' @seealso [assess_coa()], [deviance_tree()]
#'
#' @export
predictive_deviance <- function(seminr_model,
                                focal_construct,
                                technique = predict_DA,
                                noFolds = 10,
                                reps = 1,
                                cores = NULL,
                                seed = 123,
                                predict_model = NULL) {

  if (is.null(predict_model)) {
    if (!validate_for_prediction(seminr_model, "predictive_deviance")) {
      return(NULL)
    }
    set.seed(seed)
    predict_model <- predict_pls(seminr_model,
                                 technique = technique,
                                 noFolds = noFolds,
                                 reps = reps,
                                 cores = cores)
  }

  fitted <- predict_model$composites$composite_in_sample[, focal_construct]
  predicted <- predict_model$composites$composite_out_of_sample[, focal_construct]
  actual_star <- predict_model$composites$actuals_star[, focal_construct]

  IS_MSE <- mean((actual_star - fitted)^2)
  OOS_MSE <- mean((actual_star - predicted)^2)
  overfit_ratio <- (OOS_MSE - IS_MSE) / IS_MSE

  PD <- fitted - predicted
  pd_data <- cbind(as.data.frame(seminr_model$construct_scores), PD = PD)

  result <- list(
    PD = PD,
    pd_data = pd_data,
    IS_MSE = IS_MSE,
    OOS_MSE = OOS_MSE,
    overfit_ratio = overfit_ratio,
    fitted_score = fitted,
    predicted_score = predicted
  )

  class(result) <- c("coa_deviance", class(result))
  result
}

# =============================================================================
# STEP 2: DEVIANCE TREE
# =============================================================================

#' Build Deviance Tree
#'
#' Grows a decision tree on predictive deviance scores to identify groups
#' of cases with extreme (deviant) PD values.
#'
#' @param pd_result A `coa_deviance` object from `predictive_deviance()`.
#' @param deviance_bounds Two-element numeric vector of quantile bounds
#'   (default `c(0.025, 0.975)`).
#'
#' @return An object of class `coa_dtree` containing the rpart tree,
#'   deviant groups, group roots, unique deviants, and sorted PD values.
#'
#' @seealso [assess_coa()], [predictive_deviance()], [group_rules()], [competes()]
#'
#' @export
deviance_tree <- function(pd_result, deviance_bounds = c(0.025, 0.975)) {

  # Grow a full (unpruned) regression tree on PD scores:
  # cp = 0 disables complexity pruning to capture all splits
  # minsplit = 2 prevents trivial single-case leaf nodes
  tree <- rpart(PD ~ ., data = pd_result$pd_data, minsplit = 2, cp = 0)

  dev_interval <- stats::quantile(pd_result$PD, probs = deviance_bounds)

  nodes <- extract_nodes(tree$frame, dev_interval)

  sorted_PD <- sort(nodes$leaves$yval, decreasing = TRUE)

  deviants <- cases(tree, nodes$is_deviant_leaf)

  deviant_groups <- lapply(nodes$dev_parent_leaves,
                           function(group) cases(tree, nodes$names %in% group))

  if (length(deviant_groups) > 0) {
    if (length(deviant_groups) > 26) {
      warning("More than 26 deviant groups found (", length(deviant_groups),
              "); only the first 26 are labeled A-Z.", call. = FALSE)
      deviant_groups <- deviant_groups[1:26]
    }
    group_roots <- as.integer(names(deviant_groups))
    names(group_roots) <- names(deviant_groups) <- LETTERS[seq_along(deviant_groups)]
  } else {
    group_roots <- integer(0)
  }

  unique_deviants <- setdiff(deviants, unlist(deviant_groups))

  dtree <- list(
    tree = tree,
    sorted_PD = sorted_PD,
    deviant_groups = deviant_groups,
    group_roots = group_roots,
    unique_deviants = unique_deviants,
    deviant_nodes = nodes$deviants
  )

  class(dtree) <- c("coa_dtree", class(dtree))
  dtree
}

# =============================================================================
# STEP 3: UNSTABLE PARAMETERS
# =============================================================================

#' Assess Parameter Instability
#'
#' Identifies unstable model parameters by removing each deviant group
#' and re-estimating the model, then computing the difference in specified
#' parameters.
#'
#' @param seminr_model An estimated SEMinR model.
#' @param deviant_groups Named list of deviant case index vectors
#'   (from `deviance_tree()$deviant_groups`).
#' @param params Character vector of model parameters to diff
#'   (default `"path_coef"`).
#'
#' @return An object of class `coa_unstable`: a named list (one entry per group)
#'   each containing `cases` and `param_diffs`.
#'
#' @seealso [assess_coa()], [deviance_tree()]
#'
#' @export
unstable_params <- function(seminr_model,
                            deviant_groups,
                            params = "path_coef") {

  if (length(deviant_groups) == 0) {
    result <- list()
    class(result) <- c("coa_unstable", class(result))
    return(result)
  }

  result <- lapply(deviant_groups, function(group_cases) {
    diffs <- param_diffs(group_cases, seminr_model, params)
    list(cases = group_cases, param_diffs = diffs)
  })

  class(result) <- c("coa_unstable", class(result))
  result
}

# =============================================================================
# TREE EXTRACTION HELPERS
# =============================================================================

#' Trace the path from root (node 1) to a given node in an rpart binary tree.
#'
#' In rpart's node numbering, a node's parent is node_id %/% 2:
#' even nodes (2k) are left children, odd nodes (2k+1) are right children.
#' This function recursively walks up to the root and returns the full path.
#' @noRd
path_to <- function(node_id) {
  if (node_id[1] != 1)
    c(Recall(if (node_id %% 2 == 0L) node_id / 2 else (node_id - 1) / 2), node_id)
  else node_id
}

#' @noRd
cases <- function(tree, logical_frame) {
  which(tree$where %in% which(logical_frame))
}

#' @noRd
main_ancestors <- function(parent_ids) {
  ids <- as.integer(parent_ids)
  ancestor_ids <- unique(vapply(ids, function(id) {
    min(ids[ids %in% path_to(id)])
  }, integer(1)))
  as.character(ancestor_ids)
}

#' Extract deviant node information from an rpart tree frame.
#'
#' Classifies each node as leaf vs internal, deviant vs non-deviant (based on
#' whether its mean PD falls outside dev_interval), then:
#' 1. Identifies deviant leaves (individual terminal nodes with extreme PD)
#' 2. Identifies deviant parents (internal nodes whose subtrees are deviant)
#' 3. Finds the highest ancestor of each deviant parent to form groups
#' 4. Maps each group ancestor to its descendant leaf nodes
#' @noRd
extract_nodes <- function(frame, dev_interval) {
  # Classify every node in the tree

  is_leaf <- frame$var == "<leaf>"
  is_deviant <- frame$yval < dev_interval[1] | frame$yval > dev_interval[2]
  is_deviant_leaf <- is_deviant & is_leaf
  is_deviant_parent <- is_deviant & !is_leaf

  names <- row.names(frame)
  leaves <- frame[is_leaf, ]
  deviants <- frame[is_deviant, ]

  leaf_ids <- row.names(leaves)
  dev_parent_ids <- row.names(frame[is_deviant_parent, ])

  # Find the highest (most ancestral) deviant parent for each group,
  # then collect all leaves under each ancestor to form deviant groups
  if (length(dev_parent_ids) == 0) {
    dev_ancestor_ids <- character(0)
    dev_parent_leaves <- list()
  } else {
    dev_ancestor_ids <- main_ancestors(dev_parent_ids)
    dev_parent_leaves <- leaves_from_nodes(dev_ancestor_ids, leaf_ids)
  }

  list(
    is_leaf = is_leaf,
    is_deviant_leaf = is_deviant_leaf,
    names = names,
    leaves = leaves,
    deviants = deviants,
    dev_ancestor_ids = dev_ancestor_ids,
    dev_parent_leaves = dev_parent_leaves
  )
}

#' Find all leaf nodes that descend from each parent node.
#'
#' For each parent, traces root-to-leaf paths for every leaf and checks
#' whether the parent appears on that path (i.e., the leaf is a descendant).
#' @noRd
leaves_from_nodes <- function(parent_ids, leaf_ids) {
  # Pre-compute root-to-leaf paths for all leaves
  node_paths <- lapply(as.integer(leaf_ids), path_to)
  names(node_paths) <- leaf_ids

  # For each parent, find leaves whose path passes through it
  paths_list <- lapply(parent_ids, function(node_id) {
    nid <- as.integer(node_id)
    matching <- which(vapply(node_paths, function(p) nid %in% p, logical(1)))
    names(matching)
  })
  names(paths_list) <- parent_ids
  paths_list
}

#' Re-estimate model without specified cases and compute parameter differences.
#' @noRd
param_diffs <- function(remove_cases, pls_model, params = "path_coef") {
  reduced_data <- pls_model$data[-remove_cases, ]
  # suppressMessages hides the "Generating the seminr model" console output
  reduced_model <- suppressMessages(
    seminr::rerun(pls_model, data = reduced_data)
  )

  diffs <- lapply(params, function(param) {
    reduced_model[[param]] - pls_model[[param]]
  })
  names(diffs) <- params
  diffs
}

# =============================================================================
# RULES EXTRACTION
# =============================================================================

#' Extract Split Rules for a Deviant Group
#'
#' Returns the decision tree split criteria that define a deviant group.
#'
#' @param group_name Single uppercase letter identifying the group (e.g., "A").
#' @param coa_result A `coa_analysis` object from `assess_coa()`.
#'
#' @return A data frame with columns `construct`, `gte`, and `lt` describing
#'   the range of construct scores that define the group.
#'
#' @seealso [assess_coa()], [competes()]
#'
#' @export
group_rules <- function(group_name, coa_result) {
  dtree <- coa_result$deviance_tree
  group_root <- dtree$group_roots[[group_name]]
  node_path <- path_to(group_root)[-1]

  tree_info <- tree_split_index(dtree$tree)

  splits <- do.call(rbind, lapply(node_path, function(nid) {
    get_split_at_node(nid, dtree$tree, tree_info)
  }))

  splits_df <- data.frame(
    var   = splits$var,
    sign  = ifelse(splits$ncat > 0, ">=", "< "),
    value = splits$index,
    stringsAsFactors = FALSE
  )

  consolidate_all_rules(splits_df)
}

#' Report Competing Split Criteria
#'
#' Shows all competing split variables at a given tree node, ranked by
#' improvement. Useful for understanding which constructs nearly split
#' instead of the chosen variable.
#'
#' @param node_id Integer node ID from the rpart tree.
#' @param dtree A `coa_dtree` object (from `coa_result$deviance_tree`).
#'
#' @return A data frame with columns `criterion`, `sign`, `value`, and `improve`.
#'
#' @seealso [group_rules()], [assess_coa()]
#'
#' @export
competes <- function(node_id, dtree) {
  if (node_id == 1) stop("No splits before root (node 1) of tree")

  tree <- dtree$tree
  tree_info <- tree_split_index(tree)
  is_odd <- node_id %% 2 == 1

  search_node <- if (is_odd) node_id - 1 else node_id
  frame_row <- match(search_node / 2, tree_info$node_ids)

  start <- tree_info$split_index[frame_row]
  end <- start + tree$frame$ncompete[frame_row]
  splits <- as.data.frame(tree$splits[start:end, ])

  if (is_odd) splits$ncat <- splits$ncat * -1

  data.frame(
    criterion = row.names(splits),
    sign      = ifelse(splits$ncat > 0, ">=", "< "),
    value     = splits$index,
    improve   = splits$improve,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Pre-compute tree frame metadata for split lookups
#' @noRd
tree_split_index <- function(tree) {
  frame <- tree$frame
  is_leaf_vec <- frame$var == "<leaf>"
  list(
    split_index = cumsum(c(1, frame$ncompete + frame$nsurrogate + !is_leaf_vec)),
    node_ids = as.integer(row.names(frame))
  )
}

#' Get the primary split at a given node
#' @noRd
get_split_at_node <- function(node_id, tree, tree_info) {
  is_odd <- node_id %% 2 == 1
  search_node <- if (is_odd) node_id - 1 else node_id
  frame_row <- match(search_node, tree_info$node_ids)
  split_row <- tree_info$split_index[frame_row - 1]

  split <- data.frame(
    var = rownames(tree$splits)[split_row],
    ncat = tree$splits[split_row, "ncat"],
    index = tree$splits[split_row, "index"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (is_odd) split$ncat <- split$ncat * -1
  split
}

#' @noRd
consolidate_all_rules <- function(splits_df) {
  vars <- unique(splits_df$var)
  criteria <- as.data.frame(do.call(rbind, lapply(vars, consolidate_rule, splits_df)))
  rownames(criteria) <- NULL
  criteria
}

#' @noRd
consolidate_rule <- function(construct, splits_df) {
  split_rows <- splits_df[splits_df$var == construct, ]

  gt_rule <- split_rows[split_rows$sign == ">=", ]
  if (nrow(gt_rule) > 1) gt_rule <- gt_rule[gt_rule$value == max(gt_rule$value), ]
  lt_rule <- split_rows[split_rows$sign == "< ", ]
  if (nrow(lt_rule) > 1) lt_rule <- lt_rule[lt_rule$value == min(lt_rule$value), ]

  data.frame(
    construct = construct,
    gte = if (nrow(gt_rule) == 0) NA_real_ else gt_rule$value,
    lt  = if (nrow(lt_rule) == 0) NA_real_ else lt_rule$value,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' @export
print.coa_analysis <- function(x, ...) {
  cat("Composite Overfit Analysis (COA)\n")
  cat("================================\n")
  cat("Focal construct:", x$focal_construct, "\n")
  cat("Deviance bounds:", x$deviance_bounds[1], "-", x$deviance_bounds[2], "\n")
  cat("Observations:", nrow(x$pls_model$data), "\n\n")

  pd <- x$predictive_deviance
  cat("Prediction Metrics:\n")
  cat("  In-sample MSE: ", format(pd$IS_MSE, digits = 4), "\n")
  cat("  Out-of-sample MSE:", format(pd$OOS_MSE, digits = 4), "\n")
  cat("  Overfit ratio:    ", format(pd$overfit_ratio, digits = 4), "\n\n")

  dt <- x$deviance_tree
  n_groups <- length(dt$deviant_groups)
  n_unique <- length(dt$unique_deviants)
  cat("Deviant Cases:\n")
  cat("  Groups:", n_groups)
  if (n_groups > 0) {
    group_sizes <- vapply(dt$deviant_groups, length, integer(1))
    cat(" (", paste(names(dt$deviant_groups), "=", group_sizes, collapse = ", "), ")")
  }
  cat("\n")
  cat("  Unique deviants:", n_unique, "\n")

  invisible(x)
}

#' @export
summary.coa_analysis <- function(object, ...) {
  pd <- object$predictive_deviance
  dt <- object$deviance_tree

  metrics <- data.frame(
    Metric = c("IS_MSE", "OOS_MSE", "Overfit_Ratio"),
    Value = c(pd$IS_MSE, pd$OOS_MSE, pd$overfit_ratio)
  )

  group_summary <- NULL
  if (length(dt$deviant_groups) > 0) {
    group_summary <- data.frame(
      Group = names(dt$deviant_groups),
      Size = vapply(dt$deviant_groups, length, integer(1)),
      Root_Node = dt$group_roots,
      row.names = NULL
    )
  }

  unstable_summary <- NULL
  if (length(object$unstable) > 0) {
    max_diffs <- vapply(object$unstable, function(g) {
      if ("path_coef" %in% names(g$param_diffs)) {
        max(abs(g$param_diffs$path_coef), na.rm = TRUE)
      } else {
        NA_real_
      }
    }, numeric(1))
    unstable_summary <- data.frame(
      Group = names(object$unstable),
      Max_Path_Diff = max_diffs,
      row.names = NULL
    )
  }

  result <- list(
    focal_construct = object$focal_construct,
    deviance_bounds = object$deviance_bounds,
    n_observations = nrow(object$pls_model$data),
    metrics = metrics,
    group_summary = group_summary,
    unstable_summary = unstable_summary,
    n_unique_deviants = length(dt$unique_deviants)
  )

  class(result) <- c("summary.coa_analysis", class(result))
  result
}

#' @export
print.summary.coa_analysis <- function(x, ...) {
  cat("Composite Overfit Analysis (COA) Summary\n")
  cat("=========================================\n")
  cat("Focal construct:", x$focal_construct, "\n")
  cat("Observations:", x$n_observations, "\n")
  cat("Deviance bounds:", x$deviance_bounds[1], "-", x$deviance_bounds[2], "\n\n")

  cat("Prediction Metrics:\n")
  print(x$metrics, row.names = FALSE)
  cat("\n")

  if (!is.null(x$group_summary)) {
    cat("Deviant Groups:\n")
    print(x$group_summary, row.names = FALSE)
    cat("\n")
  }

  cat("Unique deviants:", x$n_unique_deviants, "\n\n")

  if (!is.null(x$unstable_summary)) {
    cat("Parameter Instability (max |path_coef diff|):\n")
    print(x$unstable_summary, row.names = FALSE)
  }

  invisible(x)
}

# =============================================================================
# PLOT METHODS
# =============================================================================

#' Plot COA Results
#'
#' @param x A `coa_analysis` object from `assess_coa()`.
#' @param type One of `"pd"` (predictive deviance scatter), `"groups"`
#'   (average construct scores per group), or `"tree"` (rpart tree diagram).
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @export
plot.coa_analysis <- function(x, type = c("pd", "groups", "tree"), ...) {
  type <- match.arg(type)

  switch(type,
    pd     = plot_pd(x, ...),
    groups = plot_group_scores(x, ...),
    tree   = plot_tree(x, ...)
  )
}

#' @noRd
plot_pd <- function(coa, ...) {
  pd <- coa$predictive_deviance
  dt <- coa$deviance_tree
  n <- length(pd$PD)

  # Assign group membership to each case
  group_id <- rep(0L, n)
  for (i in seq_along(dt$deviant_groups)) {
    group_id[dt$deviant_groups[[i]]] <- i + 1L
  }
  group_id[dt$unique_deviants] <- 1L

  pd_df <- data.frame(pd = pd$PD, group = group_id)
  pd_df <- pd_df[order(pd_df$pd), ]
  pd_df$order <- seq_len(n)

  dev_quantiles <- stats::quantile(pd$PD, probs = coa$deviance_bounds)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  plot(NULL, NULL,
       xlim = c(1, n), ylim = range(pd$PD),
       frame.plot = FALSE,
       xlab = "Cases", ylab = "Predictive Deviance",
       main = paste("Predictive Deviance:", coa$focal_construct),
       ...)
  abline(h = dev_quantiles, col = "darkgray", lty = 2)

  non_dev <- pd_df$group == 0
  points(pd_df$order[non_dev], pd_df$pd[non_dev], pch = 19, col = "lightgray", cex = 0.6)

  if (any(!non_dev)) {
    n_group_colors <- max(0, max(pd_df$group) - 1)
    cols <- if (n_group_colors > 0) {
      c("black", palette.colors(n = n_group_colors, palette = "Set1"))
    } else {
      "black"
    }
    is_unique <- pd_df$group == 1
    is_group <- !non_dev & !is_unique

    if (any(is_unique)) {
      points(pd_df$order[is_unique], pd_df$pd[is_unique], pch = 8, col = "black", cex = 0.8)
    }
    if (any(is_group)) {
      gcols <- cols[pd_df$group[is_group]]
      points(pd_df$order[is_group], pd_df$pd[is_group], pch = 19, col = gcols, cex = 0.8)

      group_labels <- LETTERS[pd_df$group[is_group] - 1]
      text(pd_df$order[is_group], pd_df$pd[is_group],
           labels = group_labels, pos = 3, cex = 0.6, col = gcols)
    }

    # Legend
    group_names <- names(dt$deviant_groups)
    if (length(group_names) > 0) {
      legend_labels <- group_names
      legend_cols <- cols[2:(length(group_names) + 1)]
      legend_pch <- rep(19, length(group_names))
      if (length(dt$unique_deviants) > 0) {
        legend_labels <- c(legend_labels, "Unique")
        legend_cols <- c(legend_cols, "black")
        legend_pch <- c(legend_pch, 8)
      }
      legend("topleft", legend = legend_labels, col = legend_cols,
             pch = legend_pch, cex = 0.7, bty = "n")
    }
  }
}

#' @noRd
plot_group_scores <- function(coa, remove = NULL, ...) {
  dt <- coa$deviance_tree
  if (length(dt$deviant_groups) == 0) {
    message("No deviant groups to plot.")
    return(invisible(NULL))
  }

  group_scores <- vapply(names(dt$deviant_groups), function(gname) {
    scores <- coa$pls_model$construct_scores[dt$deviant_groups[[gname]], , drop = FALSE]
    colMeans(scores)
  }, numeric(ncol(coa$pls_model$construct_scores)))

  if (!is.null(remove)) {
    group_scores <- group_scores[!rownames(group_scores) %in% remove, , drop = FALSE]
  }

  # Move focal construct to end
  constructs <- rownames(group_scores)
  outcome_idx <- which(constructs == coa$focal_construct)
  if (length(outcome_idx) == 1 && outcome_idx < length(constructs)) {
    reorder <- c(seq_along(constructs)[-outcome_idx], outcome_idx)
    group_scores <- group_scores[reorder, , drop = FALSE]
  }

  num_constructs <- nrow(group_scores)
  ylim_range <- max(3.5, max(abs(group_scores), na.rm = TRUE) + 0.5)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  cols <- palette.colors(n = ncol(group_scores), palette = "Set1")

  matplot(group_scores, type = "b", pch = colnames(group_scores),
          ylim = c(-ylim_range, ylim_range),
          lty = "solid", col = cols,
          xaxt = "n", bty = "n",
          ylab = "Average Construct Score",
          main = "Deviant Group Construct Scores", ...)

  # 50% zone shading
  polygon(x = c(0, num_constructs + 1, num_constructs + 1, 0),
          y = c(1.15, 1.15, -1.15, -1.15),
          col = adjustcolor("gray", alpha.f = 0.1), border = NA)
  polygon(x = c(0, num_constructs + 1, num_constructs + 1, 0),
          y = c(0.67, 0.67, -0.67, -0.67),
          col = adjustcolor("gray", alpha.f = 0.15), border = NA)

  # Separator before outcome
  if (length(outcome_idx) == 1) {
    abline(v = num_constructs - 0.5, col = "gray80", lwd = 2)
  }

  # Redraw lines on top of shading
  matlines(group_scores, type = "b", pch = colnames(group_scores),
           lty = "solid", col = cols)

  construct_labels <- rownames(group_scores)
  font_vec <- rep(1, num_constructs)
  if (length(outcome_idx) == 1) font_vec[num_constructs] <- 2
  axis(1, at = seq_len(num_constructs), labels = construct_labels, font = font_vec)

  legend("topright", legend = colnames(group_scores), col = cols,
         lty = 1, pch = colnames(group_scores), cex = 0.7, bty = "n")
}

#' @noRd
plot_tree <- function(coa, ...) {
  tree <- coa$deviance_tree$tree
  plot(tree, uniform = TRUE, main = paste("Deviance Tree:", coa$focal_construct), ...)
  text(tree, use.n = TRUE, cex = 0.7)
}
