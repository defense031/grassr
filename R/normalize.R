# Input normalization: coerce user-supplied data into two aligned, 0/1 integer
# rater vectors regardless of shape (matrix, wide, long, paired).

# Rules for picking the positive (=1) level from a 2-level categorical variable.
# Order of preference:
#   1. user-supplied `positive` argument
#   2. "yes" (case-insensitive)
#   3. one of "1", "true", "positive", "case", "pos", "present" (case-insensitive)
#   4. first-encountered level in the data
positive_rules <- c("1", "true", "positive", "case", "pos", "present")

coerce_binary <- function(x, positive = NULL, name = "rating") {
  if (is.logical(x)) {
    return(list(values = as.integer(x), positive_level = "TRUE"))
  }

  if (is.numeric(x)) {
    vals <- x[!is.na(x)]
    if (!all(vals %in% c(0, 1))) {
      stop("Column '", name, "' is numeric but contains values other than 0/1. ",
           "Coerce to factor or character first, or pass only binary integers.",
           call. = FALSE)
    }
    return(list(values = as.integer(x), positive_level = "1"))
  }

  if (is.factor(x)) {
    x <- droplevels(x)
    lv <- levels(x)
  } else if (is.character(x)) {
    lv <- unique(x[!is.na(x)])
  } else {
    stop("Column '", name, "' has unsupported type: ", paste(class(x), collapse = "/"),
         ". Use logical, numeric 0/1, factor, or character.", call. = FALSE)
  }

  if (length(lv) == 1) {
    # One level observed: treat it as whichever label matches the positive rule
    # if possible; else as negative.
    chosen_pos <- if (!is.null(positive) && positive %in% lv) positive else NA_character_
    vals <- ifelse(is.na(x), NA_integer_,
                   as.integer(!is.na(chosen_pos) & x == chosen_pos))
    return(list(values = vals,
                positive_level = if (is.na(chosen_pos)) lv else chosen_pos))
  }

  if (length(lv) > 2) {
    stop("Column '", name, "' has ", length(lv), " levels: ",
         paste(shQuote(lv), collapse = ", "),
         ". grass handles binary ratings only; filter or recode to two levels.",
         call. = FALSE)
  }

  chosen_pos <- pick_positive(lv, positive, name = name)
  vals <- ifelse(is.na(x), NA_integer_, as.integer(x == chosen_pos))
  list(values = vals, positive_level = chosen_pos)
}

pick_positive <- function(levels, positive = NULL, name = "rating") {
  if (!is.null(positive)) {
    if (!positive %in% levels) {
      stop("`positive = ", shQuote(positive),
           "` is not one of the observed levels for column '", name, "': ",
           paste(shQuote(levels), collapse = ", "), ".", call. = FALSE)
    }
    return(positive)
  }

  lc <- tolower(levels)

  if ("yes" %in% lc) {
    chosen <- levels[which(lc == "yes")[1]]
    msg_once(paste0("coerce_yes_", name),
             sprintf("grass: coercing '%s' to binary (positive = %s, negative = %s). Override with `positive =`.",
                     name, shQuote(chosen), shQuote(levels[levels != chosen][1])))
    return(chosen)
  }

  for (rule in positive_rules) {
    if (rule %in% lc) {
      chosen <- levels[which(lc == rule)[1]]
      msg_once(paste0("coerce_rule_", name),
               sprintf("grass: coercing '%s' to binary (positive = %s, negative = %s). Override with `positive =`.",
                       name, shQuote(chosen), shQuote(levels[levels != chosen][1])))
      return(chosen)
    }
  }

  chosen <- levels[1]
  # First-encountered fallback is the risky one: silently picking the clinical
  # negative as positive inverts prevalence and scrambles reference lookup.
  # Escalate to a warning so automated pipelines surface it.
  warn_msg <- sprintf(
    "grass: no standard positive-class keyword matched the levels of '%s' (%s, %s). Defaulting positive = %s by first-encountered order. Pass `positive =` to be explicit.",
    name, shQuote(levels[1]), shQuote(levels[2]), shQuote(chosen))
  if (is.null(.grass_env$msg_seen[[paste0("warn_first_", name)]])) {
    warning(warn_msg, call. = FALSE)
    .grass_env$msg_seen[[paste0("warn_first_", name)]] <- TRUE
  }
  chosen
}

# Find a column in `data` that is clearly non-binary (and so is unlikely to
# be a rater column). Returns NULL if more than one such column exists or
# none do -- the hint is only useful when unambiguous.
guess_id_col <- function(data) {
  is_binary <- function(x) {
    if (is.logical(x)) return(TRUE)
    if (is.numeric(x)) {
      u <- unique(x[!is.na(x)])
      return(length(u) <= 2 && all(u %in% c(0, 1)))
    }
    if (is.factor(x) || is.character(x)) {
      u <- unique(as.character(x[!is.na(x)]))
      return(length(u) <= 2)
    }
    FALSE
  }
  non_binary <- names(data)[!vapply(data, is_binary, logical(1))]
  if (length(non_binary) == 1L) non_binary else NULL
}

# Drop pairwise-NA rows.
drop_pairwise_na <- function(r1, r2) {
  keep <- !is.na(r1) & !is.na(r2)
  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    frac <- n_dropped / length(r1)
    if (frac > 0.5) {
      stop("More than 50% of rating pairs contain NA (", n_dropped, " of ",
           length(r1), "). Refusing to proceed.", call. = FALSE)
    }
    warning("Dropping ", n_dropped, " of ", length(r1),
            " rating pairs with NA.", call. = FALSE)
  }
  list(r1 = r1[keep], r2 = r2[keep], n_dropped = n_dropped)
}

# ---- Format dispatchers ------------------------------------------------

normalize_input <- function(data, format = c("wide", "matrix", "long", "paired"),
                            positive = NULL, ...) {
  format <- match.arg(format)
  switch(format,
    matrix = normalize_matrix(data),
    wide   = normalize_wide(data, positive = positive, ...),
    long   = normalize_long(data, positive = positive, ...),
    paired = normalize_paired(data, positive = positive, ...)
  )
}

# format = "matrix" always means a 2x2 count table. Raw rater vectors must go
# through "paired" or "wide".
normalize_matrix <- function(data) {
  if (!is.matrix(data) && !is.table(data) && !is.data.frame(data)) {
    stop("For format = \"matrix\" supply a 2x2 count matrix or table.",
         call. = FALSE)
  }
  m <- as.matrix(data)
  if (!all(dim(m) == c(2, 2))) {
    stop("For format = \"matrix\" supply a 2x2 count matrix. Got ",
         paste(dim(m), collapse = "x"), ".", call. = FALSE)
  }
  storage.mode(m) <- "integer"
  if (any(is.na(m)) || any(m < 0)) {
    stop("2x2 count matrix must contain non-negative integers.", call. = FALSE)
  }
  # Pre-built table: skip to build_table()-ready format.
  list(table = m,
       positive_level = if (!is.null(dimnames(m)))
         rownames(m)[2] %||% "1" else "1",
       coerced = FALSE,
       n_dropped = 0L)
}

normalize_wide <- function(data, positive = NULL, rater_cols = NULL,
                           id_col = NULL, response = NULL, ...) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("For format = \"wide\" supply a data.frame or matrix with two rater columns.",
         call. = FALSE)
  }
  data <- as.data.frame(data)

  # `response` is a stats-modelling alias for `rater_cols`. If both are
  # supplied and disagree, error; if only `response` is supplied, use it.
  if (!is.null(response)) {
    if (!is.null(rater_cols) && !identical(sort(response), sort(rater_cols))) {
      stop("`response =` and `rater_cols =` were both supplied and disagree. ",
           "Pick one.", call. = FALSE)
    }
    rater_cols <- response
  }

  if (!is.null(id_col)) {
    miss <- setdiff(id_col, names(data))
    if (length(miss)) {
      stop("id_col not found in data: ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    }
    data <- data[, setdiff(names(data), id_col), drop = FALSE]
  }

  if (!is.null(rater_cols)) {
    if (length(rater_cols) != 2) {
      stop("`rater_cols` / `response` must name exactly two columns.",
           call. = FALSE)
    }
    miss <- setdiff(rater_cols, names(data))
    if (length(miss)) {
      stop("Column(s) not found in data: ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    }
    cols <- data[, rater_cols, drop = FALSE]
  } else {
    if (ncol(data) != 2) {
      col_list <- paste(names(data), collapse = ", ")
      # Deliberately avoid guessing which columns are raters vs. identifiers
      # for `rater_cols`. But we can softly point at an ID candidate when one
      # column is clearly non-binary -- a cheap hint for a cheap-to-ignore knob.
      id_hint <- guess_id_col(data)
      id_suggestion <- if (!is.null(id_hint))
        sprintf(" Column %s looks like an identifier; you could try `id_col = %s`.",
                shQuote(id_hint), shQuote(id_hint))
      else ""
      stop("For format = \"wide\" the data must have exactly two rater columns. Got ",
           ncol(data), " columns: ", col_list,
           ". Use `rater_cols = c(...)` (or `response = c(...)`) to name two rater ",
           "columns, or `id_col = \"...\"` to drop a subject identifier first.",
           id_suggestion,
           call. = FALSE)
    }
    cols <- data
  }

  r1_res <- coerce_binary(cols[[1]], positive = positive, name = names(cols)[1] %||% "rater1")
  r2_res <- coerce_binary(cols[[2]], positive = positive %||% r1_res$positive_level,
                          name = names(cols)[2] %||% "rater2")

  clean <- drop_pairwise_na(r1_res$values, r2_res$values)
  list(r1 = clean$r1, r2 = clean$r2,
       positive_level = r1_res$positive_level,
       coerced = TRUE, n_dropped = clean$n_dropped)
}

normalize_long <- function(data, positive = NULL,
                           subject = "subject", rater = "rater",
                           rating = "rating", response = NULL, ...) {
  if (!is.data.frame(data)) {
    stop("For format = \"long\" supply a data.frame with `subject`, `rater`, `rating` columns.",
         call. = FALSE)
  }
  # `response` is a stats-modelling alias for `rating`. Same conflict rule.
  if (!is.null(response)) {
    if (!identical(rating, "rating") && !identical(rating, response)) {
      stop("`response =` and `rating =` were both supplied and disagree. Pick one.",
           call. = FALSE)
    }
    rating <- response
  }
  need <- c(subject, rater, rating)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Long-format data is missing columns: ", paste(miss, collapse = ", "),
         ". Rename or pass `subject =`, `rater =`, `response =`.", call. = FALSE)
  }

  raters <- unique(as.character(data[[rater]]))
  raters <- raters[!is.na(raters)]
  if (length(raters) != 2) {
    stop("Long-format data must contain exactly 2 distinct raters. Got ",
         length(raters), ": ", paste(shQuote(raters), collapse = ", "),
         ".", call. = FALSE)
  }

  # Detect duplicate subject x rater rows (ambiguous).
  key <- paste0(data[[subject]], "\u0001", data[[rater]])
  if (anyDuplicated(key)) {
    stop("Long-format data contains duplicated subject x rater combinations. ",
         "Deduplicate or aggregate before calling grass.", call. = FALSE)
  }

  # Reshape to wide using base R.
  df <- data.frame(subject = data[[subject]],
                   rater   = as.character(data[[rater]]),
                   rating  = data[[rating]],
                   stringsAsFactors = FALSE)
  wide_df <- reshape(df, idvar = "subject", timevar = "rater",
                     direction = "wide", sep = "_")
  wide_cols <- setdiff(names(wide_df), "subject")
  if (length(wide_cols) != 2) {
    stop("Internal reshape error: expected 2 rater columns after reshape, got ",
         length(wide_cols), ".", call. = FALSE)
  }

  r1_res <- coerce_binary(wide_df[[wide_cols[1]]], positive = positive,
                          name = wide_cols[1])
  r2_res <- coerce_binary(wide_df[[wide_cols[2]]],
                          positive = positive %||% r1_res$positive_level,
                          name = wide_cols[2])
  clean <- drop_pairwise_na(r1_res$values, r2_res$values)
  list(r1 = clean$r1, r2 = clean$r2,
       positive_level = r1_res$positive_level,
       coerced = TRUE, n_dropped = clean$n_dropped)
}

normalize_paired <- function(data, positive = NULL, ...) {
  # Accept list of two vectors, 2-column matrix, or 2-column data.frame.
  if (is.list(data) && !is.data.frame(data) && length(data) == 2) {
    r1 <- data[[1]]; r2 <- data[[2]]
    n1 <- names(data)[1] %||% "rater1"
    n2 <- names(data)[2] %||% "rater2"
  } else if (is.data.frame(data) && ncol(data) == 2) {
    r1 <- data[[1]]; r2 <- data[[2]]
    n1 <- names(data)[1]; n2 <- names(data)[2]
  } else if (is.matrix(data) && ncol(data) == 2) {
    r1 <- data[, 1]; r2 <- data[, 2]
    n1 <- colnames(data)[1] %||% "rater1"
    n2 <- colnames(data)[2] %||% "rater2"
  } else {
    stop("For format = \"paired\" supply a list of two vectors, a 2-column matrix, or a 2-column data.frame.",
         call. = FALSE)
  }

  if (length(r1) != length(r2)) {
    stop("Paired rater vectors must have equal length. Got ",
         length(r1), " and ", length(r2), ".", call. = FALSE)
  }

  r1_res <- coerce_binary(r1, positive = positive, name = n1)
  r2_res <- coerce_binary(r2, positive = positive %||% r1_res$positive_level, name = n2)
  clean <- drop_pairwise_na(r1_res$values, r2_res$values)
  list(r1 = clean$r1, r2 = clean$r2,
       positive_level = r1_res$positive_level,
       coerced = TRUE, n_dropped = clean$n_dropped)
}
