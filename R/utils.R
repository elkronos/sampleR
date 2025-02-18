# utils.R
# Utility functions for the sampling package

# Set the random seed for reproducibility.
set_seed <- function(seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
}

# Return the current random state (for debugging purposes).
get_random_state <- function(seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  return(.Random.seed)
}

# Setup logging. (This is a basic placeholder.)
setup_logging <- function(level = "INFO") {
  message(paste("Logging level set to", level))
}

# Validate that the input is a non-empty data.frame and that it contains all required columns.
validate_dataframe <- function(data, required_columns = NULL) {
  if (!is.data.frame(data)) {
    stop("Input must be a data.frame.")
  }
  if (nrow(data) == 0) {
    stop("The dataset is empty.")
  }
  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("The dataset is missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  invisible(TRUE)
}

# Apply a maximum rows limit to a data.frame, if specified.
apply_max_rows <- function(data, max_rows = NULL) {
  if (!is.null(max_rows) && is.numeric(max_rows)) {
    return(head(data, max_rows))
  }
  return(data)
}
