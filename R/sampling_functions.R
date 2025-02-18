# sampling_functions.R
# (Assumes that utils.R has already been sourced.)

# 1. SIMPLE RANDOM SAMPLING
simple_random_sampling <- function(data, sample_size, max_rows = NULL, replace = FALSE, seed = NULL) {
  validate_dataframe(data)
  if (!is.null(seed)) set.seed(seed)
  
  if (sample_size > nrow(data) && !replace)
    stop("Sample size cannot be larger than the population when sampling without replacement.")
  
  sampled_data <- data[sample(seq_len(nrow(data)), size = sample_size, replace = replace), , drop = FALSE]
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(sampled_data)
}

# 2. STRATIFIED SAMPLING
stratified_sampling <- function(data, strata_column, sample_size, max_rows = NULL, replace = FALSE, equal_samples = FALSE, seed = NULL) {
  validate_dataframe(data, required_columns = strata_column)
  if (!is.null(seed)) set.seed(seed)
  
  total_population <- nrow(data)
  sampled_data <- dplyr::group_by(data, .data[[strata_column]]) %>%
    dplyr::group_modify(function(x, y) {
      stratum_n <- nrow(x)
      if (equal_samples) {
        n_sample <- sample_size
      } else {
        n_sample <- max(1, round(sample_size * (stratum_n / total_population)))
      }
      if (n_sample > stratum_n && !replace)
        stop(paste("Sample size cannot be larger than the number of rows in stratum", unique(x[[strata_column]])))
      dplyr::slice_sample(x, n = n_sample, replace = replace)
    }) %>% dplyr::ungroup()
  
  sampled_data <- as.data.frame(sampled_data)
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(sampled_data)
}

# 3. SYSTEMATIC SAMPLING
systematic_sampling <- function(data, interval, start = NULL, max_rows = NULL, sort_column = NULL, seed = NULL) {
  validate_dataframe(data)
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.null(sort_column)) {
    validate_dataframe(data, required_columns = sort_column)
    data <- data[order(data[[sort_column]]), ]
  }
  
  if (is.null(start)) {
    start <- sample(1:interval, 1)
  }
  
  indices <- seq(from = start, to = nrow(data), by = interval)
  sampled_data <- data[indices, , drop = FALSE]
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(sampled_data)
}

# 4. CLUSTER SAMPLING
cluster_sampling <- function(data, cluster_column, num_clusters, max_rows = NULL, balanced = FALSE, seed = NULL) {
  validate_dataframe(data, required_columns = cluster_column)
  if (!is.null(seed)) set.seed(seed)
  
  clusters <- unique(data[[cluster_column]])
  if (num_clusters > length(clusters)) {
    stop("Number of clusters requested exceeds available clusters.")
  }
  
  sampled_clusters <- sample(clusters, num_clusters, replace = FALSE)
  sampled_data <- data[data[[cluster_column]] %in% sampled_clusters, ]
  
  if (balanced) {
    sampled_data <- dplyr::group_by(sampled_data, .data[[cluster_column]]) %>%
      dplyr::slice_sample(n = min(as.vector(table(sampled_data[[cluster_column]])))) %>%
      dplyr::ungroup()
  }
  
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(as.data.frame(sampled_data))
}

# 5. MULTI-STAGE SAMPLING
multi_stage_sampling <- function(data, cluster_column, num_clusters, stage_two_sample_size, max_rows = NULL, proportional_stage_two = FALSE, replace = FALSE, seed = NULL) {
  validate_dataframe(data, required_columns = cluster_column)
  if (!is.null(seed)) set.seed(seed)
  
  clusters <- unique(data[[cluster_column]])
  if (num_clusters > length(clusters)) {
    stop("Number of clusters requested exceeds available clusters.")
  }
  
  sampled_clusters <- sample(clusters, num_clusters, replace = FALSE)
  stage_two_sample <- lapply(sampled_clusters, function(cluster) {
    cluster_data <- data[data[[cluster_column]] == cluster, ]
    if (proportional_stage_two) {
      n_sample <- max(1, round(stage_two_sample_size * (nrow(cluster_data) / nrow(data))))
    } else {
      n_sample <- stage_two_sample_size
    }
    if (n_sample > nrow(cluster_data) && !replace)
      stop(paste("Sample size requested for cluster", cluster, "exceeds available rows."))
    cluster_data[sample(seq_len(nrow(cluster_data)), size = n_sample, replace = replace), , drop = FALSE]
  })
  
  sampled_data <- do.call(rbind, stage_two_sample)
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(as.data.frame(sampled_data))
}

# 6. WEIGHTED SAMPLING
weighted_sampling <- function(data, weights_column, sample_size, max_rows = NULL, replace = FALSE, normalization = NULL, seed = NULL) {
  validate_dataframe(data, required_columns = weights_column)
  if (!is.null(seed)) set.seed(seed)
  
  data_copy <- data
  if (!is.null(normalization)) {
    if (normalization == "min-max") {
      w_min <- min(data_copy[[weights_column]], na.rm = TRUE)
      w_max <- max(data_copy[[weights_column]], na.rm = TRUE)
      data_copy[[weights_column]] <- (data_copy[[weights_column]] - w_min) / (w_max - w_min)
      data_copy[[weights_column]] <- data_copy[[weights_column]] + .Machine$double.eps
    } else if (normalization == "z-score") {
      w_mean <- mean(data_copy[[weights_column]], na.rm = TRUE)
      w_sd <- sd(data_copy[[weights_column]], na.rm = TRUE)
      data_copy[[weights_column]] <- (data_copy[[weights_column]] - w_mean) / w_sd
      # Shift z-scores so that the minimum becomes > 0.
      data_copy[[weights_column]] <- data_copy[[weights_column]] - min(data_copy[[weights_column]], na.rm = TRUE) + .Machine$double.eps
    } else {
      stop("Normalization method must be either 'min-max' or 'z-score'.")
    }
  }
  
  if (any(data_copy[[weights_column]] <= 0, na.rm = TRUE)) {
    stop("All weights must be positive after normalization.")
  }
  
  sampled_indices <- sample(seq_len(nrow(data_copy)), size = sample_size, replace = replace, prob = data_copy[[weights_column]])
  sampled_data <- data_copy[sampled_indices, , drop = FALSE]
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(sampled_data)
}

# 7. RESERVOIR SAMPLING
reservoir_sampling <- function(data_stream, sample_size, max_rows = NULL, handle_infinite_stream = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.data.frame(data_stream)) {
    stream <- split(data_stream, seq_len(nrow(data_stream)))
  } else if (is.list(data_stream)) {
    stream <- data_stream
  } else {
    stop("data_stream must be a data.frame or a list.")
  }
  
  reservoir <- vector("list", sample_size)
  count <- 0
  for (item in stream) {
    count <- count + 1
    if (count <= sample_size) {
      reservoir[[count]] <- item
    } else {
      j <- sample.int(count, 1)
      if (j <= sample_size) {
        reservoir[[j]] <- item
      }
    }
    if (handle_infinite_stream && count > 10000) break
  }
  
  if (!is.null(max_rows)) {
    reservoir <- reservoir[seq_len(min(length(reservoir), max_rows))]
  }
  return(reservoir)
}

# 8. BOOTSTRAP SAMPLING
bootstrap_sampling <- function(data, num_samples, sample_size, max_rows = NULL, method = "simple", seed = NULL) {
  validate_dataframe(data)
  if (!is.null(seed)) set.seed(seed)
  
  bootstrap_samples <- vector("list", num_samples)
  
  if (!(method %in% c("simple", "block"))) {
    stop("Method must be either 'simple' or 'block'.")
  }
  
  if (method == "simple") {
    for (i in seq_len(num_samples)) {
      bootstrap_samples[[i]] <- data[sample(seq_len(nrow(data)), size = sample_size, replace = TRUE), , drop = FALSE]
    }
  } else if (method == "block") {
    block_starts <- seq(1, nrow(data), by = sample_size)
    blocks <- lapply(block_starts, function(start) {
      end <- min(start + sample_size - 1, nrow(data))
      data[start:end, , drop = FALSE]
    })
    if (length(blocks) == 0) stop("No blocks created; check sample_size and data length.")
    for (i in seq_len(num_samples)) {
      chosen_block <- sample(blocks, 1)[[1]]
      bootstrap_samples[[i]] <- chosen_block[sample(seq_len(nrow(chosen_block)), size = sample_size, replace = TRUE), , drop = FALSE]
    }
  }
  
  if (!is.null(max_rows)) {
    combined <- unique(do.call(rbind, bootstrap_samples))
    combined <- head(combined, max_rows)
    n_rows <- nrow(combined)
    bootstrap_samples <- split(combined, ceiling(seq_len(n_rows) / sample_size))
  }
  
  return(bootstrap_samples)
}

# 9. TEMPORAL SAMPLING (Updated with "hours" option)
temporal_sampling <- function(data, time_column, start_time, end_time, interval, sample_size, max_rows = NULL, time_zone = NULL, unit = "days", seed = NULL) {
  validate_dataframe(data, required_columns = time_column)
  if (!is.null(seed)) set.seed(seed)
  
  data_copy <- data
  if (!is.null(time_zone)) {
    data_copy[[time_column]] <- lubridate::with_tz(lubridate::ymd_hms(data_copy[[time_column]]), tzone = time_zone)
  } else {
    data_copy[[time_column]] <- lubridate::ymd_hms(data_copy[[time_column]])
  }
  
  start_time_dt <- lubridate::ymd_hms(start_time)
  end_time_dt   <- lubridate::ymd_hms(end_time)
  
  # Filter to complete intervals only.
  data_copy <- dplyr::filter(data_copy, .data[[time_column]] >= start_time_dt, .data[[time_column]] < end_time_dt)
  
  delta <- switch(unit,
                  "hours" = lubridate::dhours(interval),
                  "days" = lubridate::ddays(interval),
                  "weeks" = lubridate::dweeks(interval),
                  "months" = lubridate::dmonths(interval),
                  stop("Unsupported unit. Use 'hours', 'days', 'weeks', or 'months'."))
  
  # Only include complete intervals: from start_time_dt to (end_time_dt - delta).
  time_seq <- seq(from = start_time_dt, to = end_time_dt - delta, by = delta)
  
  sampled_data_list <- lapply(time_seq, function(t0) {
    t1 <- t0 + delta
    interval_data <- dplyr::filter(data_copy, .data[[time_column]] >= t0, .data[[time_column]] < t1)
    if (nrow(interval_data) == 0) return(NULL)
    if (nrow(interval_data) <= sample_size) {
      return(interval_data)
    } else {
      return(dplyr::slice_sample(interval_data, n = sample_size))
    }
  })
  
  sampled_data <- do.call(rbind, Filter(Negate(is.null), sampled_data_list))
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(as.data.frame(sampled_data))
}

# 10. SPATIAL SAMPLING (Updated to use st_intersects)
spatial_sampling <- function(data, latitude_column, longitude_column, region, sample_size, max_rows = NULL, complex_region = FALSE, seed = NULL) {
  validate_dataframe(data, required_columns = c(latitude_column, longitude_column))
  if (!is.null(seed)) set.seed(seed)
  
  data_sf <- sf::st_as_sf(data, coords = c(longitude_column, latitude_column), crs = 4326, remove = FALSE)
  
  if (complex_region && is.list(region)) {
    region <- do.call(sf::st_union, region)
  }
  
  # Use st_intersects to include points on the boundary.
  within_idx <- lengths(sf::st_intersects(data_sf, region)) > 0
  region_data <- data[within_idx, , drop = FALSE]
  
  if (sample_size > nrow(region_data)) {
    stop("Sample size cannot be larger than the number of points within the region.")
  }
  
  sampled_data <- region_data[sample(seq_len(nrow(region_data)), size = sample_size, replace = FALSE), , drop = FALSE]
  sampled_data <- apply_max_rows(sampled_data, max_rows)
  return(sampled_data)
}
