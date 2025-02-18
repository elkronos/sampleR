# File: tests/testthat/UAT.R

library(testthat)
library(dplyr)
library(sf)
library(lubridate)
library(here)

here::i_am("tests/testthat/UAT.R")

# Source the utility and sampling functions.
source(here::here("R", "utils.R"))
source(here::here("R", "sampling_functions.R"))

context("Sampling Functions UAT")

# Create a sample data frame for testing.
df <- data.frame(
  id    = 1:100,
  value = rnorm(100),
  group = rep(letters[1:4], each = 25),
  time  = format(rep(seq.POSIXt(from = as.POSIXct("2020-01-01 00:00:00", tz = "UTC"),
                                by = "hour", length.out = 25), 4),
                 "%Y-%m-%d %H:%M:%S"),
  lat   = runif(100, -90, 90),
  lon   = runif(100, -180, 180)
)

# Use a global polygon that covers the entire globe with a slight inset.
global_poly <- st_sfc(
  st_polygon(list(cbind(
    c(-179, 179, 179, -179, -179),
    c(-89, -89, 89, 89, -89)
  ))),
  crs = 4326
)

#############################################
# 1. SIMPLE RANDOM SAMPLING
#############################################
test_that("simple_random_sampling: basic functionality and parameters", {
  res <- simple_random_sampling(df, sample_size = 10, seed = 123)
  expect_equal(nrow(res), 10)
  
  res2 <- simple_random_sampling(df, sample_size = 10, seed = 123)
  expect_equal(res, res2)
  
  expect_error(simple_random_sampling(df, sample_size = 101, replace = FALSE))
  
  res3 <- simple_random_sampling(df, sample_size = 20, max_rows = 5, seed = 123)
  expect_equal(nrow(res3), 5)
})

#############################################
# 2. STRATIFIED SAMPLING
#############################################
test_that("stratified_sampling: grouping, equal_samples, and error conditions", {
  res <- stratified_sampling(df, strata_column = "group", sample_size = 40, equal_samples = FALSE, seed = 123)
  counts <- table(res$group)
  expect_true(all(counts >= 1))
  
  res_eq <- stratified_sampling(df, strata_column = "group", sample_size = 5, equal_samples = TRUE, seed = 123)
  counts_eq <- table(res_eq$group)
  expect_equal(as.vector(counts_eq), rep(5, 4))
  
  expect_error(stratified_sampling(df, strata_column = "nonexistent", sample_size = 10))
  
  res_max <- stratified_sampling(df, strata_column = "group", sample_size = 40, max_rows = 10, seed = 123)
  expect_equal(nrow(res_max), 10)
})

#############################################
# 3. SYSTEMATIC SAMPLING
#############################################
test_that("systematic_sampling: interval, start, sort_column, and max_rows", {
  res <- systematic_sampling(df, interval = 10, seed = 123)
  expect_true(nrow(res) >= 1)
  
  res_sorted <- systematic_sampling(df, interval = 10, sort_column = "value", seed = 123)
  expect_equal(order(res_sorted$value), seq_len(nrow(res_sorted)))
  
  res_max <- systematic_sampling(df, interval = 10, max_rows = 3, seed = 123)
  expect_equal(nrow(res_max), 3)
})

#############################################
# 4. CLUSTER SAMPLING
#############################################
test_that("cluster_sampling: selecting clusters, balanced option, and error cases", {
  res <- cluster_sampling(df, cluster_column = "group", num_clusters = 2, seed = 123)
  expect_equal(length(unique(res$group)), 2)
  
  res_bal <- cluster_sampling(df, cluster_column = "group", num_clusters = 4, balanced = TRUE, seed = 123)
  counts_bal <- table(res_bal$group)
  expect_equal(length(unique(counts_bal)), 1)
  
  expect_error(cluster_sampling(df, cluster_column = "group", num_clusters = 10, seed = 123))
  
  res_max <- cluster_sampling(df, cluster_column = "group", num_clusters = 2, max_rows = 5, seed = 123)
  expect_equal(nrow(res_max), 5)
})

#############################################
# 5. MULTI-STAGE SAMPLING
#############################################
test_that("multi_stage_sampling: cluster selection, stage two sampling, and parameters", {
  res <- multi_stage_sampling(df, cluster_column = "group", num_clusters = 2, stage_two_sample_size = 3, seed = 123)
  expect_equal(nrow(res), 6)
  
  res_prop <- multi_stage_sampling(df, cluster_column = "group", num_clusters = 4, stage_two_sample_size = 10, proportional_stage_two = TRUE, seed = 123)
  expect_true(nrow(res_prop) >= 4)
  
  small_df <- df[1:3, ]
  small_df$group <- c("a", "b", "c")
  expect_error(multi_stage_sampling(small_df, cluster_column = "group", num_clusters = 3, stage_two_sample_size = 2, replace = FALSE, seed = 123))
  
  res_max <- multi_stage_sampling(df, cluster_column = "group", num_clusters = 2, stage_two_sample_size = 10, max_rows = 3, seed = 123)
  expect_equal(nrow(res_max), 3)
})

#############################################
# 6. WEIGHTED SAMPLING
#############################################
test_that("weighted_sampling: weights, normalization, and error conditions", {
  df$weight <- runif(nrow(df), 0.1, 10)
  
  res <- weighted_sampling(df, weights_column = "weight", sample_size = 10, seed = 123)
  expect_equal(nrow(res), 10)
  
  res_mm <- weighted_sampling(df, weights_column = "weight", sample_size = 10, normalization = "min-max", seed = 123)
  expect_equal(nrow(res_mm), 10)
  
  res_z <- weighted_sampling(df, weights_column = "weight", sample_size = 10, normalization = "z-score", seed = 123)
  expect_equal(nrow(res_z), 10)
  
  df_bad <- df; df_bad$bad_weight <- -abs(df_bad$weight)
  expect_error(weighted_sampling(df_bad, weights_column = "bad_weight", sample_size = 10, seed = 123))
  
  res_max <- weighted_sampling(df, weights_column = "weight", sample_size = 20, max_rows = 5, seed = 123)
  expect_equal(nrow(res_max), 5)
})

#############################################
# 7. RESERVOIR SAMPLING
#############################################
test_that("reservoir_sampling: basic functionality and max_rows", {
  res <- reservoir_sampling(df, sample_size = 10, seed = 123)
  expect_equal(length(res), 10)
  
  df_list <- split(df, seq_len(nrow(df)))
  res_list <- reservoir_sampling(df_list, sample_size = 5, seed = 123)
  expect_equal(length(res_list), 5)
  
  res_max <- reservoir_sampling(df, sample_size = 10, max_rows = 3, seed = 123)
  expect_equal(length(res_max), 3)
})

#############################################
# 8. BOOTSTRAP SAMPLING
#############################################
test_that("bootstrap_sampling: simple and block methods, and error conditions", {
  res_simple <- bootstrap_sampling(df, num_samples = 5, sample_size = 10, method = "simple", seed = 123)
  expect_equal(length(res_simple), 5)
  expect_true(all(sapply(res_simple, nrow) == 10))
  
  res_block <- bootstrap_sampling(df, num_samples = 3, sample_size = 10, method = "block", seed = 123)
  expect_equal(length(res_block), 3)
  expect_true(all(sapply(res_block, nrow) == 10))
  
  expect_error(bootstrap_sampling(df, num_samples = 5, sample_size = 10, method = "invalid", seed = 123))
  
  res_max <- bootstrap_sampling(df, num_samples = 5, sample_size = 10, max_rows = 20, method = "simple", seed = 123)
  combined <- unique(do.call(rbind, res_max))
  expect_lte(nrow(combined), 20)
})

#############################################
# 9. TEMPORAL SAMPLING
#############################################
test_that("temporal_sampling: time window, interval, and max_rows", {
  # Use a window from 00:00 to 18:00 to generate three 6-hour intervals.
  start_time <- "2020-01-01 00:00:00"
  end_time   <- "2020-01-01 18:00:00"
  
  res <- temporal_sampling(df, time_column = "time", start_time = start_time, end_time = end_time,
                           interval = 6, sample_size = 2, unit = "hours", seed = 123)
  time_vals <- lubridate::ymd_hms(res$time)
  expect_true(all(time_vals >= lubridate::ymd_hms(start_time) & time_vals < lubridate::ymd_hms(end_time)))
  
  res_max <- temporal_sampling(df, time_column = "time", start_time = start_time, end_time = end_time,
                               interval = 6, sample_size = 2, max_rows = 3, unit = "hours", seed = 123)
  expect_equal(nrow(res_max), 3)
})

#############################################
# 10. SPATIAL SAMPLING
#############################################
test_that("spatial_sampling: region filtering, complex_region, and max_rows", {
  res <- spatial_sampling(df, latitude_column = "lat", longitude_column = "lon",
                          region = global_poly, sample_size = 5, seed = 123)
  expect_equal(nrow(res), 5)
  
  empty_poly <- st_sfc(
    st_polygon(list(cbind(
      c(100, 101, 101, 100, 100),
      c(100, 100, 101, 101, 100)
    ))),
    crs = 4326
  )
  expect_error(spatial_sampling(df, latitude_column = "lat", longitude_column = "lon",
                                region = empty_poly, sample_size = 1, seed = 123))
  
  res_max <- spatial_sampling(df, latitude_column = "lat", longitude_column = "lon",
                              region = global_poly, sample_size = 10, max_rows = 5, seed = 123)
  expect_equal(nrow(res_max), 5)
})
