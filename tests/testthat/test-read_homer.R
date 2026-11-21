test_that("read_homer_known returns a data frame with expected structure", {
  known <- read_homer_known(example_dir)

  expect_s3_class(known, "data.frame")
  expect_gt(nrow(known), 0)

  expected_cols <- c(
    "rank", "motif_id", "consensus", "motif_length",
    "pvalue_str", "ln_pvalue", "neg_log10p", "qvalue",
    "n_target", "pct_target", "n_background", "pct_background",
    "fold_enrichment", "n_target_total", "n_background_total"
  )
  expect_true(all(expected_cols %in% names(known)))
})

test_that("read_homer_known columns have correct types", {
  known <- read_homer_known(example_dir)

  expect_type(known$rank,             "integer")
  expect_type(known$motif_id,         "character")
  expect_type(known$consensus,        "character")
  expect_type(known$ln_pvalue,        "double")
  expect_type(known$neg_log10p,       "double")
  expect_type(known$fold_enrichment,  "double")
  expect_type(known$n_target_total,   "integer")
})

test_that("read_homer_known neg_log10p is positive", {
  known <- read_homer_known(example_dir)
  expect_true(all(known$neg_log10p[!is.na(known$neg_log10p)] >= 0))
})

test_that("read_homer_known fold_enrichment = pct_target / pct_background", {
  known <- read_homer_known(example_dir)
  valid <- !is.na(known$fold_enrichment) & known$pct_background > 0
  expect_equal(
    known$fold_enrichment[valid],
    known$pct_target[valid] / known$pct_background[valid],
    tolerance = 1e-6
  )
})

test_that("read_homer_known rank is 1-based sequential", {
  known <- read_homer_known(example_dir)
  expect_equal(known$rank, seq_len(nrow(known)))
})

test_that("read_homer_known errors on missing file", {
  expect_error(read_homer_known(tempdir()), "knownResults.txt not found")
})

test_that("read_homer_denovo returns a data frame with expected structure", {
  denovo <- read_homer_denovo(example_dir)

  expect_s3_class(denovo, "data.frame")
  expect_gt(nrow(denovo), 0)

  expected_cols <- c(
    "rank", "motif_num", "name", "consensus", "motif_length",
    "best_match", "best_match_name", "best_match_score",
    "pvalue_str", "ln_pvalue", "neg_log10p",
    "n_target", "pct_target", "n_background", "pct_background",
    "fold_enrichment", "n_target_total", "n_background_total",
    "pos_std_target", "pos_std_bg"
  )
  expect_true(all(expected_cols %in% names(denovo)))
})

test_that("read_homer_denovo rank is 1-based sequential", {
  denovo <- read_homer_denovo(example_dir)
  expect_equal(denovo$rank, seq_len(nrow(denovo)))
})

test_that("read_homer_denovo motif_num points to the matching motif file", {
  denovo <- read_homer_denovo(example_dir)
  expect_type(denovo$motif_num, "integer")
  file_consensus <- vapply(denovo$motif_num, function(n) {
    h <- readLines(file.path(
      example_dir, "homerResults", sprintf("motif%d.motif", n)
    ), n = 1)
    sub("^>", "", strsplit(h, "\t")[[1]][1])
  }, character(1))
  expect_equal(denovo$consensus, file_consensus)
})

test_that("read_homer_denovo keeps motif_num when a motif file is missing", {
  tmp <- file.path(tempfile(), "homer")
  dir.create(tmp, recursive = TRUE)
  file.copy(list.files(example_dir, full.names = TRUE), tmp, recursive = TRUE)
  file.remove(file.path(tmp, "homerResults", "motif3.motif"))

  denovo <- read_homer_denovo(tmp)
  expect_false(3L %in% denovo$motif_num)
  # rank stays sequential, but motif_num skips the gap
  expect_equal(denovo$rank, seq_len(nrow(denovo)))
  expect_equal(denovo$motif_num[3], 4L)
  hdr <- readLines(file.path(tmp, "homerResults", "motif4.motif"), n = 1)
  hdr_consensus <- sub("^>", "", strsplit(hdr, "\t")[[1]][1])
  expect_equal(denovo$consensus[3], hdr_consensus)
})

test_that("read_homer_denovo neg_log10p is positive", {
  denovo <- read_homer_denovo(example_dir)
  valid  <- !is.na(denovo$neg_log10p)
  expect_true(all(denovo$neg_log10p[valid] > 0))
})

test_that("read_homer_denovo returns NULL for a directory with no motif files", {
  tmp <- tempdir()
  expect_null(read_homer_denovo(tmp))
})

test_that("read_homer_output returns a homer_output object", {
  homer <- read_homer_output(example_dir)

  expect_s3_class(homer, "homer_output")
  expect_named(homer, c("known", "denovo", "params", "input",
                         "known_db", "denovo_db", "dir"))
})

test_that("read_homer_output $known and $denovo are data frames", {
  homer <- read_homer_output(example_dir)

  expect_s3_class(homer$known,  "data.frame")
  expect_s3_class(homer$denovo, "data.frame")
})

test_that("read_homer_output $dir matches the input path", {
  homer <- read_homer_output(example_dir)
  expect_equal(homer$dir, example_dir)
})

test_that("print.homer_output produces output", {
  homer <- read_homer_output(example_dir)
  expect_output(print(homer), "HOMER output")
})
