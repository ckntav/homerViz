test_that("save_known_html writes an HTML file", {
  homer <- read_homer_output(example_dir)
  out   <- save_known_html(homer, output_dir = tempdir())

  expect_true(file.exists(out))
  expect_true(grepl("\\.html$", out))

  content <- readLines(out, warn = FALSE)
  expect_true(any(grepl("<!DOCTYPE html>", content, fixed = TRUE)))
})

test_that("save_denovo_html writes an HTML file", {
  homer <- read_homer_output(example_dir)
  out   <- save_denovo_html(homer, output_dir = tempdir())

  expect_true(file.exists(out))
  expect_true(grepl("\\.html$", out))

  content <- readLines(out, warn = FALSE)
  expect_true(any(grepl("<!DOCTYPE html>", content, fixed = TRUE)))
})

test_that("save_homer_html writes both HTML files", {
  homer <- read_homer_output(example_dir)
  out   <- save_homer_html(homer, output_dir = tempdir())

  expect_type(out, "list")
  expect_true(file.exists(out$known))
  expect_true(file.exists(out$denovo))
})

test_that("save_known_html returns the output path invisibly", {
  homer <- read_homer_output(example_dir)
  out   <- save_known_html(homer, output_dir = tempdir())
  expect_type(out, "character")
  expect_length(out, 1)
})
