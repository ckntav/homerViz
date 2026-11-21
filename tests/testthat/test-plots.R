test_that("plot_known_dotplot returns a ggplot", {
  homer <- read_homer_output(example_dir)
  p <- plot_known_dotplot(homer)
  expect_s3_class(p, "ggplot")
})

test_that("plot_known_dotplot accepts color_scheme argument", {
  homer <- read_homer_output(example_dir)
  for (cs in c("homerViz", "red", "blue", "purple", "green", "viridis")) {
    p <- plot_known_dotplot(homer, top_n = 5, color_scheme = cs)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_denovo_dotplot returns a ggplot", {
  homer <- read_homer_output(example_dir)
  p <- plot_denovo_dotplot(homer)
  expect_s3_class(p, "ggplot")
})

test_that("plot_known_scatter returns a ggplot", {
  homer <- read_homer_output(example_dir)
  p <- plot_known_scatter(homer)
  expect_s3_class(p, "ggplot")
})

test_that("plot_known_scatter label_top_n = 0 suppresses labels", {
  homer <- read_homer_output(example_dir)
  p <- plot_known_scatter(homer, label_top_n = 0)
  expect_s3_class(p, "ggplot")
})

test_that("plot_denovo_scatter returns a ggplot", {
  homer <- read_homer_output(example_dir)
  p <- plot_denovo_scatter(homer)
  expect_s3_class(p, "ggplot")
})

test_that("plot_known_logos returns a patchwork", {
  homer <- read_homer_output(example_dir)
  p <- plot_known_logos(homer, top_n = 3)
  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("plot_denovo_logos returns a patchwork", {
  homer <- read_homer_output(example_dir)
  p <- plot_denovo_logos(homer, top_n = 3)
  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("plot_known_dotplot returns NULL with warning when no significant motifs", {
  homer <- read_homer_output(example_dir)
  expect_warning(
    result <- plot_known_dotplot(homer, p_cutoff = 1e-10000),
    "No significant motifs"
  )
  expect_null(result)
})

test_that("saveViz (re-exported from gVenn) writes a PDF of a homerViz plot", {
  homer <- read_homer_output(example_dir)
  p <- plot_known_dotplot(homer)
  out <- tempfile("saveViz_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  saveViz(p, output_dir = out, output_file = "dotplot", format = "pdf",
          with_date = FALSE, verbose = FALSE)
  expect_true(file.exists(file.path(out, "dotplot.pdf")))
})
