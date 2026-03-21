# homerViz

![Maintainer](https://img.shields.io/badge/maintainer-ckntav-blue)
![Last commit](https://img.shields.io/github/last-commit/ckntav/homerViz)

**Visualization tools for HOMER motif enrichment results**

**homerViz** is an R package that parses
[HOMER](http://homer.ucsd.edu/homer/) `findMotifsGenome.pl` output directories
into tidy data frames and generates publication-ready visualizations.

## Features

- **Parsing** — reads `knownResults.txt` and `homerResults/motif*.motif` files into tidy data frames with consistent column names.
- **Dot plots** — top enriched motifs ranked by significance, coloured and
  sized by fold enrichment.
- **Scatter plots** — % target vs % background, coloured by significance and
  sized by fold enrichment.
- **Logo panels** — motif logo grids with embedded enrichment statistics.
- **HTML reports** — interactive, sortable table reports with motif logos for both known and de novo results.

## Installation

To install the development version from GitHub, use:

``` r
# install.packages("pak")  # if not already installed
pak::pak("ckntav/homerViz")

# or, alternatively:
# install.packages("devtools")  # if not already installed
devtools::install_github("ckntav/homerViz")
```

## Quick start

This brief example demonstrates the visualization capabilities of homerViz using MED1 ChIP-seq data from the A549 cell line, referenced against the HOCOMOCO motif database (v14).

> command line: findMotifsGenome.pl A549_MED1_Ctrl.stdchr.chr1.bed hg38 A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001 -p 8 -size 200 -mask -mknown H14CORE_homer_format_0.0001.motif -mcheck H14CORE_homer_format_0.0001.motif


### 1. Read HOMER output files

```r
library(homerViz)

# Read a HOMER output directory
example_dir <- system.file(
    "extdata",
    "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
    package = "homerViz"
)
homer <- read_homer_output("path/to/homer_output/")
```

### 2. Visualize top motif

#### Dot plots
```r
plot_known_dotplot(homer)
plot_denovo_dotplot(homer)
```

#### Scatter plots
```r
plot_known_scatter(homer)
plot_denovo_scatter(homer)
```

#### Logo
```r
plot_known_logos(homer)
plot_denovo_logos(homer)
```

### 3. Interactive HTML reports

```r
save_homer_html(homer)
```

Example reports generated from the bundled dataset are available online:

- [Known motif enrichment report](https://ckntav.github.io/homerViz/example/knownResults_homerViz.html)
- [De novo motif enrichment report](https://ckntav.github.io/homerViz/example/homerResults_homerViz.html)


> **Note:** The HTML reports load Bootstrap, DataTables, and jQuery from
> public CDNs and require an internet connection to render correctly.
