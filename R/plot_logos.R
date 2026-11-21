# ============================================================================
# Logo plot functions
# ============================================================================

#' Plot top known motif logos
#'
#' Displays the top N known motif logos as a grid of ggplot panels, sorted by
#' significance. Logos are read from the pre-generated SVG files in the HOMER
#' output directory.
#'
#' @param results A `homer_output` object from [read_homer_output()], or
#'   its `$known` data frame. `homer_dir` is only needed for a data frame.
#' @param homer_dir Path to the HOMER output directory.
#' @param top_n Integer. Number of top motifs to display. Default 5.
#' @param ncol Integer. Number of columns in the grid. Default 1.
#' @param show Character vector. Stats to display to the right of the logo.
#'   Any combination of `"pvalue"`, `"fold_enrichment"`, `"pct_target"`,
#'   `"pct_background"`. Default: all four stats.
#'
#' @return A `patchwork` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_known_logos(homer, top_n = 3)
#' plot_known_logos(homer, top_n = 3, show = c("pvalue", "fold_enrichment"))
plot_known_logos <- function(
    results,
    homer_dir = NULL,
    top_n = 5,
    ncol = 1,
    show = c("pvalue", "fold_enrichment", "pct_target", "pct_background")
) {
    show <- match.arg(show, several.ok = TRUE)
    if (inherits(results, "homer_output")) {
        if (is.null(homer_dir)) homer_dir <- results$dir
        results <- results$known
    }
    if (is.null(homer_dir)) stop("homer_dir must be provided.")
    if (is.null(results) || nrow(results) == 0) stop("No known motifs found.")

    df <- results[order(results$ln_pvalue), ]
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(0.05), ]
    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = 0.05")
        return(NULL)
    }
    df <- head(df, top_n)

    panels <- lapply(seq_len(nrow(df)), function(i) {
        svg_path <- file.path(
            homer_dir, "knownResults",
            paste0("known", df$rank[i], ".logo.svg")
        )
        label_top <- df$motif_id[i]
        label_stats <- .build_stat_label(df[i, ], show)
        .logo_panel(svg_path, label_top, label_stats)
    })

    patchwork::wrap_plots(panels, ncol = ncol)
}


#' Plot top de novo motif logos
#'
#' Displays the top N de novo motif logos as a grid of ggplot panels, sorted by
#' significance. Logos are read from the pre-generated SVG files in the HOMER
#' output directory.
#'
#' @param denovo A `homer_output` object from [read_homer_output()], or
#'   its `$denovo` data frame. `homer_dir` is only needed for a data frame.
#' @param homer_dir Path to the HOMER output directory.
#' @param top_n Integer. Number of top motifs to display. Default 5.
#' @param ncol Integer. Number of columns in the grid. Default 1.
#' @param show Character vector. Stats to display to the right of the logo.
#'   Any combination of `"pvalue"`, `"fold_enrichment"`, `"pct_target"`,
#'   `"pct_background"`. Default: all four stats.
#'
#' @return A `patchwork` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_denovo_logos(homer, top_n = 3)
#' plot_denovo_logos(homer, top_n = 3, show = c("pvalue", "fold_enrichment"))
plot_denovo_logos <- function(
    denovo,
    homer_dir = NULL,
    top_n = 5,
    ncol = 1,
    show = c("pvalue", "fold_enrichment", "pct_target", "pct_background")
) {
    show <- match.arg(show, several.ok = TRUE)
    if (inherits(denovo, "homer_output")) {
        if (is.null(homer_dir)) homer_dir <- denovo$dir
        denovo <- denovo$denovo
    }
    if (is.null(homer_dir)) stop("homer_dir must be provided.")
    if (is.null(denovo) || nrow(denovo) == 0) stop("No de novo motifs found.")

    df <- denovo[order(denovo$ln_pvalue), ]
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(0.05), ]
    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = 0.05")
        return(NULL)
    }
    df <- head(df, top_n)

    panels <- lapply(seq_len(nrow(df)), function(i) {
        svg_path <- file.path(
            homer_dir, "homerResults",
            paste0("motif", df$motif_num[i], ".logo.svg")
        )
        label_top <- if (!is.na(df$consensus[i]) && nzchar(df$consensus[i])) {
            df$consensus[i]
        } else {
            paste0("De novo ", df$rank[i])
        }
        label_stats <- .build_stat_label(df[i, ], show)
        .logo_panel(svg_path, label_top, label_stats)
    })

    patchwork::wrap_plots(panels, ncol = ncol)
}


# ============================================================================
# Internal helpers
# ============================================================================

.build_stat_label <- function(row, show) {
    parts <- character(0)
    if ("pvalue" %in% show) {
        parts <- c(parts, paste0("p = ", row$pvalue_str))
    }
    if ("fold_enrichment" %in% show && !is.na(row$fold_enrichment)) {
        parts <- c(parts, paste0(
            "Fold enrichment = ", round(row$fold_enrichment, 1)
        ))
    }
    if ("pct_target" %in% show && !is.na(row$pct_target)) {
        parts <- c(parts, paste0("% Target = ", round(row$pct_target, 1), "%"))
    }
    if ("pct_background" %in% show && !is.na(row$pct_background)) {
        parts <- c(parts, paste0(
            "% Background = ", round(row$pct_background, 1), "%"
        ))
    }
    if (length(parts) == 0) {
        return(NULL)
    }
    paste(parts, collapse = "\n")
}

.logo_panel <- function(svg_path, label_top, label_stats = NULL) {
    if (file.exists(svg_path)) {
        svg_str <- paste(
            readLines(svg_path, warn = FALSE, encoding = "UTF-8"),
            collapse = "\n"
        )
        svg_str <- gsub('stroke="[^"]*"', 'stroke="none"', svg_str)
        svg_str <- gsub('stroke-width="[^"]*"', 'stroke-width="0"', svg_str)
        img <- rsvg::rsvg_nativeraster(
            charToRaw(enc2utf8(svg_str)), width = 1800
        )
        grob <- grid::rasterGrob(img, interpolate = TRUE)
        logo_p <- ggplot2::ggplot() +
            ggplot2::annotation_custom(grob,
                xmin = -Inf, xmax = Inf,
                ymin = -Inf, ymax = Inf
            ) +
            ggplot2::xlim(0, 1) +
            ggplot2::ylim(0, 1)
    } else {
        logo_p <- ggplot2::ggplot() +
            ggplot2::xlim(0, 1) +
            ggplot2::ylim(0, 1) +
            ggplot2::annotate("text",
                x = 0.5, y = 0.5,
                label = "logo unavailable", color = "grey60", size = 3
            )
    }

    logo_p <- logo_p +
        ggplot2::labs(title = label_top) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = 9, face = "bold",
                hjust = 0.5, color = "#134e4a"
            ),
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.margin = ggplot2::margin(6, 4, 6, 6)
        )

    if (is.null(label_stats)) {
        return(logo_p)
    }

    stats_p <- ggplot2::ggplot() +
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +
        ggplot2::annotate("text",
            x = 0.1, y = 0.5, label = label_stats,
            hjust = 0, vjust = 0.5, size = 2.5, color = "#64748b",
            lineheight = 1.4
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.margin     = ggplot2::margin(6, 6, 6, 0)
        )

    patchwork::wrap_plots(logo_p, stats_p, ncol = 2, widths = c(3, 1))
}
