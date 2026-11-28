# ============================================================================
# Enrichment dot plot functions
# ============================================================================

#' Dot plot of top known motif enrichment
#'
#' Plots the top N significant known motifs as a dot plot. The x-axis shows
#' -log10(p), the y-axis shows motif names ordered by significance, and dot
#' size and colour encode fold enrichment.
#'
#' @param results A `homer_output` object from [read_homer_output()], or
#'   its `$known` data frame.
#' @param p_cutoff Numeric. P-value threshold for inclusion. Default `0.05`.
#' @param top_n Integer. Number of top motifs to display. Default `10`.
#' @param subtitle Character. Subtitle shown below the plot title. Default
#'   `NULL` (auto-derived from `homer_output`: BED file name + database).
#' @param color_scheme Character. One of `"homerViz"`, `"red"`, `"blue"`,
#'   `"green"`, `"purple"`, or `"viridis"`. Default `"homerViz"`.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_known_dotplot(homer)
#' plot_known_dotplot(homer, top_n = 5, color_scheme = "viridis")
plot_known_dotplot <- function(
    results,
    p_cutoff = 0.05,
    top_n = 10,
    subtitle = NULL,
    color_scheme = c("homerViz", "red", "green", "blue", "purple", "viridis")
) {
    color_scheme <- match.arg(color_scheme)
    .check_p_cutoff(p_cutoff)
    if (inherits(results, "homer_output")) {
        if (is.null(subtitle)) {
            subtitle <- .build_subtitle(results$input, results$known_db)
        }
        results <- results$known
    }
    if (is.null(results) || nrow(results) == 0) stop("No known motifs found.")

    df <- results[order(results$ln_pvalue), ]
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(p_cutoff), ]
    df <- .drop_na_fe(df, "known motif")

    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = ", p_cutoff)
        return(NULL)
    }

    df <- head(df, top_n)
    uid <- make.unique(df$motif_id)
    df$label <- factor(uid, levels = rev(uid))

    .enrichment_dotplot(
        df           = df,
        y_var        = "label",
        title        = paste0("Known motif enrichment (top ", nrow(df), ")"),
        subtitle     = subtitle,
        color_scheme = color_scheme
    )
}


#' Dot plot of top de novo motif enrichment
#'
#' Plots the top N significant de novo motifs as a dot plot. The x-axis shows
#' -log10(p), the y-axis shows consensus sequences ordered by significance, and
#' dot size and colour encode fold enrichment.
#'
#' @param denovo A `homer_output` object from [read_homer_output()], or
#'   its `$denovo` data frame.
#' @param p_cutoff Numeric. P-value threshold for inclusion. Default `0.05`.
#' @param top_n Integer. Number of top motifs to display. Default `10`.
#' @param subtitle Character. Subtitle shown below the plot title. Default
#'   `NULL` (auto-derived from `homer_output`: BED file name + database).
#' @param color_scheme Character. One of `"homerViz"`, `"red"`, `"blue"`,
#'   `"green"`, `"purple"`, or `"viridis"`. Default `"homerViz"`.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_denovo_dotplot(homer)
#' plot_denovo_dotplot(homer, top_n = 5, color_scheme = "purple")
plot_denovo_dotplot <- function(
    denovo,
    p_cutoff = 0.05,
    top_n = 10,
    subtitle = NULL,
    color_scheme = c("homerViz", "red", "green", "blue", "purple", "viridis")
) {
    color_scheme <- match.arg(color_scheme)
    .check_p_cutoff(p_cutoff)
    if (inherits(denovo, "homer_output")) {
        if (is.null(subtitle)) {
            subtitle <- .build_subtitle(denovo$input, denovo$denovo_db)
        }
        denovo <- denovo$denovo
    }
    if (is.null(denovo) || nrow(denovo) == 0) stop("No de novo motifs found.")

    df <- denovo[order(denovo$ln_pvalue), ]
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(p_cutoff), ]
    df <- .drop_na_fe(df, "de novo motif")

    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = ", p_cutoff)
        return(NULL)
    }

    df <- head(df, top_n)
    df$label <- ifelse(
        !is.na(df$consensus) & nzchar(df$consensus),
        df$consensus,
        paste0("De novo ", df$rank)
    )
    df$label <- make.unique(df$label)
    df$label <- factor(df$label, levels = rev(df$label))

    .enrichment_dotplot(
        df           = df,
        y_var        = "label",
        title        = bquote(
            bolditalic("De novo") ~
                bold(.(paste0("motif (top ", nrow(df), ")")))
        ),
        subtitle     = subtitle,
        color_scheme = color_scheme
    )
}


# ============================================================================
# Scatter plot functions
# ============================================================================

#' Scatter plot of known motif enrichment
#'
#' Plots all significant known motifs as a scatter plot. The x-axis shows
#' % background, the y-axis shows % target, dot size encodes fold enrichment,
#' and dot colour encodes -log10(p-value).
#'
#' @param results A `homer_output` object from [read_homer_output()], or
#'   its `$known` data frame.
#' @param p_cutoff Numeric. P-value threshold for inclusion. Default `0.05`.
#' @param subtitle Character. Subtitle shown below the plot title. Default
#'   `NULL` (auto-derived from `homer_output`: BED file name + database).
#' @param color_scheme Character. One of `"homerViz"`, `"red"`, `"blue"`,
#'   `"green"`, `"purple"`, or `"viridis"`. Default `"homerViz"`.
#' @param label_top_n Integer or `NULL`. Number of top motifs to label.
#'   Default `NULL` = auto-label high-confidence motifs (p <= 1e-50).
#'   Set to `0` to suppress all labels.
#' @param size_by_fe Logical. If `TRUE`, dot size encodes fold enrichment.
#'   If `FALSE` (default), all dots are drawn at the same fixed size.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_known_scatter(homer)
#' plot_known_scatter(homer, label_top_n = 5)
#' plot_known_scatter(homer, size_by_fe = FALSE)
plot_known_scatter <- function(
    results,
    p_cutoff = 0.05,
    subtitle = NULL,
    color_scheme = c("homerViz", "red", "green", "blue", "purple", "viridis"),
    label_top_n = NULL,
    size_by_fe = FALSE
) {
    color_scheme <- match.arg(color_scheme)
    .check_p_cutoff(p_cutoff)
    if (inherits(results, "homer_output")) {
        if (is.null(subtitle)) {
            subtitle <- .build_subtitle(results$input, results$known_db)
        }
        results <- results$known
    }
    if (is.null(results) || nrow(results) == 0) stop("No known motifs found.")

    df <- results[
        !is.na(results$ln_pvalue) & results$ln_pvalue <= log(p_cutoff),
    ]
    if (size_by_fe) df <- .drop_na_fe(df, "known motif")

    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = ", p_cutoff)
        return(NULL)
    }

    df$label <- df$motif_id
    label_df <- .resolve_label_df(df, label_top_n)

    .enrichment_scatter(
        df           = df,
        label_df     = label_df,
        title        = "Known motif enrichment",
        subtitle     = subtitle,
        color_scheme = color_scheme,
        size_by_fe   = size_by_fe
    )
}


#' Scatter plot of de novo motif enrichment
#'
#' Plots all significant de novo motifs as a scatter plot. The x-axis shows
#' % background, the y-axis shows % target, dot size encodes fold enrichment,
#' and dot colour encodes -log10(p-value).
#'
#' @param denovo A `homer_output` object from [read_homer_output()], or
#'   its `$denovo` data frame.
#' @param p_cutoff Numeric. P-value threshold for inclusion. Default `0.05`.
#' @param subtitle Character. Subtitle shown below the plot title. Default
#'   `NULL` (auto-derived from `homer_output`: BED file name + database).
#' @param color_scheme Character. One of `"homerViz"`, `"red"`, `"blue"`,
#'   `"green"`, `"purple"`, or `"viridis"`. Default `"homerViz"`.
#' @param label_top_n Integer or `NULL`. Number of top motifs to label.
#'   Default `NULL` = auto-label high-confidence motifs (p <= 1e-50).
#'   Set to `0` to suppress all labels.
#' @param size_by_fe Logical. If `TRUE`, dot size encodes fold enrichment.
#'   If `FALSE` (default), all dots are drawn at the same fixed size.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' plot_denovo_scatter(homer)
#' plot_denovo_scatter(homer, label_top_n = 3)
#' plot_denovo_scatter(homer, size_by_fe = FALSE)
plot_denovo_scatter <- function(
    denovo,
    p_cutoff = 0.05,
    subtitle = NULL,
    color_scheme = c("homerViz", "red", "green", "blue", "purple", "viridis"),
    label_top_n = NULL,
    size_by_fe = FALSE
) {
    color_scheme <- match.arg(color_scheme)
    .check_p_cutoff(p_cutoff)
    if (inherits(denovo, "homer_output")) {
        if (is.null(subtitle)) {
            subtitle <- .build_subtitle(denovo$input, denovo$denovo_db)
        }
        denovo <- denovo$denovo
    }
    if (is.null(denovo) || nrow(denovo) == 0) stop("No de novo motifs found.")

    df <- denovo[!is.na(denovo$ln_pvalue) & denovo$ln_pvalue <= log(p_cutoff), ]
    if (size_by_fe) df <- .drop_na_fe(df, "de novo motif")

    if (nrow(df) == 0) {
        warning("No significant motifs found at p_cutoff = ", p_cutoff)
        return(NULL)
    }

    df$label <- ifelse(
        !is.na(df$consensus) & nzchar(df$consensus),
        df$consensus,
        paste0("De novo ", df$rank)
    )
    label_df <- .resolve_label_df(df, label_top_n)

    .enrichment_scatter(
        df           = df,
        label_df     = label_df,
        title        = bquote(bolditalic("De novo") ~ bold("motif")),
        subtitle     = subtitle,
        color_scheme = color_scheme,
        size_by_fe   = size_by_fe
    )
}


# ============================================================================
# Internal helpers
# ============================================================================

# Drop rows with NA fold_enrichment and print a message if any are removed.
.drop_na_fe <- function(df, type = "motif") {
    na_rows <- is.na(df$fold_enrichment)
    if (any(na_rows)) {
        n <- sum(na_rows)
        ranks <- df$rank[na_rows]
        message(
            n, " ", type, ngettext(n, "", "s"),
            " removed (fold enrichment NA, zero background occurrences), rank",
            ngettext(length(ranks), "", "s"), ": ",
            paste(ranks, collapse = ", "), "."
        )
        df <- df[!na_rows, ]
    }
    df
}

# Validate p_cutoff is a single numeric in [0, 1].
# 0 is allowed: log(0) = -Inf means nothing passes, returning NULL with a
# warning.
.check_p_cutoff <- function(p_cutoff) {
    if (!is.numeric(p_cutoff) || length(p_cutoff) != 1 || is.na(p_cutoff) ||
        p_cutoff < 0 || p_cutoff > 1) {
        stop(
            "`p_cutoff` must be a single numeric value in [0, 1]. Got: ",
            p_cutoff
        )
    }
}

# Build a two-line subtitle from BED file name and motif database name.
# Returns NULL if both inputs are NULL/empty.
.build_subtitle <- function(input, db) {
    parts <- c(
        if (!is.null(input) && nzchar(input))
            paste0("Input: ", input) else NULL,
        if (!is.null(db) && nzchar(db)) paste0("Motif database: ", db) else NULL
    )
    if (length(parts) == 0) {
        return(NULL)
    }
    paste(parts, collapse = "\n")
}

# Resolve which rows of df to label based on label_top_n:
#   NULL  -> auto: rows where neg_log10p >= 50 (p <= 1e-50), capped at 5
#   0     -> no labels (empty data frame)
#   N > 0 -> top N rows by significance (ascending ln_pvalue)
.resolve_label_df <- function(df, label_top_n) {
    if (!is.null(label_top_n)) {
        label_top_n <- as.integer(label_top_n)
        if (is.na(label_top_n) || label_top_n < 0) {
            stop("`label_top_n` must be a non-negative integer or NULL.")
        }
        if (label_top_n == 0) {
            return(df[integer(0), ])
        }
        return(head(df[order(df$ln_pvalue), ], label_top_n))
    }
    auto <- df[!is.na(df$neg_log10p) & df$neg_log10p >= 50, ]
    head(auto[order(auto$ln_pvalue), ], 5)
}

.color_scale_continuous <- function(color_scheme, legend_name) {
    palettes <- list(
        homerViz = c(low = "#e2e8f0", high = "#134e4a"),
        red = c(low = "#FEE5D9", high = "#A50F15"),
        blue = c(low = "#DEEBF7", high = "#08519C"),
        purple = c(low = "#EFEDF5", high = "#6A51A3"),
        green = c(low = "#E5F5E0", high = "#006D2C")
    )

    if (color_scheme == "viridis") {
        ggplot2::scale_color_viridis_c(
            name      = legend_name,
            option    = "viridis",
            direction = -1,
            begin     = 0.1,
            end       = 0.95
        )
    } else {
        pal <- palettes[[color_scheme]]
        ggplot2::scale_color_gradient(
            low = pal[["low"]], high = pal[["high"]],
            name = legend_name
        )
    }
}


.homer_theme <- function(plot_margin = ggplot2::margin(15, 15, 15, 15)) {
    ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            plot.title        = ggplot2::element_text(
                face = "bold", color = "#134e4a"
            ),
            plot.subtitle     = ggplot2::element_text(
                size = 9, color = "#475569"
            ),
            axis.text         = ggplot2::element_text(
                size = 9, color = "#475569"
            ),
            axis.title        = ggplot2::element_text(
                size = 10, face = "bold"
            ),
            panel.grid.major  = ggplot2::element_line(
                color = "#e2e8f0", linewidth = 0.4
            ),
            panel.grid.minor  = ggplot2::element_blank(),
            legend.position   = "right",
            legend.title      = ggplot2::element_text(
                size = 9.5, face = "bold"
            ),
            legend.text       = ggplot2::element_text(size = 8.5),
            legend.key.height = ggplot2::unit(0.8, "cm"),
            legend.key.width  = ggplot2::unit(0.4, "cm"),
            plot.margin       = plot_margin,
            panel.background  = ggplot2::element_rect(
                fill = "white", color = NA
            ),
            plot.background   = ggplot2::element_rect(
                fill = "white", color = NA
            )
        )
}


.enrichment_dotplot <- function(df, y_var, title, subtitle, color_scheme) {
    if (nrow(df) == 0) {
        warning("No motifs to plot.")
        return(NULL)
    }
    x_hi_raw <- max(df$neg_log10p, na.rm = TRUE)
    x_hi <- if (is.finite(x_hi_raw)) x_hi_raw * 1.1 else 1

    ggplot2::ggplot(df, ggplot2::aes(
        x     = .data[["neg_log10p"]],
        y     = .data[[y_var]],
        size  = .data[["fold_enrichment"]],
        color = .data[["fold_enrichment"]]
    )) +
        ggplot2::geom_point(alpha = 0.85) +
        .color_scale_continuous(color_scheme, "Fold enrichment") +
        ggplot2::scale_size_continuous(
            name   = "Fold enrichment",
            range  = c(1, 10),
            breaks = scales::pretty_breaks(n = 4)
        ) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.02, 0.1)),
            breaks = scales::pretty_breaks(n = 6)
        ) +
        ggplot2::labs(
            x        = expression(bold(-log)[bold(10)](bold(p))),
            y        = NULL,
            title    = title,
            subtitle = subtitle
        ) +
        .homer_theme(plot_margin = ggplot2::margin(15, 30, 15, 15)) +
        ggplot2::theme(
            axis.text.y        = ggplot2::element_text(
                size = 10, color = "#1e293b"
            ),
            axis.title.x       = ggplot2::element_text(
                margin = ggplot2::margin(t = 8)
            ),
            panel.grid.major.y = ggplot2::element_line(
                color = "#f1f5f9", linewidth = 0.4
            )
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend("Fold enrichment"),
            size  = ggplot2::guide_legend("Fold enrichment")
        ) +
        ggplot2::coord_cartesian(xlim = c(0, x_hi), clip = "off")
}


.enrichment_scatter <- function(
    df, label_df, title, subtitle, color_scheme,
    size_by_fe = TRUE
) {
    base_aes <- if (size_by_fe) {
        ggplot2::aes(
            x     = .data[["pct_background"]],
            y     = .data[["pct_target"]],
            size  = .data[["fold_enrichment"]],
            color = .data[["neg_log10p"]],
            label = .data[["label"]]
        )
    } else {
        ggplot2::aes(
            x     = .data[["pct_background"]],
            y     = .data[["pct_target"]],
            color = .data[["neg_log10p"]],
            label = .data[["label"]]
        )
    }

    p <- ggplot2::ggplot(df, base_aes) +
        ggplot2::geom_abline(
            slope = 1, intercept = 0, linetype = "dashed",
            color = "#94a3b8", linewidth = 0.4
        )

    if (size_by_fe) {
        p <- p +
            ggplot2::geom_point(alpha = 0.85) +
            ggplot2::scale_size_continuous(
                name   = "Fold enrichment",
                range  = c(1, 10),
                breaks = scales::pretty_breaks(n = 4)
            ) +
            ggplot2::guides(
                size = ggplot2::guide_legend("Fold enrichment")
            )
    } else {
        p <- p + ggplot2::geom_point(alpha = 0.85, size = 1.5)
    }

    p <- p + .color_scale_continuous(
        color_scheme,
        expression(bold(-log)[bold(10)](bold(p)))
    )

    p <- p +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.05)),
            breaks = scales::pretty_breaks(n = 5)
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.05)),
            breaks = scales::pretty_breaks(n = 5)
        ) +
        ggplot2::labs(
            x        = "% background sequences",
            y        = "% target sequences",
            title    = title,
            subtitle = subtitle
        ) +
        .homer_theme()

    p <- p + ggplot2::coord_cartesian(clip = "off")

    # Add labels for selected motifs
    if (!is.null(label_df) && nrow(label_df) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data               = label_df,
            ggplot2::aes(label = .data[["label"]]),
            size               = 2.5,
            color              = "#1e293b",
            box.padding        = 0.3,
            max.overlaps       = Inf,
            show.legend        = FALSE,
            segment.color      = "#94a3b8",
            segment.size       = 0.4,
            segment.linetype   = "solid",
            min.segment.length = 0
        )
    }

    p
}
