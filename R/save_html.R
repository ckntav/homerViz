# ============================================================================
# HTML report generation
# ============================================================================

#' Save a beautiful HTML report of known motif enrichment results
#'
#' Generates a self-styled HTML file with an interactive, sortable DataTable
#' showing known motif enrichment results. Motif logos are embedded inline from
#' the HOMER output directory. The page requires an internet connection for
#' Bootstrap/DataTables CDN assets.
#'
#' @param results A `homer_output` object from [read_homer_output()], or
#'   its `$known` data frame. `homer_dir` is only needed for a data frame.
#' @param homer_dir Path to the HOMER output directory (for finding SVG logos).
#' @param output_file Path to the output HTML file. Defaults to
#'   `knownResults_homerViz.html` inside `output_dir`.
#' @param output_dir Directory to write the HTML file. Defaults to `homer_dir`.
#' @param title Character. Report title shown in the page header.
#'
#' @return The output file path, invisibly.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' out <- save_known_html(homer, output_dir = tempdir())
#' file.exists(out)
save_known_html <- function(
    results,
    homer_dir = NULL,
    output_file = NULL,
    output_dir = NULL,
    title = "Known motif enrichment results"
) {
    params <- NULL
    has_denovo <- FALSE
    if (inherits(results, "homer_output")) {
        if (is.null(homer_dir)) homer_dir <- results$dir
        params <- results$params
        has_denovo <- !is.null(results$denovo)
        results <- results$known
    }
    if (is.null(homer_dir)) stop("homer_dir must be provided.")
    if (is.null(results)) stop("No known motifs found.")
    if (is.null(output_dir)) output_dir <- homer_dir
    if (!dir.exists(output_dir)) {
        stop("output_dir does not exist: ", output_dir)
    }
    # Compute de novo link filename from existing homerResults HTML
    denovo_ex <- sort(Sys.glob(file.path(homer_dir, "*homerResults*.html")))
    denovo_ex <- denovo_ex[!grepl("_homerViz\\.html$", denovo_ex)]
    denovo_stem <- if (length(denovo_ex) > 0) tools::file_path_sans_ext(basename(denovo_ex[1])) else "homerResults"
    denovo_link <- paste0(denovo_stem, "_homerViz.html")
    if (is.null(output_file)) {
        existing <- sort(Sys.glob(file.path(homer_dir, "*knownResults*.html")))
        existing <- existing[!grepl("_homerViz\\.html$", existing)]
        stem <- if (length(existing) > 0) tools::file_path_sans_ext(basename(existing[1])) else "knownResults"
        output_file <- file.path(output_dir, paste0(stem, "_homerViz.html"))
    }

    # Extract input name for subtitle
    p_parsed <- .parse_homer_params(params)
    input_subtitle <- if (!is.null(p_parsed) && !is.na(p_parsed$input %||% NA)) {
        p_parsed$input
    } else {
        NULL
    }
    dir_name <- basename(normalizePath(homer_dir, mustWork = FALSE))

    n_target <- results$n_target_total[1]
    n_bg <- results$n_background_total[1]
    n_total <- nrow(results)
    n_hc <- sum(results$neg_log10p >= 50, na.rm = TRUE)

    df <- results[order(results$ln_pvalue), ] # ascending = most negative = most significant
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(0.01), ]
    if (nrow(df) == 0) {
        warning("No significant known motifs (p <= 0.01) to include in the report.")
        return(invisible(NULL))
    }
    max_logp <- max(df$neg_log10p, na.rm = TRUE) # max of -log10(p), positive

    # ---- Build table rows ----
    rows <- paste(vapply(seq_len(nrow(df)), function(i) {
        # rank is the original file row number -> knownResults/known{N}.logo.svg
        svg_html <- .read_svg(file.path(
            homer_dir, "knownResults",
            paste0("known", df$rank[i], ".logo.svg")
        ))

        neg_log10p <- df$neg_log10p[i] # positive: -log10(p), pre-computed
        bar_pct <- round(min(neg_log10p / max_logp * 100, 100), 1)
        bar_color <- .sig_color(neg_log10p, max_logp)
        fe_badge <- .fold_badge(df$fold_enrichment[i])

        pval_str <- if (!is.na(df$pvalue_str[i])) df$pvalue_str[i] else "\u2014"
        pval_order <- if (!is.na(df$ln_pvalue[i])) df$ln_pvalue[i] else 0
        qval_order <- if (is.na(df$qvalue[i])) "1" else as.character(df$qvalue[i])
        qval_str <- if (is.na(df$qvalue[i])) "\u2014" else sprintf("%.4g", df$qvalue[i])

        bg_bar_pct <- round(min(df$pct_background[i], 100), 1)
        motif_len <- df$motif_length[i]
        motif_id_html <- sprintf('<span class="best-tf">%s</span>', .he(df$motif_id[i]))
        ln_p_str <- if (!is.na(df$ln_pvalue[i])) sprintf("%.2f", df$ln_pvalue[i]) else "\u2014"
        is_hc <- !is.na(neg_log10p) && neg_log10p >= 50
        tr_open <- if (is_hc) '<tr class="hc-row">' else "<tr>"
        rank_html <- if (is_hc) sprintf('%d<span class="hc-star">*</span>', i) else as.character(i)

        sprintf(
            '%s
        <td class="text-center rank-cell" data-order="%d">%s</td>
        <td class="logo-cell">%s</td>
        <td class="motif-id-cell" title="%s">%s</td>
        <td class="text-center mono">%d</td>
        <td class="mono text-muted">%s</td>
        <td class="text-end mono text-nowrap" data-order="%s">%s</td>
        <td class="text-end mono">%s</td>
        <td data-order="%s" class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill" style="width:%s%%;background:%s"></div>
            </div>
            <span class="bar-label">%s</span>
        </td>
        <td class="text-end mono" data-order="%s">%s</td>
        <td class="text-center">%s</td>
        <td class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill bg-tgt" style="width:%s%%"></div>
            </div>
            <span class="bar-label">%s%%</span>
        </td>
        <td class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill bg-bg" style="width:%s%%"></div>
            </div>
            <span class="bar-label">%s%%</span>
        </td>
        <td class="text-end mono" data-order="%s">%s</td>
        <td class="text-end mono" data-order="%s">%s</td>
            </tr>',
            tr_open,
            i, rank_html,
            svg_html,
            .he(df$motif_id[i]), motif_id_html,
            motif_len,
            .he(df$consensus[i]),
            pval_order, pval_str,
            ln_p_str,
            if (is.na(neg_log10p)) "0" else sprintf("%.4f", neg_log10p),
            bar_pct, bar_color,
            .fmt_logp(neg_log10p),
            qval_order, qval_str,
            fe_badge,
            min(df$pct_target[i], 100),
            round(df$pct_target[i], 1),
            bg_bar_pct,
            round(df$pct_background[i], 1),
            as.integer(df$n_target[i]),     format(as.integer(df$n_target[i]), big.mark = ","),
            as.integer(df$n_background[i]), format(as.integer(df$n_background[i]), big.mark = ",")
        )
    }, character(1)), collapse = "\n")

    # ---- Assemble page ----
    stats_html <- sprintf(
        "High confidence motifs (p&nbsp;&le;&nbsp;1e&#x2212;50): <strong>%s / %s</strong><br>Significant motifs (p&nbsp;&le;&nbsp;0.01): <strong>%s / %s</strong><br>Target sequences: <strong>%s</strong><br>Background sequences: <strong>%s</strong>",
        format(n_hc, big.mark = ","), format(n_total, big.mark = ","),
        format(nrow(df), big.mark = ","), format(n_total, big.mark = ","),
        format(n_target, big.mark = ","), format(n_bg, big.mark = ",")
    )

    thead <- '
    <thead>
        <tr>
        <th>Rank</th>
        <th>Motif</th>
        <th>Motif ID</th>
        <th class="stat-header">Motif<br>length</th>
        <th>Consensus</th>
        <th class="stat-header" style="white-space:nowrap">p-value</th>
        <th class="stat-header">ln(p)</th>
        <th class="stat-header">&#8722;log&#8321;&#8320;(p)</th>
        <th class="stat-header" style="white-space:nowrap">q-value</th>
        <th class="stat-header">Fold<br>enrichment</th>
        <th class="col-count stat-header">%<br>Target</th>
        <th class="col-count stat-header">%<br>Bckgrd</th>
        <th class="col-count stat-header">#<br>Target</th>
        <th class="col-count stat-header">#<br>Bckgrd</th>
        </tr>
    </thead>'

    html <- .html_page(
        title           = title,
        subtitle        = input_subtitle,
        stats_html      = stats_html,
        thead           = thead,
        tbody           = rows,
        link_other      = "<em>De novo</em> motif results",
        link_other_file = denovo_link,
        show_link_other = has_denovo,
        table_id        = "knownTable",
        extra_html      = .params_html(params, type = "known", dir_name = dir_name),
        col_defs        = '[
        { "orderable": false, "targets": [1] },
        { "type": "num", "targets": [5, 6, 7, 8] },
        { "className": "col-xs", "targets": [5, 6, 7, 8, 9, 10, 11, 12, 13] },
        { "visible": false, "targets": [6, 8, 12, 13] },
        { "width":  "45px", "targets": [0] },
        { "width": "300px", "targets": [1] },
        { "width": "160px", "targets": [2] },
        { "width":  "60px", "targets": [3] },
        { "width": "110px", "targets": [4] },
        { "width":  "70px", "targets": [5] },
        { "width":  "70px", "targets": [6] },
        { "width":  "70px", "targets": [7] },
        { "width":  "70px", "targets": [8] },
        { "width":  "70px", "targets": [9] },
        { "width":  "85px", "targets": [10] },
        { "width":  "85px", "targets": [11] },
        { "width":  "85px", "targets": [12] },
        { "width":  "85px", "targets": [13] }
    ]',
        cmd_footer      = params,
        show_csv        = FALSE,
        order_col       = 7
    )

    writeLines(html, output_file, useBytes = TRUE)
    message("Saved: ", output_file)
    invisible(output_file)
}


#' Save a beautiful HTML report of de novo motif results
#'
#' Generates a self-styled HTML file with an interactive, sortable DataTable
#' showing de novo motif enrichment results. Logos are embedded inline from
#' the HOMER output directory. Best-match TF names are extracted from the
#' existing `homerResults.html` if present.
#'
#' @param denovo A `homer_output` object from [read_homer_output()], or
#'   its `$denovo` data frame.
#' @param homer_dir Path to the HOMER output directory.
#' @param output_file Output HTML path. Defaults to
#'   `homerResults_homerViz.html` inside `output_dir`.
#' @param output_dir Directory to write the HTML file. Defaults to `homer_dir`.
#' @param title Character. Report title.
#'
#' @return The output file path, invisibly.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' out <- save_denovo_html(homer, output_dir = tempdir())
#' file.exists(out)
save_denovo_html <- function(
    denovo,
    homer_dir = NULL,
    output_file = NULL,
    output_dir = NULL,
    title = "<em>De novo</em> motif results"
) {
    params <- NULL
    has_known <- FALSE
    if (inherits(denovo, "homer_output")) {
        if (is.null(homer_dir)) homer_dir <- denovo$dir
        params <- denovo$params
        has_known <- !is.null(denovo$known)
        denovo <- denovo$denovo # data frame
    }
    if (is.null(homer_dir)) stop("homer_dir must be provided.")
    if (is.null(denovo) || !is.data.frame(denovo) || nrow(denovo) == 0) {
        stop("No de novo motifs found. Run read_homer_output() first.")
    }
    if (is.null(output_dir)) output_dir <- homer_dir
    if (!dir.exists(output_dir)) {
        stop("output_dir does not exist: ", output_dir)
    }

    existing <- sort(Sys.glob(file.path(homer_dir, "*homerResults*.html")))
    existing <- existing[!grepl("_homerViz\\.html$", existing)]
    if (is.null(output_file)) {
        stem <- if (length(existing) > 0) tools::file_path_sans_ext(basename(existing[1])) else "homerResults"
        output_file <- file.path(output_dir, paste0(stem, "_homerViz.html"))
    }

    known_ex <- sort(Sys.glob(file.path(homer_dir, "*knownResults*.html")))
    known_ex <- known_ex[!grepl("_homerViz\\.html$", known_ex)]
    known_stem <- if (length(known_ex) > 0) tools::file_path_sans_ext(basename(known_ex[1])) else "knownResults"
    known_link <- paste0(known_stem, "_homerViz.html")

    p_parsed_d <- .parse_homer_params(params)
    input_subtitle <- if (!is.null(p_parsed_d) && !is.na(p_parsed_d$input %||% NA)) {
        p_parsed_d$input
    } else {
        NULL
    }
    dir_name <- basename(normalizePath(homer_dir, mustWork = FALSE))

    n_total <- nrow(denovo)
    n_hc <- sum(denovo$neg_log10p >= 50, na.rm = TRUE)
    n_target <- denovo$n_target_total[1]
    n_bg_seq <- denovo$n_background_total[1]

    # Sort and (optionally) subset - ascending = most negative ln_pvalue = most significant
    df <- denovo[order(denovo$ln_pvalue), ]
    df <- df[!is.na(df$ln_pvalue) & df$ln_pvalue <= log(0.01), ]
    if (nrow(df) == 0) {
        warning("No significant de novo motifs (p <= 0.01) to include in the report.")
        return(invisible(NULL))
    }
    max_logp <- max(df$neg_log10p, na.rm = TRUE) # max of -log10(p), positive

    # ---- Build rows ----
    rows <- paste(vapply(seq_len(nrow(df)), function(i) {
        neg_log10p <- df$neg_log10p[i] # positive: -log10(p), pre-computed
        ln_p_str <- if (!is.na(df$ln_pvalue[i])) sprintf("%.2f", df$ln_pvalue[i]) else "\u2014"
        pval_str <- if (!is.na(df$pvalue_str[i])) df$pvalue_str[i] else "\u2014"
        pval_order <- if (!is.na(df$ln_pvalue[i])) df$ln_pvalue[i] else 0
        motif_len <- df$motif_length[i]
        consensus <- if (!is.na(df$consensus[i])) df$consensus[i] else ""
        pct_tgt <- df$pct_target[i]
        pct_bg <- df$pct_background[i]
        fe <- df$fold_enrichment[i]
        best_name <- df$best_match_name[i]
        best_score <- df$best_match_score[i]
        n_tgt <- df$n_target[i]
        n_bg_row <- df$n_background[i] # HOMER-normalized; may be decimal
        tgt_s <- df$pos_std_target[i]
        bg_s <- df$pos_std_bg[i]

        svg_html <- .read_svg(file.path(
            homer_dir, "homerResults",
            paste0("motif", df$motif_num[i], ".logo.svg")
        ))

        best_html <- if (!is.na(best_name) && nzchar(best_name)) {
            score_line <- if (!is.na(best_score)) {
                sprintf(
                    '<br><small class="text-muted">score: %s</small>',
                    .he(sprintf("%.3f", best_score))
                )
            } else {
                ""
            }
            sprintf('<span class="best-tf">%s</span>%s', .he(best_name), score_line)
        } else {
            "<span class='text-muted small'>\u2014</span>"
        }

        std_html <- if (!is.na(tgt_s)) {
            bg_part <- if (!is.na(bg_s)) bg_s else "\u2014"
            sprintf(
                '<span class="mono">%s</span><br><small class="text-muted mono">bg: %s</small>',
                .he(tgt_s), .he(bg_part)
            )
        } else {
            "<span class='text-muted small'>\u2014</span>"
        }

        is_hc <- !is.na(neg_log10p) && neg_log10p >= 50
        tr_open <- if (is_hc) '<tr class="hc-row">' else "<tr>"
        rank_html <- if (is_hc) sprintf('%d<span class="hc-star">*</span>', i) else as.character(i)
        bar_pct <- if (!is.na(neg_log10p)) round(min(neg_log10p / max_logp * 100, 100), 1) else 0

        bar_color <- if (!is.na(neg_log10p)) .sig_color(neg_log10p, max_logp) else "#ccc"
        tgt_bar_pct <- if (!is.na(pct_tgt)) round(min(pct_tgt, 100), 1) else 0
        bg_bar_pct <- if (!is.na(pct_bg)) round(min(pct_bg, 100), 1) else 0
        tgt_label <- if (!is.na(pct_tgt)) paste0(round(pct_tgt, 1), "%") else "\u2014"
        bg_label <- if (!is.na(pct_bg)) paste0(round(pct_bg, 1), "%") else "\u2014"
        # n_background is HOMER-normalized (can be decimal); round for display
        n_tgt_str <- if (!is.na(n_tgt)) format(as.integer(round(n_tgt)), big.mark = ",") else "\u2014"
        n_bg_str <- if (!is.na(n_bg_row)) format(round(n_bg_row, 1), big.mark = ",") else "\u2014"

        sprintf(
            '%s
        <td class="text-center rank-cell" data-order="%d">%s</td>
        <td class="logo-cell">%s</td>
        <td>%s</td>
        <td class="text-center mono">%d</td>
        <td class="mono text-muted">%s</td>
        <td class="text-end mono text-nowrap" data-order="%s">%s</td>
        <td class="text-end mono">%s</td>
        <td data-order="%s" class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill" style="width:%s%%;background:%s"></div>
            </div>
            <span class="bar-label">%s</span>
        </td>
        <td class="text-center">%s</td>
        <td class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill bg-tgt" style="width:%s%%"></div>
            </div>
            <span class="bar-label">%s</span>
        </td>
        <td class="bar-cell">
            <div class="bar-wrap">
            <div class="bar-fill bg-bg" style="width:%s%%"></div>
            </div>
            <span class="bar-label">%s</span>
        </td>
        <td class="text-end mono" data-order="%s">%s</td>
        <td class="text-end mono" data-order="%s">%s</td>
        <td class="text-end mono text-nowrap small">%s</td>
            </tr>',
            tr_open,
            i, rank_html,
            svg_html,
            best_html,
            motif_len,
            .he(consensus),
            pval_order, pval_str,
            ln_p_str,
            if (is.na(neg_log10p)) "0" else sprintf("%.4f", neg_log10p),
            bar_pct, bar_color,
            .fmt_logp(neg_log10p),
            .fold_badge(fe),
            tgt_bar_pct, tgt_label,
            bg_bar_pct, bg_label,
            if (is.na(n_tgt)) 0 else n_tgt, n_tgt_str,
            if (is.na(n_bg_row)) 0 else n_bg_row, n_bg_str,
            std_html
        )
    }, character(1)), collapse = "\n")

    stats_html <- sprintf(
        "High confidence motifs (p&nbsp;&le;&nbsp;1e&#x2212;50): <strong>%s</strong><br>Significant motifs (p&nbsp;&le;&nbsp;0.01): <strong>%s</strong><br>Target sequences: <strong>%s</strong><br>Background sequences: <strong>%s</strong>",
        format(n_hc, big.mark = ","),
        format(nrow(df), big.mark = ","),
        format(as.integer(n_target), big.mark = ","),
        format(as.integer(n_bg_seq), big.mark = ",")
    )

    thead <- '
    <thead>
        <tr>
        <th>Rank</th>
        <th>De novo motif</th>
        <th>Best match</th>
        <th class="stat-header">Motif<br>length</th>
        <th>Consensus</th>
        <th class="stat-header" style="white-space:nowrap">p-value</th>
        <th class="stat-header">ln(p)</th>
        <th class="stat-header">&#8722;log&#8321;&#8320;(p)</th>
        <th class="stat-header">Fold<br>enrichment</th>
        <th class="stat-header">% Target</th>
        <th class="stat-header">%&nbsp;Bckgrd</th>
        <th class="col-count stat-header">#<br>Target</th>
        <th class="col-count stat-header">#<br>Bckgrd</th>
        <th class="stat-header">Position<br>STD</th>
        </tr>
    </thead>'

    html <- .html_page(
        title = title,
        subtitle = input_subtitle,
        stats_html = stats_html,
        thead = thead,
        tbody = rows,
        link_other = "Known motif enrichment results",
        link_other_file = known_link,
        show_link_other = has_known,
        table_id = "denovoTable",
        extra_html = .params_html(params, type = "denovo", dir_name = dir_name),
        cmd_footer = params,
        col_defs = '[
        { "orderable": false, "targets": [1] },
        { "type": "num", "targets": [5, 6, 7, 11, 12] },
        { "className": "col-xs", "targets": [5, 6, 7, 8, 9, 10, 11, 12] },
        { "width":  "45px", "targets": [0] },
        { "width": "300px", "targets": [1] },
        { "width": "160px", "targets": [2] },
        { "width":  "60px", "targets": [3] },
        { "width": "110px", "targets": [4] },
        { "width":  "80px", "targets": [5] },
        { "width":  "70px", "targets": [6] },
        { "width":  "90px", "targets": [7] },
        { "width":  "80px", "targets": [8] },
        { "width":  "80px", "targets": [9] },
        { "width":  "80px", "targets": [10] },
        { "width":  "85px", "targets": [11] },
        { "width":  "85px", "targets": [12] },
        { "width":  "90px", "targets": [13] },
        { "visible": false, "targets": [6, 11, 12, 13] }
    ]',
        order_col = 7,
        show_csv = FALSE
    )

    writeLines(html, output_file, useBytes = TRUE)
    message("Saved: ", output_file)
    invisible(output_file)
}


#' Generate both HTML reports for a HOMER output directory
#'
#' Convenience wrapper that calls [save_known_html()] and [save_denovo_html()].
#'
#' @param homer A `homer_output` object from [read_homer_output()].
#' @param output_dir Directory to save the HTML files. Defaults to `homer$dir`.
#'
#' @return Named list of output file paths, invisibly.
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' out <- save_homer_html(homer, output_dir = tempdir())
#' file.exists(out$known)
save_homer_html <- function(homer, output_dir = NULL) {
    if (!inherits(homer, "homer_output")) stop("`homer` must be a homer_output object.")
    if (is.null(output_dir)) output_dir <- homer$dir

    out <- list()

    if (!is.null(homer$known)) {
        out$known <- save_known_html(homer, output_dir = output_dir)
    }
    if (!is.null(homer$denovo)) {
        out$denovo <- save_denovo_html(homer, output_dir = output_dir)
    }
    if (length(out) == 0) message("Nothing to save: no known or de novo results found.")
    invisible(out)
}


# ============================================================================
# Internal HTML helpers
# ============================================================================

.he <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
}

.read_svg <- function(path) {
    if (!file.exists(path)) {
        return('<span class="text-muted small">&ndash;</span>')
    }
    svg <- paste(readLines(path, warn = FALSE), collapse = "\n")
    # No explicit dimensions; CSS constrains max-height and max-width so the browser
    # scales proportionally: shorter motifs are narrower, nothing exceeds 50px tall.
    svg <- sub(
        '<svg width="([0-9]+)" height="([0-9]+)">',
        '<svg style="display:block;height:auto;width:auto;max-height:50px;max-width:100%" viewBox="0 0 \\1 \\2" preserveAspectRatio="xMinYMid meet">',
        svg
    )
    svg
}

.fmt_logp <- function(x) {
    if (is.na(x)) {
        return("&ndash;")
    }
    if (x > 999) {
        return(sprintf("%.0f", x))
    }
    sprintf("%.1f", x)
}

.sig_color <- function(log10p, max_log10p) {
    frac <- min(log10p / max(max_log10p, 1), 1)
    if (frac > 0.75) {
        "rgba(215,154,54,1.0)"
    } else if (frac > 0.50) {
        "rgba(215,154,54,0.7)"
    } else if (frac > 0.25) {
        "rgba(215,154,54,0.4)"
    } else {
        "rgba(215,154,54,0.2)"
    }
}

.fold_badge <- function(fe) {
    if (is.na(fe) || is.null(fe)) {
        return("<span class='badge bg-secondary'>&ndash;</span>")
    }
    fe <- round(fe, 1)
    cls <- if (fe >= 5) {
        "badge-high"
    } else if (fe >= 2) {
        "badge-mid"
    } else {
        "badge-low"
    }
    sprintf('<span class="badge %s">%.1f&times;</span>', cls, fe)
}

.html_page <- function(
    title, stats_html, thead, tbody, link_other,
    link_other_file, table_id, extra_html = NULL,
    col_defs = NULL, cmd_footer = NULL,
    show_csv = TRUE, show_link_other = TRUE,
    order_col = 5, subtitle = NULL
) {
    col_defs_js <- if (!is.null(col_defs)) {
        col_defs
    } else {
        '[{ "orderable": false, "targets": [1] }, { "type": "num", "targets": [4] }]'
    }

    subtitle_html <- if (!is.null(subtitle) && nzchar(as.character(subtitle))) {
        sprintf('<div class="page-subtitle">%s</div>', .he(as.character(subtitle)))
    } else {
        ""
    }

    footer_cmd <- if (!is.null(cmd_footer)) {
        sprintf(
            '<div class="cmd-footer"><span class="cmd-footer-label">Command line:</span> <code>findMotifsGenome.pl %s</code></div>',
            .he(cmd_footer)
        )
    } else {
        ""
    }

    cmd_html <- extra_html %||% ""

    link_html <- if (isTRUE(show_link_other)) {
        sprintf(
            '<a href="%s" class="header-btn">&rarr;&nbsp;%s</a>',
            .he(link_other_file), link_other
        )
    } else {
        ""
    }

    dt_dom <- paste0(
        "'<\"row align-items-center mb-3\"<\"col-auto\"f><\"col-auto ms-auto\"B>>",
        "<\"row\"<\"col-12\"t>>",
        "<\"row mt-2 align-items-center\"<\"col-auto\"l><\"col-auto ms-auto\"p>>'"
    )

    colvis_btn <- '{ extend: "colvis", text: "show/hide columns", className: "btn btn-sm btn-colvis" }'
    csv_btn <- '{ extend: "csv", text: "&#x21E9; CSV", className: "btn btn-sm btn-outline-secondary" }'
    btns <- if (show_csv) paste(colvis_btn, csv_btn, sep = ",\n      ") else colvis_btn
    dt_buttons_js <- sprintf(",\n    buttons: [\n      %s\n    ]", btns)

    css_block <- '
    :root { --accent: #3b82f6; --accent-dark: #1d4ed8; }
    body { background: #f1f5f9; font-family: "Inter", system-ui, sans-serif; scrollbar-gutter: stable; }
    .page-header {
        background: linear-gradient(135deg, #134e4a 0%, #0d3b38 100%);
        color: #fff; padding: 2rem 2.5rem 1.5rem;
        margin-bottom: 1.5rem; position: relative;
    }
    .page-header h1 { font-size: 1.6rem; font-weight: 700; margin-bottom: 0.1rem; }
    .page-subtitle { font-size: 1.1rem; font-weight: 500; color: #fff; margin-bottom: 0.4rem; }
    .page-header .stats { opacity: 0.75; font-size: 0.9rem; }
    .header-btn {
        position: absolute; bottom: 1.25rem; right: 2.5rem;
        font-size: 0.82rem; padding: 0.35rem 0.9rem;
        border: 1.5px solid rgba(255,255,255,0.45); border-radius: 6px;
        color: rgba(255,255,255,0.85); text-decoration: none; transition: all .15s;
    }
    .header-btn:hover { background: rgba(255,255,255,0.12); color: #fff; border-color: rgba(255,255,255,0.7); }
    .card { border: none; border-radius: 12px; box-shadow: 0 1px 4px rgba(0,0,0,.1); }
    .card-body { padding: 1.5rem; }
    table.dataTable { table-layout: fixed; width: auto !important; }
    table.dataTable thead th {
        background: #134e4a; color: #e2e8f0;
        font-size: 0.75rem; text-transform: uppercase;
        letter-spacing: 0.05em; border: none !important;
        white-space: normal; max-width: 120px; line-height: 1.2; vertical-align: middle;
    }
    table.dataTable thead th.stat-header { font-size: 0.6rem; }
    table.dataTable tbody tr:hover { background: #eff6ff !important; }
    table.dataTable tbody tr.hc-row { background: #f0fdf4 !important; }
    table.dataTable tbody tr.hc-row td:first-child { border-left: 3px solid #16a34a; }
    table.dataTable tbody tr.hc-row:hover { background: #dcfce7 !important; }
    table.dataTable tbody td { vertical-align: middle; padding: 0.3rem 0.5rem; font-size: 0.75rem; }
    .rank-cell { color: #94a3b8; font-size: 0.75rem; width: 40px; }
    .hc-star { color: #16a34a; font-weight: 700; margin-left: 1px; }
    .col-count { }
    .logo-cell { min-width: 300px; height: 58px; padding: 3px 5px !important; overflow: hidden; }
    .tf-name { font-weight: 600; font-size: 1rem; color: #1e293b; }
    .mono { }
    .bar-wrap { height: 8px; border-radius: 4px; background: #e2e8f0; overflow: hidden; display: inline-block; width: 55px; vertical-align: middle; }
    .bar-fill { height: 100%; border-radius: 4px; transition: width .3s; }
    .bg-tgt  { background: rgba(215,154,54,1.0) !important; }
    .bg-bg   { background: rgba(215,154,54,0.35) !important; }
    .bar-label { font-size: 0.75rem; color: #64748b; margin-left: 6px; }
    .bar-cell .bar-wrap { display: block; margin-bottom: 3px; }
    .bar-cell .bar-label { display: block; margin-left: 0; }
    .badge-high { background: #D79A36; color: #fff; }
    .badge-mid  { background: rgba(215, 154, 54, 0.5); color: #1e293b; }
    .badge-low  { background: #94a3b8; }
    .badge-high, .badge-mid, .badge-low { font-size: 0.8rem; padding: 0.35em 0.6em; }
    .best-tf { font-weight: 600; font-size: 0.9rem; color: #134e4a; }
    .btn-colvis, .btn-colvis:focus { background: #134e4a; color: #fff; border: none; border-radius: 6px; font-size: 0.75rem; padding: 0.3rem 0.8rem; }
    .btn-colvis:hover { background: #0d3b38; color: #fff; }
    div.dt-button-collection { right: 0 !important; left: auto !important; }
    div.dt-button-collection .dt-button.active { background: #134e4a !important; color: #fff !important; border-color: #134e4a !important; }
    div.dt-button-collection .dt-button.active:hover { background: #0d3b38 !important; }
    .dataTables_wrapper .dataTables_filter { overflow: visible; padding-top: 3px; }
    .dataTables_wrapper .dataTables_filter input { border-radius: 8px; border: 1px solid #cbd5e1; padding: 0.35rem 0.75rem; font-size: 0.9rem; }
    .dataTables_wrapper .dataTables_length select { border-radius: 8px; border: 1px solid #cbd5e1; padding: 0.25rem 1.75rem 0.25rem 0.5rem; min-width: 4.5rem; }
    .dataTables_wrapper .dataTables_paginate .paginate_button.current,
    .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover { background: #134e4a !important; color: #fff !important; border-color: #134e4a !important; }
    .dataTables_wrapper .dataTables_paginate .paginate_button:not(.disabled):not(.current):hover { background: #f0fdf4 !important; color: #134e4a !important; border-color: #d1fae5 !important; }
    .dataTables_wrapper .pagination .page-item.active .page-link { background-color: #134e4a !important; border-color: #134e4a !important; color: #fff !important; }
    .dataTables_wrapper .pagination .page-link { color: #475569; }
    .dataTables_wrapper .pagination .page-link:hover { color: #134e4a; }
    .dataTables_wrapper .pagination .page-link:focus { box-shadow: 0 0 0 0.2rem rgba(19,78,74,0.25); }
    .dataTables_wrapper .dataTables_filter input:focus,
    .dataTables_wrapper .dataTables_length select:focus { border-color: #134e4a !important; box-shadow: 0 0 0 0.2rem rgba(19,78,74,0.25) !important; outline: none; }
    dt.buttons-html5 { border-radius: 8px !important; }
    .motif-id-cell { max-width: 180px; overflow: hidden; font-size: 0.75rem; }
    .powered-by { text-align: right; font-size: 0.75rem; color: #94a3b8; margin-top: 0.5rem; line-height: 1.6; }
    .powered-by a { color: #94a3b8; text-decoration: none; display: inline-flex; align-items: center; gap: 0.3rem; }
    .powered-by a:hover { color: #134e4a; }
    .powered-by svg { vertical-align: middle; }
    .cmd-footer { margin-top: 0.75rem; font-size: 0.75rem; color: #94a3b8; word-break: break-all; }
    .cmd-footer-label { font-weight: 600; margin-right: 0.25rem; }
    .cmd-footer code { color: #64748b; background: none; font-size: inherit; font-family: "SFMono-Regular", Consolas, monospace; }
    .params-block { font-size: 0.82rem; margin-top: 0.5rem; opacity: 0.9; }
    .params-grid { display: flex; flex-direction: column; gap: 0.15rem; }
    .params-grid > span { white-space: nowrap; }
    .param-label { color: rgba(255,255,255,0.6); font-weight: 600; }
    .params-block code { color: #e2e8f0; background: none; font-size: inherit; font-family: inherit; }'

    sprintf(
        '<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>%s homerViz</title>
    <link rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css"
    crossorigin="anonymous">
    <link rel="stylesheet"
    href="https://cdn.datatables.net/1.13.8/css/dataTables.bootstrap5.min.css">
    <link rel="stylesheet"
    href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.bootstrap5.min.css">
    <style>%s</style>
</head>
<body>

<div class="page-header">
    <h1>%s</h1>
    %s
    <div class="stats">%s</div>
    %s
    %s
</div>

<div class="container-fluid px-4">
    <div class="card">
    <div class="card-body">
        <div style="overflow-x: auto">
        <table id="%s" class="table table-hover table-bordered">
            %s
            <tbody>%s</tbody>
        </table>
        </div>
        %s
        <div class="powered-by">
        <div>Generated on %s</div>
        <div>with <a href="https://github.com/ckntav/homerViz" target="_blank" rel="noopener">homerViz
            <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
            <path d="M12 0C5.37 0 0 5.37 0 12c0 5.31 3.435 9.795 8.205 11.385.6.105.825-.255.825-.57 0-.285-.015-1.23-.015-2.235-3.015.555-3.795-.735-4.035-1.41-.135-.345-.72-1.41-1.23-1.695-.42-.225-1.02-.78-.015-.795.945-.015 1.62.87 1.845 1.23 1.08 1.815 2.805 1.305 3.495.99.105-.78.42-1.305.765-1.605-2.67-.3-5.46-1.335-5.46-5.925 0-1.305.465-2.385 1.23-3.225-.12-.3-.54-1.53.12-3.18 0 0 1.005-.315 3.3 1.23.96-.27 1.98-.405 3-.405s2.04.135 3 .405c2.295-1.56 3.3-1.23 3.3-1.23.66 1.65.24 2.88.12 3.18.765.84 1.23 1.905 1.23 3.225 0 4.605-2.805 5.625-5.475 5.925.435.375.81 1.095.81 2.22 0 1.605-.015 2.895-.015 3.3 0 .315.225.69.825.57A12.02 12.02 0 0 0 24 12c0-6.63-5.37-12-12-12z"/>
            </svg>
        </a></div>
        </div>
    </div>
    </div>
</div>

<script src="https://code.jquery.com/jquery-3.7.0.min.js" crossorigin="anonymous"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js" crossorigin="anonymous"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/dataTables.bootstrap5.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.bootstrap5.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js" crossorigin="anonymous"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.colVis.min.js"></script>
<script>
$(document).ready(function () {
    $("#%s").DataTable({
    pageLength: 10,
    lengthMenu: [[10, 25, 50, 100], [10, 25, 50, 100]],
    autoWidth: false,
    order: [[%d, "desc"]],
    language: {
        search: "Search motif:",
        lengthMenu: "Show _MENU_ motifs",
        info: "Showing _START_ to _END_ of _TOTAL_ motifs",
        infoEmpty: "No motifs to show",
        infoFiltered: "(filtered from _MAX_ total motifs)"
    },
    columnDefs: %s,
    dom: %s%s
    });
});
</script>
</body>
</html>',
        .he(gsub("<[^>]+>", "", title)), # <title> (plain text, strip tags)
        css_block, # <style>
        title, # <h1> (raw HTML, allows <em> etc.)
        subtitle_html, # subtitle (or "")
        stats_html, # stats line
        cmd_html, # params block (or "")
        link_html, # link (conditional)
        table_id, # table id
        thead, # <thead>
        tbody, # <tbody>
        footer_cmd, # cmd footer (or "")
        format(Sys.time(), "%d %B %Y %H:%M"), # generation date
        table_id, # DataTables init
        order_col, # DataTables sort column
        col_defs_js, # DataTables columnDefs
        dt_dom, # DataTables dom string
        dt_buttons_js # DataTables buttons (or "")
    )
}

# Parse a findMotifsGenome.pl command string into a named list of arguments.
.parse_homer_params <- function(cmd_str) {
    if (is.null(cmd_str) || !nzchar(cmd_str)) {
        return(NULL)
    }
    args <- strsplit(trimws(sub("^findMotifsGenome\\.pl\\s+", "", cmd_str)), "\\s+")[[1]]

    get_flag <- function(flag) {
        idx <- which(args == flag)
        if (length(idx) == 0) {
            return(NULL)
        }
        nxt <- idx[1] + 1
        if (nxt > length(args)) {
            return(TRUE)
        }
        # Values can start with - for ranges like -100,50; flags are -[alpha]
        if (grepl("^-[0-9]", args[nxt]) || !startsWith(args[nxt], "-")) args[nxt] else TRUE
    }
    has_flag <- function(flag) any(args == flag)

    list(
        input   = if (length(args) >= 1) basename(args[1]) else NA,
        genome  = if (length(args) >= 2) args[2] else NA,
        size    = get_flag("-size"),
        mknown  = get_flag("-mknown"),
        mcheck  = get_flag("-mcheck"),
        p       = get_flag("-p"),
        len     = get_flag("-len"),
        S       = get_flag("-S"),
        mask    = has_flag("-mask"),
        nomotif = has_flag("-nomotif"),
        noknown = has_flag("-noknown")
    )
}

# Build the params HTML block for the page header.
# type: "known" or "denovo"
.params_html <- function(cmd_str, type = "known", dir_name = NULL) {
    p <- .parse_homer_params(cmd_str) %||% list()
    if (length(p) == 0 && is.null(dir_name)) {
        return(NULL)
    }

    item <- function(label, val) {
        if (is.null(val) || isTRUE(val) || (length(val) == 1 && is.na(val))) {
            return(NULL)
        }
        sprintf(
            '<span><span class="param-label">%s</span> <code>%s</code></span>',
            label, .he(as.character(val))
        )
    }

    .db_label <- function(val) {
        if (is.null(val) || isTRUE(val) || (length(val) == 1 && is.na(val))) {
            return("Homer motif database (default)")
        }
        bn <- basename(val)
        if (grepl("all\\.motifs|known\\.motifs", bn, ignore.case = TRUE)) {
            paste0(bn, " (Homer motif database)")
        } else {
            bn
        }
    }

    if (type == "known") {
        items <- list(
            item("Output folder:", dir_name),
            item("Motif database:", .db_label(p$mknown)),
            item("Genome:", p$genome)
        )
    } else {
        items <- list(
            item("Output folder:", dir_name),
            item("Motif database:", .db_label(p$mcheck)),
            item("Genome:", p$genome)
        )
    }

    items <- Filter(Negate(is.null), items)
    if (length(items) == 0) {
        return(NULL)
    }
    sprintf(
        '<div class="params-block"><div class="params-grid">%s</div></div>',
        paste(items, collapse = "\n")
    )
}
