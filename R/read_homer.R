#' Read HOMER known motif enrichment results
#'
#' Parses the `knownResults.txt` file produced by HOMER `findMotifsGenome.pl`.
#'
#' @param dir Path to the HOMER output directory, or directly to a
#'   `knownResults.txt` file.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{rank}{Integer. Original row number in knownResults.txt,
#'       used to locate `knownResults/known{N}.logo.svg`.}
#'     \item{motif_id}{Full HOMER motif identifier
#'       (e.g. `TYY1.H14CORE.0.PSM.A(Zf)/...`)}
#'     \item{consensus}{Consensus sequence}
#'     \item{motif_length}{Integer. Length of the consensus sequence.}
#'     \item{pvalue_str}{Display p-value string (e.g. `"1e-200"`; corrected when
#'       HOMER writes `"0"` due to float underflow)}
#'     \item{ln_pvalue}{ln(p-value) as stored by HOMER, a negative number}
#'     \item{neg_log10p}{-log10(p-value), positive, computed as
#'       -ln_pvalue / log(10)}
#'     \item{qvalue}{Benjamini-corrected q-value}
#'     \item{n_target}{Number of target sequences with motif}
#'     \item{pct_target}{Percent of target sequences with motif}
#'     \item{n_background}{Number of background sequences with motif}
#'     \item{pct_background}{Percent of background sequences with motif}
#'     \item{fold_enrichment}{pct_target / pct_background
#'       (NA when pct_background = 0)}
#'     \item{n_target_total}{Total number of target sequences}
#'     \item{n_background_total}{Total number of background sequences}
#'   }
#'
#' @noRd
read_homer_known <- function(dir) {
    # Allow passing either the file or the directory
    path <- if (dir.exists(dir)) .find_known_file(dir) else dir
    if (!file.exists(path)) {
        stop("knownResults.txt not found at: ", path)
    }

    # Read header separately to extract N sequences
    header_line <- readLines(path, n = 1)
    col_names <- strsplit(header_line, "\t")[[1]]

    n_target_total <- .extract_n_sequences(col_names[6])
    n_bg_total <- .extract_n_sequences(col_names[8])

    # Read data, forcing the p-value column (col 3) to character so we capture
    # the raw string before any numeric conversion (preserves "1e-322", etc.)
    # colClasses is built dynamically to handle HOMER versions with extra
    # columns.
    n_cols <- length(
        strsplit(readLines(path, n = 2, warn = FALSE)[2], "\t")[[1]]
    )
    col_cls <- rep(NA_character_, n_cols)
    col_cls[3] <- "character" # keep raw p-value string
    raw <- utils::read.table(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment = "",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = col_cls
    )

    # Standardise the first 9 column names regardless of HOMER version.
    # HOMER's 4th column is ln(p-value), a negative number; exposed as
    # ln_pvalue.
    if (ncol(raw) < 9) {
        stop(
            "knownResults.txt has fewer than 9 columns (found ",
            ncol(raw), ")."
        )
    }
    colnames(raw)[seq_len(9)] <- c(
        "motif_id", "consensus", "pvalue", "logp_homer", "qvalue",
        "n_target", "pct_target", "n_background", "pct_background"
    )

    # Original file rank (1-based row number), used for SVG logo lookup
    # (HOMER names logos knownResults/known{N}.logo.svg where N = file row
    # order)
    raw$rank <- seq_len(nrow(raw))

    # Capture raw p-value string from file before numeric conversion
    raw$pvalue_str <- raw$pvalue

    # Clean percentage strings and coerce to numeric
    raw$pct_target <- as.numeric(sub("%$", "", raw$pct_target))
    raw$pct_background <- as.numeric(sub("%$", "", raw$pct_background))
    raw$n_target <- as.numeric(raw$n_target)
    raw$n_background <- as.numeric(raw$n_background)
    raw$logp_homer <- as.numeric(raw$logp_homer)
    raw$qvalue <- as.numeric(raw$qvalue)

    # Derived columns
    # logp_homer is ln(p), a negative number, stored as-is in ln_pvalue
    raw$ln_pvalue <- raw$logp_homer
    raw$neg_log10p <- -raw$ln_pvalue / log(10) # positive: -log10(p)
    raw$motif_length <- nchar(as.character(raw$consensus))
    raw$fold_enrichment <- ifelse(raw$pct_background > 0,
        raw$pct_target / raw$pct_background,
        NA_real_
    )
    raw$n_target_total <- n_target_total
    raw$n_background_total <- n_bg_total

    # Correct pvalue_str: HOMER writes "0" when the numeric value underflows.
    # Reconstruct a readable string from neg_log10p in that case.
    raw$pvalue_str <- ifelse(
        is.na(raw$pvalue_str) | raw$pvalue_str == "0",
        sprintf("1e-%d", round(raw$neg_log10p)),
        raw$pvalue_str
    )

    # Reorder columns (logp_homer is internal only; not exported)
    raw[, c(
        "rank", "motif_id", "consensus", "motif_length",
        "pvalue_str", "ln_pvalue", "neg_log10p", "qvalue",
        "n_target", "pct_target", "n_background", "pct_background",
        "fold_enrichment", "n_target_total", "n_background_total"
    )]
}


#' Read a HOMER output directory
#'
#' Reads the known and de novo motif enrichment results produced by HOMER
#' `findMotifsGenome.pl`, along with the command used to generate them.
#' Missing parts (e.g. a run with `-nomotif` or `-noknown`) are skipped.
#'
#' @param dir Path to the HOMER output directory (containing `knownResults.txt`
#'   and/or `homerResults/motif*.motif` files).
#'
#' @return A `homer_output` object (list) with elements:
#'   \describe{
#'     \item{known}{Data frame of known motif enrichment results, one row per
#'       motif in `knownResults.txt` order, or NULL. Columns: rank, motif_id,
#'       consensus, motif_length, pvalue_str, ln_pvalue, neg_log10p, qvalue,
#'       n_target, pct_target, n_background, pct_background,
#'       fold_enrichment, n_target_total, n_background_total.}
#'     \item{denovo}{Data frame of de novo motif enrichment results, sorted by
#'       significance, or NULL. Columns: rank, motif_num, name, consensus,
#'       motif_length, best_match, best_match_name, best_match_score,
#'       pvalue_str, ln_pvalue, neg_log10p, n_target, pct_target,
#'       n_background, pct_background, fold_enrichment, n_target_total,
#'       n_background_total, pos_std_target, pos_std_bg.}
#'     \item{params}{Character string with the HOMER command, or NULL}
#'     \item{input}{Basename of the input BED file parsed from the HOMER
#'       command, or NULL if unavailable.}
#'     \item{known_db}{Name of the known motif database, or NULL if HOMER
#'       was run with `-noknown`.}
#'     \item{denovo_db}{Name of the database used for de novo best matches,
#'       or NULL if HOMER was run with `-nomotif`.}
#'     \item{dir}{The input directory path}
#'   }
#'
#' @details
#' In both data frames, `pvalue_str` is the p-value as written by HOMER
#' (reconstructed from `ln_pvalue` when HOMER writes `"0"` because of
#' floating-point underflow), `ln_pvalue` is the natural log of the p-value
#' (a negative number), `neg_log10p` is -log10(p), and `fold_enrichment` is
#' `pct_target / pct_background` (NA when `pct_background` is 0). `rank` is
#' the position by significance; for de novo motifs, `motif_num` is the
#' original HOMER motif number N (`homerResults/motif{N}.motif`).
#'
#' @export
#' @examples
#' example_dir <- system.file("extdata",
#'     "A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001",
#'     package = "homerViz"
#' )
#' homer <- read_homer_output(example_dir)
#' homer
#' head(homer$known)
read_homer_output <- function(dir) {
    if (!dir.exists(dir)) stop("Directory not found: ", dir)
    message("Reading: ", basename(dir))

    known_file <- .find_known_file(dir)
    if (file.exists(known_file)) {
        message("  * Known motifs:  ", basename(known_file))
        known <- tryCatch(read_homer_known(dir), error = function(e) {
            warning(e$message)
            NULL
        })
    } else {
        message("  * Known motifs:  (not found, skipping)")
        known <- NULL
    }

    denovo_motif_files <- Sys.glob(
        file.path(dir, "homerResults", "motif*.motif")
    )
    n_denovo <- sum(grepl("motif[0-9]+\\.motif$", basename(denovo_motif_files)))
    if (length(denovo_motif_files) > 0) {
        message("  * De novo motifs: ", n_denovo, " motifs in homerResults/")
        denovo <- tryCatch(read_homer_denovo(dir), error = function(e) {
            warning(e$message)
            NULL
        })
    } else {
        message(
            "  * De novo motifs: ",
            "(homerResults/motif*.motif not found, skipping)"
        )
        denovo <- NULL
    }

    params <- .read_homer_params(dir)
    input_bed <- .parse_input_bed(params)
    known_db <- .parse_known_db(params)
    denovo_db <- .parse_denovo_db(params)

    structure(
        list(
            known     = known,
            denovo    = denovo,
            params    = params,
            input     = input_bed,
            known_db  = known_db,
            denovo_db = denovo_db,
            dir       = dir
        ),
        class = "homer_output"
    )
}

#' @export
print.homer_output <- function(x, ...) {
    cat("HOMER output\n")
    cat("  Directory:", x$dir, "\n")
    if (!is.null(x$input)) {
        cat("  Input BED:", x$input, "\n")
    }
    if (!is.null(x$known)) {
        cat("  Known motifs:", nrow(x$known), "motifs\n")
        cat("  Known motif database:", x$known_db %||% "none (-noknown)", "\n")
    } else {
        cat("  Known motifs: none\n")
    }
    if (!is.null(x$denovo)) {
        cat("  De novo motifs:", nrow(x$denovo), "motifs\n")
        cat(
            "  De novo best-match database:",
            x$denovo_db %||% "none (-nomotif)", "\n"
        )
    } else {
        cat("  De novo motifs: none\n")
    }
    if (!is.null(x$params)) {
        cat("  Command:", x$params, "\n")
    }
    invisible(x)
}


#' Read HOMER de novo motif enrichment results
#'
#' Parses `homerResults/motif*.motif` files and (if present) `homerResults.html`
#' to build a tidy data frame of de novo motif enrichment results, sorted by
#' significance. Analogous to `read_homer_known()` for de novo results.
#'
#' @param dir Path to the HOMER output directory.
#'
#' @return A data frame with columns: rank, motif_num, name, consensus,
#'   motif_length, best_match, best_match_name, best_match_score,
#'   pvalue_str, ln_pvalue, neg_log10p, n_target, pct_target,
#'   n_background, pct_background, fold_enrichment, n_target_total,
#'   n_background_total, pos_std_target, pos_std_bg. Returns `NULL` if no
#'   motif files are found. `rank` is the position after sorting by
#'   significance; `motif_num` is the original HOMER file number N
#'   (`homerResults/motif{N}.motif`), used to locate
#'   `homerResults/motif{N}.logo.svg`.
#'
#' @noRd
read_homer_denovo <- function(dir) {
    motif_files <- Sys.glob(file.path(dir, "homerResults", "motif*.motif"))
    nums_str <- sub(
        ".*motif([0-9]+)\\.motif$", "\\1", basename(motif_files)
    )
    motif_nums <- rep(NA_integer_, length(nums_str))
    valid_idx <- grepl("^[0-9]+$", nums_str)
    motif_nums[valid_idx] <- as.integer(nums_str[valid_idx])
    valid <- !is.na(motif_nums)
    motif_files <- motif_files[valid]
    motif_nums <- motif_nums[valid]
    if (length(motif_files) == 0) {
        return(NULL)
    }

    # Sort by motif number before metadata parsing
    sort_idx <- order(motif_nums)
    motif_files <- motif_files[sort_idx]
    motif_nums <- motif_nums[sort_idx]

    # Parse metadata from each motif file header
    motif_metas <- lapply(motif_files, function(f) {
        h <- readLines(f, n = 1, warn = FALSE)
        tryCatch(
            .parse_motif_header(h),
            error = function(e) {
                list(
                    consensus = NA_character_, name = NA_character_,
                    ln_pvalue = NA_real_, pvalue_str = NA_character_,
                    target_count = NA_real_, target_pct = NA_real_,
                    bg_count = NA_real_, bg_pct = NA_real_
                )
            }
        )
    })

    # Parse best_match, pos_std, totals from homerResults.html (if present)
    homer_html_files <- Sys.glob(file.path(dir, "*homerResults*.html"))
    homer_html_files <- sort(
        homer_html_files[!grepl("_homerViz\\.html$", homer_html_files)]
    )
    if (length(homer_html_files) > 1) {
        warning(
            "Multiple homerResults HTML files found; using: ",
            basename(homer_html_files[1])
        )
    }
    homer_html <- if (length(homer_html_files) > 0) {
        homer_html_files[1]
    } else {
        file.path(dir, "homerResults.html")
    }
    homer_data <- .parse_homer_html(homer_html)
    totals <- homer_data[["__totals__"]]
    n_target <- totals$n_target %||% NA_integer_
    n_bg <- totals$n_bg %||% NA_integer_

    # Sort by significance (ln_pvalue is negative; ascending = most significant
    # first). Use NA for missing values so they sort last rather than
    # incorrectly to p=1.
    lp_all <- vapply(motif_metas, function(m) m$ln_pvalue %||% NA_real_, numeric(1))
    ord <- order(lp_all, na.last = TRUE)
    motif_metas <- motif_metas[ord]
    motif_nums <- motif_nums[ord]

    # Build one row per motif
    rows <- lapply(seq_along(motif_metas), function(i) {
        meta <- motif_metas[[i]]
        orig_num <- motif_nums[i]
        hd <- homer_data[[as.character(orig_num)]] %||% list()

        # Best match: parse name and score
        best <- hd$best_match %||% ""
        score_m <- regmatches(best, regexpr("\\(([0-9.]+)\\)$", best))
        score_str <- if (length(score_m) > 0) {
            sub("^\\((.+)\\)$", "\\1", score_m)
        } else {
            NA_character_
        }
        best_name <- if (nzchar(best)) {
            trimws(sub("\\([0-9.]+\\)$", "", best))
        } else {
            NA_character_
        }

        # Position STD: split "48.9bp (64.1bp)" into target and bg parts.
        # bg part is only present when parentheses exist; fall back to NA
        # otherwise.
        pos_std <- hd$pos_std %||% NA_character_
        tgt_std <- if (!is.na(pos_std)) {
            sub("^([0-9.]+bp).*", "\\1", pos_std)
        } else {
            NA_character_
        }
        bg_std <- if (!is.na(pos_std) && grepl("\\([0-9.]+bp\\)", pos_std)) {
            sub(".*\\(([0-9.]+bp)\\).*", "\\1", pos_std)
        } else {
            NA_character_
        }

        pct_tgt <- meta$target_pct %||% NA_real_
        pct_bg <- meta$bg_pct %||% NA_real_
        fe <- if (!is.na(pct_tgt) && !is.na(pct_bg) && pct_bg > 0) {
            pct_tgt / pct_bg
        } else {
            NA_real_
        }

        consensus <- if (
            !is.null(meta$consensus) && !is.na(meta$consensus)
        ) {
            as.character(meta$consensus)
        } else {
            NA_character_
        }

        data.frame(
            rank = i,
            motif_num = orig_num,
            name = as.character(meta$name %||% NA_character_),
            consensus = consensus,
            motif_length = if (!is.na(consensus)) {
                nchar(consensus)
            } else {
                NA_integer_
            },
            best_match = if (!is.na(best_name)) {
                paste0(
                    best_name,
                    if (!is.na(score_str)) paste0("(", score_str, ")") else ""
                )
            } else {
                NA_character_
            },
            best_match_name = best_name,
            best_match_score = if (!is.na(score_str)) as.numeric(score_str) else NA_real_,
            pvalue_str = meta$pvalue_str %||% NA_character_,
            ln_pvalue = meta$ln_pvalue %||% NA_real_,
            neg_log10p = -(meta$ln_pvalue %||% NA_real_) / log(10),
            # Column names match the known df for consistency
            n_target = meta$target_count %||% NA_real_,
            pct_target = pct_tgt,
            n_background = meta$bg_count %||% NA_real_,
            pct_background = pct_bg,
            fold_enrichment = fe,
            n_target_total = n_target,
            n_background_total = n_bg,
            pos_std_target = tgt_std,
            pos_std_bg = bg_std
        )
    })

    do.call(rbind, rows)
}

# ---- Internal file-finding helpers ------------------------------------------

# Read motifFindingParameters.txt and return the raw argument string, or NULL.
# The file contains a single line: "cmd = <bed> <genome> <outdir> [flags...]"
.read_homer_params <- function(dir) {
    path <- file.path(dir, "motifFindingParameters.txt")
    if (!file.exists(path)) {
        return(NULL)
    }
    line <- readLines(path, n = 1, warn = FALSE)
    trimws(sub("^cmd\\s*=\\s*", "", line))
}

# Extract the basename of the input BED file from the HOMER argument string.
# The BED file is always the first positional argument.
.parse_input_bed <- function(params) {
    if (is.null(params) || !nzchar(params)) {
        return(NULL)
    }
    parts <- strsplit(params, "\\s+")[[1]]
    if (length(parts) < 1) {
        return(NULL)
    }
    basename(parts[1])
}

# Extract the known motif database name from the HOMER command string.
# Returns NULL if -noknown is present (no known analysis was run).
# Returns the filename without extension if -mknown <path> is present.
# Returns "Homer motif database (default)" otherwise.
.parse_known_db <- function(params) {
    if (is.null(params) || !nzchar(params)) {
        return("Homer motif database (default)")
    }
    if (grepl("-noknown", params)) {
        return(NULL)
    }
    m <- regmatches(params, regexpr("-mknown\\s+(\\S+)", params))
    if (length(m) == 0 || !nzchar(m)) {
        return("Homer motif database (default)")
    }
    path <- sub("-mknown\\s+", "", m)
    tools::file_path_sans_ext(basename(path))
}

# Extract the de novo best-match database name from the HOMER command string.
# Returns NULL if -nomotif is present (no de novo analysis was run).
# Returns the filename without extension if -mcheck <path> is present.
# Returns "Homer motif database (default)" otherwise.
.parse_denovo_db <- function(params) {
    if (is.null(params) || !nzchar(params)) {
        return("Homer motif database (default)")
    }
    if (grepl("-nomotif", params)) {
        return(NULL)
    }
    m <- regmatches(params, regexpr("-mcheck\\s+(\\S+)", params))
    if (length(m) == 0 || !nzchar(m)) {
        return("Homer motif database (default)")
    }
    path <- sub("-mcheck\\s+", "", m)
    tools::file_path_sans_ext(basename(path))
}

# Find knownResults.txt in a directory, supporting custom filename prefixes.
# HOMER always ends the file in "knownResults.txt"; it may have a prefix like
# "20251208_sample_knownResults.txt".
.find_known_file <- function(dir) {
    # Prefer the canonical name
    standard <- file.path(dir, "knownResults.txt")
    if (file.exists(standard)) {
        return(standard)
    }

    # Fall back to any file matching *knownResults.txt
    hits <- list.files(dir, pattern = "knownResults\\.txt$", full.names = TRUE)
    if (length(hits) == 1L) {
        return(hits[1])
    }
    if (length(hits) > 1L) {
        warning(
            "Multiple knownResults.txt files found; using: ",
            basename(hits[1])
        )
        return(hits[1])
    }

    # Not found - return the standard path so the caller gets a clear error
    standard
}
