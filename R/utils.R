#' @importFrom utils head
#' @importFrom rlang .data
NULL

# Null-coalescing operator (package-wide, single definition)
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ---- Internal: parse HOMER homerResults.html for per-motif metadata ----------
# Returns a named list keyed by motif number (as character).
# Each element: list(pvalue_str, ln_pvalue, best_match, pos_std)
# Special key "__totals__": list(n_target, n_bg)
.parse_homer_html <- function(homer_html) {
    if (!file.exists(homer_html)) {
        return(list())
    }

    lines <- readLines(homer_html, warn = FALSE)
    out <- list()
    n_target <- NA_integer_
    n_bg <- NA_integer_

    for (line in lines) {
        if (grepl("Total target sequences", line, ignore.case = TRUE)) {
            m <- regmatches(line, regexpr(
                "Total target sequences\\s*=\\s*([0-9]+)",
                line, ignore.case = TRUE
            ))
            if (length(m) > 0) n_target <- as.integer(sub(".*=\\s*", "", m))
        }
        if (grepl("Total background sequences", line, ignore.case = TRUE)) {
            m <- regmatches(line, regexpr(
                "Total background sequences\\s*=\\s*([0-9]+)",
                line, ignore.case = TRUE
            ))
            if (length(m) > 0) n_bg <- as.integer(sub(".*=\\s*", "", m))
        }

        ms <- regmatches(line, gregexpr(
            "homerResults/motif([0-9]+)\\.info\\.html", line
        ))[[1]]
        if (length(ms) == 0) next
        ns <- sub("homerResults/motif([0-9]+)\\.info\\.html", "\\1", ms)

        # p-value string and ln(p-value): <TD>{pval}</TD><TD>{ln_p}</TD>
        pval_m <- regmatches(line, regexpr(
            "<TD>([0-9e+\\-\\.]+)</TD><TD>(-?[0-9e+\\.]+)</TD>", line
        ))
        pval_str <- if (length(pval_m) > 0) {
            sub("<TD>([0-9e+\\-\\.]+)</TD>.*", "\\1", pval_m)
        } else {
            NA_character_
        }
        ln_pval <- if (length(pval_m) > 0) {
            as.numeric(sub(
                ".*<TD>(-?[0-9e+\\.]+)</TD>$", "\\1", pval_m
            ))
        } else {
            NA_real_
        }

        # Best-match text before <BR/>
        bm <- if (grepl("<TD>[^<]+<BR/>", line)) {
            trimws(sub(".*<TD>([^<]+)<BR/>.*", "\\1", line))
        } else {
            ""
        }

        # STD(Bg STD): e.g. "48.9bp (64.1bp)"
        std_m <- regmatches(line, regexpr(
            "<TD>([0-9.]+bp \\([0-9.]+bp\\))</TD>", line
        ))
        pos_std <- if (length(std_m) > 0) {
            sub("<TD>(.*)</TD>", "\\1", std_m)
        } else {
            NA_character_
        }

        for (n in ns) {
            if (!n %in% names(out)) {
                out[[n]] <- list(
                    pvalue_str = pval_str, ln_pvalue = ln_pval,
                    best_match = bm, pos_std = pos_std
                )
            }
        }
    }
    out[["__totals__"]] <- list(n_target = n_target, n_bg = n_bg)
    out
}


# ---- Internal: parse a HOMER motif header line -------------------------------
# @param line Character string of the header line (starting with ">")
# @return Named list with motif metadata
.parse_motif_header <- function(line) {
    # Remove leading ">"
    line <- sub("^>", "", line)
    parts <- strsplit(line, "\t")[[1]]

    if (length(parts) < 6) {
        stop(
            "Malformed motif header: expected at least 6 tab-separated",
            " fields, got ", length(parts), ". Line: ", line
        )
    }

    consensus <- parts[1]
    name <- parts[2]
    ln_pvalue <- as.numeric(parts[4]) # HOMER stores ln(p-value) here
    # parts[3] is the score threshold; parts[5] is always "0"
    stats_str <- parts[6]

    # Parse "T:count(pct%),B:count(pct%),P:pvalue"
    t_match <- regmatches(stats_str, regexpr(
        "T:([0-9.]+)\\(([0-9.]+)%\\)", stats_str
    ))
    b_match <- regmatches(stats_str, regexpr(
        "B:([0-9.]+)\\(([0-9.]+)%\\)", stats_str
    ))
    p_match <- regmatches(stats_str, regexpr("P:([0-9e+\\-]+)", stats_str))

    t_vals <- as.numeric(regmatches(t_match, gregexpr("[0-9.]+", t_match))[[1]])
    b_vals <- as.numeric(regmatches(b_match, gregexpr("[0-9.]+", b_match))[[1]])
    pval_str_raw <- if (length(p_match) > 0) {
        sub("P:", "", p_match)
    } else {
        NA_character_
    }

    # Store ln(p) as given by HOMER in parts[4] - avoids underflow for
    # extreme p-values like e^-2584 where as.numeric() returns 0.
    # HOMER stores ln(p) as a negative number, e.g. -2584

    list(
        consensus    = consensus,
        name         = name,
        target_count = t_vals[1],
        target_pct   = t_vals[2],
        bg_count     = b_vals[1],
        bg_pct       = b_vals[2],
        pvalue_str   = pval_str_raw,
        ln_pvalue    = if (!is.na(ln_pvalue)) ln_pvalue else NA_real_
    )
}


# ---- Internal: extract N (total sequences) from a HOMER column name ---------
# @param colname Column name like "# of Target Sequences with Motif(of 4101)"
# @return Integer N, or NA if not found
.extract_n_sequences <- function(colname) {
    m <- regmatches(colname, regexpr("\\(of ([0-9]+)\\)", colname))
    if (length(m) == 0) {
        return(NA_integer_)
    }
    as.integer(sub("\\(of ([0-9]+)\\)", "\\1", m))
}
