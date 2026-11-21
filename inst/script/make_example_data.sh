#!/usr/bin/env bash
# =============================================================================
# How the example data in inst/extdata was generated
# =============================================================================
#
# Files described here:
#   inst/extdata/A549_MED1_Ctrl.stdchr.chr1.bed
#   inst/extdata/H14CORE_homer_format_0.0001.motif
#       (kept in the source repository; excluded from the package tarball)
#   inst/extdata/A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001/
#
# ---- 1. Peaks: A549 MED1 ChIP-seq, control condition ------------------------
#
# Source: MED1 ChIP-seq in A549 lung carcinoma cells, control (EtOH vehicle)
# condition, two biological replicates, from:
#   Tav C, et al. (2023). Glucocorticoid stimulation induces regionalized
#   gene responses within topologically associating domains.
#   Frontiers in Genetics 14:1237092. doi:10.3389/fgene.2023.1237092
#
# Raw data: GEO series GSE226487.
#
# Processing of the raw reads (as described in the GEO records): reads were
# analysed with the ChIP-seq pipeline developed at the McGill University and
# Genome Quebec Innovation Centre. High-quality reads were aligned to the
# human reference genome (hg38) with BWA, PCR duplicates were removed with
# Picard MarkDuplicates, and narrow peaks were called with MACS2 callpeak,
# using the corresponding sequenced input DNA as background control.
#
# The MACS2 peaks of both replicates are provided in the article's
# Supplementary Table 1 (tabs "MED1_CTRL_rep1" and "MED1_CTRL_rep2";
# 1-based coordinates). Download it, and the ENCODE hg38 blacklist
# (ENCFF356LFX; Amemiya et al. 2019, Sci Rep 9:9354):

curl -L -o Tav2023_Table_1.xlsx \
    "https://public-pages-files-2025.frontiersin.org/articles/1237092/file/Table_1.XLSX/1237092_supplementary-materials_tables_1_xlsx/1"
curl -L -o ENCFF356LFX.bed.gz \
    "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"

# Each replicate was restricted to the standard chromosomes ("stdchr") and
# peaks overlapping the ENCODE blacklist were removed. Peaks shared by both
# replicates (gVenn overlap group "11") were kept, then subset to chr1 to keep
# the example small (460 peaks, genome build hg38). Requires the R packages
# readxl, GenomicRanges, GenomeInfoDb, rtracklayer and gVenn.

Rscript - <<'EOF'
suppressPackageStartupMessages({
    library(readxl)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(rtracklayer)
    library(gVenn)
})

read_peaks <- function(sheet) {
    x <- read_excel("Tav2023_Table_1.xlsx", sheet = sheet)
    GRanges(x$seqnames, IRanges(x$start, x$end))
}

blacklist <- import("ENCFF356LFX.bed.gz")

clean_peaks <- function(gr) {
    gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
    gr[!overlapsAny(gr, blacklist)]
}

peaks <- list(
    MED1_CTRL_rep1 = clean_peaks(read_peaks("MED1_CTRL_rep1")),
    MED1_CTRL_rep2 = clean_peaks(read_peaks("MED1_CTRL_rep2"))
)

# Peaks present in both replicates
groups <- extractOverlaps(computeOverlaps(peaks))
common <- groups[["group_11"]]

# Keep chr1 only and write a BED file (0-based start)
chr1 <- sort(common[seqnames(common) == "chr1"])
write.table(
    data.frame(seqnames(chr1), start(chr1) - 1, end(chr1), ".", 0, "."),
    "A549_MED1_Ctrl.stdchr.chr1.bed",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
EOF

# ---- 2. Motif database: HOCOMOCO v14 core collection ------------------------
#
# H14CORE_homer_format_0.0001.motif is the HOCOMOCO v14 H14CORE collection
# (1,595 motifs) in HOMER format, with motif score thresholds corresponding
# to p < 0.0001. It was downloaded from the HOCOMOCO v14 website:
#   https://hocomoco14.autosome.org/downloads_v14
#
# ---- 3. HOMER motif enrichment ----------------------------------------------
#
# HOMER version: v4.11.1
#
# Run from the directory containing the BED and motif files. The hg38 genome
# must be installed in HOMER (configureHomer.pl -install hg38).

findMotifsGenome.pl A549_MED1_Ctrl.stdchr.chr1.bed \
    hg38 \
    A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001 \
    -p 8 \
    -size 200 \
    -mask \
    -mknown H14CORE_homer_format_0.0001.motif \
    -mcheck H14CORE_homer_format_0.0001.motif

# The exact command is also stored by HOMER in
# A549_MED1_Ctrl_chr1_vs_HOCOMOCOv14_0p0001/motifFindingParameters.txt
#
# ---- 4. Files kept in the package -------------------------------------------
#
# To stay under the Bioconductor tarball size limit, only the HOMER output
# files read by homerViz are included in the built package (see
# .Rbuildignore):
#   knownResults.txt               known motif enrichment table
#   knownResults/known{N}.logo.svg known motif logos
#   homerResults.html              de novo summary (best matches, STD)
#   homerResults/motif{N}.motif    de novo motif matrices and statistics
#   homerResults/motif{N}.logo.svg de novo motif logos
#   motifFindingParameters.txt     HOMER command line
#
# Excluded from the build: knownResults.html, knownResults/known{N}.motif,
# homerResults/motif{N}.info.html, homerResults/motif{N}.similar*,
# homerResults/motif{N}RV.*, homerMotifs.*, seq.autonorm.tsv, and the
# *_homerViz.html reports generated by save_homer_html().
