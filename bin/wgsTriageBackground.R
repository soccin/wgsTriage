#!/usr/bin/env Rscript
##
## Build the background QC distribution from previously mapped projects.
##
## Usage:
##   ./bin/wgsTriageBackground.R [QCDATA_ROOT] [--out DIR]
##
## Default root is ./QCData, default output directory is ./background.
##
## Writes three files:
##   backgroundSamples.tsv   every historical sample, with its gate verdict
##   backgroundStats.tsv     robust reference ranges, from clean samples only
##   backgroundFlagged.tsv   historical samples that fail the gates
##
## The reference ranges are computed only from samples that pass the gates.
## The archive is known to contain defective cohorts, and including them would
## widen the reference range enough to admit the next bad cohort. Robust
## statistics alone are not sufficient protection when the contaminated
## fraction is large.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Repo root, derived from this script's own location in bin/ and expressed
## relative to the working directory. Deliberately not here::here(): that
## resolves from the working directory, which for this tool is the data
## directory rather than the repo, and it fails silently when it guesses wrong.
##
scriptPath <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |>
    str_remove("^--file=")
repoRoot <- if (length(scriptPath) > 0) {
    path_rel(path_dir(path_dir(path_real(scriptPath))))
} else {
    "."
}

source(path(repoRoot, "R", "qcLib.R"))

args <- commandArgs(trailingOnly = TRUE)
outFlag <- which(args == "--out")
outDir <- if (length(outFlag) > 0) args[outFlag + 1] else path(repoRoot, "data", "background")
positional <- args[!args %in% c("--out", outDir)]
qcRoot <- if (length(positional) > 0) positional[1] else "QCData"

if (!dir_exists(qcRoot)) {
    stop(glue("QCData root not found: {qcRoot}"))
}
dir_create(outDir)

cat(glue("Scanning {qcRoot} for Picard metrics ...\n\n"))

picard <- collectPicardSamples(qcRoot)
if (is.null(picard) || nrow(picard) == 0) {
    stop(glue("No Picard metrics files found under {qcRoot}"))
}

## Project label for a multiqc data directory, matched to the label derived
## from Picard paths so the two sources join.
projectFromMultiqcPath <- function(path) {
    str_remove(path, "/sbam/multiqc/.*$") |>
        str_remove("/Map$") |>
        path_file()
}

mqcFiles <- dir_ls(qcRoot, recurse = TRUE, glob = "*multiqc_samtools_stats.txt", fail = FALSE)

samtools <- if (length(mqcFiles) > 0) {
    tibble(path = as.character(mqcFiles)) |>
        mutate(project = projectFromMultiqcPath(path),
               parsed = map(path, readMultiqcSamtools)) |>
        filter(!map_lgl(parsed, is.null)) |>
        select(project, parsed) |>
        unnest(parsed) |>
        distinct(project, sample, .keep_all = TRUE)
} else {
    NULL
}

background <- if (is.null(samtools)) {
    picard
} else {
    left_join(picard, samtools, by = c("project", "sample"))
}

nSamtools <- if ("supplementaryRate" %in% names(background)) sum(!is.na(background$supplementaryRate)) else 0

cat(sprintf("Parsed %d samples from %d projects\n", nrow(background), n_distinct(background$project)))
cat(sprintf("Samtools metrics available for %d samples\n\n", nSamtools))

##
## Gate every historical sample. This is the same evaluation the report applies
## to a current cohort, so the background carries verdicts on the same terms.
##
## Sample names are not unique across the archive: 100 of them occur in more
## than one project, usually because a cohort was remapped under a new name.
## The gate functions key on a single column, so pass them a project-qualified
## id. Keying on the bare sample name merges verdicts across unrelated runs and
## silently assigns one project's failure to another project's sample.
##
background <- background |> mutate(uid = str_c(project, "::", sample))

gateResults <- evaluateGates(background |> select(-sample) |> rename(sample = uid))
verdicts <- sampleVerdict(gateResults) |> rename(uid = sample)

background <- background |>
    left_join(verdicts, by = "uid") |>
    mutate(hasCore = !is.na(pctChimeras) & !is.na(meanCoverage),
           referenceSample = nFail == 0 & hasCore)

##
## Robust reference ranges from clean samples only.
##
metricCols <- c(GATES$metric, "meanCoverage", "insertSizeAverage",
                "pctImproperPairs", "interChromRate", "pctExcDupe")
metricCols <- metricCols[metricCols %in% names(background)]

referenceStats <- background |>
    filter(referenceSample) |>
    select(sampleType, all_of(metricCols)) |>
    pivot_longer(-sampleType, names_to = "metric", values_to = "value") |>
    filter(!is.na(value)) |>
    summarise(
        n      = n(),
        median = median(value),
        mad    = mad(value),
        q01    = quantile(value, 0.01),
        q05    = quantile(value, 0.05),
        q25    = quantile(value, 0.25),
        q75    = quantile(value, 0.75),
        q95    = quantile(value, 0.95),
        q99    = quantile(value, 0.99),
        min    = min(value),
        max    = max(value),
        .by = metric)

## Coverage is the one metric where tumor and normal genuinely differ, so it
## also gets a per-class reference. Everything else measures data integrity and
## should not vary with sample class.
coverageStats <- background |>
    filter(referenceSample, !is.na(meanCoverage)) |>
    summarise(
        n      = n(),
        median = median(meanCoverage),
        mad    = mad(meanCoverage),
        q05    = quantile(meanCoverage, 0.05),
        q25    = quantile(meanCoverage, 0.25),
        q75    = quantile(meanCoverage, 0.75),
        q95    = quantile(meanCoverage, 0.95),
        .by = sampleType) |>
    mutate(metric = "meanCoverage") |>
    relocate(metric, sampleType)

flagged <- background |>
    filter(verdict %in% c("FAIL", "WARN")) |>
    select(project, sample, sampleType, verdict, verdictReason, nFail, nWarn,
           failedMetrics, warnedMetrics,
           any_of(c("pctChimeras", "pctSoftclip", "alignedFrac", "pctReadUsed",
                    "supplementaryRate", "meanCoverage", "pctExcTotal"))) |>
    arrange(desc(nFail), desc(pctChimeras))

## Per-metric coverage of the background. A reference range built from three
## samples is not a reference range, and the report must be able to say so
## rather than quoting it as though it were solid.
metricCoverage <- background |>
    select(all_of(metricCols)) |>
    summarise(across(everything(), \(x) sum(!is.na(x)))) |>
    pivot_longer(everything(), names_to = "metric", values_to = "nAvailable") |>
    mutate(nTotal = nrow(background),
           pctAvailable = round(nAvailable / nTotal * 100, 1))

write_tsv(background, path(outDir, "backgroundSamples.tsv"))
write_tsv(referenceStats, path(outDir, "backgroundStats.tsv"))
write_tsv(coverageStats, path(outDir, "backgroundCoverageStats.tsv"))
write_tsv(flagged, path(outDir, "backgroundFlagged.tsv"))
write_tsv(metricCoverage, path(outDir, "backgroundMetricCoverage.tsv"))

##
## Console summary.
##
nRef <- sum(background$referenceSample)
nFailSamples <- sum(background$verdict == "FAIL")
nWarnSamples <- sum(background$verdict == "WARN")
rule <- strrep("=", 74)

cat(rule, "\n")
cat("  BACKGROUND QC IMPORT\n")
cat(rule, "\n\n")
cat(sprintf("  Samples parsed      %5d\n", nrow(background)))
cat(sprintf("  Projects            %5d\n", n_distinct(background$project)))
cat(sprintf("  Reference samples   %5d   clean, used to set the ranges\n", nRef))
cat(sprintf("  Failing gates       %5d\n", nFailSamples))
cat(sprintf("  Warning only        %5d\n\n", nWarnSamples))

if (nFailSamples > 0) {
    cat(sprintf("  %d historical samples fail the gates and are excluded from the\n", nFailSamples))
    cat("  reference ranges. Review backgroundFlagged.tsv: these are either real\n")
    cat("  defects that were processed anyway, or evidence a threshold is wrong.\n\n")

    cat("  Worst historical samples\n")
    cat("  ", strrep("-", 70), "\n", sep = "")
    flagged |>
        filter(verdict == "FAIL") |>
        head(12) |>
        mutate(line = sprintf("  %-24s %-14s chimeras %6.2f%%  softclip %6.2f%%",
                              str_trunc(sample, 24), str_trunc(project, 14),
                              pctChimeras, pctSoftclip)) |>
        pull(line) |>
        walk(\(x) cat(x, "\n"))
    if (nFailSamples > 12) cat(sprintf("  ... and %d more\n", nFailSamples - 12))
    cat("\n")
}

cat("  Reference ranges, clean samples only\n")
cat("  ", strrep("-", 70), "\n", sep = "")
referenceStats |>
    left_join(metricCoverage |> select(metric, pctAvailable), by = "metric") |>
    mutate(line = sprintf("  %-19s n=%-4d median %9.3f   q05-q95 %8.3f - %8.3f  [%s%% of archive]",
                          metric, n, median, q05, q95, format(pctAvailable, nsmall = 1))) |>
    pull(line) |>
    walk(\(x) cat(x, "\n"))

missingMetrics <- metricCoverage |> filter(nAvailable == 0)
if (nrow(missingMetrics) > 0) {
    cat("\n  No background available for: ", str_c(missingMetrics$metric, collapse = ", "), "\n", sep = "")
    cat("  These gates run on fixed thresholds only, with no historical context.\n")
}

cat("\n")
cat(sprintf("  Written to %s/\n", outDir))
cat(rule, "\n")
