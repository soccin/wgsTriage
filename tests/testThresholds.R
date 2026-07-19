#!/usr/bin/env Rscript
##
## Smoke test for the filter thresholds.
##
## A filter that blocks multi-day compute must not be able to silently stop
## filtering. This runs bin/wgsTriage.R end to end against a two-sample fixture
## and asserts the verdicts, the completeness accounting and the exit status.
##
## This is deliberately one test and not a suite. It covers the regression that
## would otherwise be invisible until it costs another three days of compute:
## thresholds that quietly pass everything.
##
## Fixtures in fixtures/miniCohort/ are two real samples from Proj_16840_N,
## trimmed to the header, the metrics block and the PAIR row, and renamed:
##   DEFECT_N01  worst sample in the cohort, 17.3% chimeras, fails 5 thresholds
##   CLEAN_N01   passes all 7
## Both classify as N, so no tumor/normal pair forms. That is expected here.
##
## Usage:
##   Rscript tests/testThresholds.R
##
## Exits 0 if every assertion holds, 1 otherwise.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Repo root from this script's own location, matching the bootstrap in bin/.
## tests/ sits one level below the root exactly as bin/ does. Deliberately not
## here::here(): it resolves from the working directory and fails silently,
## picking up a stale copy of qcLib.R from outside the repo.
##
scriptPath <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |>
    str_remove("^--file=")
repoRoot <- if (length(scriptPath) > 0) {
    path_rel(path_dir(path_dir(path_real(scriptPath))))
} else {
    "."
}

fixtureDir <- path(repoRoot, "tests", "fixtures", "miniCohort")
backgroundDir <- path(repoRoot, "data", "background")
outDir <- path(tempdir(), "wgsTriageTest")

failures <- character()

check <- function(label, ok) {
    ok <- isTRUE(ok)
    if (!ok) failures <<- c(failures, label)
    cat(glue("  {if (ok) 'ok  ' else 'FAIL'}  {label}"), "\n", sep = "")
}

##
## Run the tool. Output is captured to a file rather than piped: piping Rscript
## into another process raises SIGPIPE and kills R mid-write, which looks
## exactly like a bug in the script under test.
##
runLog <- path(outDir, "run.log")
dir_create(outDir)

status <- system2("Rscript",
                  c(path(repoRoot, "bin", "wgsTriage.R"),
                    fixtureDir,
                    "--background", backgroundDir,
                    "--out", outDir,
                    "--project", "miniCohort"),
                  stdout = runLog, stderr = runLog)

cat("\n")
check("exit status is 0", status == 0)

samplesFile <- path(outDir, "preflightQC_samples.tsv")
if (!file_exists(samplesFile)) {
    cat("\n  No sample report written. Tool output follows:\n\n", sep = "")
    walk(read_lines(runLog), \(x) cat("  ", x, "\n", sep = ""))
    quit(save = "no", status = 1)
}

samples <- read_tsv(samplesFile, show_col_types = FALSE, progress = FALSE)
report <- read_lines(path(outDir, "preflightQC.txt"), progress = FALSE)

verdictOf <- function(name) {
    samples |> filter(sample == name) |> pull(verdict)
}

check("DEFECT_N01 returns FAIL", identical(verdictOf("DEFECT_N01"), "FAIL"))
check("CLEAN_N01 returns PASS", identical(verdictOf("CLEAN_N01"), "PASS"))

## The denominator is every sample directory under out/metrics, not every
## sample that parsed. Reporting on a subset without saying so is the failure
## this accounting exists to prevent, so assert the phrasing, not just the count.
check("cohort completeness reports 2 of 2",
      any(str_detect(report, fixed("2 of 2 samples checked"))))

## The background supplies the reference ranges the report compares against.
## Without it the tool degrades to fixed thresholds and says so, which is
## correct behaviour but leaves the out-of-range path untested.
check("background reference ranges were loaded",
      any(str_detect(report, "Reference ranges from [0-9]+ previously mapped samples")))

cat("\n")
if (length(failures) > 0) {
    cat(glue("{length(failures)} assertion(s) failed:"), "\n", sep = "")
    walk(failures, \(x) cat("  ", x, "\n", sep = ""))
    cat("\nTool output is in ", runLog, "\n", sep = "")
    quit(save = "no", status = 1)
}

cat("All assertions passed.\n")
quit(save = "no", status = 0)
