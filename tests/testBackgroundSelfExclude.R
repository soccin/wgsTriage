#!/usr/bin/env Rscript
##
## Test that the background importer does not import its own test fixtures.
##
## Every wgsTriage checkout carries tests/fixtures/miniCohort with real metric
## file names. The 2026-09-27 background build scanned a working tree holding
## ten clones and imported the fixtures as a project called miniCohort. This
## builds a small archive in a temp directory that mixes real-looking cohorts
## with checkout and fixture copies, runs bin/wgsTriageBackground.R over it, and
## asserts which files were kept and which were excluded.
##
## Archive layout (all metric files are copies of the miniCohort fixtures):
##   CohortA/                            kept, a plain cohort
##   wgsTriage_tests/                    kept, name contains the words only
##   Users/X/Map/wgsTriage/tests/...     excluded, clone and fixture rules
##   Users/X/Map/wgsTriage/QCData/...    excluded, inside a clone below <QCDir>
##   scratch/tests/fixtures/...          excluded, fixture rule only
##   ProjB/wgsTriage/out/metrics/...     excluded, clone rule only
##
## A second run points <QCDir> at the QCData directory inside the clone, the
## conventional ./QCData location, and asserts that it is imported: a wgsTriage
## directory above <QCDir> must not exclude the archive itself.
##
## Usage:
##   Rscript tests/testBackgroundSelfExclude.R
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
testDir <- path(tempdir(), "selfExcludeTest")
archive <- path(testDir, "archive")

failures <- character()

check <- function(label, ok) {
    ok <- isTRUE(ok)
    if (!ok) failures <<- c(failures, label)
    cat(glue("  {if (ok) 'ok  ' else 'FAIL'}  {label}"), "\n", sep = "")
}

## Copy one fixture file to a path relative to the archive.
place <- function(fixtureRel, archiveRel) {
    dest <- path(archive, archiveRel)
    dir_create(path_dir(dest))
    file_copy(path(fixtureDir, fixtureRel), dest, overwrite = TRUE)
}

cleanAsm <- "out/metrics/CLEAN_N01/CLEAN_N01.asm.txt"
cleanWgs <- "out/metrics/CLEAN_N01/CLEAN_N01.wgs.txt"
defectAsm <- "out/metrics/DEFECT_N01/DEFECT_N01.asm.txt"
defectWgs <- "out/metrics/DEFECT_N01/DEFECT_N01.wgs.txt"
multiqc <- "sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt"

if (dir_exists(testDir)) dir_delete(testDir)

place(cleanAsm, path("CohortA", cleanAsm))
place(cleanWgs, path("CohortA", cleanWgs))
place(multiqc, path("CohortA", multiqc))
place(defectAsm, path("wgsTriage_tests", defectAsm))
place(defectWgs, path("wgsTriage_tests", defectWgs))

clone <- path("Users", "X", "Map", "wgsTriage")
dir_copy(fixtureDir, path(archive, clone, "tests", "fixtures", "miniCohort"))
place(cleanAsm, path(clone, "QCData", "CohortQ", cleanAsm))
place(cleanWgs, path(clone, "QCData", "CohortQ", cleanWgs))
place(cleanAsm, path("scratch", "tests", "fixtures", "miniCohort", cleanAsm))
place(cleanWgs, path("ProjB", "wgsTriage", cleanWgs))

## 5 fixture files in the clone, 2 in its QCData, 1 in scratch, 1 in ProjB.
expectedExcluded <- 9

##
## Output is captured to a file rather than piped: piping Rscript into another
## process raises SIGPIPE and kills R mid-write.
##
runBackground <- function(qcDir, outDir) {
    dir_create(outDir)
    runLog <- path(outDir, "run.log")
    status <- system2("Rscript",
                      c(path(repoRoot, "bin", "wgsTriageBackground.R"),
                        qcDir, "--out", outDir),
                      stdout = runLog, stderr = runLog)
    list(status = status, log = read_lines(runLog, progress = FALSE), runLog = runLog)
}

readOut <- function(outDir, name) {
    read_tsv(path(outDir, name), show_col_types = FALSE, progress = FALSE)
}

## ---------------------------------------------------------------------------
## Run 1: archive holding cohorts, clones and fixture copies
## ---------------------------------------------------------------------------

cat("\n  Archive with checkouts and fixtures\n")
out1 <- path(testDir, "out1")
run1 <- runBackground(archive, out1)
check("exit status is 0", run1$status == 0)

if (!file_exists(path(out1, "backgroundImportAudit.tsv"))) {
    cat("\n  No audit written. Tool output follows:\n\n", sep = "")
    walk(run1$log, \(x) cat("  ", x, "\n", sep = ""))
    quit(save = "no", status = 1)
}

samples1 <- readOut(out1, "backgroundSamples.tsv")
audit1 <- readOut(out1, "backgroundImportAudit.tsv")
excluded1 <- audit1 |> filter(disposition == "excludedSelf")
kept1 <- audit1 |> filter(disposition != "excludedSelf")

check("no sample has project miniCohort",
      !any(samples1$project == "miniCohort"))
check("projects are exactly CohortA and wgsTriage_tests",
      setequal(samples1$project, c("CohortA", "wgsTriage_tests")))
check(glue("{expectedExcluded} files recorded as excludedSelf"),
      nrow(excluded1) == expectedExcluded)
check("excluded files are marked not parsed",
      all(!excluded1$parsed))
check("no kept file lies under a wgsTriage directory",
      !any(str_detect(kept1$path, "/wgsTriage/")))
check("no kept file lies under tests/fixtures",
      !any(str_detect(kept1$path, "/tests/fixtures/")))
check("console summary reports the exclusion",
      any(str_detect(run1$log, "excluded \\(self\\)")))
check("excluded files are not counted as unreadable",
      !any(str_detect(run1$log, "unreadable files")))

## ---------------------------------------------------------------------------
## Run 2: <QCDir> is the QCData directory inside a clone
## ---------------------------------------------------------------------------

cat("\n  Archive at <clone>/QCData\n")
out2 <- path(testDir, "out2")
run2 <- runBackground(path(archive, clone, "QCData"), out2)
check("exit status is 0", run2$status == 0)

samples2 <- if (file_exists(path(out2, "backgroundSamples.tsv"))) {
    readOut(out2, "backgroundSamples.tsv")
} else {
    tibble(project = character(), sample = character())
}
audit2 <- if (file_exists(path(out2, "backgroundImportAudit.tsv"))) {
    readOut(out2, "backgroundImportAudit.tsv")
} else {
    tibble(disposition = character())
}

check("CohortQ CLEAN_N01 is imported",
      any(samples2$project == "CohortQ" & samples2$sample == "CLEAN_N01"))
check("nothing is excluded",
      nrow(audit2) > 0 && !any(audit2$disposition == "excludedSelf"))

cat("\n")
if (length(failures) > 0) {
    cat(glue("{length(failures)} assertion(s) failed:"), "\n", sep = "")
    walk(failures, \(x) cat("  ", x, "\n", sep = ""))
    cat("\nTool output is in ", run1$runLog, " and ", run2$runLog, "\n", sep = "")
    quit(save = "no", status = 1)
}

cat("All assertions passed.\n")
quit(save = "no", status = 0)
