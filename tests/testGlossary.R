#!/usr/bin/env Rscript
##
## Staleness test for docs/GLOSSARY.md.
##
## The glossary exists so that every number the tool emits can be traced back to
## a named field. A column added to an output without an entry in the glossary
## makes the file quietly wrong, and a hand-written document is accurate exactly
## once unless something checks it.
##
## This runs bin/wgsTriage.R against the two-sample fixture, reads the headers of
## both output tables, and asserts that every column name appears in
## docs/GLOSSARY.md. Add a column, add an entry.
##
## Usage:
##   Rscript tests/testGlossary.R
##
## Exits 0 if every column is documented, 1 otherwise.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Repo root from this script's own location, matching the bootstrap in bin/ and
## in testThresholds.R. Deliberately not here::here(): it resolves from the
## working directory and picks up whatever copy of the repo that happens to be.
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
glossaryFile <- path(repoRoot, "docs", "GLOSSARY.md")
outDir <- path(tempdir(), "wgsTriageGlossaryTest")

failures <- character()

check <- function(label, ok) {
    ok <- isTRUE(ok)
    if (!ok) failures <<- c(failures, label)
    cat(glue("  {if (ok) 'ok  ' else 'FAIL'}  {label}"), "\n", sep = "")
}

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

if (!file_exists(glossaryFile)) {
    cat("\n  docs/GLOSSARY.md does not exist.\n")
    quit(save = "no", status = 1)
}

samplesFile <- path(outDir, "wgsTriage_samples.tsv")
pairsFile <- path(outDir, "wgsTriage_pairs.tsv")

if (!file_exists(samplesFile) || !file_exists(pairsFile)) {
    cat("\n  Output tables were not written. Tool output follows:\n\n", sep = "")
    walk(read_lines(runLog), \(x) cat("  ", x, "\n", sep = ""))
    quit(save = "no", status = 1)
}

glossary <- read_file(glossaryFile)

##
## Only the header is needed, and the fixture cohort forms no pairs, so the pairs
## table is legitimately empty of rows. Reading with n_max = 0 gets the column
## names without asking readr to guess types from nothing.
##
columnsOf <- function(path) {
    read_tsv(path, n_max = 0, show_col_types = FALSE, progress = FALSE) |> names()
}

sampleColumns <- columnsOf(samplesFile)
pairColumns <- columnsOf(pairsFile)

check("sample table has columns", length(sampleColumns) > 0)
check("pair table has columns", length(pairColumns) > 0)

##
## Match on the backticked name rather than the bare string. `patient` is a word
## that occurs in prose throughout the file, so a bare search would report it as
## documented no matter what the file said about it.
##
documented <- function(column) str_detect(glossary, fixed(str_c("`", column, "`")))

undocumented <- function(columns) columns[!map_lgl(columns, documented)]

missingSample <- undocumented(sampleColumns)
missingPair <- undocumented(pairColumns)

check(glue("all {length(sampleColumns)} wgsTriage_samples.tsv columns are in docs/GLOSSARY.md"),
      length(missingSample) == 0)
if (length(missingSample) > 0) {
    walk(missingSample, \(x) cat("          undocumented: ", x, "\n", sep = ""))
}

check(glue("all {length(pairColumns)} wgsTriage_pairs.tsv columns are in docs/GLOSSARY.md"),
      length(missingPair) == 0)
if (length(missingPair) > 0) {
    walk(missingPair, \(x) cat("          undocumented: ", x, "\n", sep = ""))
}

cat("\n")
if (length(failures) > 0) {
    cat(glue("{length(failures)} assertion(s) failed:"), "\n", sep = "")
    walk(failures, \(x) cat("  ", x, "\n", sep = ""))
    cat("\nEvery column of both output tables needs an entry in docs/GLOSSARY.md.\n")
    cat("Tool output is in ", runLog, "\n", sep = "")
    quit(save = "no", status = 1)
}

cat("All assertions passed.\n")
quit(save = "no", status = 0)
