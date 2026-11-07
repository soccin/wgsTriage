#!/usr/bin/env Rscript
##
## Test for -m/--manifest.
##
## Two parts. The first calls selectNormals(), formPairs() and evaluatePairs()
## from R/qcLib.R directly on a synthetic cohort, one patient per rule, because
## the miniCohort fixture has two samples and cannot hold a patient with two
## normals and a tumor. The second runs bin/wgsTriage.R end to end on the
## fixture with manifests written here, and checks the class override, the pair
## verdicts, the notes on samples present on one side only, and the input
## errors.
##
## Usage:
##   Rscript tests/testManifest.R
##
## Exits 0 if every assertion holds, 1 otherwise.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Repo root from this script's own location, matching the bootstrap in bin/ and
## in the other tests. Deliberately not here::here(): it resolves from the
## working directory and can pick up a stale copy of qcLib.R.
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

fixtureDir <- path(repoRoot, "tests", "fixtures", "miniCohort")
backgroundDir <- path(repoRoot, "data", "background")
outRoot <- path(tempdir(), "wgsTriageManifestTest")
dir_create(outRoot)

failures <- character()

check <- function(label, ok) {
    ok <- isTRUE(ok)
    if (!ok) failures <<- c(failures, label)
    cat(glue("  {if (ok) 'ok  ' else 'FAIL'}  {label}"), "\n", sep = "")
}

##
## Part 1. Selection rules on a synthetic cohort.
##
## A  two normals, PASS 30x and WARN 50x       PASS wins on verdict over coverage
## B  two normals, PASS 30x and PASS 40x       coverage decides
## C  two normals, PASS 30x and PASS 30x       manifest order decides
## D  normal listed but not in the run         pair with it, INCOMPLETE
## E  FAIL normal in the run, another absent   absent ranks above FAIL
## F  tumor, no normal                         no pair
## G  normal, no tumor                         no pair
## H  one normal, FAIL                         used anyway, pair FAIL
##
cat("\n  Selection rules\n")

synth <- tribble(
    ~sample, ~patient, ~sampleType, ~verdict, ~meanCoverage,
    "A_N1",  "A",      "N",         "PASS",   30,
    "A_N2",  "A",      "N",         "WARN",   50,
    "A_T1",  "A",      "T",         "PASS",   90,
    "A_T2",  "A",      "T",         "PASS",   90,
    "B_N1",  "B",      "N",         "PASS",   30,
    "B_N2",  "B",      "N",         "PASS",   40,
    "B_T1",  "B",      "T",         "PASS",   90,
    "C_N1",  "C",      "N",         "PASS",   30,
    "C_N2",  "C",      "N",         "PASS",   30,
    "C_T1",  "C",      "T",         "PASS",   90,
    "D_T1",  "D",      "T",         "PASS",   90,
    "E_N1",  "E",      "N",         "FAIL",   30,
    "E_T1",  "E",      "T",         "PASS",   90,
    "F_T1",  "F",      "T",         "PASS",   90,
    "G_N1",  "G",      "N",         "PASS",   30,
    "H_N1",  "H",      "N",         "FAIL",   30,
    "H_T1",  "H",      "T",         "PASS",   90)

synthManifest <- synth |>
    select(sample, patient, sampleType) |>
    bind_rows(tribble(~sample, ~patient, ~sampleType,
                      "D_N1",  "D",      "N",
                      "E_N2",  "E",      "N")) |>
    mutate(manifestOrder = row_number())

selected <- selectNormals(synth, synthManifest)

normalFor <- function(tumor) selected |> filter(sampleTumor == tumor) |> pull(sampleNormal)
choiceFor <- function(tumor) selected |> filter(sampleTumor == tumor) |> pull(normalChoice) |> as.character()

check("one pair per tumor that has a normal (7)", nrow(selected) == 7)
check("A: PASS normal beats WARN normal with higher coverage",
      identical(normalFor("A_T1"), "A_N1") && identical(choiceFor("A_T1"), "best of 2 by verdict"))
check("A: both tumors share one normal", identical(normalFor("A_T2"), "A_N1"))
check("B: coverage decides between two PASS normals",
      identical(normalFor("B_T1"), "B_N2") && identical(choiceFor("B_T1"), "best of 2 by coverage"))
check("C: manifest order decides a full tie",
      identical(normalFor("C_T1"), "C_N1") && identical(choiceFor("C_T1"), "best of 2 by manifest order"))
check("D: normal absent from the run is still paired",
      identical(normalFor("D_T1"), "D_N1") && identical(choiceFor("D_T1"), "not in this Map run"))
check("E: absent normal ranks above a FAIL normal", identical(normalFor("E_T1"), "E_N2"))
check("F: tumor with no normal forms no pair", length(normalFor("F_T1")) == 0)
check("G: normal with no tumor forms no pair", !"G_N1" %in% selected$sampleNormal)
check("H: a lone FAIL normal is used", identical(choiceFor("H_T1"), "only normal"))

## Without a manifest every tumor x normal combination is formed. Kept as an
## assertion so the difference between the two modes is stated somewhere.
check("formPairs gives 2 x 2 = 4 pairs for patient A",
      nrow(formPairs(synth) |> filter(patient == "A")) == 4)

## No insertSizeAverage column here, which is also the no-samtools path.
evaluated <- evaluatePairs(synth, synth |> select(sample, verdict), selected)
pairOf <- function(tumor) evaluated |> filter(sampleTumor == tumor)

check("D: pair with an absent normal is INCOMPLETE",
      identical(pairOf("D_T1")$pairVerdict, "INCOMPLETE") &&
          identical(as.character(pairOf("D_T1")$pairReason), "normal not in this Map run"))
check("D: absent normal is still class N", identical(pairOf("D_T1")$sampleTypeNormal, "N"))
check("H: pair with a FAIL normal is FAIL", identical(pairOf("H_T1")$pairVerdict, "FAIL"))
check("A: two PASS samples with no insert size pass", identical(pairOf("A_T1")$pairVerdict, "PASS"))

##
## Part 2. End to end on the fixture. Both fixture samples are named as
## normals, so every pair below exists only because the manifest says so.
##
## Output is captured to a file rather than piped: piping Rscript into another
## process raises SIGPIPE and kills R mid-write.
##
runTool <- function(name, extra) {
    outDir <- path(outRoot, name)
    dir_create(outDir)
    runLog <- path(outDir, "run.log")
    status <- system2("Rscript",
                      c(path(repoRoot, "bin", "wgsTriage.R"),
                        fixtureDir,
                        "--background", backgroundDir,
                        "--out", outDir,
                        "--project", "miniCohort",
                        extra),
                      stdout = runLog, stderr = runLog)
    list(status = status, outDir = outDir, log = read_lines(runLog))
}

readOut <- function(run, file) {
    read_tsv(path(run$outDir, file), show_col_types = FALSE, progress = FALSE)
}

## Run 1: DEFECT_N01 reclassified as the tumor of CLEAN_N01, plus a second
## normal for the same patient that is not in the run. CLEAN_N01 passes, so it
## outranks the absent normal. sarek_patient is there to show extra columns are
## ignored.
cat("\n  Manifest run 1: class override, choice over an absent normal\n")

manifest1 <- path(outRoot, "manifest1.csv")
write_lines(c("sample,patient,status,sarek_patient",
              "CLEAN_N01,P1,0,P1",
              "DEFECT_N01,P1,1,P1",
              "GHOST_N02,P1,0,P1+N2"),
            manifest1)

run1 <- runTool("run1", c("-m", manifest1))
check("exit status is 0", run1$status == 0)

samples1 <- readOut(run1, "wgsTriage_samples.tsv")
pairs1 <- readOut(run1, "wgsTriage_pairs.tsv")
report1 <- read_lines(path(run1$outDir, "wgsTriage.txt"))
html1 <- read_file(path(run1$outDir, "wgsTriage.html"))

defect <- samples1 |> filter(sample == "DEFECT_N01")
check("DEFECT_N01 is class T from the manifest, despite its name", identical(defect$sampleType, "T"))
check("DEFECT_N01 takes patient P1 from the manifest", identical(defect$patient, "P1"))
check("DEFECT_N01 is held to the tumor coverage floor",
      isTRUE(defect$coverageFloor == COVERAGE_WARN[["T"]]))
check("exactly one pair", nrow(pairs1) == 1)
check("pair is DEFECT_N01 with CLEAN_N01",
      identical(pairs1$sampleTumor, "DEFECT_N01") && identical(pairs1$sampleNormal, "CLEAN_N01"))
check("pair fails because the tumor failed",
      identical(pairs1$pairVerdict, "FAIL") && identical(pairs1$pairReason, "tumor failed QC"))
check("present PASS normal chosen over the absent one",
      identical(pairs1$normalChoice, "best of 2 by verdict"))
check("console notes the manifest sample not in the run",
      any(str_detect(report1, fixed("1 manifest sample(s) are not in this Map run"))) &&
          any(str_detect(report1, fixed("GHOST_N02"))))
check("console lists GHOST_N02 nowhere as unused or unpaired",
      !any(str_detect(report1, "^  (Unpaired|Normals not used)")))
check("HTML states the pairing source is the manifest",
      str_detect(html1, fixed("Class and patient come from the --manifest file")))
check("HTML shows why the normal was chosen",
      str_detect(html1, fixed('<span class="r">best of 2 by verdict</span>')))

## Run 2: CLEAN_N01 as a tumor whose only normal is not in the run, and
## DEFECT_N01 left out of the manifest so it falls back to its name.
cat("\n  Manifest run 2: normal not in the run, sample not in the manifest\n")

manifest2 <- path(outRoot, "manifest2.csv")
write_lines(c("sample,patient,status",
              "CLEAN_N01,P2,1",
              "GHOST_N02,P2,0"),
            manifest2)

run2 <- runTool("run2", c("--manifest", manifest2))
check("exit status is 0", run2$status == 0)

samples2 <- readOut(run2, "wgsTriage_samples.tsv")
pairs2 <- readOut(run2, "wgsTriage_pairs.tsv")
report2 <- read_lines(path(run2$outDir, "wgsTriage.txt"))

check("DEFECT_N01 keeps its name-based class N",
      identical(samples2 |> filter(sample == "DEFECT_N01") |> pull(sampleType), "N"))
check("pair is CLEAN_N01 with the absent GHOST_N02",
      nrow(pairs2) == 1 && identical(pairs2$sampleNormal, "GHOST_N02"))
check("pair is INCOMPLETE with the reason stated",
      identical(pairs2$pairVerdict, "INCOMPLETE") &&
          identical(pairs2$pairReason, "normal not in this Map run"))
check("absent normal has no verdict", is.na(pairs2$verdictNormal))
check("console counts the pair as not assessed",
      any(str_detect(report2, fixed("1 of 1 pairs could not be assessed."))))
check("console notes the sample not in the manifest",
      any(str_detect(report2, fixed("1 sample(s) are not in the manifest"))))
check("console lists DEFECT_N01 as unpaired",
      any(str_detect(report2, fixed("Unpaired samples (1):"))))

## Input errors. Each must stop before writing a report.
cat("\n  Input errors\n")

pairingDummy <- path(outRoot, "pairing.tsv")
write_lines(c("NORMAL_ID\tTUMOR_ID", "CLEAN_N01\tDEFECT_N01"), pairingDummy)
run3 <- runTool("run3", c("-m", manifest1, "--pairing", pairingDummy))
check("--manifest with --pairing is refused",
      run3$status != 0 && any(str_detect(run3$log, fixed("cannot be used together"))))

manifestBad <- path(outRoot, "manifestBad.csv")
write_lines(c("sample,patient,status", "CLEAN_N01,P1,N", "DEFECT_N01,P1,T"), manifestBad)
run4 <- runTool("run4", c("-m", manifestBad))
check("status other than 0 or 1 is refused",
      run4$status != 0 && any(str_detect(run4$log, fixed("status other than 1 (tumor) or 0 (normal)"))))

manifestElsewhere <- path(outRoot, "manifestElsewhere.csv")
write_lines(c("sample,patient,status", "OTHER_T01,P9,1", "OTHER_N01,P9,0"), manifestElsewhere)
run5 <- runTool("run5", c("-m", manifestElsewhere))
check("manifest with no sample in the run is refused",
      run5$status != 0 && any(str_detect(run5$log, fixed("No sample in manifest"))))

cat("\n")
if (length(failures) > 0) {
    cat(glue("{length(failures)} assertion(s) failed:"), "\n", sep = "")
    walk(failures, \(x) cat("  ", x, "\n", sep = ""))
    cat("\nTool output is under ", outRoot, "\n", sep = "")
    quit(save = "no", status = 1)
}

cat("All assertions passed.\n")
quit(save = "no", status = 0)
