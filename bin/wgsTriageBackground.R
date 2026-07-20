#!/usr/bin/env Rscript
##
## Build the background QC distribution from previously mapped projects.
##
## Run with --help for usage. Writes six files; three are name-free aggregates
## that are committed, three carry sample names and are gitignored.
##
## The reference ranges are computed only from samples that pass the filter thresholds.
## The archive is known to contain defective cohorts, and including them would
## widen the reference range enough to admit the next bad cohort. Robust
## statistics alone are not sufficient protection when the contaminated
## fraction is large.
##
## Import contract: this script takes whatever the archive happens to look like
## and reads as much of it as can be read. Directory layout varies between
## cohorts and even within them, files go missing, and the same cohort turns up
## copied to a second location. None of that is an error. Anything the importer
## cannot use, cannot match, or had to choose between is counted and reported
## rather than dropped in silence, because a background that quietly shrank is
## indistinguishable from a background that was always small.
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
defaultOut <- path(repoRoot, "data", "background")

usage <- function() {
    glue("
wgsTriageBackground.R {WGSTRIAGE_VERSION} -- build the background QC
distribution that gives wgsTriage.R its historical reference ranges.

Usage:
  Rscript bin/wgsTriageBackground.R <QCDir> [--out <OutDir>]
  Rscript bin/wgsTriageBackground.R --help

Arguments:
  <QCDir>           Archive of previously mapped projects. Searched recursively
                    for Picard <sample>.asm.txt and <sample>.wgs.txt, and for
                    multiqc_samtools_stats.txt. Layout is not assumed: any
                    arrangement of wrapper directories is accepted, and samples
                    missing some of their files are still imported. Read only;
                    nothing is modified.
                    Default: ./QCData

Options:
  --out <OutDir>    Directory for the six output files.
                    Default: <repoRoot>/data/background
  -h, --help        Show this message and exit.

Writes into <OutDir>, overwriting all six on every run:
  backgroundStats.tsv           reference ranges          no names   committed
  backgroundCoverageStats.tsv   coverage by sample class  no names   committed
  backgroundMetricCoverage.tsv  per-metric availability   no names   committed
  backgroundSamples.tsv         every sample and verdict  NAMES      gitignored
  backgroundFlagged.tsv         samples below threshold   NAMES      gitignored
  backgroundImportAudit.tsv     every file, and its fate  NAMES      gitignored

Only backgroundStats.tsv is read back by wgsTriage.R. Expect 1 to 2 minutes
over a roughly 1.4 GB archive.
")
}

if (any(args %in% c("-h", "--help"))) {
    cat(usage(), "\n", sep = "")
    quit(save = "no", status = 0)
}

##
## Argument parsing walks the vector by position. Matching option values by
## string equality, as this did before, misparses the command line whenever the
## --out directory and the QC directory are spelled the same way.
##
parseArgs <- function(args) {
    out <- NULL
    positional <- character(0)
    i <- 1
    while (i <= length(args)) {
        a <- args[i]
        if (a == "--out") {
            if (i == length(args)) stop("--out requires a directory argument. See --help.", call. = FALSE)
            out <- args[i + 1]
            i <- i + 2
        } else if (str_starts(a, "--out=")) {
            out <- str_remove(a, "^--out=")
            i <- i + 1
        } else if (str_starts(a, "-")) {
            stop(glue("Unknown option: {a}. See --help."), call. = FALSE)
        } else {
            positional <- c(positional, a)
            i <- i + 1
        }
    }
    list(out = out, positional = positional)
}

parsed <- parseArgs(args)
outDir <- parsed$out %||% defaultOut
qcRoot <- if (length(parsed$positional) > 0) parsed$positional[1] else "QCData"

if (!dir_exists(qcRoot)) {
    stop(glue("
QCData root not found: {qcRoot}

  Usage: Rscript bin/wgsTriageBackground.R <QCDir> [--out <OutDir>]

  <QCDir> is the archive of previously mapped projects to scan.
  Run with --help for the full description and the list of outputs.
"), call. = FALSE)
}
dir_create(outDir)

##
## Warnings are collected rather than emitted through warning(), which at the
## end of an Rscript run prints after the summary and is routinely missed. They
## are printed as their own block below the summary.
##
importWarnings <- character(0)
warnImport <- function(...) {
    importWarnings <<- c(importWarnings, glue(..., .envir = parent.frame()))
}

## ---------------------------------------------------------------------------
## Discovery
## ---------------------------------------------------------------------------

##
## Directory names that are pipeline scaffolding rather than cohort identity,
## matched case insensitively against one path component at a time.
##
## The project label has to come out the same for a sample's alignment metrics,
## its coverage metrics and its multiqc table, because that label is half the
## join key. The previous implementation derived it with a fixed sequence of
## path edits, which assumed the wrapper directories always appeared in the same
## order and at the same depth. They do not. The same cohort appears as
## <proj>/out/metrics, <proj>/Map/out/metrics, <proj>/Normal03/out/metrics and
## <proj>/results/r_002/mapping/metrics, its multiqc lives under
## <proj>/Set1/sbam/multiqc/multiqc_data, and an rsync -R copy buries the whole
## thing under a replica of the source absolute path. Walking up and discarding
## scaffolding is stable under all of those; string surgery on the tail is not.
##
SCAFFOLD_DIR <- str_c(
    "^(",
    str_c(c("out", "output", "metrics", "picard", "qc", "stats",
            "sbam(\\.[0-9]+)?", "bam", "bams",
            "map", "smap", "mapping", "align", "alignment", "alignments",
            "multiqc", "multiqc_data",
            "result", "results", "r_?[0-9]+",
            "module[0-9]*",
            "set[0-9]+", "batch[0-9]+", "b[0-9]+", "grp[0-9]+", "group[0-9]+",
            "normal[0-9]+", "tumor[0-9]+", "tumour[0-9]+",
            "down[a-z0-9_]*"),
          collapse = "|"),
    ")$")

## File name patterns. The canonical archive spells these .asm.txt and .wgs.txt;
## the alternatives are Picard's own default output names, accepted so a
## differently configured pipeline does not silently contribute nothing.
ASM_PATTERN <- "\\.(asm|alignment_summary_metrics)(\\.txt)?$"
WGS_PATTERN <- "\\.(wgs|wgs_metrics|raw_wgs_metrics)(\\.txt)?$"
MQC_PATTERN <- "multiqc_samtools_stats\\.txt$"

## Sample name is the file name with the metric extension and any pipeline stage
## suffixes removed. Stage suffixes stack (.md.recal), so strip to a fixed point
## rather than once. Sample names themselves contain dots (CTCL.AM_CTCL26), so
## only the known stage tokens are removed, never everything after the first dot.
STAGE_SUFFIX <- "\\.(smap|md|recal|sorted|dedup|markdup|bqsr|final)$"

rawSampleFromFile <- function(paths) {
    path_file(paths) |>
        str_remove(ASM_PATTERN) |>
        str_remove(WGS_PATTERN)
}

sampleFromFile <- function(paths) {
    base <- rawSampleFromFile(paths)
    repeat {
        stripped <- str_remove(base, STAGE_SUFFIX)
        if (identical(stripped, base)) break
        base <- stripped
    }
    base
}

##
## Project label: walk up from the file, discarding scaffolding directories and
## the per-sample directory, and take the first real name that remains.
##
projectFromPath <- function(filePath, sample, rawSample, root) {
    parts <- path_rel(filePath, root) |> str_split_1("/")
    fileName <- parts[length(parts)]
    dirs <- parts[-length(parts)]
    while (length(dirs) > 0) {
        last <- dirs[length(dirs)]
        ## Only exact matches count as the per-sample directory. A prefix rule
        ## would eat a cohort directory whose name happens to start the sample
        ## names it contains, which is common (UMich holding Umich10_T).
        ## fileName covers the layout where the directory repeats the file name
        ## in full, extension included.
        isScaffold <- str_detect(str_to_lower(last), SCAFFOLD_DIR) ||
            last == sample || last == rawSample || last == fileName
        if (!isScaffold) break
        dirs <- dirs[-length(dirs)]
    }
    if (length(dirs) == 0) "root" else dirs[length(dirs)]
}

##
## One traversal of the archive, then classify by file name. Three separate
## globbed traversals cost three walks of a 1.4 GB tree and, more importantly,
## made it impossible to report on files that matched nothing.
##
cat(glue("Scanning {qcRoot} for Picard metrics ...\n\n"))

allFiles <- tryCatch(
    dir_ls(qcRoot, recurse = TRUE, type = c("file", "symlink"), fail = FALSE) |> as.character(),
    error = function(e) {
        warnImport("Directory scan hit an error and may be incomplete: {conditionMessage(e)}")
        character(0)
    })

asmFiles <- str_subset(allFiles, ASM_PATTERN)
wgsFiles <- str_subset(allFiles, WGS_PATTERN)
mqcFiles <- str_subset(allFiles, MQC_PATTERN)

## ---------------------------------------------------------------------------
## Parsing
## ---------------------------------------------------------------------------

##
## Read a set of metrics files into one row per file, keeping a record for every
## file including the ones that could not be read. A file that fails to parse is
## evidence about the archive and belongs in the audit, not in a dropped row.
##
collectMetricFiles <- function(paths, reader, kind, root) {
    empty <- tibble(kind = character(), path = character(), project = character(),
                    sample = character(), mtime = as.POSIXct(character()),
                    parsed = logical(), disposition = character(), detail = character())
    if (length(paths) == 0) return(list(data = tibble(), audit = empty))

    sample <- sampleFromFile(paths)
    rawSample <- rawSampleFromFile(paths)
    project <- pmap_chr(list(paths, sample, rawSample),
                        \(p, s, r) projectFromPath(p, s, r, root))
    info <- file_info(paths)

    ## Errors and warnings are both captured per file and neither is allowed to
    ## reach the console directly. A malformed file produces both, and letting
    ## them print interleaves a wall of vroom and file() noise ahead of the
    ## summary while saying nothing about which file caused it. The audit says
    ## which file, which is the part worth having.
    attempt <- map(paths, safely(quietly(reader)))
    errMsg <- map_chr(attempt, \(x) if (is.null(x$error)) NA_character_ else conditionMessage(x$error))
    result <- map(attempt, \(x) if (is.null(x$result)) NULL else x$result$result)
    warnMsg <- map_chr(attempt, \(x) {
        if (is.null(x$result) || length(x$result$warnings) == 0) return(NA_character_)
        str_c(unique(x$result$warnings), collapse = "; ")
    })
    usable <- map_lgl(result, \(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0)

    ## Whether a file could be read and what then became of it are separate
    ## facts. Folding them into one column makes a file that parsed perfectly
    ## and was then superseded by a newer copy indistinguishable from a file
    ## that could not be read at all.
    audit <- tibble(
        kind = kind,
        path = paths,
        project = project,
        sample = sample,
        mtime = info$modification_time,
        parsed = usable & is.na(errMsg),
        disposition = case_when(!is.na(errMsg) ~ "parseError",
                                !usable       ~ "noMetricsFound",
                                .default      = "imported"),
        ## Squished because R's parser warnings are multi-line, and a newline
        ## inside a field silently splits one audit row into several when the
        ## file is read back as TSV.
        detail = str_squish(coalesce(errMsg, warnMsg, "")))

    ## multiqc tables name their own samples, one file holding many; Picard files
    ## name one sample in the file name. Carry the file-derived name separately
    ## and let a name from inside the file win, so both shapes come out with a
    ## single unambiguous sample column.
    data <- if (any(usable)) {
        out <- tibble(kind = kind, path = paths, project = project,
                      fileSample = sample, mtime = info$modification_time,
                      parsed = result) |>
            filter(usable) |>
            unnest(parsed)
        if ("sample" %in% names(out)) out |> select(-fileSample) else out |> rename(sample = fileSample)
    } else {
        tibble()
    }

    list(data = data, audit = audit)
}

##
## Read a multiqc samtools stats table.
##
## Deliberately not qcLib's readMultiqcSamtools, which this script used before.
## That reader keeps only rows suffixed .recal and discards the whole file if any
## one of the six wanted columns is absent. Across this archive 512 of 611
## samtools rows carry .md and no .recal counterpart at all, so requiring .recal
## threw away 84% of the samtools evidence and left supplementaryRate, the
## tightest metric in the whole filter set, resting on 98 samples.
##
## .recal is still preferred where it exists: it is the final BAM state and the
## same object Picard measured. .md is the same library one step earlier and is
## a sound fallback for these metrics. Which one was used is recorded per sample
## so a reader can tell the two apart, and mixing is reported.
##
readMultiqcSamtoolsAny <- function(path) {
    dat <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    if (!"Sample" %in% names(dat) || nrow(dat) == 0) return(NULL)

    ## Absent columns become an all-NA column of the right length rather than
    ## sinking the file. A table missing insert size still carries a usable
    ## supplementary alignment rate.
    column <- function(nm) {
        if (nm %in% names(dat)) suppressWarnings(as.numeric(dat[[nm]])) else rep(NA_real_, nrow(dat))
    }
    total <- column("raw_total_sequences")

    tibble(
        rawSample         = as.character(dat$Sample),
        samtoolsStage     = case_when(str_detect(rawSample, "\\.recal$") ~ "recal",
                                      str_detect(rawSample, "\\.md$")    ~ "md",
                                      .default                           = "raw"),
        sample            = str_remove(rawSample, "\\.(md|recal)$"),
        supplementaryRate = column("supplementary_alignments") / total * 100,
        interChromRate    = column("pairs_on_different_chromosomes") / total * 100,
        insertSizeAverage = column("insert_size_average"),
        pctProperlyPaired = column("reads_properly_paired_percent"),
        pctMapped         = column("reads_mapped_percent")) |>
        mutate(stageRank = match(samtoolsStage, c("recal", "md", "raw"))) |>
        arrange(stageRank) |>
        distinct(sample, .keep_all = TRUE) |>
        select(-stageRank, -rawSample)
}

asmCollected <- collectMetricFiles(asmFiles, readAsmMetrics, "asm", qcRoot)
wgsCollected <- collectMetricFiles(wgsFiles, readWgsMetrics, "wgs", qcRoot)
mqcCollected <- collectMetricFiles(mqcFiles, readMultiqcSamtoolsAny, "multiqc", qcRoot)

## ---------------------------------------------------------------------------
## Duplicate resolution
## ---------------------------------------------------------------------------

##
## The same sample can be present more than once under one project label, most
## often because a cohort was remapped into a second directory that survives
## alongside the first. Newest file wins, and every collision is recorded so the
## choice can be checked. Ties break on path so repeated runs agree.
##
resolveDuplicates <- function(dat, kind) {
    if (nrow(dat) == 0) return(list(data = dat, dups = tibble()))
    repeated <- dat |>
        summarise(n = n(), .by = c(project, sample)) |>
        filter(n > 1)
    dups <- dat |>
        semi_join(repeated, by = c("project", "sample")) |>
        select(kind, project, sample, path, mtime) |>
        arrange(project, sample, desc(mtime), path) |>
        mutate(kept = row_number() == 1, .by = c(project, sample))
    kept <- dat |>
        arrange(desc(mtime), path) |>
        distinct(project, sample, .keep_all = TRUE)
    if (nrow(repeated) > 0) {
        warnImport("{nrow(repeated)} {kind} sample(s) had more than one file under the same project; kept the newest. See backgroundImportAudit.tsv.")
    }
    list(data = kept, dups = dups)
}

asmResolved <- resolveDuplicates(asmCollected$data, "asm")
wgsResolved <- resolveDuplicates(wgsCollected$data, "wgs")
mqcResolved <- resolveDuplicates(mqcCollected$data, "multiqc")

##
## Per-kind columns are named apart before joining so the only shared columns
## are the join keys and nothing gets silently suffixed .x / .y. A kind that
## matched no file at all yields a zero-row frame carrying just the join keys,
## which the joins then absorb; dropping columns from a frame that has none is
## an error, and an archive holding only coverage metrics is a normal thing for
## this tool to be pointed at.
##
labelKind <- function(dat, pathName, mtimeName) {
    if (nrow(dat) == 0) return(tibble(project = character(), sample = character()))
    dat |> select(-kind) |> rename("{pathName}" := path, "{mtimeName}" := mtime)
}

asm <- labelKind(asmResolved$data, "asmPath", "asmMtime")
wgs <- labelKind(wgsResolved$data, "wgsPath", "wgsMtime")
mqc <- labelKind(mqcResolved$data, "mqcPath", "mqcMtime")

## ---------------------------------------------------------------------------
## Joining
## ---------------------------------------------------------------------------

##
## Join on (project, sample) and keep whatever finds no partner as its own row.
##
## An earlier version of this fix also made a second pass that matched leftovers
## on sample name alone, on the theory that a differing project label meant the
## same run had been filed in two places. Checking the Picard headers of the
## samples that pass caught shows that is not what it means here.
## NK_KHYG1_CL_D, for one, has alignment metrics taken from
## sbam/bam/NK_KHYG1_CL_D.smap.bam against human_g1k_v37_decoy in November 2025
## and coverage metrics taken from a markduplicates .md.cram against b37 in
## February 2026. Those are two separate mappings of one biological sample, and
## merging them would invent a sample whose chimera rate and whose coverage were
## never measured on the same alignment.
##
## Nothing is lost by leaving them apart. Reference ranges are built per metric
## from whatever is not NA, so both rows contribute exactly what they measured.
## Name collisions across differing project labels are reported instead, because
## a genuine split of one run across two labels looks the same from here and is
## worth a human deciding.
##
joinOnProjectAndSample <- function(left, right) {
    ## Typed rather than bare tibble(), so the caller can read crossProject$sample
    ## without tripping an uninitialised-column warning when one side is empty.
    noCross <- tibble(sample = character(), projectLeft = character(),
                      projectRight = character())
    if (nrow(left) == 0) return(list(data = right, crossProject = noCross))
    if (nrow(right) == 0) return(list(data = left, crossProject = noCross))

    keys <- c("project", "sample")
    matched <- inner_join(left, right, by = keys)
    leftOnly <- anti_join(left, right, by = keys)
    rightOnly <- anti_join(right, left, by = keys)

    ## Same sample name on both sides, but filed under different projects.
    crossProject <- inner_join(
        leftOnly |> distinct(sample, projectLeft = project),
        rightOnly |> distinct(sample, projectRight = project),
        by = "sample")

    list(data = bind_rows(matched, leftOnly, rightOnly), crossProject = crossProject)
}

picardJoin <- joinOnProjectAndSample(asm, wgs)
picard <- picardJoin$data

##
## Only a total absence is fatal. An archive holding coverage metrics but no
## alignment metrics, or nothing but multiqc tables, still has a background
## worth building out of what it has, and refusing it would be the same
## all-or-nothing behaviour this rewrite exists to remove.
##
if (nrow(picard) == 0 && nrow(mqc) == 0) {
    found <- glue("{length(asmFiles)} alignment, {length(wgsFiles)} coverage, {length(mqcFiles)} multiqc")
    reason <- if (length(asmFiles) + length(wgsFiles) + length(mqcFiles) == 0) {
        "No candidate files were found at all."
    } else {
        glue("Candidate files were found ({found}) but none of them could be read; see backgroundImportAudit.tsv.")
    }
    stop(glue("
No usable QC metrics found under {qcRoot}

  {reason}

  Expected, in any arrangement of directories:
    <sample>.asm.txt   Picard CollectAlignmentSummaryMetrics
    <sample>.wgs.txt   Picard CollectWgsMetrics
    multiqc_samtools_stats.txt

  Run with --help for the full description.
"), call. = FALSE)
}

samtoolsJoin <- joinOnProjectAndSample(picard, mqc)
background <- samtoolsJoin$data

crossProject <- bind_rows(
    picardJoin$crossProject |> mutate(pairing = "alignment/coverage"),
    samtoolsJoin$crossProject |> mutate(pairing = "picard/samtools"))

if (nrow(crossProject) > 0) {
    names <- crossProject |> distinct(sample) |> pull(sample)
    shown <- str_c(head(names, 6), collapse = ", ")
    more <- if (length(names) > 6) glue(" and {length(names) - 6} more") else ""
    warnImport("{length(names)} sample name(s) appear under more than one project label and were left as separate rows, not merged: {shown}{more}. Check whether these are one run filed twice or genuinely separate mappings.")
}

##
## A sample that only ever appeared in a multiqc table is still a sample. It
## carries real samtools metrics and contributes to those reference ranges; it
## simply has no Picard files. Dropping it would understate the one background
## that is already thinnest.
##
background <- background |>
    mutate(sampleType = classifySampleType(sample),
           patient = patientStem(sample)) |>
    relocate(project, sample, sampleType, patient)

## The join is constructed to produce one row per (project, sample). Check it
## rather than assume it, because a silent duplicate here would double-weight a
## sample in every reference range downstream.
keyDupes <- background |> summarise(n = n(), .by = c(project, sample)) |> filter(n > 1)
if (nrow(keyDupes) > 0) {
    warnImport("{nrow(keyDupes)} (project, sample) key(s) are duplicated after joining; reference ranges will double-count them.")
}

##
## Mixing BAM stages is a real, and small, inconsistency. Measured on the 99
## samples in this archive that carry both rows, recal runs below md by at most
## 0.055 points of supplementary rate, 0.52 of properly paired percent and 1.3
## of insert size. That is around 3% of the 1.0% supplementary threshold and
## does not move a clean sample anywhere near it, but the reference range does
## blend two stages and a reader is entitled to know. The alternative, keeping
## recal only, is worse than mixed: every cohort with an elevated supplementary
## rate in this archive exists as md rows alone, so a recal-only background
## contains no failing sample at all and the threshold ends up validated against
## evidence that excludes the thing it is meant to catch.
##
if ("samtoolsStage" %in% names(background)) {
    stageMix <- background |> filter(!is.na(samtoolsStage)) |> distinct(samtoolsStage) |> nrow()
    if (stageMix > 1) {
        warnImport("samtools metrics come from more than one BAM stage; recal was preferred where present and md used otherwise. Per-sample stage is in the samtoolsStage column of backgroundSamples.tsv.")
    }
}

## ---------------------------------------------------------------------------
## Threshold evaluation
## ---------------------------------------------------------------------------

##
## Every threshold metric is materialised, as an all-NA column where the archive
## has no such data. Without this a metric that is absent everywhere is not
## evaluated at all and therefore never reported as MISSING, which reads as
## though it had been checked and passed.
##
ensureNumericColumns <- function(dat, cols) {
    for (m in setdiff(cols, names(dat))) dat <- dat |> mutate("{m}" := NA_real_)
    dat
}

background <- background |>
    ensureNumericColumns(c(THRESHOLDS$metric, "meanCoverage", "insertSizeAverage",
                           "pctImproperPairs", "interChromRate", "pctExcDupe",
                           "pctSoftclip", "alignedFrac", "pctReadUsed", "pctExcTotal"))

##
## Apply the filter thresholds to every historical sample. This is the same
## evaluation the report applies to a current cohort, so the background carries
## verdicts on the same terms.
##
## The threshold functions key on a single column, so pass them a
## project-qualified id. Keying on the bare sample name merges verdicts across
## unrelated runs and silently assigns one project's failure to another
## project's sample.
##
background <- background |> mutate(uid = str_c(project, "::", sample))

thresholdResults <- evaluateThresholds(background |> select(-sample) |> rename(sample = uid))
verdicts <- sampleVerdict(thresholdResults) |> rename(uid = sample)

evaluable <- thresholdResults |>
    summarise(nEvaluable = sum(status != "MISSING"), .by = sample) |>
    rename(uid = sample)

##
## A reference sample is one that failed nothing it could be checked on.
##
## This previously also required both Picard files to be present, which sounds
## conservative and is not. It means a sample missing one of its two files
## contributes to no reference range at all, not even the ranges built from the
## file it does have. On an archive where the alignment and coverage metrics
## were collected into different subtrees that rule reduced 401 usable coverage
## samples to 2, and a coverage range built from 2 samples is not a range. The
## gate is what a sample failed, not which files it happens to have; a sample
## with only coverage metrics is still checked against every coverage threshold.
## referenceTier records which samples carry the full pair so the distinction
## stays visible in the outputs.
##
background <- background |>
    left_join(verdicts, by = "uid") |>
    left_join(evaluable, by = "uid") |>
    mutate(nEvaluable = coalesce(nEvaluable, 0L),
           hasCore = !is.na(pctChimeras) & !is.na(meanCoverage),
           referenceSample = nFail == 0 & nEvaluable > 0,
           referenceSample = coalesce(referenceSample, FALSE),
           referenceTier = if_else(hasCore, "full", "partial"))

##
## Robust reference ranges from clean samples only.
##
metricCols <- c(THRESHOLDS$metric, "meanCoverage", "insertSizeAverage",
                "pctImproperPairs", "interChromRate", "pctExcDupe")
metricCols <- metricCols[metricCols %in% names(background)]

emptyStats <- tibble(metric = character(), n = integer(), median = numeric(),
                     mad = numeric(), q01 = numeric(), q05 = numeric(),
                     q25 = numeric(), q75 = numeric(), q95 = numeric(),
                     q99 = numeric(), min = numeric(), max = numeric())

referenceLong <- if (length(metricCols) == 0) {
    tibble(sampleType = character(), metric = character(), value = numeric())
} else {
    background |>
        filter(referenceSample) |>
        select(sampleType, all_of(metricCols)) |>
        pivot_longer(-sampleType, names_to = "metric", values_to = "value") |>
        filter(!is.na(value))
}

referenceStats <- if (nrow(referenceLong) == 0) emptyStats else {
    referenceLong |>
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
}

## Coverage is the one metric where tumor and normal genuinely differ, so it
## also gets a per-class reference. Everything else measures data integrity and
## should not vary with sample class.
coverageInput <- background |> filter(referenceSample, !is.na(meanCoverage))
coverageStats <- if (nrow(coverageInput) == 0) {
    tibble(metric = character(), sampleType = character(), n = integer(),
           median = numeric(), mad = numeric(), q05 = numeric(),
           q25 = numeric(), q75 = numeric(), q95 = numeric())
} else {
    coverageInput |>
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
}

##
## Selected on the counts rather than on the verdict. A sample carrying warnings
## is reported as INCOMPLETE the moment any one metric is missing, and since
## pctReadUsed is absent from every archived file that describes most of the
## archive. Filtering on the verdict therefore hid every warned sample in the
## one file whose purpose is to list them.
##
flagged <- background |>
    filter(nFail > 0 | nWarn > 0) |>
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

##
## Import audit: one row per file the scan considered, whether or not it made it
## into the background, plus the duplicate and rescue decisions.
##
fileAudit <- bind_rows(asmCollected$audit, wgsCollected$audit, mqcCollected$audit)

duplicateAudit <- bind_rows(asmResolved$dups, wgsResolved$dups, mqcResolved$dups)
if (nrow(duplicateAudit) > 0) {
    ## A multiqc file holds many samples and can be the superseded copy for one
    ## of them while remaining the kept copy for the rest, so collapse to one
    ## verdict per file before folding the outcome back into the audit.
    dupStatus <- duplicateAudit |>
        summarise(anySuperseded = any(!kept), .by = c(kind, path)) |>
        transmute(kind, path,
                  dupDisposition = if_else(anySuperseded, "supersededDuplicate", NA_character_),
                  dupDetail = "more than one file for this project and sample")
    fileAudit <- fileAudit |>
        left_join(dupStatus, by = c("kind", "path")) |>
        mutate(disposition = coalesce(dupDisposition, disposition),
               detail = if_else(is.na(dupDetail), detail, dupDetail)) |>
        select(-dupDisposition, -dupDetail)
}

write_tsv(background, path(outDir, "backgroundSamples.tsv"))
write_tsv(referenceStats, path(outDir, "backgroundStats.tsv"))
write_tsv(coverageStats, path(outDir, "backgroundCoverageStats.tsv"))
write_tsv(flagged, path(outDir, "backgroundFlagged.tsv"))
write_tsv(metricCoverage, path(outDir, "backgroundMetricCoverage.tsv"))
write_tsv(fileAudit, path(outDir, "backgroundImportAudit.tsv"))

##
## Console summary.
##
## Counted on nFail and nWarn rather than on the verdict, for the reason given
## above the flagged table: the verdict collapses to INCOMPLETE whenever any
## metric is missing, which would report zero warnings across the whole archive.
nRef <- sum(background$referenceSample)
nRefFull <- sum(background$referenceSample & background$referenceTier == "full")
nFailSamples <- sum(background$nFail > 0, na.rm = TRUE)
nWarnSamples <- sum(background$nFail == 0 & background$nWarn > 0, na.rm = TRUE)
nParseFail <- sum(!fileAudit$parsed)
nSuperseded <- sum(fileAudit$disposition == "supersededDuplicate")
nCrossProject <- n_distinct(crossProject$sample)
rule <- strrep("=", 74)

## An archive can be missing a whole class of file, in which case that column
## never gets created. Absent means no sample has it, not an error.
presentIn <- function(dat, column) {
    if (column %in% names(dat)) !is.na(dat[[column]]) else rep(FALSE, nrow(dat))
}
hasAsm <- presentIn(background, "asmPath")
hasWgs <- presentIn(background, "wgsPath")
hasMqc <- presentIn(background, "mqcPath")

cat(rule, "\n")
cat("  BACKGROUND QC IMPORT\n")
cat(rule, "\n\n")
cat(sprintf("  Samples parsed      %5d\n", nrow(background)))
cat(sprintf("  Projects            %5d\n", n_distinct(background$project)))
cat(sprintf("  Reference samples   %5d   clean, used to set the ranges\n", nRef))
cat(sprintf("                      %5d   of those carry both Picard files\n", nRefFull))
cat(sprintf("  Below threshold     %5d\n", nFailSamples))
cat(sprintf("  Warning only        %5d\n\n", nWarnSamples))

cat("  Input files\n")
cat("  ", strrep("-", 70), "\n", sep = "")
parsedCount <- function(which) sum(fileAudit$kind == which & fileAudit$parsed)
cat(sprintf("  %-22s %5d found   %5d parsed\n", "alignment (asm)",
            length(asmFiles), parsedCount("asm")))
cat(sprintf("  %-22s %5d found   %5d parsed\n", "coverage (wgs)",
            length(wgsFiles), parsedCount("wgs")))
cat(sprintf("  %-22s %5d found   %5d parsed\n", "multiqc samtools",
            length(mqcFiles), parsedCount("multiqc")))
cat("\n")
cat(sprintf("  %-22s %5d\n", "samples with both", sum(hasAsm & hasWgs)))
cat(sprintf("  %-22s %5d\n", "alignment only", sum(hasAsm & !hasWgs)))
cat(sprintf("  %-22s %5d\n", "coverage only", sum(!hasAsm & hasWgs)))
cat(sprintf("  %-22s %5d\n", "samtools metrics", sum(hasMqc)))
if (nParseFail > 0)  cat(sprintf("  %-22s %5d   listed in backgroundImportAudit.tsv\n", "unreadable files", nParseFail))
if (nSuperseded > 0) cat(sprintf("  %-22s %5d   older copies, newest kept\n", "superseded files", nSuperseded))
if (nCrossProject > 0) cat(sprintf("  %-22s %5d   kept separate, see warnings\n", "name in 2+ projects", nCrossProject))
cat("\n")

if ("samtoolsStage" %in% names(background)) {
    stages <- background |> filter(!is.na(samtoolsStage)) |> count(samtoolsStage)
    if (nrow(stages) > 1) {
        cat("  samtools metrics were read from a mix of BAM stages:\n")
        stages |>
            mutate(line = sprintf("    %-8s %5d sample(s)", samtoolsStage, n)) |>
            pull(line) |>
            walk(\(x) cat(x, "\n"))
        cat("\n")
    }
}

if (nFailSamples > 0) {
    cat(sprintf("  %d historical samples fall below the filter thresholds and are\n", nFailSamples))
    cat("  excluded from the reference ranges. Review backgroundFlagged.tsv: these\n")
    cat("  are either real defects that were processed anyway, or evidence that a\n")
    cat("  threshold is wrong.\n\n")

    ## Names the metrics that actually failed rather than printing two fixed
    ## columns. Which metrics exist depends on which files the archive holds, and
    ## the fixed version printed a column of NA the moment a cohort arrived
    ## without its alignment metrics.
    cat("  Worst historical samples\n")
    cat("  ", strrep("-", 70), "\n", sep = "")
    flagged |>
        filter(nFail > 0) |>
        head(12) |>
        mutate(line = sprintf("  %-24s %-14s %d failed: %s",
                              str_trunc(sample, 24), str_trunc(project, 14),
                              nFail, str_trunc(failedMetrics, 28))) |>
        pull(line) |>
        walk(\(x) cat(x, "\n"))
    if (nFailSamples > 12) cat(sprintf("  ... and %d more\n", nFailSamples - 12))
    cat("\n")
}

cat("  Reference ranges, clean samples only\n")
cat("  ", strrep("-", 70), "\n", sep = "")
if (nrow(referenceStats) == 0) {
    cat("  None. No sample passed every threshold it could be evaluated on.\n")
} else {
    referenceStats |>
        left_join(metricCoverage |> select(metric, pctAvailable), by = "metric") |>
        mutate(line = sprintf("  %-19s n=%-4d median %9.3f   q05-q95 %8.3f - %8.3f  [%s%% of archive]",
                              metric, n, median, q05, q95, format(pctAvailable, nsmall = 1))) |>
        pull(line) |>
        walk(\(x) cat(x, "\n"))
}

missingMetrics <- metricCoverage |> filter(nAvailable == 0)
if (nrow(missingMetrics) > 0) {
    cat("\n  No background available for: ", str_c(missingMetrics$metric, collapse = ", "), "\n", sep = "")
    cat("  These filters run on fixed values only, with no historical context.\n")
}

if (length(importWarnings) > 0) {
    cat("\n  Warnings\n")
    cat("  ", strrep("-", 70), "\n", sep = "")
    walk(importWarnings, \(w) cat("  ", str_wrap(w, width = 70, exdent = 4), "\n", sep = ""))
}

cat("\n")
cat(sprintf("  Written to %s/\n", outDir))
cat(rule, "\n")
