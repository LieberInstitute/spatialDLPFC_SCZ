################# NOTE #########################
# Code adapted from:                           #
# https://github.com/LieberInstitute/ranger_metrics/tree/master/code/spaceranger
##################33333#########################



library("readxl")
library("here")
library("sessioninfo")

## Rscript -e 'sgejobs::job_single("01_read_summaries_spaceranger", create_shell = TRUE, queue = "bluejay")'

## Output dirs
dir_rdata <-
    here(
        "processed-data",
        "spaceranger",
        "01_read_experimental_summaries"
    )
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Locate summaries
## NOTE (boyiguo1): Customize for SCZ project
summaries <-
  list.files(
    here("raw-data", "experiment_info"),
    pattern = "*.xlsx",
    full.names = TRUE,
    recursive = TRUE
  )


## Figure out number of rows
summ_rows <- lapply(summaries, function(x) {
    print(x)
    info <- read_excel(x, sheet = 1)
    sample_col <- grep("sample \\#", tolower(colnames(info)))
    stopifnot(length(sample_col) == 1)
    sum(!is.na(info[, sample_col])) - sum(tolower(info[, sample_col]) == "sample")
})

## Read tables
summ_data <- mapply(function(x, max_row) {
    print(x)
    info <- read_excel(x, sheet = 1, n_max = max_row)
    ## We don't need the 2 dilution columns (cDNA and library)
    ## We'll drop the empty title columns
    info <-
        info[, !grepl("dilution|^\\.+", tolower(colnames(info)))]

    ## Estimated Read Pairs are manually computed differently for each study

    ## Add which sheet file the data comes from
    info[["sheet_file"]] <- x

    ## Fix the cDNA Input ng
    colnames(info)[colnames(info) == "cDNA Input"] <-
        "cDNA Input ng"

    ## Fix "Est. Read Pairs", "Est Read Pairs (million)", "Est Read Pairs (50000 per spot)"
    colnames(info)[colnames(info) %in% c(
        "Est. Read Pairs",
        "Est Read Pairs (million)",
        "Est Read Pairs (50000 per spot)"
    )] <- "Est Read Pairs"

    ## Fix "Ave Frag Size [bp]" and "Ave frag length library"
    colnames(info)[colnames(info) %in% c("Ave Frag Size [bp]", "Ave frag length library")] <-
        "Ave frag length"

    ## Fix "Slide SN #"
    colnames(info)[colnames(info) == "Slide SN #"] <- "Slide #"

    ## Fix "Br_Region"
    colnames(info)[colnames(info) == "Br_Region"] <- "Tissue"

    ## Fix "BrNumbr" and "Br####"
    colnames(info)[colnames(info) %in% c("BrNumbr", "Br####")] <-
        "Brain"

    ## We'll keep the second Agilent [pg / ul] column
    agilent_cols <-
        grep("agilent|final conc", tolower(colnames(info)))
    if (length(agilent_cols) == 2) {
        agi_check <- info[, agilent_cols[2]] >= info[, agilent_cols[1]]
        stopifnot(all(agi_check[!is.na(agi_check)]))
        colnames(info)[agilent_cols[2]] <- "Agilent [pg/ul]"
        info <- info[, -agilent_cols[1]]
    } else if (length(agilent_cols) != 1) {
        stop("Unexpected number of alignent columns", call. = FALSE)
    }
    return(info)
}, summaries, summ_rows, SIMPLIFY = FALSE)

## Get column names, aim to keep the same ones
summ_cols <- lapply(summ_data, colnames)
# unique(unlist(summ_cols))
sort(table(unlist(summ_cols)))


## Function for manually checking things
which_file <- function(x) {
    res <- sapply(summ_cols, function(y) {
        any(x %in% y)
    })
    which(res)
}
# which_file("Final Conc. [pg/ul]")
# which_file("KAPA (pM)")

## Keep only columns present in all sheets
keep_col <- table(unlist(summ_cols))
keep_col <- names(keep_col[keep_col == length(summaries)])

## Merge all tables
spaceranger_summary <-
    do.call(rbind, lapply(summ_data, function(x) {
        x[, keep_col]
    }))

## Some rows are weird
n_NAs <- apply(spaceranger_summary, 1, function(x) {
    sum(is.na(x))
})
names(n_NAs) <- seq_len(length(n_NAs))
## We manually explored things and just those with 17 or 18 NAs are to be
## discarded
spaceranger_summary <- spaceranger_summary[n_NAs < 17, ]


## Fix "Ct" with a "failed" value
spaceranger_summary$Ct[spaceranger_summary$Ct == "failed"] <- 0

## Export for later use
write.csv(
    spaceranger_summary,
    file.path(dir_rdata, "presequencing_summary.csv"),
    row.names = FALSE
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
