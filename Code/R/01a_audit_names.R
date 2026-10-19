
# 01a_audit_names.R
# Run once, before anything else. Reads every file, standardises its district
# names against that state's coordinates file, and reports what failed.
# Nothing is written.

source("/.../.../Code/R/00_config.R")
source("/.../.../Code/R/01_district_names.R")

refs <- load_all_references(STATES, file_of)
for (s in STATES)
  cat(s, ": ", length(refs[[s]]), " districts in the coordinates file\n", sep = "")
cat("\n")

# District names as they appear in each file. The wide files carry them as
# column headers after the leading time columns; the campaign file in a column.
hdr <- function(path, drop) {
  h <- read.csv(path, nrows = 1, check.names = FALSE, stringsAsFactors = FALSE)
  x <- trimws(setdiff(names(h), drop))
  x[nzchar(x) & !grepl("^X\\.?[0-9]*$", x)]
}

named <- list()
for (s in STATES) {
  add <- function(tag, path, drop) {
    if (is.na(path) || !file.exists(path)) return(invisible(NULL))
    named[[paste(s, tag)]] <<- list(state = s, names = hdr(path, drop))
  }
  add("cases",  file_of(s, "cases"),  c("Year", "year", "Month", "month"))
  add("births", file_of(s, "births"), c("Year", "year", "Month", "month"))
  add("popn",   file_of(s, "popn"),   c("Year", "year"))
  add("mcv1",   file_of(s, "mcv1"),   c("Year", "year"))
  add("mcv2",   file_of(s, "mcv2"),   c("Year", "year"))

  fc <- file_of(s, "campaign")
  if (!is.na(fc) && file.exists(fc)) {
    d <- read.csv(fc, check.names = FALSE, stringsAsFactors = FALSE)
    dc <- intersect(c("District", "district"), names(d))[1]
    if (!is.na(dc))
      named[[paste(s, "campaign")]] <- list(state = s, names = trimws(d[[dc]]))
  } else {
    cat(s, "has no campaign file\n\n")
  }
}

res <- audit_district_names(named, refs)

# ---- alias lines for anything still unmatched --------------------------------
unmatched <- unique(unlist(lapply(names(named), function(nm) {
  x   <- named[[nm]]$names
  ref <- refs[[named[[nm]]$state]]
  std <- suppressWarnings(standardise_districts(x, ref, source = nm))
  x[is.na(std)]
})))

if (length(unmatched)) {
  cat("\n---- add to DISTRICT_ALIASES in 01_district_names.R ----\n")
  for (u in unmatched)
    cat(sprintf('  "%s" = "",   # seen as: %s\n', normalise_name(u), u))
  cat("\nfill in the reference name on the right of each line\n")
} else {
  cat("\nevery district name matched\n")
}

