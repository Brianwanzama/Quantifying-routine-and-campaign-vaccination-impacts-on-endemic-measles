
# 01_district_names.R
# standardise_districts(), called wherever a file is read.
# Names are resolved within a state, so two states sharing a district name do
# not collide.


# ---- the reference set, per state --------------------------------------------
load_reference_districts <- function(path) {
  h <- read.csv(path, nrows = 1, check.names = FALSE, stringsAsFactors = FALSE)
  ref <- setdiff(names(h), c("District", "district", "", "X"))
  ref <- trimws(ref)
  ref[nzchar(ref)]
}

load_all_references <- function(states, file_of_fn) {
  out <- lapply(states, function(s) load_reference_districts(file_of_fn(s, "coords")))
  names(out) <- states
  out
}

# Left is any name seen in a file, normalised; right is the reference name.
DISTRICT_ALIASES <- c(
  # Uttar Pradesh, renamed districts
  "prayagraj"               = "Allahabad",
  "ayodhya"                 = "Faizabad",
  "faizabadayodhya"         = "Faizabad",
  "jyotibaphulenagar"       = "Amroha",
  "mahamayanagar"           = "Hathras",
  "kanshiramnagar"          = "Kasganj",
  "kashiramnagar"           = "Kasganj",
  "csmnagar"                = "Amethi",   # Chhatrapati Shahuji Maharaj Nagar
  "bhadohi"                 = "Sant Ravi Das Nagar",
  "badohi"                  = "Sant Ravi Das Nagar",
  "santravidasnagar"        = "Sant Ravi Das Nagar",
  "santravidasnagarbhadohi" = "Sant Ravi Das Nagar",

  # abbreviations
  "kheri"                   = "Lakhimpur Kheri",
  "lakhimpurkheri"          = "Lakhimpur Kheri",
  "mahrajganj"              = "Maharajganj",
  "maharajgunj"             = "Maharajganj",
  "maunathbhanjan"          = "Mau",

  # spelling
  "badaun"                  = "Budaun",
  "bagpat"                  = "Baghpat",
  "bulandshahar"            = "Bulandshahr",
  "ferozabad"               = "Firozabad",
  "srawasti"                = "Shravasti",
  "shrawasti"               = "Shravasti",
  "khushinagar"             = "Kushinagar",
  "raebareli"               = "Rae Bareli",
  "barbanki"                = "Barabanki",
  "siddharthanagar"         = "Siddharth Nagar",
  "gautambudhnagar"         = "Gautam Buddha Nagar",
  "gautambudhanagar"        = "Gautam Buddha Nagar",
  "sonebhadra"              = "Sonbhadra",
  "chitrakut"               = "Chitrakoot",
  "unnav"                   = "Unnao",

  # West Bengal. The reference uses the older transliterations, so most names
  # resolve by normalisation alone. Only the modern spellings, where a file
  # uses one, are mapped back to the reference form. Mapping the other way
  # would send a name that already matches to one that does not exist.
  "darjeeling"              = "Darjiling",
  "howrah"                  = "Haora",
  "hooghly"                 = "Hugli",
  "coochbehar"              = "Koch Bihar",
  "malda"                   = "Maldah",
  "purulia"                 = "Puruliya",
  "north24parganas"         = "North Twenty Four Parganas",
  "south24parganas"         = "South Twenty Four Parganas",
  "north24pgs"              = "North Twenty Four Parganas",
  "south24pgs"              = "South Twenty Four Parganas",
  "medinipureast"           = "Purba Medinipur",
  "medinipurwest"           = "Paschim Medinipur",
  # the shapefile spells this with an h that the reference does not
  "pashchimmedinipur"       = "Paschim Medinipur"
)

normalise_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
}

# ---- the function called at every read ---------------------------------------
# reference is that state's names. Returns NA where nothing matched.
standardise_districts <- function(x, reference, source = "", strict = FALSE) {
  raw  <- trimws(as.character(x))
  norm <- normalise_name(raw)

  ref_norm <- normalise_name(reference)
  ali_norm <- normalise_name(names(DISTRICT_ALIASES))

  out <- rep(NA_character_, length(raw))

  hit <- match(raw, reference)                       # exact
  out[!is.na(hit)] <- reference[hit[!is.na(hit)]]

  todo <- is.na(out)                                 # alias list
  if (any(todo)) {
    hit <- match(norm[todo], ali_norm)
    got <- !is.na(hit)
    out[todo][got] <- unname(DISTRICT_ALIASES[hit[got]])
  }

  todo <- is.na(out)                                 # normalised
  if (any(todo)) {
    hit <- match(norm[todo], ref_norm)
    got <- !is.na(hit)
    out[todo][got] <- reference[hit[got]]
  }

  # An alias may resolve to a name that is not in this state's reference, which
  # happens when a name belongs to the other state. Treated as unmatched.
  bad <- !is.na(out) & !(out %in% reference)
  if (any(bad)) out[bad] <- NA_character_

  miss <- unique(raw[is.na(out)])
  if (length(miss)) {
    msg <- paste0("unmatched district names",
                  if (nzchar(source)) paste0(" in ", source) else "", ":\n  ",
                  paste(miss, collapse = "\n  "))
    if (strict) stop(msg) else warning(msg, call. = FALSE)
  }
  out
}

# ---- the audit, run once before anything else --------------------------------
audit_district_names <- function(named, references) {
  cov <- list()
  for (nm in names(named)) {
    state <- named[[nm]]$state
    x     <- named[[nm]]$names
    ref   <- references[[state]]
    std   <- suppressWarnings(standardise_districts(x, ref, source = nm))

    miss   <- unique(x[is.na(std)])
    got    <- sort(unique(std[!is.na(std)]))
    absent <- setdiff(ref, got)

    cat("=== ", nm, " ===\n", sep = "")
    cat("  in file ", length(unique(x)), " | matched ", length(got),
        " | unmatched ", length(miss), "\n", sep = "")
    if (length(miss)) { cat("  unmatched:\n"); for (m in miss) cat("    ", m, "\n") }
    if (length(absent)) {
      cat("  absent from this file (", length(absent), "):\n", sep = "")
      for (a in absent) cat("    ", a, "\n")
    }
    cat("\n")
    cov[[nm]] <- paste(state, got, sep = "__")
  }

  common <- Reduce(intersect, cov)
  all_ref <- unlist(lapply(names(references), function(s)
    paste(s, references[[s]], sep = "__")))
  cat("districts present in every file: ", length(common), " of ",
      length(all_ref), "\n", sep = "")
  if (length(common) < length(all_ref)) {
    cat("missing from at least one file:\n")
    for (d in setdiff(all_ref, common)) cat("  ", d, "\n")
  }
  invisible(list(coverage = cov, common = common))
}
