
# Run this after 19,21 and 29 and after compiling S1_Appendix.tex.
# Output : Derived/supporting_information/

source("/.../.../.../Code/R/00_config.R")

suppressPackageStartupMessages({ library(grid); library(png) })

SI_DIR <- file.path(OUT_DIR, "supporting_information")
dir.create(SI_DIR, showWarnings = FALSE, recursive = TRUE)

# where the compiled appendix lives; change if you compile it elsewhere
SI_APPENDIX_PDF <- file.path(ROOT, "manuscript", "S1_Appendix.pdf")

# ---- figures: label, title and caption around the image -----------------------
si_figures <- list(
  list(id = "S1 Fig",
       png = file.path(FIG_DIR, "fig_cf_campaign.png"),
       title = "Infections averted per 1,000 population by the measles-rubella campaign, by district.",
       text = paste("The twenty Uttar Pradesh districts with the largest effect, April 2013",
                    "to March 2020. Points are medians and bars 95% intervals across 1,000",
                    "bootstrap replicates.")),
  list(id = "S2 Fig",
       png = file.path(FIG_DIR, "fig_nfhs_gain_vs_population.png"),
       title = "Effect of the static covariates on forecast error against district population.",
       text = "")
)

# A4 portrait, 2 cm margins; the image is scaled to the space left by the text
draw_si_figure <- function(f, out) {
  W <- 8.27; H <- 11.69; m <- 0.79
  img <- png::readPNG(f$png)
  asp <- dim(img)[1] / dim(img)[2]                       # height / width

  grDevices::cairo_pdf(out, width = W, height = H)   # open first, so text can be measured
  on.exit(dev.off())

  head_txt <- paste0(f$id, ". ", f$title)
  wrap <- function(s, n) paste(strwrap(s, n), collapse = "\n")
  head_g <- textGrob(wrap(head_txt, 85), x = 0, just = c("left", "top"),
                     gp = gpar(fontsize = 11, fontface = "bold"))
  cap_g  <- textGrob(wrap(f$text, 95), x = 0, just = c("left", "top"),
                     gp = gpar(fontsize = 10))
  head_h <- convertHeight(grobHeight(head_g), "inches", valueOnly = TRUE)
  cap_h  <- if (nzchar(f$text)) convertHeight(grobHeight(cap_g), "inches", valueOnly = TRUE) else 0

  gap   <- 0.2
  img_w <- W - 2 * m
  img_h <- min(img_w * asp, H - 2 * m - head_h - cap_h - 3 * gap)
  img_w <- img_h / asp

  grid.newpage()
  y <- H - m
  pushViewport(viewport(x = m / W, y = y / H, width = (W - 2 * m) / W,
                        height = head_h / H, just = c("left", "top")))
  grid.draw(head_g); popViewport()
  y <- y - head_h - gap
  grid.raster(img, x = m / W, y = y / H, width = img_w / W, height = img_h / H,
              just = c("left", "top"))
  y <- y - img_h - gap
  if (nzchar(f$text)) {
    pushViewport(viewport(x = m / W, y = y / H, width = (W - 2 * m) / W,
                          height = cap_h / H, just = c("left", "top")))
    grid.draw(cap_g); popViewport()
  }
  invisible(out)
}

made <- character(0); missing <- character(0)

for (f in si_figures) {
  out <- file.path(SI_DIR, paste0(gsub(" ", "_", f$id), ".pdf"))
  if (file.exists(f$png)) { draw_si_figure(f, out); made <- c(made, basename(out)) }
  else missing <- c(missing, paste(f$id, "<-", f$png))
}

# ---- tables: already labelled and captioned by script 29 ------------------------
for (i in 1:6) {
  src <- file.path(TAB_DIR, "supplement", sprintf("S%d_Table.pdf", i))
  if (file.exists(src)) {
    file.copy(src, file.path(SI_DIR, basename(src)), overwrite = TRUE)
    made <- c(made, basename(src))
  } else missing <- c(missing, paste0("S", i, " Table <- ", src))
}

# ---- appendix ---------------------------------------------------------------------
if (file.exists(SI_APPENDIX_PDF)) {
  file.copy(SI_APPENDIX_PDF, file.path(SI_DIR, "S1_Appendix.pdf"), overwrite = TRUE)
  made <- c(made, "S1_Appendix.pdf")
} else missing <- c(missing, paste("S1 Appendix <-", SI_APPENDIX_PDF))

cat("\nWritten to", SI_DIR, ":\n"); cat(paste0("  ", sort(made), collapse = "\n"), "\n")
if (length(missing)) {
  cat("\nNot found (run the script that makes it, or fix the path):\n")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
}
