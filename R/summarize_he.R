library(tidyverse)

files <- c(
  "./data/raw/He_50_WGULF.txt",
  "./data/raw/He_50_EGULF.txt"
)

summ <- map_dfr(files, function(f){
  x <- read_tsv(f, show_col_types = FALSE)
  tibble(
    file = f,
    n_windows = nrow(x),
    mean_He = mean(x$He, na.rm = TRUE),
    median_He = median(x$He, na.rm = TRUE),
    pct_He_gt_0.5 = mean(x$He > 0.5, na.rm = TRUE) * 100,
    pct_He_gt_0.7 = mean(x$He > 0.7, na.rm = TRUE) * 100,
    mean_numInds = mean(x$numInds, na.rm = TRUE),
    median_numInds = median(x$numInds, na.rm = TRUE),
    mean_failed_geno = mean(x$numFailGeno, na.rm = TRUE),
    mean_invalid_phase = mean(x$numInvalidPhase, na.rm = TRUE)
  )
})

write_tsv(summ, "./data/derived/he_file_summary.tsv")
print(summ)
