library(tidyverse)

runs <- c(
  "./data/raw/red_drum_11_75"
)

summ <- map_dfr(runs, function(prefix){
  panel_file <- paste0(prefix, "panel.txt")
  primer_file <- paste0(prefix, "primers.txt")
  pos_file <- paste0(prefix, "pos.txt")

  panel <- read_tsv(panel_file, show_col_types = FALSE,
                    col_types = cols(Pos = col_character()))
  primers <- read_tsv(primer_file, show_col_types = FALSE)
  pos <- read_tsv(pos_file, show_col_types = FALSE)

  amplicons <- primers %>%
    group_by(Locus) %>%
    summarize(
      left_start = start[orient == "left"],
      right_start = start[orient == "right"],
      right_len = len[orient == "right"],
      amplicon_length = right_start + right_len - left_start,
      .groups = "drop"
    )

  tibble(
    run = prefix,
    n_loci = nrow(panel),
    mean_He = mean(panel$summaryHe, na.rm = TRUE),
    min_He = min(panel$summaryHe, na.rm = TRUE),
    median_He = median(panel$summaryHe, na.rm = TRUE),
    max_He = max(panel$summaryHe, na.rm = TRUE),
    mean_snps_per_locus = mean(str_count(panel$Pos, ",") + 1),
    median_snps_per_locus = median(str_count(panel$Pos, ",") + 1),
    n_primers = nrow(primers),
    n_primer_pairs = nrow(primers) / 2,
    mean_primer_tm = mean(primers$tm, na.rm = TRUE),
    mean_amplicon_snp_count = nrow(pos) / nrow(panel),
    mean_amplicon_length = mean(amplicons$amplicon_length),
    median_amplicon_length = median(amplicons$amplicon_length),
    min_amplicon_length = min(amplicons$amplicon_length),
    max_amplicon_length = max(amplicons$amplicon_length)
  )
})

write_tsv(summ, "./data/derived/panel_run_summary.tsv")
print(summ)
