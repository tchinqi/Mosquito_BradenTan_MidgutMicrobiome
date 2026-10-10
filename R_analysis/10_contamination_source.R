# =============================================================================
# 10_contamination_source.R - where did the reads in the negative controls
#                             come from?
#
# Predictions if the blanks were contaminated by NEIGHBOURING SAMPLES
# (well-to-well transfer / tag misassignment) rather than by reagents:
#   A. Plate specificity: each blank resembles its OWN plate's samples more
#      than the samples on other plates
#   B. OTU origin: nearly all blank reads belong to OTUs that also occur in
#      samples on the same plate; few reads from OTUs found only in controls
#   C. Wolbachia OTUs: the blank carries the Wolbachia OTU(s) typical of the
#      mosquitoes on its plate (Wolbachia cannot come from reagents)
#   D. Distance decay (only if the metadata has well positions): samples
#      physically closer to H12 resemble the blank more (Minich et al., 2019)
#
# Outputs
#   output/stats/10_*.tsv
#   output/stats/10_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))
plate_of <- function(id) str_extract(id, "(AG|ZX|BT)_P0[12]")

meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")
ctrl    <- read_otu_matrix("controls")
bio_rel  <- bio  / rowSums(bio)
ctrl_rel <- ctrl / rowSums(ctrl)

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Lab_Taxon), by = "Sample_ID") |>
  mutate(Plate = plate_of(Sample_ID),
         Host  = case_when(Origin == "Wild" ~ "Wild",
                           Lab_Taxon == "Ae. albopictus" ~ "Lab Ae. albopictus",
                           Lab_Taxon == "Culex" ~ "Lab Culex"))
plates <- sort(unique(sd$Plate))
blank_plate <- setNames(plate_of(rownames(ctrl)), rownames(ctrl))
stopifnot(!anyNA(sd$Plate), !anyNA(blank_plate))

# ---- A. Plate specificity ---------------------------------------------------
# Bray-Curtis between each blank and (i) the mean OTU profile of every plate,
# (ii) every individual sample; own plate vs other plates (Wilcoxon).
bc <- function(a, b) sum(abs(a - b)) / sum(a + b)
plate_means <- t(sapply(plates, function(p) colMeans(bio_rel[sd$Plate == p, , drop = FALSE])))

A_means <- expand_grid(Control = rownames(ctrl), Plate = plates) |>
  mutate(own = Plate == blank_plate[Control],
         bray_curtis = map2_dbl(Control, Plate, ~ bc(ctrl_rel[.x, ], plate_means[.y, ])))
save_table(A_means, "10_A_blank_vs_plate_means", paths$stats)

A_samples <- map_dfr(rownames(ctrl), function(cid) {
  d   <- apply(bio_rel, 1, function(s) bc(ctrl_rel[cid, ], s))
  own <- sd$Plate == blank_plate[cid]
  w   <- wilcox.test(d[own], d[!own], alternative = "less", exact = FALSE)
  closest <- names(sort(d))[1:10]
  tibble(Control = cid, Plate = blank_plate[cid],
         median_bc_own = median(d[own]), median_bc_other = median(d[!own]),
         p_own_more_similar = w$p.value,
         closest10_on_own_plate = sum(sd$Plate[match(closest, sd$Sample_ID)] == blank_plate[cid]))
})
save_table(A_samples, "10_A_blank_vs_individual_samples", paths$stats)

# ---- B. OTU origin ----------------------------------------------------------
in_any_sample <- colSums(bio > 0) > 0
B <- map_dfr(rownames(ctrl), function(cid) {
  r   <- ctrl[cid, ]
  own <- colSums(bio[sd$Plate == blank_plate[cid], , drop = FALSE] > 0) > 0
  tibble(Control = cid, Plate = blank_plate[cid], reads = sum(r), otus = sum(r > 0),
         pct_reads_otus_on_own_plate   = 100 * sum(r[own]) / sum(r),
         pct_reads_otus_in_any_sample  = 100 * sum(r[in_any_sample]) / sum(r),
         pct_reads_otus_only_in_controls = 100 * sum(r[!in_any_sample]) / sum(r),
         otus_only_in_controls = sum(r > 0 & !in_any_sample))
})
save_table(B, "10_B_blank_otu_origin", paths$stats)

# ---- C. Wolbachia OTU types -------------------------------------------------
wol <- otu_map$OTU[otu_map$Genus == "Wolbachia"]
wol_reads_bio <- colSums(bio[, wol, drop = FALSE])
wol_main <- names(sort(wol_reads_bio, decreasing = TRUE))[seq_len(min(5, length(wol)))]
# share of each Wolbachia OTU within all Wolbachia reads, per host group and per blank
wol_share <- function(m) {
  w <- m[, wol_main, drop = FALSE]
  tot <- sum(m[, wol, drop = FALSE])
  if (tot == 0) return(setNames(rep(NA_real_, length(wol_main)), wol_main))
  100 * colSums(w) / tot
}
C <- bind_rows(
  map_dfr(c("Wild", "Lab Ae. albopictus", "Lab Culex"), function(h) {
    as_tibble_row(as.list(wol_share(bio[sd$Host == h, , drop = FALSE]))) |> mutate(Source = h, .before = 1)
  }),
  map_dfr(rownames(ctrl), function(cid) {
    as_tibble_row(as.list(wol_share(ctrl[cid, , drop = FALSE]))) |>
      mutate(Source = paste0("Blank ", blank_plate[cid]), .before = 1)
  })
)
# which host groups are on each plate (context for C)
plate_hosts <- sd |> count(Plate, Host) |>
  group_by(Plate) |> summarise(hosts = paste(sprintf("%s %d", Host, n), collapse = ", "), .groups = "drop")
save_table(C, "10_C_wolbachia_otu_shares", paths$stats)
save_table(plate_hosts, "10_C_hosts_per_plate", paths$stats)

# ---- D. Distance decay (only if well positions exist) ------------------------
well_col <- names(meta)[str_detect(str_to_lower(names(meta)), "^well$|well_?pos|position|^well_id$")]
D_lines <- "  skipped: no well-position column found in the metadata (needs e.g. 'Well' = A01...H12)"
if (length(well_col) >= 1) {
  wc <- well_col[1]
  pos <- meta |> transmute(Sample_ID, well = toupper(as.character(.data[[wc]]))) |>
    mutate(row = match(str_sub(well, 1, 1), LETTERS[1:8]),
           col = as.integer(str_extract(well, "\\d+")))
  D <- map_dfr(rownames(ctrl), function(cid) {
    s   <- sd |> filter(Plate == blank_plate[cid]) |> left_join(pos, by = "Sample_ID") |>
      filter(!is.na(row), !is.na(col))
    if (nrow(s) < 10) return(NULL)
    dist <- sqrt((s$row - 8)^2 + (s$col - 12)^2)          # distance from H12
    sim  <- 1 - sapply(s$Sample_ID, function(x) bc(ctrl_rel[cid, ], bio_rel[x, ]))
    ct <- cor.test(dist, sim, method = "spearman", exact = FALSE)
    tibble(Control = cid, Plate = blank_plate[cid], n = nrow(s),
           rho_distance_vs_similarity = unname(ct$estimate), p = ct$p.value)
  })
  if (nrow(D) > 0) {
    save_table(D, "10_D_distance_decay", paths$stats)
    D_lines <- sprintf("  %-36s rho = %+.2f, %s (n = %d; negative rho = closer wells more similar)",
                       D$Control, D$rho_distance_vs_similarity, fmt_p(D$p), D$n)
  }
}

# ---- Summary ----------------------------------------------------------------
own_rank <- A_means |> group_by(Control) |>
  summarise(own_bc = bray_curtis[own], best_other = min(bray_curtis[!own]),
            own_is_most_similar = own_bc < best_other, .groups = "drop")

summary_lines <- c(
  "10_contamination_source.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "A. PLATE SPECIFICITY (Bray-Curtis on OTU profiles; lower = more similar)",
  sprintf("  %-36s own plate %.2f vs closest other plate %.2f -> %s",
          own_rank$Control, own_rank$own_bc, own_rank$best_other,
          if_else(own_rank$own_is_most_similar, "OWN PLATE MOST SIMILAR", "another plate more similar")),
  "  Individual samples (one-sided Wilcoxon: own-plate samples more similar?):",
  sprintf("  %-36s median BC own %.2f vs other %.2f, %s | %d of 10 most similar samples on own plate",
          A_samples$Control, A_samples$median_bc_own, A_samples$median_bc_other,
          fmt_p(A_samples$p_own_more_similar), A_samples$closest10_on_own_plate),
  "",
  "B. OTU ORIGIN (% of each blank's reads)",
  sprintf("  %-36s %5d reads | %5.1f%% from OTUs on own plate | %5.1f%% from OTUs in any sample | %4.1f%% from control-only OTUs (%d OTUs)",
          B$Control, as.integer(B$reads), B$pct_reads_otus_on_own_plate, B$pct_reads_otus_in_any_sample,
          B$pct_reads_otus_only_in_controls, as.integer(B$otus_only_in_controls)),
  "",
  "C. WOLBACHIA OTU TYPES (% of Wolbachia reads in each OTU)",
  paste0("  ", sprintf("%-24s", "Source"), paste(sprintf("%12s", wol_main), collapse = "")),
  sprintf("  %-24s%s", C$Source,
          apply(as.matrix(C[, wol_main]), 1, function(v) paste(sprintf("%12s", ifelse(is.na(v), "-", sprintf("%.0f%%", v))), collapse = ""))),
  "  Hosts per plate:",
  sprintf("    %s: %s", plate_hosts$Plate, plate_hosts$hosts),
  "",
  "D. DISTANCE DECAY FROM H12",
  D_lines
)
writeLines(summary_lines, file.path(paths$stats, "10_results_summary.txt"))
cat(summary_lines, sep = "\n")
