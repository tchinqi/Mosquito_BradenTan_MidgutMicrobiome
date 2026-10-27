# =============================================================================
# 12_discussion_checks.R - three targeted checks of claims in the Discussion
#
#   1. Wolbachia (4.2): is it genuinely lower in Lab mosquitoes, or only diluted
#      by Elizabethkingia? Wolbachia is recalculated as a share of the reads
#      that are NOT Elizabethkingia, and as a log-ratio to Stenotrophomonas A
#      (a shared-core genus present in almost every mosquito).
#   2. Depth and richness (3.4 / 4.4): is Kranji Marshes richer only because
#      it was sequenced more deeply? (a) richness adjusted for log10 depth
#      (residuals + quasi-Poisson model), (b) rarefied to 1,000 reads.
#      The same check is run for Wild vs Lab.
#   3. Source of negative-control reads (4.3): well-to-well transfer vs tag
#      misassignment. Tag layout (verified against the demultiplexing sheet):
#        - wells numbered down columns: idx = (col - 1) * 8 + row (A = 1)
#        - reverse tag = ((idx - 1) %% 30) + 1, same positions on every plate
#        - forward tag block = ceiling(idx / 30); 4 forward tags per plate,
#          never reused between plates
#      So every H12 control carries reverse tag 6, shared with wells F1, D5 and
#      B9 on ALL plates (and with the other five controls), and carries the
#      4th forward tag of its plate, shared only with C12-G12 on its own plate.
#      Predictions
#        well-to-well transfer -> physical neighbours of H12 resemble the
#                                 control most; similarity decays with distance
#        tag misassignment     -> samples sharing a tag with the control
#                                 resemble it more than other samples, including
#                                 samples on OTHER plates
#      Tests (Bray-Curtis on OTU relative abundance, as in R script 10)
#        3a distance decay on the control's own plate (Spearman), excluding
#           tag-sharing wells
#        3b tag-sharing samples vs plate-matched random samples (permutation)
#        3c fingerprint OTUs: control OTUs found in <= 5 biological samples,
#           and where those samples sit
#
# Needs: data/raw/metadata/demux_metadata_audit.csv (or data/raw/species/)
# Outputs
#   output/stats/12_*.tsv, output/stats/12_results_summary.txt
#   output/figures/FigS_wolbachia_excl_elizabethkingia.pdf/.png
#   output/figures/FigS_control_source.pdf/.png
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))
out <- character()
say <- function(...) { line <- paste0(...); out <<- c(out, line); message(line) }

meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")
ctrl    <- read_otu_matrix("controls")
depth   <- rowSums(bio)

sd <- tibble(Sample_ID = rownames(bio), depth = depth) |>
  left_join(select(meta, Sample_ID, Origin, Lab_Taxon, Location), by = "Sample_ID") |>
  mutate(Group = case_when(Origin == "Wild" ~ "Wild",
                           Lab_Taxon == "Ae. albopictus" ~ "Lab Ae. albopictus",
                           Lab_Taxon == "Culex" ~ "Lab Culex"),
         Group = factor(Group, levels = c("Wild", "Lab Ae. albopictus", "Lab Culex")),
         observed = rowSums(bio > 0))
stopifnot(!anyNA(sd$Group))

say("12_discussion_checks.R - results summary")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M"))

# =============================================================================
# 1. WOLBACHIA WITHOUT ELIZABETHKINGIA
# =============================================================================
genus_reads <- function(pattern) {
  cols <- otu_map$OTU[str_detect(otu_map$Genus, regex(pattern, ignore_case = TRUE))]
  stopifnot(length(cols) > 0)
  rowSums(bio[, cols, drop = FALSE])
}
wol   <- genus_reads("^Wolbachia")
eliz  <- genus_reads("^Elizabethkingia")
steno <- genus_reads("^Stenotrophomonas_A")

W <- sd |>
  mutate(wol_pct       = 100 * wol / depth,
         non_eliz      = depth - eliz,
         wol_pct_excl  = if_else(non_eliz > 0, 100 * wol / non_eliz, NA_real_),
         log_wol_steno = log10((wol + 1) / (steno + 1)),
         is_wild       = factor(Group == "Wild", levels = c(TRUE, FALSE), labels = c("Wild", "Lab")))
save_table(select(W, Sample_ID, Group, depth, wol_pct, wol_pct_excl, log_wol_steno),
           "12_1_wolbachia_per_sample", paths$stats)

say("")
say("1. WOLBACHIA: is the lower Lab value only dilution by Elizabethkingia?")
W_sum <- W |> group_by(Group) |>
  summarise(n = n(), n_excl = sum(!is.na(wol_pct_excl)),
            median_raw = median(wol_pct), median_excl = median(wol_pct_excl, na.rm = TRUE),
            mean_excl = mean(wol_pct_excl, na.rm = TRUE),
            median_log_ratio = median(log_wol_steno), .groups = "drop")
save_table(W_sum, "12_1_wolbachia_summary", paths$stats)
for (i in seq_len(nrow(W_sum))) with(W_sum[i, ], say(sprintf(
  "  %-20s n = %3d | Wolbachia %% of all reads median %5.1f | %% of non-Elizabethkingia reads median %5.1f (mean %5.1f, n = %d) | log10(Wolbachia/Stenotrophomonas A) median %+.2f",
  as.character(Group), n, median_raw, median_excl, mean_excl, n_excl, median_log_ratio)))

kw1 <- kruskal.test(wol_pct_excl ~ Group, data = W)
kw2 <- kruskal.test(log_wol_steno ~ Group, data = W)
pw1 <- pairwise.wilcox.test(W$wol_pct_excl, W$Group, p.adjust.method = "BH", exact = FALSE)
pw2 <- pairwise.wilcox.test(W$log_wol_steno, W$Group, p.adjust.method = "BH", exact = FALSE)
wl1 <- wilcox.test(wol_pct_excl ~ is_wild, data = W, exact = FALSE)
say(sprintf("  Excluding Elizabethkingia: Kruskal-Wallis H = %.2f, %s | Wild vs all Lab Wilcoxon %s",
            kw1$statistic, fmt_p(kw1$p.value), fmt_p(wl1$p.value)))
say("    pairwise (BH): Wild vs Lab Ae. alb ", fmt_p(pw1$p.value["Lab Ae. albopictus", "Wild"]),
    "; Wild vs Lab Culex ", fmt_p(pw1$p.value["Lab Culex", "Wild"]),
    "; Lab Ae. alb vs Lab Culex ", fmt_p(pw1$p.value["Lab Culex", "Lab Ae. albopictus"]))
say(sprintf("  Log-ratio Wolbachia : Stenotrophomonas A: Kruskal-Wallis H = %.2f, %s",
            kw2$statistic, fmt_p(kw2$p.value)))
say("    pairwise (BH): Wild vs Lab Ae. alb ", fmt_p(pw2$p.value["Lab Ae. albopictus", "Wild"]),
    "; Wild vs Lab Culex ", fmt_p(pw2$p.value["Lab Culex", "Wild"]),
    "; Lab Ae. alb vs Lab Culex ", fmt_p(pw2$p.value["Lab Culex", "Lab Ae. albopictus"]))
say("  Reading: if Lab is still lower after removing Elizabethkingia, Wolbachia is")
say("  genuinely reduced in Lab mosquitoes, not just diluted.")

grp_lab <- c("Wild" = "Wild", "Lab Ae. albopictus" = "Lab\nAe. albopictus", "Lab Culex" = "Lab\nCulex")
W_long <- W |> select(Group, wol_pct, wol_pct_excl) |>
  pivot_longer(-Group, names_to = "measure", values_to = "pct") |>
  filter(!is.na(pct)) |>
  mutate(measure = factor(measure, levels = c("wol_pct", "wol_pct_excl"),
                          labels = c("% of all classified reads", "% of non-Elizabethkingia reads")))
p_wol <- ggplot(W_long, aes(Group, pct)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "#f3f2ee") +
  geom_jitter(width = 0.18, height = 0, size = 0.6, alpha = 0.35, colour = "#eb6834") +
  facet_wrap(~ measure, labeller = as_labeller(c(
    "% of all classified reads"      = "A   % of all classified reads",
    "% of non-Elizabethkingia reads" = "B   % of non-Elizabethkingia reads"))) +
  scale_x_discrete(labels = grp_lab) +
  labs(x = NULL, y = expression(italic(Wolbachia)~"(%)"),
       subtitle = sprintf("Excluding Elizabethkingia: Kruskal-Wallis %s", fmt_p(kw1$p.value))) +
  theme_fyp() + theme(plot.subtitle = element_text(size = 8))

# Panel C: log10 ratio of Wolbachia to Stenotrophomonas A (pseudocount of 1 read)
p_ratio <- W |>
  mutate(panel = "C   Ratio to Stenotrophomonas A") |>
  ggplot(aes(Group, log_wol_steno)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "#f3f2ee") +
  geom_jitter(width = 0.18, height = 0, size = 0.6, alpha = 0.35, colour = "#eb6834") +
  facet_wrap(~ panel) +
  scale_x_discrete(labels = grp_lab) +
  labs(x = NULL,
       y = expression(log[10]~"("*italic(Wolbachia)*" / "*italic(Stenotrophomonas)~"A)"),
       subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw2$p.value))) +
  theme_fyp() + theme(plot.subtitle = element_text(size = 8))

fig_wol <- (p_wol | p_ratio) + plot_layout(widths = c(2, 1))
save_fig(fig_wol, "FigS_wolbachia_excl_elizabethkingia", width = 240, height = 95)

# =============================================================================
# 2. RICHNESS ADJUSTED FOR SEQUENCING DEPTH
# =============================================================================
say("")
say("2. DOES SEQUENCING DEPTH EXPLAIN THE RICHNESS DIFFERENCES?")
shannon <- diversity(bio, "shannon")
sd <- sd |> mutate(shannon = shannon, log_depth = log10(depth))

depth_checks <- function(d, group, label) {
  d <- d |> mutate(G = droplevels(factor(.data[[group]])))
  # (a) residuals of richness on log10 depth, then Kruskal-Wallis / Wilcoxon
  d$res_obs <- resid(lm(observed ~ log_depth, data = d))
  d$res_sha <- resid(lm(shannon  ~ log_depth, data = d))
  k_obs <- kruskal.test(res_obs ~ G, data = d)
  k_sha <- kruskal.test(res_sha ~ G, data = d)
  # (b) quasi-Poisson GLM: observed ~ log10 depth + group
  m1 <- glm(observed ~ log_depth + G, family = quasipoisson, data = d)
  m0 <- glm(observed ~ log_depth,     family = quasipoisson, data = d)
  f  <- anova(m0, m1, test = "F")
  # (c) rarefy to 1,000 reads (samples with >= 1,000 reads only)
  keep <- d$depth >= 1000
  set.seed(4288)
  rr   <- rrarefy(bio[d$Sample_ID[keep], , drop = FALSE], 1000)
  dr   <- d[keep, ] |> mutate(obs_r = rowSums(rr > 0), sha_r = diversity(rr, "shannon"))
  k_r1 <- kruskal.test(obs_r ~ G, data = dr)
  k_r2 <- kruskal.test(sha_r ~ G, data = dr)
  meds <- d |> group_by(G) |> summarise(n = n(), depth = median(depth), observed = median(observed),
                                         res_obs = median(res_obs), .groups = "drop") |>
    left_join(dr |> group_by(G) |> summarise(n_rare = n(), obs_rare = median(obs_r),
                                            sha_rare = median(sha_r), .groups = "drop"), by = "G")
  say("  ", label)
  for (i in seq_len(nrow(meds))) with(meds[i, ], say(sprintf(
    "    %-16s n = %3d | median depth %5.0f | observed %4.1f | depth-adjusted residual %+5.2f | rarefied 1,000 reads: observed %4.1f, Shannon %.2f (n = %d)",
    as.character(G), n, depth, observed, res_obs, obs_rare, sha_rare, n_rare)))
  say(sprintf("    Depth-adjusted (residuals): observed K-W H = %.2f, %s; Shannon H = %.2f, %s",
              k_obs$statistic, fmt_p(k_obs$p.value), k_sha$statistic, fmt_p(k_sha$p.value)))
  say(sprintf("    Quasi-Poisson observed ~ log10 depth + group: group F = %.2f, %s",
              f$F[2], fmt_p(f$`Pr(>F)`[2])))
  say(sprintf("    Rarefied to 1,000 reads: observed K-W H = %.2f, %s; Shannon H = %.2f, %s",
              k_r1$statistic, fmt_p(k_r1$p.value), k_r2$statistic, fmt_p(k_r2$p.value)))
  cf <- coef(summary(m1)); cf <- cf[str_detect(rownames(cf), "^G"), , drop = FALSE]
  for (r in rownames(cf)) say(sprintf("      GLM %s vs %s: richness ratio %.2f, %s",
                                       sub("^G", "", r), levels(d$G)[1], exp(cf[r, 1]), fmt_p(cf[r, 4])))
  if (nlevels(d$G) > 2) {
    pw <- pairwise.wilcox.test(d$res_obs, d$G, p.adjust.method = "BH", exact = FALSE)$p.value
    say("    Pairwise depth-adjusted richness (BH): ",
        paste(apply(which(!is.na(pw), arr.ind = TRUE), 1, function(ix)
          sprintf("%s vs %s %s", rownames(pw)[ix[1]], colnames(pw)[ix[2]], fmt_p(pw[ix[1], ix[2]]))),
          collapse = "; "))
  }
  meds
}
site_df <- sd |> filter(!is.na(Location)) |>
  mutate(Location = factor(Location, levels = c("Kranji Marshes", "Mandai Zoo", "Bird Paradise", "River Safari")))
save_table(depth_checks(site_df, "Location", "Wild collection sites (n = 345; reference = Kranji Marshes)"),
           "12_2_depth_adjusted_sites", paths$stats)
save_table(depth_checks(sd |> mutate(Origin = factor(Origin, levels = c("Wild", "Lab"))), "Origin",
                        "Wild vs Lab (n = 570; reference = Wild)"),
           "12_2_depth_adjusted_origin", paths$stats)
say("  Reading: if Kranji Marshes stays richest after adjustment and rarefaction,")
say("  its higher richness is not a depth artefact.")

# =============================================================================
# 3. WHERE DID THE CONTROL READS COME FROM? (wells and tags)
# =============================================================================
say("")
say("3. NEGATIVE-CONTROL READS: WELL-TO-WELL TRANSFER OR TAG MISASSIGNMENT?")
f_audit <- c(here("data", "raw", "metadata", "demux_metadata_audit.csv"),
             file.path(paths$species, "demux_metadata_audit.csv"))
f_audit <- f_audit[file.exists(f_audit)][1]
stopifnot("demux_metadata_audit.csv not found" = !is.na(f_audit))

plates <- c("AG_P01", "AG_P02", "ZX_P01", "ZX_P02", "BT_P01", "BT_P02")
wells <- read_csv(f_audit, show_col_types = FALSE) |>
  transmute(Sample_ID = paste0("sample_", SampleID),
            Plate = str_remove(Plate_ID, "_16S"), Well = str_to_upper(str_trim(Well))) |>
  filter(Plate %in% plates) |>
  mutate(row = match(str_sub(Well, 1, 1), LETTERS[1:8]),
         col = as.integer(str_remove(Well, "^[A-H]")),
         idx = (col - 1) * 8 + row,
         rtag = ((idx - 1) %% 30) + 1,          # reverse tag number (same on every plate)
         fblock = (idx - 1) %/% 30 + 1)          # forward-tag block (unique to plate)
all_ids <- c(rownames(bio), rownames(ctrl))
miss <- setdiff(all_ids, wells$Sample_ID)
if (length(miss)) stop("No well found for: ", paste(head(miss, 10), collapse = ", "))
wells <- filter(wells, Sample_ID %in% all_ids)
stopifnot(nrow(wells) == 576, !anyNA(wells$idx), all(wells$col %in% 1:12))
ctrl_w <- wells |> filter(Sample_ID %in% rownames(ctrl))
stopifnot(all(ctrl_w$Well == "H12"), all(ctrl_w$rtag == 6))
say("  Wells matched for all 576 samples; all 6 controls are H12, reverse tag 6.")

bio_rel  <- bio  / rowSums(bio)
ctrl_rel <- ctrl / rowSums(ctrl)
bc <- function(a, b) sum(abs(a - b)) / sum(a + b)
sw <- wells |> filter(Sample_ID %in% rownames(bio))
sw <- sw[match(rownames(bio), sw$Sample_ID), ]

rel_tbl <- map_dfr(seq_len(nrow(ctrl_w)), function(i) {
  cw <- ctrl_w[i, ]
  sw |> mutate(Control = cw$Sample_ID, Control_plate = cw$Plate,
               bray_curtis = apply(bio_rel, 1, function(s) bc(ctrl_rel[cw$Sample_ID, ], s)),
               same_plate = Plate == cw$Plate,
               dist = if_else(same_plate, sqrt((row - cw$row)^2 + (col - cw$col)^2), NA_real_),
               relation = case_when(
                 same_plate & fblock == cw$fblock ~ "Shares forward tag (own plate)",
                 same_plate & rtag == cw$rtag     ~ "Shares reverse tag (own plate)",
                 !same_plate & rtag == cw$rtag    ~ "Shares reverse tag (other plate)",
                 same_plate                       ~ "No shared tag (own plate)",
                 TRUE                             ~ "No shared tag (other plate)"),
               neighbour = same_plate & dist <= sqrt(2) & dist > 0)
})
save_table(select(rel_tbl, Control, Control_plate, Sample_ID, Plate, Well, rtag, fblock,
                  relation, dist, neighbour, bray_curtis),
           "12_3_control_vs_samples_by_well_and_tag", paths$stats)

# 3a. distance decay on own plate, tag-sharing wells removed
say("  3a. Distance decay on own plate (Spearman: Bray-Curtis vs distance from H12;")
say("      tag-sharing wells excluded; positive rho = closer wells more similar)")
A3 <- rel_tbl |> filter(same_plate, relation == "No shared tag (own plate)") |>
  group_by(Control, Control_plate) |>
  group_modify(function(d, k) {
    ct <- suppressWarnings(cor.test(d$dist, d$bray_curtis, method = "spearman", exact = FALSE))
    nb <- d$bray_curtis[d$neighbour]
    tibble(n = nrow(d), rho = unname(ct$estimate), p = ct$p.value,
           neighbours = paste(d$Well[d$neighbour], collapse = ","),
           bc_neighbours = if (length(nb)) mean(nb) else NA_real_,
           bc_rest = mean(d$bray_curtis[!d$neighbour]))
  }) |> ungroup()
save_table(A3, "12_3a_distance_decay", paths$stats)
for (i in seq_len(nrow(A3))) with(A3[i, ], say(sprintf(
  "    %-7s rho = %+.2f, %s (n = %d) | untagged neighbours %s: mean BC %.2f vs rest %.2f",
  Control_plate, rho, fmt_p(p), n, neighbours, bc_neighbours, bc_rest)))

# 3b. tag-sharing samples vs plate-matched random samples
say("  3b. Tag-sharing samples vs random samples from the same plates")
say("      (permutation, 9,999 draws; statistic = mean Bray-Curtis to the control;")
say("       p = fraction of random sets at least as similar as the tag-sharing set)")
perm_test <- function(d, n_perm = 9999) {
  tagged <- d |> filter(str_detect(relation, "^Shares"))
  pool   <- d |> filter(str_detect(relation, "^No shared"))
  obs    <- mean(tagged$bray_curtis)
  need   <- count(tagged, Plate)
  set.seed(4288)
  null <- replicate(n_perm, mean(unlist(map2(need$Plate, need$n, function(p, n)
    sample(pool$bray_curtis[pool$Plate == p], n)))))
  tibble(n_tagged = nrow(tagged), mean_bc_tagged = obs, mean_bc_random = mean(null),
         p = (sum(null <= obs) + 1) / (n_perm + 1))
}
B3_each <- rel_tbl |> group_by(Control, Control_plate) |>
  group_modify(~ perm_test(.x)) |> ungroup() |> mutate(set = "All tag-sharing samples")
B3_rev  <- rel_tbl |>
  filter(relation %in% c("Shares reverse tag (other plate)", "No shared tag (other plate)")) |>
  group_by(Control, Control_plate) |>
  group_modify(~ perm_test(.x)) |> ungroup() |> mutate(set = "Reverse-tag sharers on OTHER plates only")
B3 <- bind_rows(B3_each, B3_rev)
save_table(B3, "12_3b_tag_sharing_permutation", paths$stats)
for (s in unique(B3$set)) {
  say("    ", s, ":")
  d <- filter(B3, set == s)
  for (i in seq_len(nrow(d))) with(d[i, ], say(sprintf(
    "      %-7s n = %2d | mean BC tag-sharing %.2f vs random %.2f, %s",
    Control_plate, n_tagged, mean_bc_tagged, mean_bc_random, fmt_p(p))))
  say(sprintf("      Combined (Fisher's method across 6 controls): %s",
              fmt_p(pchisq(-2 * sum(log(d$p)), df = 2 * nrow(d), lower.tail = FALSE))))
}

# 3c. fingerprint OTUs
say("  3c. Fingerprint OTUs: control OTUs found in <= 5 biological samples")
prev <- colSums(bio > 0)
C3 <- map_dfr(seq_len(nrow(ctrl_w)), function(i) {
  cw   <- ctrl_w[i, ]
  otus <- colnames(ctrl)[ctrl[cw$Sample_ID, ] > 0 & prev <= 5 & prev > 0]
  map_dfr(otus, function(o) {
    hits <- rel_tbl |> filter(Control == cw$Sample_ID, Sample_ID %in% rownames(bio)[bio[, o] > 0])
    tibble(Control_plate = cw$Plate, OTU = o,
           Genus = otu_map$Genus[otu_map$OTU == o],
           reads_in_control = ctrl[cw$Sample_ID, o],
           n_samples = nrow(hits),
           samples = paste(sprintf("%s %s [%s]", hits$Plate, hits$Well,
                                   str_remove(hits$relation, " \\(.*")), collapse = "; "))
  })
})
save_table(C3, "12_3c_fingerprint_otus", paths$stats)
if (nrow(C3) == 0) say("    none") else for (i in seq_len(nrow(C3))) with(C3[i, ], say(sprintf(
  "    %-7s %-9s %-28s %4d reads | in %d sample(s): %s",
  Control_plate, OTU, str_trunc(Genus, 28), reads_in_control, n_samples, samples)))
say("  Reading: well-to-well transfer predicts 3a (rho > 0, neighbours most similar);")
say("  tag misassignment predicts 3b (tag sharers more similar, even on other plates).")

# Figure: Bray-Curtis by relation to the control
rel_plot <- rel_tbl |>
  mutate(relation = factor(relation, levels = c("Shares forward tag (own plate)",
                                                "Shares reverse tag (own plate)",
                                                "No shared tag (own plate)",
                                                "Shares reverse tag (other plate)",
                                                "No shared tag (other plate)")),
         Control_plate = factor(Control_plate, levels = plates))
p_src <- ggplot(rel_plot, aes(relation, 1 - bray_curtis)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#f3f2ee") +
  geom_jitter(data = filter(rel_plot, str_detect(relation, "^Shares")),
              width = 0.15, height = 0, size = 1, colour = "#2a78d6") +
  facet_wrap(~ Control_plate, nrow = 2) +
  scale_x_discrete(labels = c("Fwd tag\nown", "Rev tag\nown", "None\nown", "Rev tag\nother", "None\nother")) +
  labs(x = "Relation of sample to the H12 control", y = "Similarity to control (1 - Bray-Curtis)") +
  theme_fyp(9)
save_fig(p_src, "FigS_control_source", width = 190, height = 120)

writeLines(out, file.path(paths$stats, "12_results_summary.txt"))
message("Done. Summary: output/stats/12_results_summary.txt")
