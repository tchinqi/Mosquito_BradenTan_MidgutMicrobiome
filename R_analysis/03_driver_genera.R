# =============================================================================
# 03_driver_genera.R - Supplementary figure: Which genera differ between Wild and Lab?
# (updated 2026-10-06: Unclassified at genus excluded; direction from Wilcoxon ranks)
#
# For every genus present in >= 10% of the 570 biological samples:
#   - mean relative abundance + prevalence (% of mosquitoes with >= 1 read)
#   - Wilcoxon rank-sum test on relative abundance (Wild vs Lab)
#   - Fisher's exact test on prevalence (Wild vs Lab)
#   - Benjamini-Hochberg correction across genera, separately for each test
# Relative abundances are compositional: an increase in one genus lowers the
# others, so these are differences in proportion, not absolute load.
#
# Outputs
#   output/figures/Fig_driver_genera.pdf/.png
#   output/stats/03_driver_genera_wild_vs_lab.tsv   (all tested genera)
#   output/stats/03_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")

origin <- meta$Origin[match(rownames(bio), meta$Sample_ID)] |> droplevels()
stopifnot(!anyNA(origin), identical(levels(origin), c("Wild", "Lab")))

# ---- 1. Genus x sample matrix (all 68 genus labels) -------------------------
genus_of <- otu_map$Genus[match(colnames(bio), otu_map$OTU)]
stopifnot(!anyNA(genus_of))
gen_counts <- t(rowsum(t(bio), group = genus_of))      # samples x genera
gen_rel    <- gen_counts / rowSums(gen_counts)
stopifnot(ncol(gen_counts) == 68, all(abs(rowSums(gen_rel) - 1) < 1e-9))

# ---- 2. Filter: present in >= 10% of samples --------------------------------
prev_all <- colMeans(gen_counts > 0)
# "Unclassified_at_genus" is not a genus, so it is excluded from testing (Methods 2.6)
tested   <- setdiff(names(prev_all)[prev_all >= 0.10], "Unclassified_at_genus")
message(length(tested), " genera present in >= 10% of samples are tested.")

# ---- 3. Tests ---------------------------------------------------------------
is_w <- origin == "Wild"
is_l <- origin == "Lab"

res <- map_dfr(tested, function(g) {
  x <- gen_rel[is_w, g]; y <- gen_rel[is_l, g]
  pw <- mean(gen_counts[is_w, g] > 0); pl <- mean(gen_counts[is_l, g] > 0)
  ft <- fisher.test(matrix(c(sum(gen_counts[is_w, g] > 0), sum(gen_counts[is_w, g] == 0),
                             sum(gen_counts[is_l, g] > 0), sum(gen_counts[is_l, g] == 0)),
                           nrow = 2))
  wt <- wilcox.test(x, y, exact = FALSE)
  tibble(genus = g,
         mean_pct_wild = 100 * mean(x), mean_pct_lab = 100 * mean(y),
         diff_pct_points = 100 * (mean(x) - mean(y)),
         prev_pct_wild = 100 * pw, prev_pct_lab = 100 * pl,
         W = as.numeric(wt$statistic), p_abundance = wt$p.value,
         fisher_OR = as.numeric(ft$estimate), p_prevalence = ft$p.value)
}) |>
  mutate(p_abundance_BH  = p.adjust(p_abundance,  method = "BH"),
         p_prevalence_BH = p.adjust(p_prevalence, method = "BH"),
         # Direction follows the Wilcoxon RANKS (what the test compares), not the
         # mean: W above its null expectation n_wild * n_lab / 2 means Wild ranks higher.
         rank_higher = if_else(W > sum(is_w) * sum(is_l) / 2, "Wild", "Lab"),
         mean_higher = if_else(diff_pct_points > 0, "Wild", "Lab"),
         higher_in   = if_else(p_abundance_BH >= 0.05, "n.s.", rank_higher)) |>
  arrange(desc(abs(diff_pct_points)))
save_table(res, "03_driver_genera_wild_vs_lab", paths$stats)

# ---- 4. Figure: dumbbell plots of the 20 most abundant tested genera --------
top <- res |>
  mutate(overall = (mean_pct_wild * sum(is_w) + mean_pct_lab * sum(is_l)) / length(origin)) |>
  slice_max(overall, n = 20, with_ties = FALSE) |>
  arrange(diff_pct_points) |>                      # Lab-higher at bottom, Wild-higher at top
  mutate(genus = factor(genus, levels = genus))

y_lab_expr <- as.expression(lapply(levels(top$genus), function(g) {
  lab_g <- pretty_genus(g)
  if (g == "Unclassified_at_genus") lab_g else bquote(italic(.(lab_g)))
}))

star <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

long_ab <- top |>
  select(genus, Wild = mean_pct_wild, Lab = mean_pct_lab) |>
  pivot_longer(c(Wild, Lab), names_to = "Origin", values_to = "pct") |>
  mutate(Origin = factor(Origin, levels = c("Wild", "Lab")))
long_pr <- top |>
  select(genus, Wild = prev_pct_wild, Lab = prev_pct_lab) |>
  pivot_longer(c(Wild, Lab), names_to = "Origin", values_to = "pct") |>
  mutate(Origin = factor(Origin, levels = c("Wild", "Lab")))

dumbbell <- function(long, wide_lo, wide_hi, stars, xlab, xmax) {
  ggplot(long, aes(x = pct, y = genus)) +
    geom_segment(data = top, aes(x = {{ wide_lo }}, xend = {{ wide_hi }}, y = genus, yend = genus),
                 inherit.aes = FALSE, colour = "grey70", linewidth = 0.8) +
    geom_point(aes(colour = Origin), size = 2.4) +
    geom_text(data = top, aes(x = xmax, y = genus, label = {{ stars }}),
              inherit.aes = FALSE, hjust = 1, size = 3.2) +
    scale_colour_manual(values = origin_colours, name = NULL) +
    scale_y_discrete(labels = y_lab_expr) +
    scale_x_continuous(limits = c(0, xmax), breaks = function(l) scales::breaks_pretty()(c(0, min(l[2], 100))), expand = expansion(mult = c(0.01, 0.02))) +    labs(x = xlab, y = NULL) +
    theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
          legend.position = "top")
}

xmax_ab <- ceiling(max(long_ab$pct) / 5) * 5 + 5
p_ab <- dumbbell(long_ab, pmin(mean_pct_wild, mean_pct_lab), pmax(mean_pct_wild, mean_pct_lab),
                 star(p_abundance_BH), "Mean relative abundance (%)", xmax_ab)
p_pr <- dumbbell(long_pr, pmin(prev_pct_wild, prev_pct_lab), pmax(prev_pct_wild, prev_pct_lab),
                 star(p_prevalence_BH), "Prevalence (% of mosquitoes)", 118) +
  theme(axis.text.y = element_blank())

fig <- (p_ab | p_pr) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")
save_fig(fig, "Fig_driver_genera", width = 190, height = 140)

# ---- 5. Summary -------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, digits = 2, format = "g"))
line_for <- function(r) sprintf(
  "  %-26s Wild %5.2f%% vs Lab %5.2f%% (BH p %s) | prevalence %3.0f%% vs %3.0f%% (BH p %s)",
  pretty_genus(r$genus), r$mean_pct_wild, r$mean_pct_lab, fmt_p(r$p_abundance_BH),
  r$prev_pct_wild, r$prev_pct_lab, fmt_p(r$p_prevalence_BH))

summary_lines <- c(
  "03_driver_genera.R - Wild vs Lab, per genus",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Genera tested (present in >= 10%% of 570 samples, excl. Unclassified at genus): %d", length(tested)),
  "Direction = Wilcoxon ranks. Genera where mean and rank direction disagree:",
  if (any(res$rank_higher != res$mean_higher & res$higher_in != "n.s."))
    sprintf("  %s: ranks higher in %s, mean higher in %s",
            pretty_genus(res$genus[res$rank_higher != res$mean_higher & res$higher_in != "n.s."]),
            res$rank_higher[res$rank_higher != res$mean_higher & res$higher_in != "n.s."],
            res$mean_higher[res$rank_higher != res$mean_higher & res$higher_in != "n.s."])
  else "  none",
  sprintf("Higher relative abundance in Wild (BH p < 0.05): %d | in Lab: %d | n.s.: %d",
          sum(res$higher_in == "Wild"), sum(res$higher_in == "Lab"), sum(res$higher_in == "n.s.")),
  "",
  "HIGHER IN WILD (sorted by difference in mean %)",
  unlist(lapply(split(res, seq_len(nrow(res)))[res$higher_in == "Wild"], line_for)),
  "",
  "HIGHER IN LAB",
  unlist(lapply(split(res, seq_len(nrow(res)))[res$higher_in == "Lab"], line_for)),
  "",
  "NOT SIGNIFICANT (abundance)",
  unlist(lapply(split(res, seq_len(nrow(res)))[res$higher_in == "n.s."], line_for))
)
writeLines(summary_lines, file.path(paths$stats, "03_results_summary.txt"))
cat(summary_lines, sep = "\n")