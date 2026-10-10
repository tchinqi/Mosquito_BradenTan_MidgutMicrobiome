# =============================================================================
# 05_classification.R - Figure 1: Greengenes2 vs Kraken2 classification
#
#   A. Greengenes2: deepest rank reached by LC and MC OTUs at >=80% / >=97%
#   B. Kraken2 rescue of GG2-80 unassigned + d__Bacteria-only OTUs,
#      confidence 0 (default) vs 0.97 - same input sequences in both runs
#   C. MC, GG2 >=97%: % of OTUs vs % of reads that are classified
#
# Counts are the verified pipeline outputs (Table 1, Table 2, script 44):
#   GG2 nested counts: any / family+ / genus+ / species
#   Kraken2 summaries: data/raw/kraken_rescue_summary/
#   MC reads: 1,682,795 total, 1,432,996 in GG2-97 classified OTUs
# Nested counts are converted to exclusive "deepest rank reached" bins.
# NOTE: Kraken2 confidence is k-mer support, not sequence identity.
#
# Outputs
#   output/figures/Fig1_classification.pdf/.png
#   output/tables/05_classification_deepest_rank.tsv
# =============================================================================

source(here::here("R", "00_setup.R"))

# "\u2265" (>=) needs the Cairo PDF device; fall back to ">=" without it
use_cairo <- isTRUE(capabilities("cairo"))
GE <- if (use_cairo) "\u2265" else ">="

# ---- 1. Verified counts -----------------------------------------------------
nested <- tribble(
  ~panel,        ~dataset, ~setting,          ~total, ~any, ~family, ~genus, ~species,
  "Greengenes2", "LC",     paste0(GE, "80% identity"), 7642, 2772,   2622,   2263,      873,
  "Greengenes2", "LC",     paste0(GE, "97% identity"), 7642,  508,    490,    456,      170,
  "Greengenes2", "MC",     paste0(GE, "80% identity"),  867,  681,    658,    564,      233,
  "Greengenes2", "MC",     paste0(GE, "97% identity"),  867,  328,    316,    295,      106,
  "Kraken2",     "LC",     "Confidence 0",      4884,   86,     77,     75,       59,
  "Kraken2",     "LC",     "Confidence 0.97",   4884,    0,      0,      0,        0,
  "Kraken2",     "MC",     "Confidence 0",       188,   16,     15,     14,        7,
  "Kraken2",     "MC",     "Confidence 0.97",    188,    0,      0,      0,        0
)

rank_levels <- c("Species", "Genus", "Family", "Above family", "Unassigned")
rank_cols   <- c("Species" = "#1d4f91", "Genus" = "#2a78d6", "Family" = "#7aa9e6",
                 "Above family" = "#c3d7f2", "Unassigned" = "#d9d8d2")

bins <- nested |>
  mutate(Species = species, Genus = genus - species, Family = family - genus,
         `Above family` = any - family, Unassigned = total - any) |>
  select(panel, dataset, setting, total, all_of(rank_levels))

stopifnot(all(as.matrix(bins[, rank_levels]) >= 0),
          all(rowSums(bins[, rank_levels]) == bins$total))
save_table(bins, "05_classification_deepest_rank")

long <- bins |>
  pivot_longer(all_of(rank_levels), names_to = "rank", values_to = "n") |>
  mutate(rank = factor(rank, levels = rev(rank_levels)),     # Unassigned at the left
         pct  = 100 * n / total,
         bar  = paste0(dataset, ", ", setting))

# Percent assigned, printed at the end of each bar
assigned_lab <- bins |>
  mutate(bar = paste0(dataset, ", ", setting),
         label = sprintf("%s / %s assigned (%.1f%%)",
                         scales::comma(total - Unassigned), scales::comma(total),
                         100 * (total - Unassigned) / total))

stack_panel <- function(which, title) {
  d   <- filter(long, panel == which)
  lab <- filter(assigned_lab, panel == which)
  lev <- rev(unique(d$bar))                       # first row on top
  d$bar   <- factor(d$bar, levels = lev)
  lab$bar <- factor(lab$bar, levels = lev)
  ggplot(d, aes(x = pct, y = bar, fill = rank)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(data = lab, aes(x = 101, y = bar, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.7) +
    scale_fill_manual(values = rank_cols, breaks = rank_levels,
                      name = "Deepest rank reached", drop = FALSE) +
    scale_x_continuous(limits = c(0, 160), breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    guides(x = guide_axis(cap = "both")) +
    labs(x = "OTUs (%)", y = NULL, title = title) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

p_gg <- stack_panel("Greengenes2", "Greengenes2 (all OTUs)")
p_kr <- stack_panel("Kraken2",
                    paste0("Kraken2 (GG2 ", GE, "80% unassigned or domain-only OTUs)"))

# ---- 2. Panel C: OTUs vs reads for the final dataset ------------------------
otu_vs_reads <- tibble(
  measure    = factor(c("OTUs", "Reads"), levels = c("Reads", "OTUs")),
  classified = c(328, 1432996),
  total      = c(867, 1682795)
) |>
  mutate(pct = 100 * classified / total,
         label = sprintf("%s / %s (%.1f%%)", scales::comma(classified),
                         scales::comma(total), pct))

p_reads <- ggplot(otu_vs_reads, aes(y = measure)) +
  geom_col(aes(x = 100), width = 0.6, fill = rank_cols[["Unassigned"]]) +
  geom_col(aes(x = pct), width = 0.6, fill = rank_cols[["Genus"]]) +
  geom_text(aes(x = 101, label = label), hjust = 0, size = 2.7) +
  scale_x_continuous(limits = c(0, 160), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  guides(x = guide_axis(cap = "both")) +
  labs(x = paste0("Classified by Greengenes2 ", GE, "97% (blue)"), y = NULL,
       title = paste0("Final dataset (MC, GG2 ", GE, "97%): OTUs vs reads")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ---- 3. Assemble ------------------------------------------------------------
fig1 <- p_gg / p_kr / p_reads +
  plot_layout(heights = c(4, 4, 2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8.5, face = "bold"),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(3.5, "mm"))
save_fig(fig1, "Fig1_classification", width = 180, height = 150)
if (use_cairo) {   # overwrite the PDF with a Cairo version so "\u2265" renders
  ggsave(file.path(paths$figures, "Fig1_classification.pdf"), fig1,
         width = 180, height = 150, units = "mm", device = cairo_pdf)
}

cat("\nDeepest-rank counts:\n")
print(bins, width = Inf)
