# =============================================================================
# 05_classification.R - Figure 1: taxonomic classification of LC vs MC OTUs
#
#   A. Greengenes2 >=80% identity (default): deepest rank reached, LC vs MC
#      side by side - LC leaves far more OTUs unassigned even at default
#   B. Greengenes2 >=97% identity: same comparison at the stringent threshold
#   C. Greengenes2 (>=80%) vs Kraken2 on ALL OTUs: pie charts of the deepest
#      rank reached, LC and MC (counts from the "Updated Greengenes2 Kraken2
#      LC MC 28 September" spreadsheet, Summary rows)
#
# The Kraken2 rescue at confidence 0 vs 0.97 is reported as Table S2C,
# not plotted (per Priya, Oct 2026).
#
# Counts are the verified pipeline outputs (Table S2A):
#   GG2 nested counts: any / family+ / genus+ / species
# Nested counts are converted to exclusive "deepest rank reached" bins.
#
# Outputs
#   output/figures/Fig1_classification.pdf/.png
#   output/tables/05_classification_deepest_rank.tsv
#   output/tables/05_classification_gg2_vs_kraken2_all_OTUs.tsv
# =============================================================================

source(here::here("R", "00_setup.R"))

# "\u2265" (>=) needs the Cairo PDF device; fall back to ">=" without it
use_cairo <- isTRUE(capabilities("cairo"))
GE <- if (use_cairo) "\u2265" else ">="

# ---- 1. Verified Greengenes2 counts (Table S2A) ------------------------------
nested <- tribble(
  ~dataset, ~identity, ~total, ~any, ~family, ~genus, ~species,
  "LC",     80,        7642,   2772, 2622,    2263,   873,
  "LC",     97,        7642,    508,  490,     456,   170,
  "MC",     80,         867,    681,  658,     564,   233,
  "MC",     97,         867,    328,  316,     295,   106
)

rank_levels <- c("Species", "Genus", "Family", "Above family", "Unassigned")

bins <- nested |>
  mutate(Species = species, Genus = genus - species, Family = family - genus,
         `Above family` = any - family, Unassigned = total - any) |>
  select(dataset, identity, total, all_of(rank_levels))

stopifnot(all(as.matrix(bins[, rank_levels]) >= 0),
          all(rowSums(bins[, rank_levels]) == bins$total))
save_table(bins, "05_classification_deepest_rank")

long <- bins |>
  pivot_longer(all_of(rank_levels), names_to = "rank", values_to = "n") |>
  mutate(rank    = factor(rank, levels = rank_levels),
         pct     = 100 * n / total,
         dataset = factor(dataset, levels = c("LC", "MC")))

# LC and MC in clearly different colours
dataset_cols <- c("LC" = "#E69F00", "MC" = "#2a78d6")
dataset_labs <- c("LC" = "LC (0.1% filter; 7,642 OTUs)",
                  "MC" = "MC (0.5% / 20-read filter; 867 OTUs)")

ymax <- ceiling(max(long$pct) / 10) * 10 + 10

rank_axis_labs <- c("Species" = "Species", "Genus" = "Genus", "Family" = "Family",
                    "Above family" = "Above\nfamily", "Unassigned" = "Un-\nassigned")

rank_panel <- function(id, title) {
  d <- filter(long, identity == id)
  ggplot(d, aes(x = rank, y = pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.8),
              vjust = -0.4, size = 2.1) +
    scale_x_discrete(labels = rank_axis_labs) +
    scale_fill_manual(values = dataset_cols, labels = dataset_labs, name = NULL) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 100, 20),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Deepest rank assigned", y = "OTUs (%)", title = title) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 7.5, lineheight = 0.9),
          legend.position = "top")
}

p_80 <- rank_panel(80, paste0("Greengenes2, ", GE, "80% identity (default)"))
p_97 <- rank_panel(97, paste0("Greengenes2, ", GE, "97% identity"))

# ---- 2. Panel C: Greengenes2 vs Kraken2 on ALL OTUs (pie charts) ------------
# Source: "Updated Greengenes2 Kraken2 LC MC 28 September.xlsx", Summary rows
# (sheets LC and MC). Kraken2 = full classification of every OTU, no GG2
# pre-filter. Nested counts -> exclusive deepest-rank bins, as for GG2.
kraken_nested <- tribble(
  ~dataset, ~total, ~unclassified, ~below_family, ~family, ~genus, ~species,
  "LC",     7642,   4885,          934,           1823,    1351,   728,
  "MC",      867,    173,          228,            466,     342,   156
)
stopifnot(with(kraken_nested, unclassified + below_family + family == total))

kraken_bins <- kraken_nested |>
  transmute(dataset, classifier = "Kraken2", total,
            Species = species, Genus = genus - species, Family = family - genus,
            `Above family` = below_family, Unassigned = unclassified)

gg80_bins <- bins |>
  filter(identity == 80) |>
  transmute(dataset, classifier = paste0("Greengenes2 (", GE, "80%)"), total,
            Species, Genus, Family, `Above family`, Unassigned)

pie_bins <- bind_rows(gg80_bins, kraken_bins)
stopifnot(all(rowSums(pie_bins[, rank_levels]) == pie_bins$total))
save_table(pie_bins, "05_classification_gg2_vs_kraken2_all_OTUs")

rank_cols <- c("Species" = "#1d4f91", "Genus" = "#2a78d6", "Family" = "#7aa9e6",
               "Above family" = "#c3d7f2", "Unassigned" = "#d9d8d2")

pie_long <- pie_bins |>
  pivot_longer(all_of(rank_levels), names_to = "rank", values_to = "n") |>
  mutate(rank       = factor(rank, levels = rank_levels),
         pct        = 100 * n / total,
         classifier = factor(classifier, levels = unique(pie_bins$classifier)),
         dataset    = factor(paste0(dataset, " (", scales::comma(total), " OTUs)"),
                             levels = c("LC (7,642 OTUs)", "MC (867 OTUs)")))

p_pie <- ggplot(pie_long, aes(x = 1, y = pct, fill = rank)) +
  geom_col(width = 1, colour = "white", linewidth = 0.4) +
  geom_text(aes(x      = 1.2,             # push labels towards the rim
                label  = ifelse(pct >= 5, sprintf("%.0f%%", pct), ""),
                colour = ifelse(rank %in% c("Species", "Genus"), "white", "black"),
                group  = rank),            # stack labels in the same order as the slices
            position = position_stack(vjust = 0.5), size = 2.6) +
  scale_colour_identity() +
  coord_polar(theta = "y", direction = -1) +
  facet_grid(dataset ~ classifier, switch = "y") +
  scale_fill_manual(values = rank_cols, name = "Deepest rank assigned") +
  labs(title = "Greengenes2 vs Kraken2, all OTUs") +
  theme_void(base_size = 9) +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        strip.text  = element_text(face = "bold", size = 8.5, margin = margin(2, 2, 2, 2)),
        strip.text.y.left = element_text(angle = 90),
        legend.position = "right",
        legend.title = element_text(size = 8.5, face = "bold"),
        legend.text  = element_text(size = 8))

# ---- 3. Assemble ------------------------------------------------------------
top <- (p_80 | p_97) + plot_layout(guides = "collect") &
  theme(legend.position = "top")
fig1 <- top / p_pie + plot_layout(heights = c(1, 1.15))
fig1 <- fig1 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold"))

fig_w <- 200
fig_h <- 190
save_fig(fig1, "Fig1_classification", width = fig_w, height = fig_h)
if (use_cairo) {   # overwrite the PDF with a Cairo version so "\u2265" renders
  ggsave(file.path(paths$figures, "Fig1_classification.pdf"), fig1,
         width = fig_w, height = fig_h, units = "mm", device = cairo_pdf)
}

cat("\nDeepest-rank counts (% of OTUs):\n")
print(long |> select(dataset, identity, rank, n, pct) |>
        mutate(pct = round(pct, 1)), n = Inf)