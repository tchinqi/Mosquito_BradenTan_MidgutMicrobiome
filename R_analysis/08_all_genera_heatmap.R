# =============================================================================
# 08_all_genera_heatmap.R - Supplementary: every genus, every group in one view
# (each main figure's panel A shows its own slice; this combines them)
#
# Answers the feedback "show all genera, not just the top 8".
# A stacked bar cannot show 68 genera legibly, so all genera are shown as a
# heatmap (rows = genera, columns = groups):
#   A. mean relative abundance (%), log10 colour scale so rare genera are
#      visible next to dominant ones; white = not detected in that group
#   B. prevalence (% of mosquitoes with >= 1 read), same colours as Fig. 3A
# Mean relative abundance is calculated exactly as in the stacked bars
# (within-mosquito proportions averaged across mosquitoes).
# Rows are ordered by mean relative abundance across all 570 mosquitoes;
# "Unclassified at genus" is shown last.
#
# Groups: Wild, Lab | Mandai Zoo, Kranji Marshes, Bird Paradise, River Safari |
#         Lab Ae. albopictus, Lab Culex   (same as 07_core_microbiome.R)
#
# Outputs
#   output/figures/FigS_all_genera_heatmap.pdf/.png
#   output/stats/08_all_genera_mean_relative_abundance_by_group.tsv
#   output/stats/08_all_genera_prevalence_by_group.tsv
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

# ---- 1. Load ----------------------------------------------------------------
meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")

genus_of <- otu_map$Genus[match(colnames(bio), otu_map$OTU)]
gen      <- t(rowsum(t(bio), group = genus_of))     # 570 x 68
gen_rel  <- gen / rowSums(gen)
stopifnot(ncol(gen) == 68, all(abs(rowSums(gen_rel) - 1) < 1e-9))

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Location, Lab_Taxon), by = "Sample_ID")

# ---- 2. Groups (identical to script 07) -------------------------------------
group_levels <- c("Wild", "Lab", "Mandai Zoo", "Kranji Marshes", "Bird Paradise",
                  "River Safari", "Ae. albopictus", "Culex")
set_levels   <- c("Origin", "Wild sites", "Lab species")
membership <- bind_rows(
  sd |> filter(Origin == "Wild")  |> transmute(Sample_ID, Group = "Wild", Set = "Origin"),
  sd |> filter(Origin == "Lab")   |> transmute(Sample_ID, Group = "Lab",  Set = "Origin"),
  sd |> filter(!is.na(Location))  |> transmute(Sample_ID, Group = as.character(Location), Set = "Wild sites"),
  sd |> filter(!is.na(Lab_Taxon)) |> transmute(Sample_ID, Group = Lab_Taxon, Set = "Lab species")
) |>
  mutate(Group = factor(Group, levels = group_levels), Set = factor(Set, levels = set_levels))
group_n <- count(membership, Set, Group)
stopifnot(identical(as.integer(group_n$n), c(380L, 190L, 194L, 88L, 39L, 24L, 165L, 25L)))
n_of <- setNames(group_n$n, as.character(group_n$Group))

grp_short <- c("Wild" = "Wild", "Lab" = "Lab", "Mandai Zoo" = "MZ", "Kranji Marshes" = "KM",
               "Bird Paradise" = "BP", "River Safari" = "RS",
               "Ae. albopictus" = "Ae. alb.", "Culex" = "Culex")
x_lab_expr <- setNames(as.expression(lapply(group_levels, function(g) {
  if (g %in% c("Ae. albopictus", "Culex")) bquote(atop(italic(.(grp_short[[g]])), "(" * .(n_of[[g]]) * ")"))
  else bquote(atop(.(grp_short[[g]]), "(" * .(n_of[[g]]) * ")"))
})), group_levels)

# ---- 3. Mean relative abundance + prevalence per genus per group ------------
genera <- colnames(gen)
stats <- membership |>
  group_by(Set, Group) |>
  group_modify(function(d, k) {
    tibble(Genus      = genera,
           mean_rel   = unname(100 * colMeans(gen_rel[d$Sample_ID, , drop = FALSE])),
           prevalence = unname(100 * colMeans(gen[d$Sample_ID, , drop = FALSE] > 0)))
  }) |>
  ungroup()

# Row order: overall mean relative abundance (all 570), Unclassified last
overall <- sort(colMeans(gen_rel), decreasing = TRUE)
row_order <- c(setdiff(names(overall), "Unclassified_at_genus"), "Unclassified_at_genus")
stats <- stats |> mutate(Genus = factor(Genus, levels = rev(row_order)))

save_table(stats |> select(Genus, Group, mean_rel) |>
             pivot_wider(names_from = Group, values_from = mean_rel) |>
             arrange(desc(Genus)),
           "08_all_genera_mean_relative_abundance_by_group", paths$stats)
save_table(stats |> select(Genus, Group, prevalence) |>
             pivot_wider(names_from = Group, values_from = prevalence) |>
             arrange(desc(Genus)),
           "08_all_genera_prevalence_by_group", paths$stats)

g_expr <- setNames(as.expression(lapply(row_order, function(g) {
  if (g == "Unclassified_at_genus") pretty_genus(g) else bquote(italic(.(pretty_genus(g))))
})), row_order)

# ---- 4. Heatmaps ------------------------------------------------------------
heat_theme <- theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    axis.text.y = element_text(size = 6.5),
                    axis.text.x = element_text(size = 7.5),
                    panel.spacing.x = unit(2, "mm"),
                    strip.text = element_text(size = 8.5),
                    legend.position = "top",
                    legend.title = element_text(size = 8), legend.text = element_text(size = 7))

ab_breaks <- c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100)
p_ab <- stats |>
  mutate(mean_rel = if_else(mean_rel > 0, mean_rel, NA_real_)) |>
  ggplot(aes(x = Group, y = Genus, fill = mean_rel)) +
  geom_tile(colour = "grey88", linewidth = 0.25) +
  facet_grid(~ Set, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10",
                       limits = c(1e-4, 100), breaks = ab_breaks,
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       oob = scales::squish, na.value = "white",
                       name = "Mean relative\nabundance (%)") +
  scale_x_discrete(labels = function(b) x_lab_expr[b], expand = c(0, 0)) +
  scale_y_discrete(labels = function(b) g_expr[b], expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = unit(45, "mm"), barheight = unit(3, "mm"),
                                title.position = "top", title.hjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  heat_theme

prev_palette <- c("#f7f7f5", "#c3d7f2", "#7aa9e6", "#2a78d6", "#1d4f91")   # as Fig. 3A
p_pr <- stats |>
  mutate(prevalence = if_else(prevalence > 0, prevalence, NA_real_)) |>
  ggplot(aes(x = Group, y = Genus, fill = prevalence)) +
  geom_tile(colour = "grey88", linewidth = 0.25) +
  facet_grid(~ Set, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = prev_palette, limits = c(0, 100), na.value = "white",
                       name = "Prevalence\n(% of mosquitoes)") +
  scale_x_discrete(labels = function(b) x_lab_expr[b], expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = unit(35, "mm"), barheight = unit(3, "mm"),
                                title.position = "top", title.hjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  heat_theme +
  theme(axis.text.y = element_blank())

fig <- (p_ab | p_pr) + plot_annotation(tag_levels = "A")
save_fig(fig, "FigS_all_genera_heatmap", width = 200, height = 250)

# ---- 5. Console summary -----------------------------------------------------
n_detected <- stats |> group_by(Group) |> summarise(genera_detected = sum(prevalence > 0), .groups = "drop")
cat("\nGenus labels detected per group (of 68, incl. Unclassified at genus):\n")
print(n_detected)
cat("\nGenera detected in only one of the 8 groups:\n")
print(stats |> filter(prevalence > 0) |> count(Genus, name = "groups_detected") |>
        filter(groups_detected == 1) |> arrange(Genus), n = Inf)
