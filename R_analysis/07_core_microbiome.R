# =============================================================================
# 07_core_microbiome.R - Figure 3: core midgut bacteria across origin, wild
#                        collection sites and laboratory species
#
# Definition (Methods 2.11):
#   - present   = >= 1 read assigned to the genus in a mosquito. Under the MC
#                 workflow every retained consensus has >= 20 reads, so
#                 "present" means >= 20 reads.
#   - prevalence = % of mosquitoes in a group in which the genus is present
#   - core      = prevalence >= 90% within a group (Baltar et al., 2023)
#   - sensitivity thresholds: 50% (Rosso et al., 2018; Custer et al., 2023)
#                             and 70%
#   - shared core = genera that are core in EVERY group
#   - "Unclassified at genus" is not a genus and is excluded
#
# Groups (a mosquito can belong to more than one group):
#   Origin      : Wild (380), Lab (190)
#   Wild sites  : Mandai Zoo (194), Kranji Marshes (88), Bird Paradise (39),
#                 River Safari (24)  - the 345 location-assigned wild samples
#   Lab species : Ae. albopictus (165), Culex (25)
#   Wild host species are not split here (see 11_host_species.R).
#
# Figure panels
#   A. Prevalence heatmap: every genus that reaches >= 50% prevalence in at
#      least one group; cells outlined = core at 90%
#   B. Number of core genera per group across thresholds (50-100%)
#   C. Share of each group's reads that come from the shared core,
#      group-specific core and non-core genera (90% definition)
#
# Outputs
#   output/figures/Fig3_core_microbiome.pdf/.png
#   output/stats/07_*.tsv
#   output/stats/07_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

core_threshold   <- 90                       # primary definition
sens_thresholds  <- c(50, 70)                # sensitivity
curve_thresholds <- seq(50, 100, by = 10)    # panel B
is_core <- function(prev, t) prev >= t - 1e-9   # tolerance for e.g. 342/380 = 90%

# ---- 1. Load ----------------------------------------------------------------
meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")
ctrl    <- read_otu_matrix("controls")

genus_of <- function(m) otu_map$Genus[match(colnames(m), otu_map$OTU)]
gen      <- t(rowsum(t(bio), group = genus_of(bio)))     # 570 samples x 68 genus labels
gen_rel  <- gen / rowSums(gen)
stopifnot(ncol(gen) == 68, all(rowSums(gen) > 0), all(abs(rowSums(gen_rel) - 1) < 1e-9))

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Location, Lab_Taxon), by = "Sample_ID")
stopifnot(identical(sd$Sample_ID, rownames(gen)))

# ---- 2. Groups --------------------------------------------------------------
group_levels <- c("Wild", "Lab", "Mandai Zoo", "Kranji Marshes", "Bird Paradise",
                  "River Safari", "Ae. albopictus", "Culex")
set_levels   <- c("Origin", "Wild sites", "Lab species")

membership <- bind_rows(
  sd |> filter(Origin == "Wild")  |> transmute(Sample_ID, Group = "Wild", Set = "Origin"),
  sd |> filter(Origin == "Lab")   |> transmute(Sample_ID, Group = "Lab",  Set = "Origin"),
  sd |> filter(!is.na(Location))  |> transmute(Sample_ID, Group = as.character(Location), Set = "Wild sites"),
  sd |> filter(!is.na(Lab_Taxon)) |> transmute(Sample_ID, Group = Lab_Taxon, Set = "Lab species")
) |>
  mutate(Group = factor(Group, levels = group_levels),
         Set   = factor(Set,   levels = set_levels))

group_n <- count(membership, Set, Group)
stopifnot(!anyNA(membership$Group),
          identical(as.integer(group_n$n), c(380L, 190L, 194L, 88L, 39L, 24L, 165L, 25L)))
n_of     <- setNames(group_n$n, as.character(group_n$Group))
n_groups <- length(group_levels)

# Labels: short names, species in italics
grp_short <- c("Wild" = "Wild", "Lab" = "Lab", "Mandai Zoo" = "MZ", "Kranji Marshes" = "KM",
               "Bird Paradise" = "BP", "River Safari" = "RS",
               "Ae. albopictus" = "Ae. alb.", "Culex" = "Culex")
is_species <- function(g) g %in% c("Ae. albopictus", "Culex")
x_lab_expr <- setNames(as.expression(lapply(group_levels, function(g) {
  if (is_species(g)) bquote(atop(italic(.(grp_short[[g]])), "(" * .(n_of[[g]]) * ")"))
  else               bquote(atop(.(grp_short[[g]]), "(" * .(n_of[[g]]) * ")"))
})), group_levels)
row_lab_expr <- setNames(as.expression(c(lapply(group_levels, function(g) {
  if (is_species(g)) bquote(italic(.(grp_short[[g]]))) else grp_short[[g]]
}), list(quote(bold("Shared core"))))), c(group_levels, "Shared core"))
genus_expr <- function(gs) setNames(as.expression(lapply(gs, function(g) {
  if (g == "Unclassified_at_genus") pretty_genus(g) else bquote(italic(.(pretty_genus(g))))
})), gs)

# ---- 3. Prevalence and mean relative abundance per genus per group ----------
genera <- colnames(gen)
grp_stats <- membership |>
  group_by(Set, Group) |>
  group_modify(function(d, k) {
    m <- gen[d$Sample_ID, , drop = FALSE]
    r <- gen_rel[d$Sample_ID, , drop = FALSE]
    tibble(Genus      = genera,
           n          = nrow(m),
           prevalence = unname(100 * colMeans(m > 0)),
           mean_rel   = unname(100 * colMeans(r)))
  }) |>
  ungroup()

save_table(grp_stats |>
             select(Genus, Group, prevalence) |>
             pivot_wider(names_from = Group, values_from = prevalence) |>
             arrange(desc(Wild + Lab)),
           "07_genus_prevalence_by_group", paths$stats)

# ---- 4. Core membership -----------------------------------------------------
core_tab <- grp_stats |>
  filter(Genus != "Unclassified_at_genus") |>
  mutate(core_50 = is_core(prevalence, 50),
         core_70 = is_core(prevalence, 70),
         core_90 = is_core(prevalence, core_threshold))

shared_at <- function(t) {
  core_tab |>
    filter(is_core(prevalence, t)) |>
    group_by(Genus) |>
    summarise(k = n_distinct(Group), .groups = "drop") |>
    filter(k == n_groups) |>
    pull(Genus)
}
shared_core <- shared_at(core_threshold)
member90    <- core_tab |> filter(core_90) |> transmute(Group = as.character(Group), Genus)

save_table(core_tab |> filter(core_50) |> arrange(Group, desc(prevalence)),
           "07_core_membership_50_70_90", paths$stats)

# ---- 5. Core size across thresholds (Fig. 3B) -------------------------------
core_size <- map_dfr(curve_thresholds, function(t) {
  per_group <- core_tab |>
    group_by(Group) |>
    summarise(n_core = sum(is_core(prevalence, t)), .groups = "drop") |>
    mutate(Group = as.character(Group))
  bind_rows(per_group, tibble(Group = "Shared core", n_core = length(shared_at(t)))) |>
    mutate(threshold = t)
})
save_table(core_size |> pivot_wider(names_from = threshold, values_from = n_core,
                                    names_prefix = "core_at_"),
           "07_core_size_by_threshold", paths$stats)

# ---- 6. Read share of shared / group-specific / non-core (Fig. 3C) ----------
part_levels <- c("Shared core", "Group-specific core", "Not core")
read_share <- membership |>
  group_by(Set, Group) |>
  group_modify(function(d, k) {
    r    <- gen_rel[d$Sample_ID, , drop = FALSE]
    own  <- member90$Genus[member90$Group == as.character(k$Group)]
    spec <- setdiff(own, shared_core)
    sh   <- 100 * mean(rowSums(r[, shared_core, drop = FALSE]))
    sp   <- 100 * mean(rowSums(r[, spec, drop = FALSE]))
    tibble(Part = part_levels, pct = c(sh, sp, 100 - sh - sp),
           genera = c(paste(pretty_genus(shared_core), collapse = ", "),
                      paste(pretty_genus(spec), collapse = ", "), ""))
  }) |>
  ungroup() |>
  mutate(Part = factor(Part, levels = part_levels))
save_table(read_share, "07_core_read_share", paths$stats)

# ---- 7. Core genera in the negative controls --------------------------------
ctrl_gen   <- t(rowsum(t(ctrl), group = genus_of(ctrl)))       # 6 controls x genera
core_union <- sort(unique(member90$Genus))
ctrl_check <- map_dfr(core_union, function(g) {
  cr <- ctrl_gen[, g]; br <- gen[, g]
  tibble(Genus = g,
         shared_core                   = g %in% shared_core,
         controls_with_genus           = sum(cr > 0),
         control_reads_total           = sum(cr),
         control_reads_max_single      = max(cr),
         pct_of_all_control_reads      = 100 * sum(cr) / sum(ctrl_gen),
         bio_prevalence_pct            = 100 * mean(br > 0),
         bio_reads_total               = sum(br),
         bio_median_reads_when_present = median(br[br > 0]),
         control_share_of_genus_reads_pct = 100 * sum(cr) / (sum(cr) + sum(br)))
})
save_table(ctrl_check, "07_core_genera_in_negative_controls", paths$stats)
save_table(as_tibble(ctrl_gen[, shared_core, drop = FALSE], rownames = "Control") |>
             mutate(total_control_reads = as.numeric(rowSums(ctrl_gen))),
           "07_shared_core_reads_per_control", paths$stats)

# ---- 8. Figure 3 ------------------------------------------------------------
prev_palette <- c("#f7f7f5", "#c3d7f2", "#7aa9e6", "#2a78d6", "#1d4f91")
heat_theme <- theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    panel.spacing.x = unit(2.5, "mm"),
                    strip.text = element_text(size = 9),
                    legend.title = element_text(size = 8.5), legend.text = element_text(size = 8))

# A. Prevalence heatmap
heat_genera <- core_tab |>
  group_by(Genus) |>
  summarise(maxp = max(prevalence), meanp = mean(prevalence), .groups = "drop") |>
  filter(is_core(maxp, min(sens_thresholds))) |>
  arrange(desc(meanp)) |>
  pull(Genus)
ha <- core_tab |>
  filter(Genus %in% heat_genera) |>
  mutate(Genus = factor(Genus, levels = rev(heat_genera)),
         face  = if_else(core_90, "bold", "plain"))
g_expr <- genus_expr(heat_genera)

p_heat <- ggplot(ha, aes(x = Group, y = Genus, fill = prevalence)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(ha, core_90), fill = NA, colour = "black", linewidth = 0.7) +
  geom_text(aes(label = if_else(round(prevalence) >= 90 & !core_90,
                                sprintf("%.1f", prevalence),        # e.g. 89.7, not "90"
                                as.character(round(prevalence))),
                colour = prevalence > 60, fontface = face),
            size = 2.6, show.legend = FALSE) +
  facet_grid(~ Set, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = prev_palette, limits = c(0, 100),
                       name = "Prevalence\n(% of mosquitoes)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(labels = function(b) x_lab_expr[b], expand = c(0, 0)) +
  scale_y_discrete(labels = function(b) g_expr[b], expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  heat_theme

# B. Core size across thresholds
cs <- core_size |>
  mutate(Row = factor(Group, levels = rev(c(group_levels, "Shared core"))),
         threshold = factor(threshold, levels = curve_thresholds))
col_90 <- which(curve_thresholds == core_threshold)

p_size <- ggplot(cs, aes(x = threshold, y = Row, fill = n_core)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n_core, colour = n_core > 0.55 * max(n_core)), size = 2.7,
            show.legend = FALSE) +
  geom_hline(yintercept = 1.5, colour = "grey40", linewidth = 0.4) +
  annotate("rect", xmin = col_90 - 0.5, xmax = col_90 + 0.5,
           ymin = 0.5, ymax = nlevels(cs$Row) + 0.5,
           fill = NA, colour = "black", linewidth = 0.7) +
  scale_fill_gradientn(colours = prev_palette, name = "Core\ngenera",
                       breaks = function(l) unique(round(scales::breaks_pretty(4)(l)))) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(labels = function(b) paste0("≥", b, "%"), expand = c(0, 0)) +
  scale_y_discrete(labels = function(b) row_lab_expr[b], expand = c(0, 0)) +
  labs(x = "Prevalence threshold", y = NULL) +
  heat_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# C. Read share
part_cols <- c("Shared core" = "#1d4f91", "Group-specific core" = "#7aa9e6", "Not core" = "#d9d8d2")
rs <- read_share |> mutate(Group = factor(as.character(Group), levels = rev(group_levels)))

p_share <- ggplot(rs, aes(x = pct, y = Group, fill = Part)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = if_else(pct >= 8, sprintf("%.0f", pct), ""), colour = Part),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.5, show.legend = FALSE) +
  facet_grid(Set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = part_cols, name = NULL) +
  scale_colour_manual(values = c("Shared core" = "white", "Group-specific core" = "black",
                                 "Not core" = "black"), guide = "none") +
  scale_x_continuous(limits = c(0, 100.01), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_y_discrete(labels = function(b) row_lab_expr[b]) +
  labs(x = "Mean share of reads", y = NULL) +
  theme(strip.text.y = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 8), legend.key.size = unit(3.5, "mm"),
        panel.spacing.y = unit(2, "mm"))

fig3 <- p_heat / (p_size | p_share) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A")
save_fig(fig3, "Fig3_core_microbiome", width = 200, height = 235)
if (isTRUE(capabilities("cairo"))) {   # Cairo PDF so the ">=" symbol renders
  ggsave(file.path(paths$figures, "Fig3_core_microbiome.pdf"), fig3,
         width = 200, height = 235, units = "mm", device = cairo_pdf)
}

# ---- 9. Plain-text summary --------------------------------------------------
core_line <- function(g, t) {
  r <- core_tab |> filter(Group == g, is_core(prevalence, t)) |> arrange(desc(prevalence))
  if (nrow(r) == 0) return(sprintf("  %-15s (n = %3d): none", g, n_of[[g]]))
  sprintf("  %-15s (n = %3d): %d genera - %s", g, n_of[[g]], nrow(r),
          paste(sprintf("%s %.0f%%", pretty_genus(r$Genus), r$prevalence), collapse = ", "))
}

summary_lines <- c(
  "07_core_microbiome.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Core = genus present (>= 1 read; >= 20 reads under MC) in >= 90% of mosquitoes in a group.",
  "",
  sprintf("SHARED CORE (core in all %d groups)", n_groups),
  sprintf("  at %3d%%: %d genera - %s", c(50, 70, 90),
          sapply(c(50, 70, 90), function(t) length(shared_at(t))),
          sapply(c(50, 70, 90), function(t) paste(pretty_genus(shared_at(t)), collapse = ", "))),
  "",
  "CORE GENERA PER GROUP AT 90% (prevalence)",
  unlist(lapply(group_levels, core_line, t = 90)),
  "",
  "CORE GENERA PER GROUP AT 70% (sensitivity)",
  unlist(lapply(group_levels, core_line, t = 70)),
  "",
  "CORE GENERA PER GROUP AT 50% (sensitivity)",
  unlist(lapply(group_levels, core_line, t = 50)),
  "",
  "NUMBER OF CORE GENERA BY THRESHOLD",
  core_size |>
    pivot_wider(names_from = threshold, values_from = n_core) |>
    mutate(line = paste0(sprintf("  %-15s", Group),
                         apply(across(-Group), 1, function(v) paste(sprintf("%3d", v), collapse = " ")))) |>
    pull(line),
  sprintf("  %-15s%s", "threshold (%):", paste(sprintf("%3d", curve_thresholds), collapse = " ")),
  "",
  "MEAN SHARE OF READS (90% definition): shared core / group-specific core / not core",
  read_share |>
    select(Group, Part, pct) |>
    pivot_wider(names_from = Part, values_from = pct) |>
    mutate(line = sprintf("  %-15s %5.1f%% / %5.1f%% / %5.1f%%", Group,
                          `Shared core`, `Group-specific core`, `Not core`)) |>
    pull(line),
  "",
  "CORE GENERA IN NEGATIVE CONTROLS (6 controls)",
  sprintf("  %-26s in %d/6 controls, %s reads (max %s in one control; %.1f%% of all control reads) | biological: %.0f%% prevalence, median %s reads when present",
          pretty_genus(ctrl_check$Genus), as.integer(ctrl_check$controls_with_genus),
          scales::comma(ctrl_check$control_reads_total), scales::comma(ctrl_check$control_reads_max_single),
          ctrl_check$pct_of_all_control_reads, ctrl_check$bio_prevalence_pct,
          scales::comma(ctrl_check$bio_median_reads_when_present))
)
writeLines(summary_lines, file.path(paths$stats, "07_results_summary.txt"))
cat(summary_lines, sep = "\n")
