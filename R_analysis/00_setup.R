# =============================================================================
# 00_setup.R - shared setup for every analysis script
# Usage at the top of each script:  source(here::here("R", "00_setup.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(patchwork)
  library(here)
})

set.seed(4288)   # same permutations every run

# ---- Paths ------------------------------------------------------------------
paths <- list(
  plotting  = here("data", "raw", "plotting_final"),
  metadata  = here("data", "raw", "metadata", "sample_metadata_current.tsv"),
  kraken    = here("data", "raw", "kraken_rescue_summary"),
  species   = here("data", "raw", "species"),
  processed = here("data", "processed"),
  figures   = here("output", "figures"),
  tables    = here("output", "tables"),
  stats     = here("output", "stats")
)

# ---- Metadata ---------------------------------------------------------------
# Origin is standardised to Wild / Lab / Control.
# Lab_Taxon comes from the sample ID (BT_ALB / BT_CLX). Wild host species
# come from COI barcodes (11a_build_host_species.R) - never inferred from 16S.
read_metadata <- function() {
  read_tsv(paths$metadata, show_col_types = FALSE) |>
    mutate(
      Origin = case_when(
        str_detect(str_to_lower(Origin), "^wild") ~ "Wild",
        str_detect(str_to_lower(Origin), "^lab")  ~ "Lab",
        TRUE                                      ~ "Control"
      ),
      Origin = factor(Origin, levels = c("Wild", "Lab", "Control")),
      Lab_Taxon = case_when(
        Origin == "Lab" & str_detect(Sample_ID, "BT_ALB") ~ "Ae. albopictus",
        Origin == "Lab" & str_detect(Sample_ID, "BT_CLX") ~ "Culex",
        TRUE ~ NA_character_
      ),
      Location = if_else(Include_in_Location_Analysis == "YES",
                         Location_Name, NA_character_),
      Location = factor(Location, levels = c("Mandai Zoo", "Kranji Marshes",
                                             "Bird Paradise", "River Safari"))
    )
}

# ---- OTU -> genus map -------------------------------------------------------
read_otu_map <- function() {
  read_tsv(file.path(paths$plotting, "MC_GG97_328_OTU_to_genus.tsv"),
           show_col_types = FALSE) |>
    rename(OTU = Feature_ID, Taxon = GG2_97_Taxon, Genus = Genus_for_plotting)
}

# ---- OTU tables -------------------------------------------------------------
# Returns a numeric matrix: rows = samples, columns = the 328 OTUs.
# which: "biological" (570), "wild" (380), "lab" (190), "location" (345),
#        "controls" (6)
read_otu_matrix <- function(which = "biological") {
  file <- switch(which,
    biological = "MC_GG97_biological_570_samples_x_328_OTUs.tsv",
    wild       = "MC_GG97_wild_380_samples_x_328_OTUs.tsv",
    lab        = "MC_GG97_lab_190_samples_x_328_OTUs.tsv",
    location   = "MC_GG97_location_345_samples_x_328_OTUs.tsv",
    controls   = "MC_GG97_controls_6_samples_x_328_OTUs.tsv",
    stop("unknown table: ", which))
  otus <- read_otu_map()$OTU
  tab  <- read_tsv(file.path(paths$plotting, file), show_col_types = FALSE)
  ids  <- tab[[1]]
  m    <- as.matrix(tab[, otus])
  rownames(m) <- ids
  storage.mode(m) <- "double"
  m
}

# ---- Sanity check: stop early if the inputs ever change ---------------------
check_inputs <- function() {
  meta <- read_metadata()
  bio  <- read_otu_matrix("biological")
  stopifnot(
    nrow(meta) == 576,
    sum(meta$Origin == "Wild") == 380,
    sum(meta$Origin == "Lab") == 190,
    sum(meta$Origin == "Control") == 6,
    sum(!is.na(meta$Location)) == 345,
    nrow(bio) == 570, ncol(bio) == 328,
    sum(bio) == 1422797,   # after script 53 (BT_P02 control = BT_ALB_192)
    all(rownames(bio) %in% meta$Sample_ID[meta$Origin != "Control"])
  )
  message("Inputs OK: 570 biological samples x 328 OTUs, 1,422,797 reads.")
  invisible(TRUE)
}

# ---- Genus colours (fixed for EVERY figure) ---------------------------------
# Top 8 genera by overall mean relative abundance get the 8 validated
# colour-blind-safe hues, in this order (adjacent pairs pass CVD checks).
# Everything else -> "Other"; GG2 OTUs without a genus -> "Unclassified".
# Colour follows the genus, never its rank in a given plot.
genus_colours <- c(
  "Elizabethkingia"           = "#2a78d6",  # blue
  "Wolbachia"                 = "#eb6834",  # orange
  "Stenotrophomonas_A_615274" = "#1baf7a",  # aqua
  "Sphingobacterium"          = "#eda100",  # yellow
  "Pseudomonas_E_650326"      = "#e87ba4",  # magenta
  "Sphingomonas_L_486704"     = "#008300",  # green
  "Rhodococcus_C_375578"      = "#4a3aa7",  # violet
  "Comamonas_F_589250"        = "#e34948",  # red
  "Unclassified_at_genus"     = "#9b9a95",  # mid grey
  "Other"                     = "#d9d8d2"   # light grey
)
focal_genera <- names(genus_colours)[1:8]

# Collapse any genus label to the plotting categories above.
lump_genus <- function(genus) {
  case_when(
    genus %in% focal_genera            ~ genus,
    genus == "Unclassified_at_genus"   ~ "Unclassified_at_genus",
    TRUE                               ~ "Other"
  ) |> factor(levels = names(genus_colours))
}

# Display names: drop GG2 numeric suffixes ("Stenotrophomonas_A_615274" ->
# "Stenotrophomonas A"). Note the GG2 suffix convention in the Methods.
pretty_genus_base <- function(x) {
  x <- as.character(x)
  x[x == "Unclassified_at_genus"] <- "Unclassified at genus"
  x <- str_remove(x, "_\\d+$")
  str_replace_all(x, "_", " ")
}
# GG2 has two different genera each called "Pantoea A" and "Pseudomonas E"
# (different numeric IDs). Where names would clash, keep the ID:
# "Pseudomonas_E_650326" -> "Pseudomonas E (650326)".
genus_clashes <- local({
  l <- pretty_genus_base(unique(read_otu_map()$Genus))
  unique(l[duplicated(l)])
})
pretty_genus <- function(x) {
  l <- pretty_genus_base(x)
  i <- l %in% genus_clashes
  l[i] <- paste0(l[i], " (", str_extract(as.character(x)[i], "\\d+$"), ")")
  l
}
genus_labels <- setNames(pretty_genus(names(genus_colours)), names(genus_colours))

# ---- Group colours ----------------------------------------------------------
origin_colours   <- c(Wild = "#1baf7a", Lab = "#4a3aa7")
location_colours <- c("Mandai Zoo"     = "#2a78d6", "Kranji Marshes" = "#eb6834",
                      "Bird Paradise"  = "#008300", "River Safari"   = "#e87ba4",
                      "Lab"            = "#9b9a95")

# ---- Figure theme + saving --------------------------------------------------
theme_fyp <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.title     = element_text(face = "bold"),
      plot.tag         = element_text(face = "bold")
    )
}
theme_set(theme_fyp())

# Saves PDF (for the report) + 600 dpi PNG (for slides) in output/figures/
save_fig <- function(plot, name, width = 180, height = 120, units = "mm") {
  dir.create(paths$figures, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(paths$figures, paste0(name, ".pdf")), plot,
         width = width, height = height, units = units)
  ggsave(file.path(paths$figures, paste0(name, ".png")), plot,
         width = width, height = height, units = units, dpi = 600, bg = "white")
  message("Saved output/figures/", name, ".pdf/.png")
}

# ---- All-genera heatmap (panel A of Figs 2, 4, 5; FigS_all_genera) ----------
# genus_counts():     samples x OTUs  ->  samples x genus labels (68)
# genus_row_order():  fixed row order for every heatmap = mean relative
#                     abundance across all 570 mosquitoes, Unclassified last
# genus_mean_table(): mean relative abundance (%) + prevalence (%) per genus
#                     per group, means of within-mosquito proportions (exactly
#                     as in the stacked bars)
# genus_heatmap():    log10 heatmap of mean relative abundance; white = not
#                     detected in that group; genera absent from every column
#                     of the figure are dropped
genus_counts <- function(m) {
  map <- read_otu_map()
  g <- map$Genus[match(colnames(m), map$OTU)]
  stopifnot(!anyNA(g))
  t(rowsum(t(m), group = g))
}

genus_row_order <- function() {
  g <- genus_counts(read_otu_matrix("biological"))
  o <- names(sort(colMeans(g / rowSums(g)), decreasing = TRUE))
  c(setdiff(o, "Unclassified_at_genus"), "Unclassified_at_genus")
}

genus_mean_table <- function(counts, groups, facet = NULL) {
  g <- genus_counts(counts)
  r <- g / rowSums(g)
  stopifnot(length(groups) == nrow(g), all(rowSums(g) > 0))
  tibble(Sample_ID = rownames(g), Group = groups,
         Facet = if (is.null(facet)) "all" else facet) |>
    group_by(Facet, Group) |>
    group_modify(function(d, k) {
      tibble(Genus      = colnames(g),
             n          = nrow(d),
             mean_rel   = unname(100 * colMeans(r[d$Sample_ID, , drop = FALSE])),
             prevalence = unname(100 * colMeans(g[d$Sample_ID, , drop = FALSE] > 0)))
    }) |>
    ungroup()
}

genus_heatmap <- function(tab, x_labels = waiver(), row_order = genus_row_order(),
                          drop_absent = TRUE,
                          breaks = c(1e-3, 1e-2, 0.1, 1, 10, 100),
                          legend_width = 40) {
  if (drop_absent) {
    keep <- tab |> group_by(Genus) |> summarise(any = any(mean_rel > 0), .groups = "drop") |>
      filter(any) |> pull(Genus)
    tab <- filter(tab, Genus %in% keep)
  }
  lev <- rev(intersect(row_order, unique(tab$Genus)))
  tab <- tab |> mutate(Genus = factor(Genus, levels = lev),
                       mean_rel = if_else(mean_rel > 0, mean_rel, NA_real_))
  g_expr <- setNames(as.expression(lapply(lev, function(g) {
    if (g == "Unclassified_at_genus") pretty_genus(g) else bquote(italic(.(pretty_genus(g))))
  })), lev)
  p <- ggplot(tab, aes(x = Group, y = Genus, fill = mean_rel)) +
    geom_tile(colour = "grey88", linewidth = 0.25) +
    scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10",
                         limits = c(1e-4, 100), breaks = breaks,
                         labels = function(b) format(b, scientific = FALSE, drop0trailing = TRUE, trim = TRUE),
                         oob = scales::squish, na.value = "white",
                         name = "Mean relative\nabundance (%)") +
    scale_x_discrete(labels = x_labels, expand = c(0, 0)) +
    scale_y_discrete(labels = function(b) g_expr[b], expand = c(0, 0)) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0,
                                  barwidth = unit(legend_width, "mm"),
                                  barheight = unit(2.5, "mm"))) +
    labs(x = NULL, y = NULL) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6.3),
          axis.text.x = element_text(size = 7.5),
          panel.spacing.x = unit(1.5, "mm"),
          strip.text = element_text(size = 8.5),
          legend.position = "top", legend.justification = "left",
          legend.title = element_text(size = 7.5), legend.text = element_text(size = 6.5))
  if (!all(tab$Facet == "all")) p <- p + facet_grid(~ Facet, scales = "free_x", space = "free_x")
  p
}

save_table <- function(df, name, dir = paths$tables) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  write_tsv(df, file.path(dir, paste0(name, ".tsv")))
  message("Saved ", file.path(basename(dir), paste0(name, ".tsv")))
}

# ---- Self-check: the duplicate-name fix must be active ----------------------
stopifnot(
  "pretty_genus() is not the updated version - replace 00_setup.R" =
    identical(pretty_genus(c("Pseudomonas_E_650326", "Pantoea_A_680069")),
              c("Pseudomonas E (650326)", "Pantoea A (680069)"))
)
