# =============================================================================
# 11a_build_host_species.R - tag COI host-species IDs onto the wild 16S samples
#
# Wild 16S plates: AG_P01, AG_P02, ZX_P01, ZX_P02 (380 wild midguts).
# COI identifications (one per specimen):
#   Agabus plate 1 : CO1_P1_Gbif-Results.xlsx (or .csv)  - GBIF/BOLD BLAST
#   Agabus plate 2 : CO1_P2_Gbif-Results.xlsx (or .csv)  - OPTIONAL, used if present
#   Zi Xin         : ZiXin_Taxa_classified_COI.xlsx, Sheet1 (Sheet2 is a subset)
# 16S sample <-> specimen: demux_metadata_audit.csv (from Atlas9 01_demux_metadata/)
#
# Match key = specimen name in upper case with all non-alphanumerics removed
# (COI files write e.g. "COI07km3-13"; 16S writes "07KM3-13").
#
# Put the input files here:
#   data/raw/species/CO1_P1_Gbif-Results.xlsx
#   data/raw/species/ZiXin_Taxa_classified_COI.xlsx
#   data/raw/metadata/demux_metadata_audit.csv
# Output:
#   data/raw/species/wild_host_species_COI.tsv   (one row per wild 16S sample)
# When Aggy's plate 2 file arrives: save it as data/raw/species/CO1_P2_Gbif-Results.xlsx
# and re-run this script.
# =============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(readxl))

sp_dir  <- paths$species
f_audit <- c(here("data", "raw", "metadata", "demux_metadata_audit.csv"),
             file.path(paths$species, "demux_metadata_audit.csv"))
f_audit <- f_audit[file.exists(f_audit)][1]
if (is.na(f_audit)) f_audit <- "missing"
find_file <- function(base) {
  f <- file.path(sp_dir, paste0(base, c(".xlsx", ".csv")))
  f[file.exists(f)][1]
}
read_any <- function(f) if (grepl("\\.csv$", f)) read_csv(f, show_col_types = FALSE) else read_excel(f)
norm <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))

f_ag1 <- find_file("CO1_P1_Gbif-Results")
f_ag2 <- find_file("CO1_P2_Gbif-Results")          # NA until Aggy sends it
f_zx  <- file.path(sp_dir, "ZiXin_Taxa_classified_COI.xlsx")
stopifnot("demux_metadata_audit.csv not found in data/raw/metadata/" = file.exists(f_audit),
          "CO1_P1_Gbif-Results not found in data/raw/species/"       = !is.na(f_ag1),
          "ZiXin_Taxa_classified_COI.xlsx not found in data/raw/species/" = file.exists(f_zx))

# ---- 1. The 380 wild 16S samples --------------------------------------------
wild <- read_csv(f_audit, show_col_types = FALSE) |>
  filter(Plate_ID %in% c("AG_P01_16S", "AG_P02_16S", "ZX_P01_16S", "ZX_P02_16S"),
         !str_detect(Specimen, regex("neg", ignore_case = TRUE))) |>
  transmute(Sample_ID = paste0("sample_", SampleID), Specimen,
            Plate = str_remove(Plate_ID, "_16S"), Well, key = norm(Specimen))
stopifnot(nrow(wild) == 380, !anyDuplicated(wild$key))

# ---- 2. COI identifications -------------------------------------------------
last_n <- function(x, n) vapply(str_split(x, "_"), function(p) p[length(p) - n], character(1))
read_gbif <- function(f, label) {
  d <- read_any(f)
  tibble(key          = norm(sub("^COI", "", sub("_all.*$", "", d$occurrenceId))),
         COI_source   = label,
         COI_call_raw = last_n(d$classification, 0),
         COI_rank     = "SPECIES",
         COI_support  = sprintf("%.2f%% identity", as.numeric(d$identity)),
         Genus_raw    = last_n(d$classification, 1))
}
zx <- read_excel(f_zx, sheet = "Sheet1")
zx_tab <- tibble(key          = norm(zx$SpecimenID),
                 COI_source   = "Zi Xin (tree placement)",
                 COI_call_raw = coalesce(as.character(zx$Species), as.character(zx$Genus)),
                 COI_rank     = as.character(zx$`Tax Rank`),
                 COI_support  = ifelse(is.na(zx$Confidence), "",
                                       sprintf("%.0f%% confidence", as.numeric(zx$Confidence))),
                 Genus_raw    = as.character(zx$Genus))

coi <- bind_rows(read_gbif(f_ag1, "Agabus P1 (GBIF/BOLD BLAST)"),
                 if (!is.na(f_ag2)) read_gbif(f_ag2, "Agabus P2 (GBIF/BOLD BLAST)"),
                 zx_tab)
stopifnot("a specimen appears in more than one COI file" = !anyDuplicated(coi$key))

# ---- 3. Join + harmonise names ----------------------------------------------
host <- wild |>
  left_join(coi, by = "key") |>
  mutate(Host_group = case_when(
    Plate == "AG_P02" & is.na(COI_call_raw)                  ~ "Pending (AG_P02 COI not yet available)",
    is.na(COI_call_raw)                                      ~ "Unidentified",
    COI_call_raw == "Aedes albopictus"                       ~ "Ae. albopictus",
    COI_call_raw == "Aedes aegypti"                          ~ "Ae. aegypti",
    Genus_raw == "Culex"                                     ~ "Culex spp.",
    Genus_raw == "Mansonia" & str_detect(COI_call_raw, "indiana") ~ "Mansonia indiana",
    Genus_raw == "Mansonia"                                  ~ "Mansonia spp.",
    Genus_raw == "Anopheles"                                 ~ "Anopheles spp.",
    TRUE                                                     ~ paste("Other:", COI_call_raw))) |>
  select(Sample_ID, Specimen, Plate, Well, Host_group, COI_call_raw, COI_rank, COI_support, COI_source) |>
  arrange(Plate, Sample_ID)

out <- file.path(sp_dir, "wild_host_species_COI.tsv")
write_tsv(host, out)

# ---- 4. Report --------------------------------------------------------------
cat("\nHost species of the 380 wild 16S samples:\n")
print(count(host, Host_group, sort = TRUE), n = Inf)
cat("\nBy plate:\n")
print(count(host, Plate, Host_group), n = Inf)
cat("\nCulex detail:\n")
print(count(filter(host, Host_group == "Culex spp."), COI_call_raw, sort = TRUE))
unmatched <- setdiff(coi$key, wild$key)
cat("\nCOI records with no wild 16S sample (other projects / not sequenced):", length(unmatched), "\n")
if (is.na(f_ag2)) {
  expected <- c("Ae. albopictus" = 219, "Mansonia indiana" = 36, "Culex spp." = 23,
                "Pending (AG_P02 COI not yet available)" = 95, "Unidentified" = 4)
  got <- table(host$Host_group)[names(expected)]
  if (!isTRUE(all(got == expected))) warning("Counts differ from the 2026-10-08 build - check the input files")
}
message("Saved ", out)
