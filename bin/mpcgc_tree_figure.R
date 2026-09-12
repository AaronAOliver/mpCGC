#!/usr/bin/env Rscript
# Render a family phylogeny with ggtree, coloured by host taxonomy or by
# sequence-similarity-network community.
#
# Tip labels in the FastTree output are protein identifiers, which carry no
# taxonomy of their own, so two joins are needed:
#
#   protein --(protein map from CGC_PROTEINS)--> genome
#   genome  --(metadata / samplesheet taxonomy)--> lineage
#
# --color-by taxonomy uses the second join and falls back to community colouring
# when no taxonomy is available. --color-by community skips it entirely.
#
# Dependencies are ggplot2 and ggtree only, so a stock bioconductor-ggtree
# container runs this without additions.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
})

# Minimal "--flag value" parsing, so the only dependencies are ggplot2 and
# ggtree and a stock bioconductor-ggtree container can run this unchanged.
args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- which(args == paste0("--", flag))
  if (length(hit) && hit[1] < length(args)) {
    nxt <- args[hit[1] + 1]
    if (!startsWith(nxt, "--")) return(nxt)
  }
  default
}

opt <- list(
  tree        = arg_value("tree"),
  annotation  = arg_value("annotation", ""),
  protein_map = arg_value("protein-map", ""),
  taxonomy    = arg_value("taxonomy", ""),
  color_by    = arg_value("color-by", "taxonomy"),
  rank        = arg_value("rank", "phylum"),
  family      = arg_value("family", ""),
  out_prefix  = arg_value("out-prefix"),
  width       = as.numeric(arg_value("width", "8")),
  height      = as.numeric(arg_value("height", "10")),
  max_groups  = as.integer(arg_value("max-groups", "12"))
)

usage <- paste(
  "Usage: mpcgc_tree_figure.R --tree <nwk> --out-prefix <prefix>",
  "         [--annotation <communities.tsv>] [--protein-map <map.tsv>]",
  "         [--taxonomy <metadata.tsv>] [--color-by taxonomy|community]",
  "         [--rank phylum|class|order|family|genus] [--family <name>]",
  "         [--width N] [--height N] [--max-groups N]", sep = "\n")

if (any(args %in% c("--help", "-h"))) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}
for (req in c("tree", "out_prefix")) {
  if (is.null(opt[[req]]) || !nzchar(opt[[req]])) {
    cat(usage, "\n")
    stop(sprintf("--%s is required", gsub("_", "-", req)), call. = FALSE)
  }
}
if (!file.exists(opt$tree)) {
  stop(sprintf("tree file not found: %s", opt$tree), call. = FALSE)
}
if (!opt$color_by %in% c("taxonomy", "community")) {
  stop("--color-by must be 'taxonomy' or 'community'", call. = FALSE)
}
if (is.na(opt$max_groups) || opt$max_groups < 1) opt$max_groups <- 12

tree <- read.tree(opt$tree)
tips <- tree$tip.label

# ---- the mpCGC 11-group palette, so tree figures match the other figures ----
GROUP_COLORS <- c(
  "Actinobacteriota"    = "#1f77b4",
  "Alphaproteobacteria" = "#ff7f0e",
  "Archaea"             = "#2ca02c",
  "Bacteroidota"        = "#d62728",
  "Chloroflexota"       = "#9467bd",
  "Cyanobacteria"       = "#8c564b",
  "Desulfobacterota"    = "#e377c2",
  "Firmicutes"          = "#7f7f7f",
  "Gammaproteobacteria" = "#bcbd22",
  "Planctomycetota"     = "#17becf",
  "Verrucomicrobiota"   = "#1a3f70",
  "other"               = "#c9c9c9"
)
FALLBACK <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#1a3f70", "#a05d56")

PHYLA <- c("Actinobacteriota", "Bacteroidota", "Chloroflexota", "Cyanobacteria",
           "Desulfobacterota", "Firmicutes", "Planctomycetota",
           "Verrucomicrobiota")
CLASSES <- c("Alphaproteobacteria", "Gammaproteobacteria")
RANK_PREFIX <- c(domain = "d__", phylum = "p__", class = "c__", order = "o__",
                 family = "f__", genus = "g__", species = "s__")

# Collapse a GTDB string onto the mpCGC lineage groups: the eight phyla that
# dominate the catalog, Proteobacteria split into Alpha and Gamma because the
# phylum is too coarse to be informative, Archaea pooled, everything else grey.
lineage_group <- function(tax) {
  parts <- trimws(strsplit(tax, ";", fixed = TRUE)[[1]])
  parts <- parts[nchar(parts) > 3]
  if (!length(parts)) return("other")
  names_only <- substring(parts, 4)
  if (startsWith(names_only[1], "Archaea")) return("Archaea")
  if (length(names_only) >= 2) {
    ph <- sub("_.*$", "", names_only[2])
    if (ph %in% PHYLA) return(ph)
  }
  if (length(names_only) >= 3 && names_only[3] %in% CLASSES) return(names_only[3])
  "other"
}

# Pull one rank verbatim, for --rank order/family/genus etc.
rank_value <- function(tax, rank) {
  pre <- RANK_PREFIX[[rank]]
  if (is.null(pre)) return("other")
  parts <- trimws(strsplit(tax, ";", fixed = TRUE)[[1]])
  hit <- parts[startsWith(parts, pre)]
  if (!length(hit) || nchar(hit[1]) <= 3) return("other")
  substring(hit[1], 4)
}

read_tsv_safe <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  df <- try(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                       quote = ""), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) return(NULL)
  df
}

# ---- resolve a group label per tip -----------------------------------------
group_of <- setNames(rep(NA_character_, length(tips)), tips)
legend_title <- "SSN community"
used <- "community"

if (identical(opt$color_by, "taxonomy")) {
  pmap <- read_tsv_safe(opt$protein_map)
  meta <- read_tsv_safe(opt$taxonomy)

  if (!is.null(pmap) && !is.null(meta) &&
      all(c("protein", "mag") %in% names(pmap))) {

    id_col <- intersect(c("Bin_id", "sample", "MAG", "genome"), names(meta))
    tax_col <- grep("taxonomy", names(meta), ignore.case = TRUE, value = TRUE)

    if (length(id_col) && length(tax_col)) {
      tax_of <- setNames(meta[[tax_col[1]]], meta[[id_col[1]]])
      mag_of <- setNames(pmap$mag, pmap$protein)
      tax_for_tip <- tax_of[mag_of[tips]]

      labels <- vapply(tax_for_tip, function(t) {
        if (is.na(t) || !nzchar(t)) return("other")
        if (identical(opt$rank, "phylum")) lineage_group(t) else rank_value(t, opt$rank)
      }, character(1), USE.NAMES = FALSE)

      if (any(labels != "other")) {
        group_of <- setNames(labels, tips)
        legend_title <- if (identical(opt$rank, "phylum")) "Lineage" else
          paste0(toupper(substring(opt$rank, 1, 1)), substring(opt$rank, 2))
        used <- "taxonomy"
      } else {
        message("[ggtree] no tip could be assigned a lineage; ",
                "falling back to community colouring")
      }
    } else {
      message("[ggtree] metadata lacks a genome column or a taxonomy column; ",
              "falling back to community colouring")
    }
  } else {
    message("[ggtree] no protein map or taxonomy table given; ",
            "falling back to community colouring")
  }
}

if (identical(used, "community")) {
  ann <- read_tsv_safe(opt$annotation)
  if (!is.null(ann) && all(c("protein", "leiden_community") %in% names(ann))) {
    comm <- setNames(paste0("C", ann$leiden_community), ann$protein)
    group_of <- setNames(comm[tips], tips)
  }
  group_of[is.na(group_of)] <- "other"
}

group_of[is.na(group_of) | !nzchar(group_of)] <- "other"

# keep the largest groups, pool the rest, so the legend stays readable
sizes <- sort(table(group_of), decreasing = TRUE)
keep <- setdiff(names(sizes), "other")
keep <- keep[seq_len(min(length(keep), opt$max_groups))]
labels <- ifelse(group_of %in% keep, group_of, "other")

levels_in_order <- c(
  intersect(names(GROUP_COLORS), keep),
  sort(setdiff(keep, names(GROUP_COLORS))),
  if (any(labels == "other")) "other" else NULL
)
meta_df <- data.frame(label = tips,
                      group = factor(labels, levels = levels_in_order),
                      stringsAsFactors = FALSE)

cols <- character(0)
spare <- FALLBACK
for (lv in levels_in_order) {
  if (lv %in% names(GROUP_COLORS)) {
    cols[lv] <- GROUP_COLORS[[lv]]
  } else if (identical(lv, "other")) {
    cols[lv] <- "#c9c9c9"
  } else {
    cols[lv] <- spare[1]
    spare <- c(spare[-1], spare[1])
  }
}

# ---- draw -------------------------------------------------------------------
show_tips <- length(tips) <= 120

p <- ggtree(tree, layout = "rectangular", linewidth = 0.3) %<+% meta_df +
  geom_tippoint(aes(colour = group), size = 1.7, na.rm = TRUE) +
  scale_colour_manual(values = cols, na.value = "#c9c9c9",
                      name = legend_title, drop = FALSE) +
  theme_tree2() +
  theme(legend.position = "right",
        legend.key.size = unit(0.9, "lines"),
        plot.title = element_text(size = 13, face = "bold"))

# tip labels are drawn past the last node, so the panel needs room on the right
# or long protein identifiers are clipped at the edge
p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.28)))

if (nzchar(opt$family)) {
  p <- p + ggtitle(sprintf("%s  (%d sequences)", opt$family, length(tips)))
}
if (show_tips) {
  p <- p + geom_tiplab(size = 1.8, offset = 0.02)
}

height <- max(opt$height, length(tips) * 0.035)

ggsave(paste0(opt$out_prefix, ".png"), p, width = opt$width,
       height = height, dpi = 300, limitsize = FALSE)

# ggplot2 >= 4 needs svglite for SVG and no longer falls back to cairo, so try
# svglite, then the cairo device, and settle for PNG alone rather than failing.
svg_path <- paste0(opt$out_prefix, ".svg")
svg_written <- FALSE

if (requireNamespace("svglite", quietly = TRUE)) {
  svg_written <- tryCatch({
    ggsave(svg_path, p, width = opt$width, height = height, limitsize = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

if (!svg_written && isTRUE(capabilities("cairo"))) {
  svg_written <- tryCatch({
    ggsave(svg_path, p, width = opt$width, height = height, limitsize = FALSE,
           device = grDevices::svg)
    TRUE
  }, error = function(e) FALSE)
}

if (!svg_written) {
  if (file.exists(svg_path)) file.remove(svg_path)
  message("[ggtree] SVG not written: neither svglite nor a working cairo SVG ",
          "device is available. The PNG was still produced.")
}

cat(sprintf("[ggtree] %s: %d tips, coloured by %s, %d groups, %s\n",
            opt$family, length(tips), used, length(levels_in_order),
            if (svg_written) "png + svg" else "png only"))
