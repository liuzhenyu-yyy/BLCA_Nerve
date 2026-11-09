# TCGA-BLCA analysis and visualization
# ===================== 00_common.R =====================
# Shared plotting settings. This script does not select values or PNI labels by significance.
PNI_COLORS <- c("Non-PNI" = "#4A77A1", "PNI" = "#BD819B")
INK <- "#243447"
MUTED <- "#667085"
FONT_FAMILY <- "Arial"

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < .001, formatC(p, digits = 2, format = "e"), sprintf("%.4f", p)))
}

theme_tcga <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = INK, size = 18),
      plot.subtitle = ggplot2::element_text(colour = MUTED, size = 11),
      axis.text = ggplot2::element_text(colour = INK), axis.title = ggplot2::element_text(colour = INK),
      panel.grid.major.y = ggplot2::element_line(colour = "#E4E7EC", linewidth = .25),
      panel.grid.minor = ggplot2::element_blank(), legend.position = "none",
      plot.margin = ggplot2::margin(9, 10, 7, 8)
    )
}

# Some slide renderers omit R vector lines with an exactly zero bounding dimension.
# Set only that invisible dimension to one EMU; preserve endpoints and statistical data.
repair_pptx_lines <- function(pptx_path) {
  tmp <- tempfile("tcga_pptx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(pptx_path, exdir = tmp)
  slide_dir <- file.path(tmp, "ppt", "slides")
  files <- list.files(slide_dir, pattern = "^slide[0-9]+[.]xml$", full.names = TRUE)
  fixed <- 0L
  for (file in files) {
    doc <- xml2::read_xml(file)
    shapes <- xml2::xml_find_all(doc, ".//p:sp", ns)
    changed <- FALSE
    for (shape in shapes) {
      ext <- xml2::xml_find_first(shape, "p:spPr/a:xfrm/a:ext", ns)
      line <- xml2::xml_find_first(shape, "p:spPr/a:ln", ns)
      paths <- xml2::xml_find_all(shape, "p:spPr/a:custGeom/a:pathLst/a:path", ns)
      if (inherits(ext, "xml_missing") || inherits(line, "xml_missing") || !length(paths)) next
      dims <- c("cx", "cy")[vapply(c("cx", "cy"), function(k) identical(xml2::xml_attr(ext, k), "0"), logical(1))]
      if (length(dims) != 1L) next
      dimension <- dims[[1]]
      coord <- if (dimension == "cx") "x" else "y"
      pathdim <- if (dimension == "cx") "w" else "h"
      valid <- all(vapply(paths, function(p) {
        points <- xml2::xml_find_all(p, ".//a:pt", ns)
        identical(xml2::xml_attr(p, pathdim), "0") && all(xml2::xml_attr(points, coord) == "0")
      }, logical(1)))
      if (!valid) stop("Unexpected zero-extent geometry in PPTX; original file preserved.")
      for (p in paths) xml2::xml_set_attr(p, pathdim, "1")
      xml2::xml_set_attr(ext, dimension, "1")
      fixed <- fixed + 1L
      changed <- TRUE
    }
    if (changed) xml2::write_xml(doc, file)
  }
  if (fixed) {
    patched <- tempfile(fileext = ".pptx")
    on.exit(unlink(patched), add = TRUE)
    zip::zipr(patched, files = list.files(tmp, all.files = TRUE, no.. = TRUE), root = tmp, include_directories = TRUE)
    listing <- utils::unzip(patched, list = TRUE)$Name
    stopifnot("[Content_Types].xml" %in% listing, "ppt/presentation.xml" %in% listing)
    if (!file.copy(patched, pptx_path, overwrite = TRUE)) stop("Could not save the PPTX compatibility patch.")
  }
  invisible(fixed)
}


# ===================== 01_molecular.R =====================
# TCGA molecular comparisons and AJCC overall-stage visualization.
# This module reads no files and writes no files. Molecular scores are supplied
# in dat; it does NOT recompute xCell, ssGSEA, or gene expression normalization.
# Shared settings are defined earlier. The 14 molecular comparisons use two-sided P values and BH q values.

make_molecular_figures <- function(dat, config, metric_specs) {
  keys <- c(
    "CD4_T", "CD8_T", "Neutrophils", "Tregs", "Exhaustion",
    "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX",
    "N1_like", "N2_like", "PMN_MDSC_like"
  )
  required <- c("sample_id", "PNI", "ajcc_stage", keys)
  stopifnot(
    all(required %in% names(dat)), !anyDuplicated(dat$sample_id),
    all(dat$PNI %in% c("Non-PNI", "PNI")),
    nrow(dat) == config$cohort_n,
    all(c("metric_key", "metric_label_en", "category", "unit") %in% names(metric_specs))
  )
  meta <- metric_specs[match(keys, metric_specs$metric_key), , drop = FALSE]
  stopifnot(!anyNA(meta$metric_key), !anyDuplicated(meta$metric_key))
  # Legacy directional metadata is intentionally ignored.
  meta$alternative <- "two.sided"
  meta$hypothesis <- "PNI versus Non-PNI: two-sided comparison"
  meta$source_column <- meta$metric_key
  dat$PNI <- factor(dat$PNI, levels = c("Non-PNI", "PNI"))
  long <- tidyr::pivot_longer(
    dat[, c("sample_id", "PNI", keys)],
    cols = dplyr::all_of(keys),
    names_to = "metric_key", values_to = "value"
  ) |>
    dplyr::left_join(meta, by = "metric_key")
  available <- dplyr::filter(long, is.finite(value))
  descriptive <- available |>
    dplyr::group_by(metric_key, PNI) |>
    dplyr::summarise(
      n = dplyr::n(), median = stats::median(value),
      q1 = stats::quantile(value, .25), q3 = stats::quantile(value, .75),
      .groups = "drop"
    )
  tests <- dplyr::bind_rows(lapply(seq_len(nrow(meta)), function(i) {
    m <- meta[i, , drop = FALSE]
    d <- available[available$metric_key == m$metric_key, , drop = FALSE]
    x <- d$value[d$PNI == "PNI"]
    y <- d$value[d$PNI == "Non-PNI"]
    if (!length(x) || !length(y)) stop("Both PNI groups need data for ", m$metric_key)
    two <- stats::wilcox.test(x, y,
      alternative = "two.sided",
      exact = FALSE, correct = TRUE
    )
    data.frame(
      metric_key = m$metric_key, source_column = m$metric_key,
      metric_label_en = m$metric_label_en, category = m$category, unit = m$unit,
      alternative = m$alternative, hypothesis = m$hypothesis,
      n_minus = length(y), n_plus = length(x),
      median_minus = stats::median(y), median_plus = stats::median(x),
      q1_minus = unname(stats::quantile(y, .25)),
      q3_minus = unname(stats::quantile(y, .75)),
      q1_plus = unname(stats::quantile(x, .25)),
      q3_plus = unname(stats::quantile(x, .75)),
      median_difference_plus_minus = stats::median(x) - stats::median(y),
      rank_biserial_plus_higher = 2 * unname(two$statistic) / (length(x) * length(y)) - 1,
      p_two_sided = two$p.value,
      method = "Wilcoxon rank-sum; asymptotic; continuity corrected",
      analysis_status = "Exploratory", row.names = NULL
    )
  }))
  tests$q_bh <- stats::p.adjust(tests$p_two_sided, method = "BH", n = nrow(tests))
  tests$multiplicity_family <- "14 molecular measures"
  rna_n <- sum(rowSums(is.finite(as.matrix(dat[, keys]))) > 0)
  boxplot_for <- function(plot_keys, title, ncol) {
    idx <- match(plot_keys, tests$metric_key)
    s <- tests[idx, , drop = FALSE]
    panel_labels <- paste0(
      s$metric_label_en, "\nTwo-sided P = ", fmt_p(s$p_two_sided),
      "; BH q = ", fmt_p(s$q_bh)
    )
    equal_n <- length(unique(s$n_minus)) == 1L && length(unique(s$n_plus)) == 1L
    if (!equal_n) {
      panel_labels <- paste0(
        panel_labels, "\nn (Non-PNI / PNI): ",
        s$n_minus, " / ", s$n_plus
      )
    }
    d <- available[available$metric_key %in% plot_keys, , drop = FALSE]
    d$panel <- factor(panel_labels[match(d$metric_key, plot_keys)], levels = panel_labels)
    x_labels <- if (equal_n) {
      c(
        "Non-PNI" = sprintf("Non-PNI\n(n = %d)", s$n_minus[1]),
        "PNI" = sprintf("PNI\n(n = %d)", s$n_plus[1])
      )
    } else {
      c("Non-PNI" = "Non-PNI", "PNI" = "PNI")
    }
    y_label <- if (all(plot_keys %in% keys[1:4])) {
      "xCell enrichment score"
    } else if (all(plot_keys %in% keys[6:11])) {
      "log2(TPM + 1)"
    } else {
      "ssGSEA score"
    }
    ggplot2::ggplot(d, ggplot2::aes(PNI, value, fill = PNI, colour = PNI)) +
      ggplot2::geom_boxplot(width = .54, alpha = .42, outlier.shape = NA, linewidth = .6) +
      ggplot2::geom_point(position = ggplot2::position_jitter(
        width = .14, height = 0,
        seed = 20261101
      ), alpha = .4, size = .9, stroke = 0) +
      ggplot2::scale_fill_manual(values = PNI_COLORS) +
      ggplot2::scale_colour_manual(values = PNI_COLORS) +
      ggplot2::scale_x_discrete(labels = x_labels) +
      ggplot2::facet_wrap(~panel, ncol = ncol, scales = "free_y") +
      ggplot2::labs(
        title = title,
        subtitle = sprintf(
          "Fixed %d-patient cohort, RNA n = %d. Two-sided exploratory comparisons.",
          nrow(dat), rna_n
        ), x = NULL, y = y_label
      ) +
      theme_tcga() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 10.5, colour = INK)
      )
  }
  plots <- list(
    immune_cells = boxplot_for(keys[1:4], "Immune-cell enrichment by PNI status", 2),
    exhaustion = boxplot_for("Exhaustion", "Exhaustion-associated expression by PNI status", 1),
    markers = boxplot_for(keys[6:11], "Immune-inhibitory marker expression", 3),
    neu = boxplot_for(keys[12:14], "Neutrophil-related transcriptional scores", 3)
  )
  effect <- tests
  effect$label <- factor(effect$metric_label_en, levels = rev(meta$metric_label_en))
  effect$direction <- ifelse(effect$rank_biserial_plus_higher >= 0, "PNI", "Non-PNI")
  plots$effect_overview <- ggplot2::ggplot(
    effect,
    ggplot2::aes(rank_biserial_plus_higher, label)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#94A3B8", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = 0, xend = rank_biserial_plus_higher,
      yend = label
    ), colour = "#CBD5E1") +
    ggplot2::geom_point(ggplot2::aes(colour = direction), size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("P = ", fmt_p(p_two_sided), "; q = ", fmt_p(q_bh))),
      nudge_y = .23, size = 2.9, colour = MUTED
    ) +
    ggplot2::scale_colour_manual(values = PNI_COLORS) +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .25)) +
    ggplot2::labs(
      title = "Comparisons across immune measures",
      subtitle = "Two-sided P values; BH q values across 14 molecular measures.",
      x = "Rank-biserial correlation (positive = higher in PNI)", y = NULL
    ) +
    theme_tcga() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 10)
    )

  # AJCC overall pathological stage is distinct from the T category in module 02.
  stage_levels <- paste("Stage", c("I", "II", "III", "IV"))
  stage_tab <- table(dat$PNI, factor(dat$ajcc_stage, levels = stage_levels))
  nonzero_stage_tab <- stage_tab[, colSums(stage_tab) > 0, drop = FALSE]
  stage_ft <- stats::fisher.test(nonzero_stage_tab)
  stage_counts <- as.data.frame(stage_tab, stringsAsFactors = FALSE)
  names(stage_counts) <- c("PNI", "stage", "n")
  stage_counts <- stage_counts |>
    dplyr::group_by(PNI) |>
    dplyr::mutate(group_n = sum(n), proportion = n / group_n) |>
    dplyr::ungroup()
  stage_counts$PNI <- factor(stage_counts$PNI, levels = c("Non-PNI", "PNI"))
  stage_counts$stage <- factor(stage_counts$stage, levels = stage_levels)
  stage_n <- rowSums(stage_tab)
  stage_test <- data.frame(
    method = stage_ft$method, n = sum(stage_tab),
    p_value = stage_ft$p.value, alternative = "Omnibus association across stages",
    unknown_excluded = nrow(dat) - sum(stage_tab)
  )
  plots$stage <- ggplot2::ggplot(stage_counts, ggplot2::aes(PNI, proportion, fill = stage)) +
    ggplot2::geom_col(width = .58, colour = "white", linewidth = .3) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(proportion >= .05, paste0(n, "\n", round(100 * proportion), "%"), ""),
        colour = ifelse(stage %in% c("Stage III", "Stage IV"), "white", INK)
      ),
      position = ggplot2::position_stack(vjust = .5), size = 4
    ) +
    ggplot2::scale_colour_identity(guide = "none") +
    ggplot2::scale_x_discrete(labels = stats::setNames(
      sprintf("%s\n(n = %d)", names(stage_n), stage_n), names(stage_n)
    )) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, .04))
    ) +
    ggplot2::scale_fill_manual(values = stats::setNames(
      c("#D8EFF2", "#86C5C8", "#3E8B91", "#174E57"), stage_levels
    )) +
    ggplot2::labs(
      title = "AJCC pathological stage by PNI status",
      subtitle = paste0("Omnibus Fisher P = ", fmt_p(stage_ft$p.value)),
      x = NULL, y = "Proportion within PNI group", fill = "Pathological stage"
    ) +
    theme_tcga() +
    ggplot2::theme(legend.position = "right")
  # Stable figure order matches the previous comprehensive analysis.
  plots <- plots[c("effect_overview", "stage", "immune_cells", "exhaustion", "markers", "neu")]
  list(
    plots = plots,
    tables = list(
      molecular_tests = tests, molecular_descriptive = descriptive,
      ajcc_stage_counts = stage_counts, ajcc_stage_test = stage_test
    )
  )
}


# ===================== 02_pathology.R =====================
# TCGA PNI / T-stage / nodal / LVI figures, reconstructed from patient-level data.
# This module reads no files and writes no files. Unknown categories are retained
# in alluvial and T-distribution figures, but excluded from feature-specific tests.

make_pathology_figures <- function(dat, config) {
  required <- c("sample_id", "PNI", "t_stage_raw", "n_stage_raw", "lvi_raw")
  stopifnot(
    all(required %in% names(dat)), !anyDuplicated(dat$sample_id),
    nrow(dat) == config$cohort_n,
    all(dat$PNI %in% c("Non-PNI", "PNI"))
  )
  d <- dat[, required]
  # Validate source codes explicitly instead of treating an unexpected code as negative.
  t <- toupper(trimws(as.character(d$t_stage_raw)))
  n <- toupper(trimws(as.character(d$n_stage_raw)))
  lvi <- toupper(trimws(as.character(d$lvi_raw)))
  t[is.na(t) | t == ""] <- "ND"
  n[is.na(n) | n == ""] <- "NX"
  lvi[is.na(lvi) | lvi == ""] <- "ND"
  stopifnot(
    all(t %in% c("T2", "T3", "T4", "ND")),
    all(n %in% c("N0", "N1", "N2", "N3", "NX")),
    all(lvi %in% c("NO", "YES", "ND"))
  )
  d$PNI <- factor(d$PNI, levels = c("Non-PNI", "PNI"))
  d$T_stage <- factor(ifelse(t == "ND", "T unknown", t),
    levels = c("T unknown", "T2", "T3", "T4")
  )
  d$Nodal <- factor(ifelse(n %in% c("N1", "N2", "N3"), "N+", n),
    levels = c("N0", "NX", "N+")
  )
  d$LVI <- factor(ifelse(lvi == "NO", "No", ifelse(lvi == "YES", "Yes", "ND")),
    levels = c("No", "ND", "Yes")
  )
  group_n <- table(d$PNI)
  cohort_subtitle <- sprintf(
    "TCGA-BLCA: fixed cohort n = %d (Non-PNI %d, PNI %d)",
    nrow(d), group_n[["Non-PNI"]], group_n[["PNI"]]
  )
  paths <- dplyr::count(d, PNI, T_stage, Nodal, LVI, name = "n", .drop = FALSE)
  stopifnot(sum(paths$n) == nrow(d))
  unknown_t_n <- sum(d$T_stage == "T unknown")
  unknown_n_n <- sum(d$Nodal == "NX")
  unknown_lvi_n <- sum(d$LVI == "ND")

  # Fisher 2x2 orientation: rows PNI / Non-PNI; columns positive / negative.
  bar_rows <- test_rows <- list()
  for (feature in c("Nodal", "LVI")) {
    known <- if (feature == "Nodal") d$Nodal != "NX" else d$LVI != "ND"
    positive <- if (feature == "Nodal") d$Nodal == "N+" else d$LVI == "Yes"
    tab <- table(
      factor(d$PNI[known], levels = c("PNI", "Non-PNI")),
      factor(positive[known], levels = c(TRUE, FALSE))
    )
    two <- stats::fisher.test(tab, alternative = "two.sided")
    test_rows[[feature]] <- data.frame(
      feature = feature, n = sum(known),
      odds_ratio = unname(two$estimate), ci_lower_two_sided = two$conf.int[1],
      ci_upper_two_sided = two$conf.int[2], p_two_sided = two$p.value, method = "Fisher exact",
      alternative = "two.sided",
      analysis_status = "Exploratory"
    )
    for (g in levels(d$PNI)) {
      denominator <- sum(d$PNI == g & known)
      numerator <- sum(d$PNI == g & known & positive)
      bar_rows[[paste(feature, g)]] <- data.frame(
        feature = feature, group = g,
        positive = numerator, n_known = denominator,
        n_unknown = sum(d$PNI == g & !known), proportion = numerator / denominator,
        p_two_sided = two$p.value
      )
    }
  }
  bars <- dplyr::bind_rows(bar_rows)
  bars$group <- factor(bars$group, levels = c("Non-PNI", "PNI"))
  bars$feature <- factor(bars$feature, levels = c("Nodal", "LVI"))
  bars$label <- sprintf("%d/%d\n%.1f%%", bars$positive, bars$n_known, 100 * bars$proportion)
  tests <- dplyr::bind_rows(test_rows)
  tests$q_bh <- stats::p.adjust(tests$p_two_sided, method = "BH", n = nrow(tests))
  tests$multiplicity_family <- "Nodal and LVI comparisons (2 tests)"
  bars$q_bh <- tests$q_bh[match(as.character(bars$feature), tests$feature)]
  p_labels <- tests
  p_labels$feature <- factor(p_labels$feature, levels = c("Nodal", "LVI"))
  p_labels$label <- paste0(
    "Two-sided P = ", fmt_p(p_labels$p_two_sided),
    "\nBH q = ", fmt_p(p_labels$q_bh)
  )
  bar_plot <- ggplot2::ggplot(bars, ggplot2::aes(feature, proportion, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = .66), width = .57) +
    ggplot2::geom_text(ggplot2::aes(label = label),
      position = ggplot2::position_dodge(width = .66), vjust = -.18, size = 3.4,
      lineheight = .95, family = "Arial", colour = INK
    ) +
    ggplot2::geom_segment(
      data = p_labels, ggplot2::aes(
        x = as.numeric(feature) - .34,
        xend = as.numeric(feature) + .34, y = .96, yend = .96
      ),
      inherit.aes = FALSE, linewidth = .35, colour = MUTED
    ) +
    ggplot2::geom_text(
      data = p_labels, ggplot2::aes(x = feature, y = 1.025, label = label),
      inherit.aes = FALSE, size = 3.6, family = "Arial", colour = INK
    ) +
    ggplot2::scale_fill_manual(values = PNI_COLORS, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = c(
      Nodal = "Lymph node metastasis",
      LVI = "Lymphovascular invasion"
    ), expand = ggplot2::expansion(add = .6)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.10), breaks = seq(0, 1, .25),
      labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(x = NULL, y = "Patients with feature", fill = NULL) +
    theme_tcga() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 6)),
      plot.margin = ggplot2::margin(7, 10, 3, 10)
    )
  bar_standalone <- bar_plot + ggplot2::labs(
    title = "Nodal and vascular invasion by PNI status",
    subtitle = sprintf(
      "TCGA-BLCA cohort n = %d. Known nodal status: n = %d. Known LVI status: n = %d.",
      nrow(d), nrow(d) - unknown_n_n, nrow(d) - unknown_lvi_n
    )
  ) +
    ggplot2::theme(legend.position = "top", legend.justification = "left")

  flow <- ggplot2::ggplot(
    dplyr::filter(paths, n > 0),
    ggplot2::aes(axis1 = PNI, axis2 = T_stage, axis3 = Nodal, axis4 = LVI, y = n)
  ) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = PNI),
      width = .26, alpha = .60,
      knot.pos = .45, colour = NA
    ) +
    ggalluvial::geom_stratum(width = .26, fill = "#F3F5F7", colour = "#8B97A8", linewidth = .6) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum, ggplot2::aes(
        label = ifelse(ggplot2::after_stat(stratum) == "T unknown", "",
          ifelse(ggplot2::after_stat(stratum) == "NX",
            paste0("NX: ", ggplot2::after_stat(count)),
            paste0(ggplot2::after_stat(stratum), "\n(n = ", ggplot2::after_stat(count), ")")
          )
        ),
        size = ifelse(ggplot2::after_stat(stratum) == "NX", 2.6, 3.2)
      ),
      lineheight = .93, family = "Arial", colour = INK
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_discrete(
      limits = c("PNI status", "T stage", "Lymph node status", "LVI status"),
      expand = ggplot2::expansion(add = c(.18, .18))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.01, .13))) +
    ggplot2::scale_fill_manual(values = PNI_COLORS, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_void(base_family = "Arial", base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = INK, size = 11,
        margin = ggplot2::margin(t = 6)
      ), legend.position = "top", legend.justification = "left",
      legend.text = ggplot2::element_text(size = 11), legend.key.size = grid::unit(.28, "cm"),
      plot.margin = ggplot2::margin(2, 10, 8, 10)
    )
  # Place the small unknown-T label outside its node without changing its width.
  if (unknown_t_n > 0) {
    nodes <- ggplot2::ggplot_build(flow)$data[[2]]
    unknown_node <- nodes[nodes$x == 2 & as.character(nodes$stratum) == "T unknown", ]
    stopifnot(nrow(unknown_node) == 1L, unknown_node$count == unknown_t_n)
    flow <- flow + ggplot2::annotate("segment",
      x = 2, xend = 2,
      y = unknown_node$ymax, yend = nrow(d) * 1.045, linewidth = .35, colour = "#8B97A8"
    ) +
      ggplot2::annotate("text",
        x = 2, y = nrow(d) * 1.081,
        label = sprintf("Unknown T (n = %d)", unknown_t_n), size = 3,
        family = "Arial", colour = INK
      )
  }
  combined <- patchwork::wrap_plots(flow, bar_plot, ncol = 1, heights = c(1.25, 1)) +
    patchwork::plot_annotation(
      title = "PNI, T stage and invasive pathological features",
      subtitle = cohort_subtitle, tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = "Arial", size = 18,
          face = "bold", colour = INK
        ), plot.subtitle = ggplot2::element_text(
          family = "Arial",
          size = 11, colour = MUTED
        ),
        plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 14, colour = INK),
        plot.margin = ggplot2::margin(10, 12, 8, 12)
      )
    )
  sankey <- flow + ggplot2::labs(
    title = "PNI, T stage, nodal status and LVI",
    subtitle = cohort_subtitle
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(family = "Arial", size = 19, face = "bold", colour = INK),
      plot.subtitle = ggplot2::element_text(
        family = "Arial", size = 11.5, colour = MUTED,
        margin = ggplot2::margin(b = 8)
      ),
      plot.margin = ggplot2::margin(10, 14, 10, 14)
    )

  # Latest requested 100%-stacked T-stage chart: PNI is left, Non-PNI right.
  # Denominators include all patients; unknown T remains a gray bar segment.
  stacked_data <- data.frame(
    group = factor(d$PNI, levels = c("PNI", "Non-PNI")),
    stage = factor(ifelse(t == "ND", "Unknown", t), levels = c("T2", "T3", "T4", "Unknown"))
  )
  t_counts <- dplyr::count(stacked_data, group, stage, name = "n", .drop = FALSE) |>
    dplyr::group_by(group) |>
    dplyr::mutate(group_n = sum(n), proportion = n / group_n, percent = 100 * proportion) |>
    dplyr::ungroup()
  known_t <- stacked_data[stacked_data$stage != "Unknown", ]
  t_tab <- table(known_t$group, factor(known_t$stage, levels = c("T2", "T3", "T4")))
  t_test <- stats::fisher.test(t_tab[, colSums(t_tab) > 0, drop = FALSE])
  t_tests <- data.frame(
    method = t_test$method, n = sum(t_tab), unknown_excluded = unknown_t_n,
    p_value = t_test$p.value, alternative = "Omnibus association; known T stages only"
  )
  fill <- c(T2 = "#CED1DF", T3 = "#DDD6E2", T4 = "#E5CDD9", Unknown = "#ECECEC")
  edge <- c(T2 = "#7B86A6", T3 = "#9D86A8", T4 = "#BD819B", Unknown = "#A4A4A4")
  legend_order <- c("T4", "T3", "T2", "Unknown")
  stacked <- ggplot2::ggplot(t_counts, ggplot2::aes(group, proportion, fill = stage, colour = stage)) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), width = .64, linewidth = .9) +
    ggplot2::annotate("segment",
      x = 1, xend = 2, y = 1.035, yend = 1.035,
      colour = "black", linewidth = .7
    ) +
    ggplot2::annotate("text",
      x = 1.5, y = 1.105, label = paste0("Fisher P = ", fmt_p(t_test$p.value)), fontface = "bold",
      size = 4, family = "Arial", colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = fill, breaks = legend_order) +
    ggplot2::scale_colour_manual(values = edge, breaks = legend_order) +
    ggplot2::guides(colour = "none", fill = ggplot2::guide_legend(
      override.aes = list(colour = unname(edge[legend_order]), linewidth = .9)
    )) +
    ggplot2::scale_x_discrete(
      labels = c(
        PNI = sprintf("PNI\n(n = %d)", group_n[["PNI"]]),
        "Non-PNI" = sprintf("non-PNI\n(n = %d)", group_n[["Non-PNI"]])
      ),
      expand = ggplot2::expansion(add = .6)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.16), breaks = c(0, .5, 1), labels = c("0", "50", "100"),
      minor_breaks = c(.25, .75), expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "T-stage distribution by PNI status",
      subtitle = sprintf("TCGA-BLCA cohort (n = %d)", nrow(d)), x = NULL, y = "Proportion (%)", fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      aspect.ratio = 1.45,
      plot.title = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
      plot.subtitle = ggplot2::element_text(size = 11.5, colour = "#525866", margin = ggplot2::margin(b = 8)),
      axis.title.y = ggplot2::element_text(size = 17, colour = "black", margin = ggplot2::margin(r = 9)),
      axis.text.x = ggplot2::element_text(size = 15, colour = "black", margin = ggplot2::margin(t = 8)),
      axis.text.y = ggplot2::element_text(size = 15, colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = .7),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = .7),
      axis.ticks.length = grid::unit(.15, "cm"), legend.position = "right",
      legend.text = ggplot2::element_text(size = 14, colour = "black"),
      legend.key.height = grid::unit(.65, "cm"), legend.key.width = grid::unit(.55, "cm"),
      legend.spacing.y = grid::unit(.2, "cm"),
      plot.margin = ggplot2::margin(12, 14, 12, 14)
    )
  list(
    plots = list(
      Sankey_and_bar = combined, Sankey = sankey, Bar = bar_standalone,
      T_stage_stacked = stacked
    ),
    tables = list(
      pathology_patients = d, pathology_paths = paths,
      invasion_bar_counts = bars, invasion_tests = tests,
      T_stage_counts = t_counts, T_stage_test = t_tests
    )
  )
}


# ===================== 03_survival_cox.R =====================
# TCGA-BLCA survival and Cox figures, rebuilt from patient-level data.
# No cached models or plot objects are required. OS is primary, PFI supplementary,
# and optional DFI exploratory. This module intentionally contains no DSS analysis.
# All fits use the full recorded follow-up; the plot horizon is a display limit.

make_survival_cox_figures <- function(dat, config) {
  for (pkg in c("survival", "dplyr", "broom", "ggplot2", "patchwork", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the required package: ", pkg)
  }
  suppressPackageStartupMessages({
    library(survival)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
  required <- c(
    "sample_id", "case_id", "PNI", "age", "sex", "t_stage_raw",
    "n_stage_raw", "grade_raw", "molecular_subtype_raw", "purity_absolute",
    "OS_event", "OS_days", "PFI_event", "PFI_days"
  )
  if (isTRUE(config$include_dfi)) required <- c(required, "DFI_event", "DFI_days")
  stopifnot(
    all(required %in% names(dat)), !anyDuplicated(dat$sample_id),
    !anyDuplicated(dat$case_id), all(dat$PNI %in% c("Non-PNI", "PNI")),
    nrow(dat) == config$cohort_n
  )
  horizon <- config$display_months
  if (is.null(horizon)) horizon <- 60
  B <- as.integer(config$bootstrap_reps)
  seed <- as.integer(config$seed)
  if (length(B) != 1L || is.na(B) || B < 2L) stop("bootstrap_reps must be at least 2.")
  if (length(seed) != 1L || is.na(seed)) stop("A numeric bootstrap seed is required.")
  days_per_month <- 30.4375
  risk_times <- seq(0, horizon, by = 12)
  line_types <- c("Non-PNI" = "solid", "PNI" = "22")
  plots <- list()
  tables <- list()
  saved_models <- list()
  d <- dat |> mutate(
    group = factor(PNI, levels = c("Non-PNI", "PNI")),
    t_stage_source = t_stage_raw, n_stage_source = n_stage_raw,
    sex_source = sex, grade_source = grade_raw, subtype_source = molecular_subtype_raw,
    t_stage_raw = toupper(trimws(as.character(t_stage_raw))),
    n_stage_raw = toupper(trimws(as.character(n_stage_raw)))
  )

  # Return warnings as explicit data rather than losing them in console output.
  fit_safe <- function(formula, data) {
    ws <- character()
    fit <- tryCatch(
      withCallingHandlers(
        coxph(formula,
          data = data, ties = "efron", x = TRUE, y = TRUE, model = TRUE,
          control = coxph.control(iter.max = 100)
        ),
        warning = function(w) {
          ws <<- c(ws, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    list(
      fit = fit, warnings = paste(unique(ws), collapse = " | "),
      error = if (inherits(fit, "error")) conditionMessage(fit) else ""
    )
  }
  make_formula <- function(terms) as.formula(paste("Surv(time_days, event) ~", paste(terms, collapse = " + ")))
  c_index <- function(fit, data) {
    lp <- as.numeric(predict(fit, newdata = data, type = "lp", reference = "zero"))
    if (any(!is.finite(lp))) stop("Nonfinite prediction")
    as.numeric(concordance(Surv(time_days, event) ~ lp,
      data = data,
      reverse = TRUE, timewt = "n"
    )$concordance)
  }
  endpoint_meta <- list(
    OS = list(title = "Overall survival", y = "Overall survival probability"),
    PFI = list(title = "Progression-free interval", y = "Progression-free probability")
  )
  if (isTRUE(config$include_dfi)) {
    endpoint_meta$DFI <- list(
      title = "Disease-free interval (exploratory)", y = "Disease-free probability"
    )
  }
  survival_stats <- list()
  inclusion <- list()
  endpoint_plot_parts <- list()

  for (ep in names(endpoint_meta)) {
    meta <- endpoint_meta[[ep]]
    z <- d |>
      transmute(sample_id, case_id, group,
        event = .data[[paste0(ep, "_event")]],
        time_days = .data[[paste0(ep, "_days")]]
      ) |>
      mutate(reason = case_when(
        is.na(event) ~ "Event status unavailable",
        !event %in% c(0, 1) ~ "Invalid event status",
        is.na(time_days) | !is.finite(time_days) ~ "Time unavailable",
        time_days < 0 ~ "Negative time", TRUE ~ "Included"
      ), included = reason == "Included")
    inclusion[[ep]] <- z |> mutate(endpoint = ep, .before = 1)
    z <- z |>
      filter(included) |>
      mutate(time_months = time_days / days_per_month)
    if (sum(z$event) == 0 || any(table(z$group) == 0)) stop(ep, ": insufficient survival data")
    fit <- survfit(Surv(time_months, event) ~ group, data = z, conf.type = "log")
    logrank <- survdiff(Surv(time_months, event) ~ group, data = z, rho = 0)
    p_lr2 <- pchisq(unname(logrank$chisq), 1, lower.tail = FALSE)
    fitted <- fit_safe(Surv(time_months, event) ~ group, z)
    if (fitted$error != "") stop(ep, ": ", fitted$error)
    cx <- broom::tidy(fitted$fit, exponentiate = TRUE, conf.int = TRUE)
    cox_usable <- fitted$warnings == "" && nrow(cx) == 1L &&
      all(is.finite(coef(fitted$fit))) && all(is.finite(cx$std.error)) && all(cx$std.error > 0)
    ph <- tryCatch(cox.zph(fitted$fit, transform = "km"), error = function(e) NULL)
    ph_p <- if (is.null(ph)) NA_real_ else unname(ph$table[1, "p"])
    revfit <- survfit(Surv(time_months, 1 - event) ~ 1, data = z, conf.type = "log")
    followup <- unname(summary(revfit)$table["median"])
    groups <- z |>
      group_by(group) |>
      summarise(n = n(), events = sum(event), max_time = max(time_months), .groups = "drop")
    legend_labels <- setNames(sprintf("%s (n = %d)", groups$group, groups$n), groups$group)
    full_curve <- broom::tidy(fit) |>
      mutate(group = factor(sub("^group=", "", strata), levels = levels(d$group)))
    initial <- data.frame(
      time = 0, n.risk = groups$n, n.event = 0, n.censor = 0,
      estimate = 1, std.error = 0, conf.high = 1, conf.low = 1,
      strata = paste0("group=", groups$group), group = groups$group
    )
    full_curve <- bind_rows(initial, full_curve) |> arrange(group, time)
    curve <- bind_rows(lapply(levels(d$group), function(g) {
      endpoint <- min(horizon, max(z$time_months[z$group == g]))
      a <- full_curve |>
        filter(group == g, time <= endpoint) |>
        arrange(time)
      if (tail(a$time, 1) < endpoint) {
        terminal <- a[nrow(a), , drop = FALSE]
        terminal$time <- endpoint
        terminal$n.event <- 0
        terminal$n.censor <- 0
        a <- bind_rows(a, terminal)
      }
      a
    }))
    ribbon <- curve |>
      group_by(group) |>
      group_modify(~ {
        a <- arrange(.x, time)
        if (nrow(a) < 2) {
          return(a)
        }
        bind_rows(a[1, , drop = FALSE], bind_rows(lapply(seq.int(2, nrow(a)), function(i) {
          bind_rows(mutate(a[i - 1, , drop = FALSE], time = a$time[i]), a[i, , drop = FALSE])
        })))
      }) |>
      ungroup()
    risk <- bind_rows(lapply(levels(d$group), function(g) {
      a <- z[z$group == g, ]
      data.frame(
        group = g, time = risk_times,
        n_risk = vapply(risk_times, function(t) sum(a$time_months >= t), integer(1))
      )
    }))
    risk$group <- factor(risk$group, levels = rev(levels(d$group)))
    rr <- summary(fit, times = risk_times, extend = TRUE)
    risk_check <- left_join(risk, data.frame(
      group = sub("^group=", "", rr$strata),
      time = rr$time, model_risk = rr$n.risk
    ), by = c("group", "time"))
    stopifnot(!anyNA(risk_check$model_risk), all(risk_check$n_risk == risk_check$model_risk))
    raw_risk <- bind_rows(lapply(levels(d$group), function(g) {
      data.frame(
        group = g, time = risk_times, raw_risk = vapply(risk_times, function(t) {
          sum(z$group == g & z$time_days >= t * days_per_month)
        }, integer(1))
      )
    }))
    raw_check <- left_join(risk, raw_risk, by = c("group", "time"))
    stopifnot(all(raw_check$n_risk == raw_check$raw_risk))
    main <- ggplot(curve, aes(time, estimate, colour = group, linetype = group)) +
      geom_ribbon(
        data = ribbon, aes(x = time, ymin = conf.low, ymax = conf.high, fill = group),
        inherit.aes = FALSE, alpha = .13, colour = NA
      ) +
      geom_step(linewidth = .9, direction = "hv") +
      geom_point(
        data = curve |> filter(n.censor > 0), shape = 3, size = 1.8,
        stroke = .65, show.legend = FALSE
      ) +
      scale_colour_manual(values = PNI_COLORS, breaks = names(PNI_COLORS), labels = legend_labels[names(PNI_COLORS)]) +
      scale_fill_manual(values = PNI_COLORS) +
      scale_linetype_manual(
        values = line_types, breaks = names(PNI_COLORS),
        labels = legend_labels[names(PNI_COLORS)]
      ) +
      scale_x_continuous(limits = c(0, horizon), breaks = risk_times, expand = expansion(add = c(.8, .8))) +
      scale_y_continuous(
        limits = c(0, 1.02), breaks = seq(0, 1, .25),
        labels = scales::label_percent(accuracy = 1), expand = expansion(mult = 0)
      ) +
      labs(
        title = paste0(meta$title, " by PNI status"),
        subtitle = sprintf(
          "TCGA-BLCA fixed cohort n = %d: %d analyzed, %d events, %d unavailable",
          config$cohort_n, nrow(z), sum(z$event), config$cohort_n - nrow(z)
        ),
        x = NULL, y = meta$y, colour = NULL, linetype = NULL
      ) +
      guides(fill = "none") +
      theme_tcga() +
      theme(
        legend.position = "top", legend.justification = "left",
        legend.key.width = grid::unit(1, "cm"), legend.margin = margin(0, 0, 1, 0),
        plot.margin = margin(8, 8, 2, 5)
      )
    risk_plot <- ggplot(risk, aes(time, group, label = n_risk, colour = group)) +
      geom_text(size = 3.7, family = "Arial", show.legend = FALSE) +
      scale_colour_manual(values = PNI_COLORS) +
      scale_x_continuous(limits = c(0, horizon), breaks = risk_times, expand = expansion(add = c(.8, .8))) +
      labs(title = "Number at risk", x = "Recorded follow-up (months)", y = NULL) +
      theme_tcga() +
      theme(
        plot.title = element_text(size = 12, face = "bold", colour = INK),
        axis.title.x = element_text(size = 12), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#E4E7EC", linewidth = .25),
        plot.margin = margin(1, 8, 3, 5)
      )
    endpoint_plot_parts[[ep]] <- list(main = main, risk = risk_plot)
    survival_stats[[ep]] <- data.frame(
      endpoint = ep, cohort_n = config$cohort_n, n = nrow(z),
      excluded_unavailable = config$cohort_n - nrow(z), events = sum(z$event),
      non_pni_n = sum(z$group == "Non-PNI"), pni_n = sum(z$group == "PNI"),
      non_pni_events = sum(z$event[z$group == "Non-PNI"]), pni_events = sum(z$event[z$group == "PNI"]),
      logrank_chisq = unname(logrank$chisq), logrank_p = p_lr2,
      cox_inference_usable = cox_usable, hr_raw_fit = cx$estimate,
      ci_low_raw_fit = cx$conf.low, ci_high_raw_fit = cx$conf.high, cox_p_raw_fit = cx$p.value,
      hr = if (cox_usable) cx$estimate else NA_real_,
      ci_low = if (cox_usable) cx$conf.low else NA_real_,
      ci_high = if (cox_usable) cx$conf.high else NA_real_,
      cox_wald_p = if (cox_usable) cx$p.value else NA_real_,
      ph_p = ph_p, median_followup_months = followup, fit_warning = fitted$warnings
    )
    tables[[paste0(ep, "_risk_table")]] <- risk
    tables[[paste0(ep, "_curve_estimates")]] <- full_curve
    saved_models[[ep]] <- list(data = z, km = fit, cox = fitted$fit, ph = ph, warnings = fitted$warnings)
  }
  tables$survival_summary <- bind_rows(survival_stats) |>
    mutate(
      logrank_q_bh = stats::p.adjust(logrank_p, method = "BH", n = length(endpoint_meta)),
      cox_q_bh = stats::p.adjust(cox_wald_p, method = "BH", n = length(endpoint_meta)),
      multiplicity_family = paste0("Available endpoints (", n(), ") within each test method")
    )
  for (ep in names(endpoint_meta)) {
    st <- tables$survival_summary |> filter(endpoint == ep)
    annotation <- sprintf(
      paste0(
        "Two-sided log-rank P = %s; BH q = %s\n",
        "Two-sided Cox P = %s; BH q = %s\n",
        "Unadjusted HR = %.2f (95%% CI %.2f to %.2f)"
      ),
      fmt_p(st$logrank_p), fmt_p(st$logrank_q_bh), fmt_p(st$cox_wald_p),
      fmt_p(st$cox_q_bh), st$hr, st$ci_low, st$ci_high
    )
    if (!isTRUE(st$cox_inference_usable)) {
      annotation <- paste0(
        "Two-sided log-rank P = ", fmt_p(st$logrank_p), "; BH q = ", fmt_p(st$logrank_q_bh),
        "\nCox estimate not reliably estimable; see diagnostics"
      )
    }
    main <- endpoint_plot_parts[[ep]]$main +
      annotate("text",
        x = horizon / 20, y = if (ep == "DFI") .22 else .17,
        label = annotation, hjust = 0, vjust = .5, size = 3.2, family = "Arial", colour = INK
      )
    plots[[ep]] <- main / endpoint_plot_parts[[ep]]$risk + plot_layout(heights = c(3.7, 1.0))
  }
  tables$endpoint_inclusion <- bind_rows(inclusion)

  # Normalize case/whitespace consistently; reject unknown nonmissing category codes.
  normalized_code <- function(x) toupper(trimws(as.character(x)))
  unknown_codes <- c(
    "", "NA", "N/A", "ND", "UNKNOWN",
    "NOT REPORTED", "[NOT AVAILABLE]", "[NOT EVALUATED]", "[UNKNOWN]"
  )
  check_codes <- function(x, allowed, field) {
    bad <- !is.na(x) & !x %in% c(allowed, unknown_codes)
    if (any(bad)) stop("Unrecognized ", field, " codes: ", paste(unique(x[bad]), collapse = ", "))
  }
  sex_code <- normalized_code(d$sex)
  grade_code <- normalized_code(d$grade_raw)
  subtype_code <- normalized_code(d$molecular_subtype_raw)
  check_codes(sex_code, c("FEMALE", "MALE"), "sex")
  check_codes(grade_code, c("LOW GRADE", "HIGH GRADE"), "grade")
  subtype_source <- c("Luminal_papillary", "Basal_squamous", "Luminal_infiltrated", "Luminal", "Neuronal")
  check_codes(subtype_code, toupper(subtype_source), "molecular subtype")
  d$sex <- c("Female", "Male")[match(sex_code, c("FEMALE", "MALE"))]
  d$grade_raw <- c("Low Grade", "High Grade")[match(grade_code, c("LOW GRADE", "HIGH GRADE"))]
  d$molecular_subtype_raw <- subtype_source[match(subtype_code, toupper(subtype_source))]
  if (!is.numeric(d$age) || any(!is.na(d$age) & (!is.finite(d$age) | d$age < 0))) {
    stop("Age must contain finite nonnegative numeric values or NA.")
  }
  if (!is.numeric(d$purity_absolute) || any(!is.na(d$purity_absolute) &
    (!is.finite(d$purity_absolute) | d$purity_absolute < 0 | d$purity_absolute > 1))) {
    stop("ABSOLUTE purity must be between 0 and 1 or NA.")
  }
  # Fixed exploratory main/sensitivity models; no P-value-based selection.
  d <- d |> mutate(
    pni_status = group,
    t_stage = factor(ifelse(t_stage_raw %in% c("T2", "T3", "T4"), t_stage_raw, NA_character_),
      levels = c("T2", "T3", "T4")
    ),
    n_stage = factor(ifelse(n_stage_raw %in% c("N0", "N1", "N2", "N3"), n_stage_raw, NA_character_),
      levels = c("N0", "N1", "N2", "N3")
    ),
    nodal_metastasis = factor(case_when(
      n_stage_raw == "N0" ~ "N0",
      n_stage_raw %in% c("N1", "N2", "N3") ~ "N-positive", TRUE ~ NA_character_
    ), levels = c("N0", "N-positive")),
    grade = factor(case_when(
      grade_raw == "Low Grade" ~ "Low", grade_raw == "High Grade" ~ "High",
      TRUE ~ NA_character_
    ), levels = c("Low", "High")),
    sex = factor(sex, levels = c("Female", "Male")),
    molecular_subtype = factor(case_when(
      molecular_subtype_raw == "Luminal_papillary" ~ "Luminal-papillary",
      molecular_subtype_raw == "Basal_squamous" ~ "Basal-squamous",
      molecular_subtype_raw == "Luminal_infiltrated" ~ "Luminal-infiltrated",
      molecular_subtype_raw == "Luminal" ~ "Luminal", molecular_subtype_raw == "Neuronal" ~ "Neuronal",
      TRUE ~ NA_character_
    ), levels = c("Luminal-papillary", "Basal-squamous", "Luminal-infiltrated", "Luminal", "Neuronal")),
    age_per10 = age / 10, purity_per10 = purity_absolute / .1
  )
  main_terms <- c("t_stage", "nodal_metastasis", "age_per10", "sex")
  full_terms <- c("t_stage", "n_stage", "grade", "age_per10", "sex", "molecular_subtype", "purity_per10")
  label_map <- c(
    t_stageT3 = "T3 versus T2", t_stageT4 = "T4 versus T2",
    "nodal_metastasisN-positive" = "N1-3 versus N0", n_stageN1 = "N1 versus N0", n_stageN2 = "N2 versus N0",
    n_stageN3 = "N3 versus N0", age_per10 = "Age per 10 years", sexMale = "Male versus female",
    gradeHigh = "High versus low grade", "molecular_subtypeBasal-squamous" = "Basal-squamous versus luminal-papillary",
    "molecular_subtypeLuminal-infiltrated" = "Luminal-infiltrated versus luminal-papillary",
    molecular_subtypeLuminal = "Luminal versus luminal-papillary",
    molecular_subtypeNeuronal = "Neuronal versus luminal-papillary",
    purity_per10 = "ABSOLUTE purity per 0.10", pni_statusPNI = "PNI versus Non-PNI"
  )
  all_models <- list()
  all_coefficients <- list()
  model_summary <- list()
  likelihood <- list()
  ph_tests <- list()
  nonlinear <- list()
  analysis_rows <- list()
  diagnostics <- list()
  factor_counts <- list()
  for (framework in c("Main", "Eight-variable sensitivity")) {
    terms0 <- if (framework == "Main") main_terms else full_terms
    for (ep in c("OS", "PFI")) {
      z <- d
      z$event <- z[[paste0(ep, "_event")]]
      z$time_days <- z[[paste0(ep, "_days")]]
      z$exclusion_reason <- vapply(seq_len(nrow(z)), function(i) {
        why <- character()
        if (is.na(z$event[i]) || !z$event[i] %in% c(0, 1)) why <- c(why, "Event status missing or invalid")
        if (is.na(z$time_days[i]) || !is.finite(z$time_days[i]) || z$time_days[i] < 0) {
          why <- c(why, "Time missing or invalid")
        }
        for (v in c(terms0, "pni_status")) if (is.na(z[[v]][i])) why <- c(why, paste0(v, " missing"))
        paste(why, collapse = "; ")
      }, character(1))
      z$included <- z$exclusion_reason == ""
      tag <- paste(gsub("[^a-zA-Z0-9]+", "_", framework), ep, sep = "_")
      analysis_rows[[tag]] <- z |>
        select(
          sample_id, case_id, pni_status, event, time_days,
          t_stage_source, n_stage_source, sex_source, grade_source, subtype_source,
          t_stage_raw, n_stage_raw, t_stage, nodal_metastasis, n_stage, grade, age, sex,
          molecular_subtype, purity_absolute, included, exclusion_reason
        ) |>
        mutate(framework = framework, endpoint = ep, .before = 1)
      z <- as.data.frame(z[z$included, ])
      rownames(z) <- z$case_id
      if (!nrow(z) || sum(z$event) == 0) stop(tag, ": no usable outcome records")
      pair <- list()
      for (model_name in c("M0", "M1")) {
        terms_now <- if (model_name == "M0") terms0 else c(terms0, "pni_status")
        fitted <- fit_safe(make_formula(terms_now), z)
        if (fitted$error != "") stop(tag, " ", model_name, ": ", fitted$error)
        fit <- fitted$fit
        pair[[model_name]] <- fit
        key <- paste(tag, model_name)
        cs <- as.data.frame(summary(fit)$coefficients)
        cf <- data.frame(framework,
          endpoint = ep, model = model_name, term = rownames(cs),
          comparison = unname(label_map[rownames(cs)]), beta = cs$coef, se = cs$`se(coef)`, z = cs$z,
          hr = exp(cs$coef), ci_lower = exp(cs$coef - qnorm(.975) * cs$`se(coef)`),
          ci_upper = exp(cs$coef + qnorm(.975) * cs$`se(coef)`), p_two_sided_wald = cs$`Pr(>|z|)`
        )
        cf$inference_usable <- all(is.finite(cs$coef)) && all(is.finite(cs$`se(coef)`)) &&
          all(cs$`se(coef)` > 0) && fitted$warnings == ""
        if (framework == "Main" && !all(cf$inference_usable)) {
          stop(tag, " ", model_name, ": unreliable main model; ", fitted$warnings)
        }
        cf$fit_warning <- fitted$warnings
        cf$estimate_status <- ifelse(cf$term == "gradeHigh" & grepl("infinite", fitted$warnings),
          "Not reliably estimable: low-grade stratum has zero events; monotone likelihood",
          ifelse(fitted$warnings == "", "Finite estimate; inspect diagnostics",
            "Finite estimate in a model with convergence warning; interpret cautiously"
          )
        )
        all_coefficients[[key]] <- cf
        ph <- tryCatch(cox.zph(fit, transform = "km", terms = TRUE, singledf = FALSE), error = function(e) e)
        ph_df <- if (inherits(ph, "error")) {
          data.frame(framework,
            endpoint = ep, model = model_name,
            term = "UNAVAILABLE", chisq = NA_real_, df = NA_real_, p = NA_real_, error = conditionMessage(ph)
          )
        } else {
          data.frame(framework,
            endpoint = ep, model = model_name, term = rownames(ph$table),
            chisq = ph$table[, "chisq"], df = ph$table[, "df"], p = ph$table[, "p"], error = ""
          )
        }
        ph_tests[[key]] <- ph_df
        pni <- cf[cf$term == "pni_statusPNI", , drop = FALSE]
        model_summary[[key]] <- data.frame(framework,
          endpoint = ep, model = model_name,
          n = nrow(z), events = sum(z$event),
          n_non_pni = sum(z$pni_status == "Non-PNI"), n_pni = sum(z$pni_status == "PNI"),
          events_non_pni = sum(z$event[z$pni_status == "Non-PNI"]), events_pni = sum(z$event[z$pni_status == "PNI"]),
          coefficients_df = length(coef(fit)), events_per_coefficient = sum(z$event) / length(coef(fit)),
          pni_hr = if (nrow(pni)) pni$hr else NA_real_, pni_ci_lower = if (nrow(pni)) pni$ci_lower else NA_real_,
          pni_ci_upper = if (nrow(pni)) pni$ci_upper else NA_real_,
          pni_p_two_sided = if (nrow(pni)) pni$p_two_sided_wald else NA_real_,
          inference_usable = all(cf$inference_usable),
          log_likelihood = as.numeric(logLik(fit)),
          apparent_harrell_c = if (all(cf$inference_usable)) c_index(fit, z) else NA_real_,
          ph_global_p = if (!inherits(ph, "error")) ph$table["GLOBAL", "p"] else NA_real_, fit_warning = fitted$warnings
        )
        diagnostics[[key]] <- data.frame(framework,
          endpoint = ep, model = model_name, n = nrow(z),
          iterations = fit$iter, nonfinite_beta = any(!is.finite(coef(fit))), warnings = fitted$warnings
        )
        spline_terms <- sub("^age_per10$", "splines::ns(age_per10, df = 3)", terms_now)
        spfit <- fit_safe(make_formula(spline_terms), z)
        if (spfit$error == "") {
          lr <- 2 * (as.numeric(logLik(spfit$fit)) - as.numeric(logLik(fit)))
          delta_df <- attr(logLik(spfit$fit), "df") - attr(logLik(fit), "df")
          nonlinear[[key]] <- data.frame(framework,
            endpoint = ep, model = model_name,
            check = "Age: natural spline df=3 versus linear", lr_chisq = lr, df = delta_df,
            p_omnibus = pchisq(lr, delta_df, lower.tail = FALSE), warnings = spfit$warnings, error = ""
          )
        } else {
          nonlinear[[key]] <- data.frame(framework,
            endpoint = ep, model = model_name,
            check = "Age: natural spline df=3 versus linear", lr_chisq = NA_real_, df = 2,
            p_omnibus = NA_real_, warnings = spfit$warnings, error = spfit$error
          )
        }
      }
      stopifnot(identical(rownames(model.frame(pair$M0)), rownames(model.frame(pair$M1))))
      lr <- 2 * (as.numeric(logLik(pair$M1)) - as.numeric(logLik(pair$M0)))
      likelihood[[tag]] <- data.frame(framework,
        endpoint = ep, comparison = "M1 versus M0: add PNI",
        n = nrow(z), events = sum(z$event), lr_chisq = lr, df = 1L,
        p_lrt_omnibus = pchisq(lr, 1, lower.tail = FALSE)
      )
      for (v in c(terms0, "pni_status")) {
        if (is.factor(z[[v]])) {
          factor_counts[[paste(tag, v)]] <- z |>
            group_by(level = .data[[v]]) |>
            summarise(n = n(), events = sum(event), .groups = "drop") |>
            mutate(framework = framework, endpoint = ep, variable = v, .before = 1)
        }
      }
      all_models[[tag]] <- list(models = pair, data = z)
    }
  }

  # Patient-paired internal bootstrap: identical resamples for M0 and M1.
  # Keep input row order and a fixed seed for reproducible resampling.
  os <- all_models$Main_OS
  c0 <- c_index(os$models$M0, os$data)
  c1 <- c_index(os$models$M1, os$data)
  set.seed(seed)
  boot <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(os$data), size = nrow(os$data), replace = TRUE)
    zb <- os$data[idx, , drop = FALSE]
    f0 <- fit_safe(make_formula(main_terms), zb)
    f1 <- fit_safe(make_formula(c(main_terms, "pni_status")), zb)
    fail <- paste(c(f0$error, f1$error)[c(f0$error, f1$error) != ""], collapse = " | ")
    ws <- paste(c(f0$warnings, f1$warnings)[c(f0$warnings, f1$warnings) != ""], collapse = " | ")
    values <- rep(NA_real_, 4L)
    if (fail == "") {
      if (any(!is.finite(c(coef(f0$fit), coef(f1$fit))))) fail <- "Nonfinite coefficient"
      if (grepl("did not converge|infinite", ws, ignore.case = TRUE)) fail <- paste("Fit convergence warning:", ws)
    }
    if (fail == "") {
      values <- tryCatch(
        c(
          c_index(f0$fit, zb), c_index(f1$fit, zb),
          c_index(f0$fit, os$data), c_index(f1$fit, os$data)
        ),
        error = function(e) {
          fail <<- conditionMessage(e)
          rep(NA_real_, 4L)
        }
      )
      if (any(!is.finite(values))) fail <- paste(fail, "Nonfinite C index")
    }
    boot[[b]] <- data.frame(
      replicate = b, n = nrow(zb), events = sum(zb$event),
      success = fail == "", failure_reason = fail, warnings = ws,
      m0_c_boot_train = values[1], m1_c_boot_train = values[2],
      m0_c_test_original = values[3], m1_c_test_original = values[4]
    )
    if (b %% 100 == 0) message("Paired OS bootstrap: ", b, "/", B)
  }
  boot <- bind_rows(boot) |> mutate(
    m0_optimism = m0_c_boot_train - m0_c_test_original,
    m1_optimism = m1_c_boot_train - m1_c_test_original,
    apparent_delta_in_bootstrap = m1_c_boot_train - m0_c_boot_train,
    delta_optimism = m1_optimism - m0_optimism
  )
  ok <- boot |> filter(success)
  if (!nrow(ok)) stop("No valid bootstrap resamples; inspect convergence.")
  opt0 <- mean(ok$m0_optimism)
  opt1 <- mean(ok$m1_optimism)
  ci_delta <- quantile(ok$apparent_delta_in_bootstrap, c(.025, .975), names = FALSE)
  pred <- data.frame(
    endpoint = "OS", framework = "Main", n = nrow(os$data), events = sum(os$data$event),
    bootstrap_requested = B, bootstrap_success = nrow(ok), bootstrap_failed = B - nrow(ok), seed = seed,
    m0_apparent_c = c0, m1_apparent_c = c1, apparent_delta_c = c1 - c0,
    apparent_delta_bootstrap_percentile_ci_lower = ci_delta[1],
    apparent_delta_bootstrap_percentile_ci_upper = ci_delta[2],
    m0_mean_optimism = opt0, m1_mean_optimism = opt1,
    m0_optimism_corrected_c = c0 - opt0, m1_optimism_corrected_c = c1 - opt1,
    optimism_corrected_delta_c = (c1 - opt1) - (c0 - opt0),
    delta_optimism_mcse = sd(ok$delta_optimism) / sqrt(nrow(ok)),
    interval_scope = "95% percentile interval is for apparent delta C, not optimism-corrected delta C.",
    validation_scope = "Internal bootstrap only; no external validation or fixed-time calibration."
  )
  tables <- c(tables, list(
    summary_models = bind_rows(model_summary), all_coefficients = bind_rows(all_coefficients),
    likelihood_comparisons = bind_rows(likelihood), prediction_summary = pred, ph_tests = bind_rows(ph_tests),
    age_nonlinearity = bind_rows(nonlinear), fit_diagnostics = bind_rows(diagnostics),
    analysis_rows = bind_rows(analysis_rows),
    factor_level_counts = bind_rows(factor_counts), bootstrap_replicates = boot
  ))

  tables$all_coefficients <- tables$all_coefficients |>
    mutate(
      p_two_sided_wald_raw_fit = p_two_sided_wald,
      p_two_sided_wald = ifelse(inference_usable, p_two_sided_wald, NA_real_)
    ) |>
    group_by(framework, endpoint, model) |>
    mutate(q_bh = stats::p.adjust(p_two_sided_wald, method = "BH", n = n())) |>
    ungroup()
  tables$summary_models <- tables$summary_models |>
    mutate(
      pni_p_two_sided_raw_fit = pni_p_two_sided,
      pni_p_two_sided = ifelse(inference_usable, pni_p_two_sided, NA_real_),
      pni_q_bh = NA_real_
    )
  m1_rows <- tables$summary_models$model == "M1"
  tables$summary_models$pni_q_bh[m1_rows] <-
    stats::p.adjust(tables$summary_models$pni_p_two_sided[m1_rows], method = "BH", n = 4L)
  tables$likelihood_comparisons <- tables$likelihood_comparisons |>
    mutate(p_lrt_omnibus_raw_fit = p_lrt_omnibus)
  for (i in seq_len(nrow(tables$likelihood_comparisons))) {
    a <- tables$likelihood_comparisons[i, ]
    usable <- tables$summary_models |>
      filter(framework == a$framework, endpoint == a$endpoint)
    if (nrow(usable) != 2L || !all(usable$inference_usable)) {
      tables$likelihood_comparisons$p_lrt_omnibus[i] <- NA_real_
    }
  }
  tables$likelihood_comparisons$q_bh <-
    stats::p.adjust(tables$likelihood_comparisons$p_lrt_omnibus, method = "BH", n = 4L)

  cf <- tables$all_coefficients |>
    filter(framework == "Main", endpoint == "OS", model == "M1") |>
    mutate(
      label = factor(comparison, levels = rev(comparison)), is_pni = term == "pni_statusPNI",
      p_label = paste0("P = ", fmt_p(p_two_sided_wald), "; q = ", fmt_p(q_bh)),
      text = sprintf("%.3f (%.3f to %.3f)   %s", hr, ci_lower, ci_upper, p_label)
    )
  os_summary <- tables$summary_models |> filter(framework == "Main", endpoint == "OS", model == "M1")
  plots$cox_os <- ggplot(cf, aes(hr, label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "#94A3B8") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = .18, colour = "#667085", linewidth = .65) +
    geom_point(aes(colour = is_pni), size = 3) +
    scale_colour_manual(values = c("FALSE" = PNI_COLORS[["Non-PNI"]], "TRUE" = PNI_COLORS[["PNI"]])) +
    scale_x_log10() +
    scale_y_discrete(labels = stats::setNames(paste0(as.character(cf$label), "\n", cf$text), as.character(cf$label))) +
    labs(
      title = "Multivariable overall survival model",
      subtitle = sprintf(
        "%d complete cases, %d deaths. T stage, nodal metastasis, age, sex, and PNI.",
        os_summary$n, os_summary$events
      ),
      x = "Adjusted hazard ratio (two-sided 95% CI)", y = NULL
    ) +
    theme_tcga() +
    theme(legend.position = "none", panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10.5))
  cs <- tables$summary_models |>
    filter(model == "M1") |>
    mutate(
      label_text = paste0(endpoint, ifelse(framework == "Main", "  Main", "  Full sensitivity"), "  (n = ", n, ")"),
      label = factor(label_text, levels = rev(label_text)),
      text = ifelse(inference_usable,
        sprintf(
          "%.3f (%.3f to %.3f)   P = %s; q = %s", pni_hr, pni_ci_lower, pni_ci_upper,
          fmt_p(pni_p_two_sided), fmt_p(pni_q_bh)
        ), "Not reliably estimable; see diagnostics"
      ),
      pni_hr = ifelse(inference_usable, pni_hr, NA_real_),
      pni_ci_lower = ifelse(inference_usable, pni_ci_lower, NA_real_),
      pni_ci_upper = ifelse(inference_usable, pni_ci_upper, NA_real_)
    )
  plots$cox_pni <- ggplot(cs, aes(pni_hr, label, colour = framework)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "#94A3B8") +
    geom_errorbarh(aes(xmin = pni_ci_lower, xmax = pni_ci_upper), height = .18, linewidth = .65, na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_colour_manual(values = c("Main" = PNI_COLORS[["PNI"]], "Eight-variable sensitivity" = "#8B929C")) +
    scale_x_log10() +
    scale_y_discrete(drop = FALSE, labels = stats::setNames(paste0(cs$label_text, "\n", cs$text), cs$label_text)) +
    labs(
      title = "PNI after adjustment for clinical covariates",
      subtitle = "Two-sided Wald P and 95% CI; BH q across four PNI model estimates.",
      x = "Adjusted PNI hazard ratio", y = NULL
    ) +
    theme_tcga() +
    theme(legend.position = "none", panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10.5))
  cd <- data.frame(
    model = factor(rep(c("Clinical model", "+ PNI"), 2), levels = c("Clinical model", "+ PNI")),
    estimate_type = rep(c("Apparent", "Optimism-corrected"), each = 2),
    C = c(pred$m0_apparent_c, pred$m1_apparent_c, pred$m0_optimism_corrected_c, pred$m1_optimism_corrected_c)
  )
  cp <- ggplot(cd, aes(model, C, group = estimate_type, colour = estimate_type)) +
    geom_line(linewidth = .7) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.4f", C)), nudge_y = .008, size = 3.8) +
    scale_colour_manual(values = c("Apparent" = PNI_COLORS[["PNI"]], "Optimism-corrected" = PNI_COLORS[["Non-PNI"]])) +
    scale_y_continuous(expand = expansion(mult = c(.10, .20))) +
    labs(title = "OS discrimination", x = NULL, y = "Harrell C", colour = NULL) +
    theme_tcga() +
    theme(legend.position = "bottom", plot.title = element_text(size = 15))
  dd <- data.frame(
    label = factor(c("Apparent delta C", "Corrected delta C"), levels = c("Corrected delta C", "Apparent delta C")),
    delta = c(pred$apparent_delta_c, pred$optimism_corrected_delta_c),
    lower = c(ci_delta[1], NA_real_), upper = c(ci_delta[2], NA_real_)
  )
  dp <- ggplot(dd, aes(delta, label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#94A3B8") +
    geom_errorbarh(data = dd[1, ], aes(xmin = lower, xmax = upper), height = .18, linewidth = .7) +
    geom_point(size = 3, colour = PNI_COLORS[["PNI"]]) +
    geom_text(aes(label = sprintf("%.4f", delta)), nudge_y = -.2, size = 3.8, colour = INK) +
    scale_x_continuous(expand = expansion(mult = c(.10, .10))) +
    labs(title = "Increment from adding PNI", x = "Change in C", y = NULL) +
    theme_tcga() +
    theme(plot.title = element_text(size = 15), axis.text.y = element_text(size = 10))
  plots$prediction <- (cp | dp) + plot_annotation(
    title = "Internal comparison of clinical models with and without PNI",
    subtitle = sprintf(
      "OS: %d patients, %d deaths. Same patients in both models. %s paired bootstrap resamples.",
      os_summary$n, os_summary$events, format(B, big.mark = ",")
    ),
    theme = theme(
      plot.title = element_text(family = "Arial", size = 20, face = "bold", colour = INK),
      plot.subtitle = element_text(family = "Arial", size = 10.5, colour = MUTED)
    )
  )
  low <- tables$factor_level_counts |>
    filter(framework == "Eight-variable sensitivity", variable == "grade", level == "Low")
  subtype_ph <- tables$ph_tests |>
    filter(framework == "Eight-variable sensitivity", model == "M1", term == "molecular_subtype", p < .05)
  ws <- unique(tables$fit_diagnostics$warnings[tables$fit_diagnostics$warnings != ""])
  saved_models$multivariable <- all_models
  list(plots = plots, tables = tables, models = saved_models)
}


tcga_script_dir <- function() {
  src <- unlist(lapply(sys.frames(), function(f) if (!is.null(f$ofile)) f$ofile else NULL), use.names = FALSE)
  if (length(src)) {
    return(dirname(normalizePath(tail(src, 1), mustWork = TRUE)))
  }
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    return(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)))
  }
  getwd()
}

run_tcga_visualization <- function(base_dir = tcga_script_dir(), args = commandArgs(trailingOnly = TRUE)) {
  opt <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (length(hit)) sub(paste0("^--", name, "="), "", tail(hit, 1)) else default
  }
  if ("--help" %in% args) {
    cat("Usage: Rscript TCGA_visualization.R --cohort=247 [--output=DIR] [--dpi=600] [--bootstrap=1000] [--no-dfi] [--no-pptx]\n")
    return(invisible(NULL))
  }
  n <- as.integer(opt("cohort", DEFAULT_COHORT_N))
  dpi <- as.integer(opt("dpi", DEFAULT_DPI))
  if (!n %in% c(247L, 228L)) stop("Choose --cohort=247 or --cohort=228.")
  if (is.na(dpi) || dpi < 72L) stop("dpi must be at least72.")
  config <- list(
    cohort_n = n, include_dfi = DEFAULT_INCLUDE_DFI && !"--no-dfi" %in% args,
    bootstrap_reps = as.integer(opt("bootstrap", DEFAULT_BOOTSTRAP_REPS)), seed = 20261101L, display_months = 60,
    export_pptx = DEFAULT_EXPORT_PPTX && !"--no-pptx" %in% args
  )
  packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "survival", "broom", "ggalluvial", "ragg", "jsonlite")
  if (config$export_pptx) packages <- c(packages, "eoffice", "officer", "xml2", "zip")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". See README.md; this script does not install packages automatically."
    )
  }
  input <- file.path(base_dir, "data", paste0("TCGA_", n, "_analysis_input.csv"))
  specs <- file.path(base_dir, "data", "metric_specs.csv")
  if (!file.exists(input) || !file.exists(specs)) {
    stop("Missing input files. Extract the complete ZIP and keep data/ beside TCGA_visualization.R.")
  }
  dat <- read.csv(input,
    check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
  )
  metric_specs <- read.csv(specs, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  stopifnot(
    nrow(dat) == n, !anyDuplicated(dat$sample_id), !anyDuplicated(dat$case_id),
    all(dat$PNI %in% c("Non-PNI", "PNI"))
  )
  output_dir <- normalizePath(opt("output", file.path(base_dir, paste0("results_", n, "_two_sided"))), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
  message("Computing fixed cohort ", n, ": molecular and pathology summaries...")
  molecular <- make_molecular_figures(dat, config, metric_specs)
  pathology <- make_pathology_figures(dat, config)
  message("Computing survival models and ", config$bootstrap_reps, " paired bootstrap resamples...")
  surv <- make_survival_cox_figures(dat, config)
  plots <- c(molecular$plots, surv$plots, pathology$plots)
  tables <- c(molecular$tables, surv$tables, pathology$tables)
  stopifnot(!anyDuplicated(names(plots)), !anyDuplicated(names(tables)))
  for (key in names(tables)) {
    utils::write.csv(tables[[key]],
      file.path(output_dir, "tables", paste0(key, ".csv")),
      row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
  }
  pptx <- file.path(output_dir, paste0("TCGA_", n, "_visualizations_eoffice.pptx"))
  for (i in seq_along(plots)) {
    key <- names(plots)[i]
    dimensions <- if (key == "T_stage_stacked") {
      c(7.4, 7.8)
    } else if (key == "Sankey_and_bar") {
      c(8.8, 7.6)
    } else if (key == "effect_overview") {
      c(8.8, 7.4)
    } else {
      c(8.8, 6.2)
    }
    ggplot2::ggsave(file.path(output_dir, "figures", paste0(sprintf("%02d_", i), key, ".png")), plots[[i]],
      device = ragg::agg_png, width = dimensions[1], height = dimensions[2], dpi = dpi, bg = "white"
    )
    if (config$export_pptx) {
      eoffice::topptx(
        figure = plots[[i]], filename = pptx, title = "",
        width = 8.8, height = 6.2, append = i > 1, devsize = FALSE, units = "in"
      )
    }
  }
  if (config$export_pptx) {
    deck <- officer::read_pptx(pptx)
    for (i in seq_along(plots)) {
      deck <- officer::on_slide(deck, i)
    }
    print(deck, target = pptx)
    repair_pptx_lines(pptx)
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  jsonlite::write_json(list(
    config = config, input_md5 = unname(tools::md5sum(input)), metric_specs_md5 = unname(tools::md5sum(specs)),
    analysis_revision = "two-sided exploratory with BH families", figures = names(plots), tables = names(tables)
  ), file.path(output_dir, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  message("Completed: ", length(plots), " figures. Output: ", output_dir)
  invisible(list(plots = plots, tables = tables, models = surv$models, config = config))
}

# source() in RStudio also leaves TCGA_results available for further ggplot editing.
TCGA_results <- run_tcga_visualization()
