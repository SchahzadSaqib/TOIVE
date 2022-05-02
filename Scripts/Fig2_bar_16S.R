# Load packages ----
pkgs <- list.files(here::here(
  "Scripts",
  "to_load"
),
full.names = TRUE
)

pkgs <- purrr::map(pkgs, source)


# read and process data ----

TOIVE_phy_16S <- readRDS(here::here(
  "Data",
  "Analysis",
  "3_TOIVE_phy_16S.rds"
))

# Convert raw counts to relative abundance
TOIVE_phy_16S_rel <- phyloseq::transform_sample_counts(
  TOIVE_phy_16S,
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  )

## Extract L.crispatus abundance in each sample
L.crispatus_prc <- data.frame(TOIVE_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.crispatus")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  TOIVE_phy_16S_rel
)$L.crispatus_prc <- L.crispatus_prc$Lactobacillus.crispatus


## Extract L.iners abundance in each sample
L.iners_prc <- data.frame(TOIVE_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.iners")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  TOIVE_phy_16S_rel
)$L.iners_prc <- L.iners_prc$Lactobacillus.iners


## Extract L.crispatus abundance in each sample
G.vag_prc <- data.frame(TOIVE_phy_16S_rel@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  TOIVE_phy_16S_rel
)$G.vag_prc <- G.vag_prc$Gardnerella.vaginalis

## Extract ASV table and elongate
ASV_table <- data.frame(TOIVE_phy_16S_rel@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(
    cols = !SampleN,
    names_to = "taxa",
    values_to = "Abundance"
  )

## Extract meta table
meta_plot <- data.frame(TOIVE_phy_16S_rel@sam_data) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~ as.factor(.x)))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot) %>%
  dplyr::filter(!CsCtrl == 0) %>%
  as.data.frame()



# define plot variables ----
path_to_fig <- here::here(
  "Results",
  "16S"
)

## Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}


# Figure 2: stacked bar plot ----

p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl, ends_with("_prc")) %>%
  dplyr::mutate(
    CsCtrl = case_when(
      (CsCtrl == 1) ~ "RPL",
      (CsCtrl == 2) ~ "Control"
    ),
    taxa = if_else(Abundance < 0.05, "Others", taxa)
  ) %>%
  dplyr::arrange(SampleN) %>%
  dplyr::mutate(
    SampleN = factor(SampleN, levels = unique(SampleN)),
    taxa = str_remove(taxa, "X."),
    taxa = str_replace(taxa, "\\.", " "),
    across(starts_with(c("id")) | ends_with("_prc"),
      .fns = ~ as.numeric(as.character(.x))
    )
  ) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(Clin_bar = 1.07)


## define species colours
mycols <- viridis::viridis(5,
  begin = 0.30,
  end = 0.90,
  option = "mako"
) %>%
  base::append(viridis::viridis(5,
    option = "rocket",
    begin = 0.3,
    end = 0.90,
    direction = -1
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = 5,
    name = "Sailboat"
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = 5,
    name = "Winter"
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = 5,
    name = "Mushroom"
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = nlevels(p1$taxa) - 25,
    name = "Cascades"
  )) %>%
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1,
    palette = "grayC",
    begin = 0.65
  ))


## Create custom colours for case/control annotation
pal <- PNWColors::pnw_palette(name = "Moth", n = 5) %>%
  .[c(1, 3)] %>%
  purrr::set_names(unique(p1$CsCtrl))



# define sorting order based on top taxa abundances in samples
p1 <- p1 %>%
  dplyr::mutate(taxa_abd = case_when(
    (L.crispatus_prc >= 0.50) ~ 1,
    (L.iners_prc >= 0.50) ~ 2,
    (G.vag_prc >= 0.50) ~ 3,
    TRUE ~ 4
  )) %>%
  dplyr::arrange(desc(CsCtrl), Sample_grp, taxa_abd, desc(Abundance)) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN))) %>%
  dplyr::mutate(Sample_grp = ifelse(Sample_grp == "Endometrium", "(A) Endometrium", "(B) Vagina"))


## Define plot theme
plot_theme(text_size = 10)
plot_colour <- "grey50"

# initiate and construct the plot
stack_bar <- ggplot2::ggplot(p1, aes(
  x = SampleN,
  y = Abundance
)) +

  # add bar for cases vs control
  ggplot2::geom_bar(aes(
    y = Clin_bar,
    fill = CsCtrl
  ),
  stat = "identity",
  width = 1,
  position = position_dodge()
  ) +

  # add manual fill
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::guides(fill = guide_legend(ncol = 1)) +
  ggplot2::facet_wrap(. ~ Sample_grp, scales = "free_x") +

  # start new scale fill
  ggnewscale::new_scale_fill() +

  # add bar for species
  ggplot2::geom_bar(aes(fill = taxa),
    stat = "identity",
    width = 1,
    position = position_stack(reverse = T)
  ) +
  ggplot2::scale_fill_manual(values = mycols, name = "Species") +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.x = element_text(
      colour = plot_colour,
      size = 14,
      face = "italic", hjust = 0
    ),
    axis.text.x = element_text(colour = plot_colour, size = 6),
    axis.title.x = element_text(hjust = -0.03)
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 7))

ggplot2::ggsave(stack_bar,
  path = path_to_fig,
  filename = "Fig2_TOIVE_stacked_bar_16S.pdf",
  device = "pdf",
  height = 14,
  width = 52,
  units = "cm",
  dpi = 320
)



# Figure 2: polar stacked plot ----

p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl, ends_with("_prc")) %>%
  dplyr::mutate(CsCtrl = case_when(
    (CsCtrl == 1) ~ "Case",
    (CsCtrl == 2) ~ "Control"
  )) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.05, "Others", taxa)) %>%
  dplyr::arrange(SampleN) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN))) %>%
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(across(starts_with(c("id")) | ends_with("_prc"),
    .fns = ~ as.numeric(as.character(.x))
  )) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(value = 1.07) %>%
  dplyr::mutate(value2 = 1.14)

# define species colours
mycols <- viridis::viridis(5,
  option = "mako",
  begin = 0.30
) %>%
  base::append(viridis::viridis(5,
    option = "rocket",
    begin = 0.3,
    end = 0.90,
    direction = -1
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = nlevels(p1$taxa) - 10,
    name = "Shuksan"
  )) %>%
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1,
    palette = "grayC",
    begin = 0.65
  ))

pal <- PNWColors::pnw_palette(name = "Moth", n = 5) %>%
  .[c(1, 5)] %>%
  purrr::set_names(unique(p1$Sample_grp))


p1 <- p1 %>%
  dplyr::arrange(Sample_grp, SampleN) %>%
  dplyr::mutate(CsCtrl = factor(CsCtrl))



# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nbreaks <- nlevels(p1$CsCtrl)
nObs <- nlevels(factor(plot_data$taxa))
to_add <- data.frame(matrix(NA, empty_bar * nbreaks * nObs, ncol(p1)))
colnames(to_add) <- colnames(p1)
to_add$CsCtrl <- rep(levels(p1$CsCtrl), each = empty_bar * nObs)

pcmpl <- p1 %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::mutate(taxa_abd = case_when(
    (L.crispatus_prc >= 0.50 | is.na(L.crispatus_prc)) ~ 4,
    (L.iners_prc >= 0.50) ~ 3,
    (G.vag_prc >= 0.50) ~ 2,
    TRUE ~ 1
  )) %>%
  dplyr::arrange(CsCtrl, Sample_grp, taxa_abd, SampleN) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(labels = rep(seq(1, nlevels(SampleN) + empty_bar * nbreaks), each = nObs)) %>%
  dplyr::mutate(SampleN = as.character(SampleN)) %>%
  dplyr::mutate(labels = factor(labels, levels = unique(labels))) %>%
  dplyr::mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  dplyr::mutate(CsCtrl = factor(CsCtrl, levels = unique(CsCtrl)))

# Get the name and the y position of each label
label_data <- pcmpl %>%
  dplyr::group_by(labels, SampleN, CsCtrl) %>%
  dplyr::summarise(tot = sum(Abundance)) %>%
  dplyr::mutate(labels = as.numeric(as.character(labels)))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$labels - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)


# prepare a data frame for base lines
base_data <- pcmpl %>%
  dplyr::group_by(CsCtrl) %>%
  dplyr::mutate(labels = as.numeric(as.character(labels))) %>%
  dplyr::summarize(start = min(labels), end = max(labels) - empty_bar) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(title = mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

# Make the plot
polar_stack <- ggplot2::ggplot(pcmpl) +

  # add bar for sample group
  ggplot2::geom_bar(aes(
    x = labels,
    y = value,
    fill = Sample_grp
  ),
  stat = "identity",
  position = "dodge"
  ) +
  ggplot2::scale_fill_manual(values = pal) +

  # start new scale
  ggnewscale::new_scale_fill() +

  # add stacked bar for taxa
  ggplot2::geom_bar(aes(
    x = labels,
    y = Abundance,
    fill = taxa
  ),
  stat = "identity",
  position = position_stack(reverse = T)
  ) +
  ggplot2::scale_fill_manual(values = mycols, name = "Species") +
  ggplot2::theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +

  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text",
    x = max(as.numeric(as.character(pcmpl$labels))) + 3,
    y = c(0, .25, .50, .75, 1.0),
    label = c("0", "25", "50", "75", "100"),
    color = "grey60",
    size = 2,
    angle = 0,
    fontface = "bold",
    hjust = 1
  ) +
  ggplot2::guides(fill = guide_legend(ncol = 2, order = 1)) +
  ggplot2::ylim(-.50, 1.33) + # control inner circle radius
  ggplot2::coord_polar() +
  ggplot2::scale_x_discrete() +

  # Add labels on top of each bar
  ggplot2::geom_text(
    data = label_data,
    aes(
      x = labels,
      y = 1.15,
      label = SampleN,
      hjust = hjust
    ),
    color = "grey60",
    size = 1.5,
    angle = label_data$angle,
    inherit.aes = FALSE
  ) +

  # Add base line information
  ggplot2::geom_segment(
    data = base_data,
    aes(
      x = start,
      y = -0.03,
      xend = end,
      yend = -0.03
    ),
    colour = "grey50",
    alpha = 0.8,
    size = 0.6,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = base_data,
    aes(
      x = title,
      y = -0.45,
      label = CsCtrl
    ),
    hjust = c(-1.2, 1.8),
    colour = "grey50",
    alpha = 0.8, size = 2,
    fontface = "bold",
    inherit.aes = FALSE
  )


# Save as pdf
ggplot2::ggsave(polar_stack,
  path = path_to_fig,
  file = "Fig2_circ_cs_vs_ctrl.pdf",
  device = "pdf",
  dpi = 320,
  width = 30,
  height = 20,
  units = "cm"
)



# Figure 2: stacked summarised bar plot ----

sample_tally <- plot_data %>%
  dplyr::distinct(SampleN, .keep_all = T) %>%
  dplyr::group_by(CsCtrl) %>%
  dplyr::count() %>%
  dplyr::mutate(CsCtrl = if_else(CsCtrl == 1, "Case", "Control"))


p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl) %>%
  dplyr::filter(!CsCtrl == 0) %>%
  dplyr::mutate(CsCtrl = if_else(CsCtrl == 1, "Case", "Control")) %>%
  dplyr::group_by(taxa, CsCtrl) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.005, "Others", taxa)) %>%
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)))


# define species colours
mycols <- viridis::viridis(5,
  option = "mako",
  begin = 0.30
) %>%
  base::append(viridis::viridis(5,
    option = "rocket",
    begin = 0.3,
    end = 0.90,
    direction = -1
  )) %>%
  base::append(PNWColors::pnw_palette(
    n = nlevels(p1$taxa) - 10,
    name = "Shuksan"
  )) %>%
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1,
    palette = "grayC",
    begin = 0.65
  )) %>%
  as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))

pal <- PNWColors::pnw_palette("Sunset", 7) %>%
  .[c(1, 6)] %>%
  purrr::set_names(levels(p1$CsCtrl))


p2 <- p1 %>%
  dplyr::inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  as.data.frame()


stack_bar <- ggplot2::ggplot(p2, aes(
  x = CsCtrl,
  y = Abundance
)) +

  # add stackedd bar
  ggplot2::geom_bar(aes(fill = taxa),
    stat = "identity",
    position = position_stack(reverse = T)
  ) +

  # add manual colour to bars
  ggplot2::scale_fill_manual(values = unique(p2$colours), name = "Species") +
  ggplot2::geom_text(aes(
    label = if_else(Abundance > 0.01,
      paste(
        round(Abundance * 100, digits = 2),
        "%"
      ),
      ""
    ),
    size = Abundance,
    colour = shades,
    group = taxa,
    angle = if_else(Abundance > 0.10, 0, 270)
  ),
  position = position_stack(0.5, reverse = T),
  show.legend = F
  ) +

  # add scale colour to text
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), name = "Species") +

  # add sample tally
  ggplot2::geom_text(
    data = sample_tally,
    aes(
      x = CsCtrl,
      y = 1.03,
      label = paste("n = ", n, sep = "")
    ),
    colour = plot_colour
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 270, colour = plot_colour)
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::coord_flip() +
  ggplot2::guides(fill = guide_legend(ncol = 4))



ggplot2::ggsave(stack_bar,
  path = path_to_fig,
  filename = "Fig2_TOIVE_stacked_bar_summ_16S.pdf",
  device = "pdf",
  height = 10,
  width = 25,
  units = "cm",
  dpi = 320
)
