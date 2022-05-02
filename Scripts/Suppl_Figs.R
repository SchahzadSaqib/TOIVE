# load packages ----
pkgs <- list.files(here::here(
  "Scripts",
  "to_load"
),
full.names = TRUE
)

pkgs <- purrr::map(pkgs, source)


# Read and clean data ----
Toive_phy_ITS <- base::readRDS(here::here(
  "Data",
  "Analysis",
  "3_Toive_phy_ITS.rds"
))

# Convert raw counts to relative abundance
Toive_phy_ITS_rel <- phyloseq::transform_sample_counts(
  Toive_phy_ITS, function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  )

## Extract C.albicans abundance in each sample
C.albicans_prc <- data.frame(Toive_phy_ITS_rel@otu_table) %>%
  dplyr::select("Candida.albicans")

## Add update C.albicans abundance to metadata
phyloseq::sample_data(
  Toive_phy_ITS_rel
)$C.albicans_prc <- C.albicans_prc$Candida.albicans

# Rearrange data for plotting ----
## Extract ASV table and elongate
ASV_table <- data.frame(Toive_phy_ITS_rel@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(
    cols = !SampleN,
    names_to = "taxa",
    values_to = "Abundance"
  )

## Extract meta table
meta_plot <- data.frame(Toive_phy_ITS_rel@sam_data) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  dplyr::mutate(across(
    .cols = everything(),
    .fns = ~ as.factor(.x)
  ))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot) %>%
  as.data.frame()


# define plot variables ----
## Define path to figures
path_to_fig <- here::here(
  "Results",
  "ITS"
)

## Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

## Define plot theme
plot_theme(text_size = 10)


# Supplementary figure ITS: stacked bar and polar ----

p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl, ends_with("_prc")) %>%
  dplyr::mutate(CsCtrl = case_when(
    (CsCtrl == 1) ~ "RPL",
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
  dplyr::mutate(Clin_bar = 1.07)


# define species colours
mycols <- viridis::viridis(5,
  begin = 0.30,
  option = "mako"
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

## Create custom colours for case/control annotation
pal <- PNWColors::pnw_palette(name = "Moth", n = 5) %>%
  .[c(1, 5)] %>%
  purrr::set_names(unique(p1$CsCtrl))


# define sorting order based on top taxa abundances in samples
p1 <- p1 %>%
  dplyr::mutate(taxa_abd = case_when(
    (C.albicans_prc >= 0.50) ~ 1,
    TRUE ~ 2
  )) %>%
  dplyr::arrange(desc(CsCtrl), taxa_abd, SampleN) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN)))


# intiate and construct the plot
stack_bar <- ggplot2::ggplot(p1, aes(
  x = SampleN,
  y = Abundance,
  fill = taxa
)) +

  # add bar for cases vs control
  ggplot2::geom_bar(aes(
    y = Clin_bar,
    fill = CsCtrl
  ),
  stat = "identity",
  position = position_dodge()
  ) +

  # add manual fill
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::facet_wrap(. ~ Sample_grp, scales = "free_x") +

  # start new scale fill
  ggnewscale::new_scale_fill() +
  ggplot2::geom_bar(aes(fill = taxa),
    stat = "identity",
    position = position_stack(reverse = T)
  ) +
  ggplot2::scale_fill_manual(values = mycols, name = "Species") +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 14, face = "italic"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(hjust = -0.03)
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 3))

ggplot2::ggsave(stack_bar,
  path = path_to_fig,
  filename = "Supp_Toive_stacked_bar_ITS.pdf",
  device = "pdf",
  height = 15,
  width = 40,
  units = "cm",
  dpi = 320
)



#### ordination plots ----

# load data
TOIVE_phy_16S <- readr::read_rds(here::here(
  "Data",
  "Analysis",
  "3_TOIVE_phy_16S.rds"
))


# modify phyloseq object for plotting
TOIVE_pcoa <- TOIVE_phy_16S %>%
  phyloseq::transform_sample_counts(function(OTU) OTU / sum(OTU)) %>%
  phyloseq::subset_samples(!Sample_grp == "Blanks") %>%
  {
    ps <- .
    df <- data.frame(ps@sam_data) %>%
      dplyr::mutate(CsCtrl = ifelse(CsCtrl == 1,
        "RPL",
        "Control"
      ))
    phyloseq::sample_data(ps) <- df
    ps
  }


# set plot parameters
plot_theme(text_size = 10)

# define variables
col <- c("CsCtrl", "Sample_grp")
grid <- c("Sample_grp", "CsCtrl")

# plot
for (step in 1:length(col)) {

  # ordination
  pcoa <- phyloseq::ordinate(TOIVE_pcoa,
    method = "PCoA",
    distance = "euclidean"
  )

  # plot object
  plo <- phyloseq::plot_ordination(TOIVE_pcoa,
    pcoa,
    color = col[step]
  ) +
    ggplot2::facet_grid(as.formula(paste(". ~", grid[step]))) +
    ggplot2::geom_point(size = 5, alpha = 0.5) +
    ggplot2::stat_ellipse(type = "euclid", level = 0.95) +
    viridis::scale_colour_viridis(
      option = "G",
      begin = 0.2,
      end = 0.8,
      discrete = T
    )

  assign(paste("plo_",
    col[step],
    sep = ""
  ),
  plo,
  envir = .GlobalEnv
  )
}


grd_plot <- cowplot::plot_grid(plo_CsCtrl,
  plo_Sample_grp,
  labels = c("A)", "B)"),
  label_colour = "grey40"
)

# save pdf
ggsave(grd_plot,
  filename = here::here(
    "Results",
    "16S",
    "PCoA_euc.pdf"
  ),
  height = 5,
  width = 10,
  dpi = 320,
  device = "pdf"
)
