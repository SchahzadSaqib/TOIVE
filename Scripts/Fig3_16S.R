# Load packages ----
pkgs <- list.files(here::here(
  "Scripts",
  "to_load"
),
full.names = TRUE
)

pkgs <- purrr::map(pkgs, source)

grp <- here::here(
  "Scripts",
  "GroupTest_mod.R"
) %>%
  purrr::map(source)


# read and process data ----

Toive_phy_16S <- readRDS(here::here(
  "Data",
  "Analysis",
  "3_Toive_phy_16S.rds"
)) %>%
  phyloseq::subset_samples(!Sample_grp == "Blanks")



# Convert raw counts to relative abundance
Toive_phy_16S_rel <- phyloseq::transform_sample_counts(
  Toive_phy_16S,
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) {
      ifelse(is.na(OTU), 0, OTU)
    }
  )


## Extract L.crispatus abundance in each sample
L.crispatus_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.crispatus")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel
)$L.crispatus_prc <- L.crispatus_prc$Lactobacillus.crispatus


## Extract L.iners abundance in each sample
L.iners_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.iners")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel
)$L.iners_prc <- L.iners_prc$Lactobacillus.iners


## Extract L.crispatus abundance in each sample
G.vag_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis")

## Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel
)$G.vag_prc <- G.vag_prc$Gardnerella.vaginalis



# Rearrange data for plotting ----
## Extract ASV table and elongate
ASV_table <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(
    cols = !SampleN,
    names_to = "taxa",
    values_to = "Abundance"
  )

## Extract meta table
meta_plot <- data.frame(Toive_phy_16S_rel@sam_data) %>%
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


# Figure 3: stacked bar + grouptest infographic -----

sample_tally <- plot_data %>%
  dplyr::distinct(SampleN, .keep_all = T) %>%
  dplyr::group_by(CsCtrl, Sample_grp) %>%
  dplyr::count() %>%
  dplyr::mutate(CsCtrl = if_else(CsCtrl == 1, "RPL", "Control"))

p1 <- plot_data %>%
  dplyr::mutate(taxa = ifelse(taxa == "uncultured.bacterium",
    "uncultured.bacteria",
    taxa
  )) %>%
  dplyr::select(SampleN:CsCtrl) %>%
  dplyr::filter(!CsCtrl == 0) %>%
  dplyr::mutate(CsCtrl = if_else(CsCtrl == 1, "RPL", "Control")) %>%
  dplyr::mutate(Sample_grp = as.character(Sample_grp)) %>%
  dplyr::group_by(taxa, CsCtrl, Sample_grp) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.005, "Others", taxa)) %>%
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  as.data.frame()


# define species colours
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
    n = 3,
    name = "Mushroom"
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


p2 <- p1 %>%
  inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  as.data.frame()


## Define plot theme
plot_theme(text_size = 13)
plot_colour <- "grey50"

stack_bar_summ <- ggplot2::ggplot(p2, aes(
  x = CsCtrl,
  y = Abundance,
  fill = taxa
)) +

  # add stackedd bar
  ggplot2::geom_bar(
    stat = "identity",
    position = position_stack(reverse = T)
  ) +

  # add manual colour to bars
  ggplot2::scale_fill_manual(values = unique(p2$colours), name = "Species") +

  # Add percentage within the bars
  ggplot2::geom_text(
    data = p2,
    aes(
      x = CsCtrl,
      y = Abundance,
      label = if_else(Abundance > 0.02,
        paste(
          round(Abundance * 100, digits = 2),
          "%"
        ),
        ""
      ),
      size = Abundance,
      group = taxa,
      angle = ifelse(Abundance > 0.10, 0, 270),
      colour = shades
    ),
    position = position_stack(0.5, reverse = T),
    inherit.aes = F,
    show.legend = F
  ) +

  # add discrete text colour based on background shades
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), name = "Species") +

  # add sample tally
  ggplot2::geom_text(
    data = sample_tally,
    aes(
      x = CsCtrl,
      y = 1.06,
      label = paste("n = ", n, sep = "")
    ),
    colour = plot_colour,
    inherit.aes = F
  ) +

  # facet wrap between cases and control
  ggplot2::facet_wrap(. ~ Sample_grp) +

  # convert y-axis scale to percentage
  scale_y_continuous(labels = scales::percent) +


  # plot parameters
  ggplot2::coord_flip() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.x = element_text(
      colour = plot_colour,
      hjust = 0.5
    ),
    axis.text.x = element_text(angle = 270)
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 5))


gg_objs <- list()
gg_objs <- append(gg_objs, list(stack_bar_summ))


## define plot variables
tidyvar <- sym("CsCtrl")
var <- "CsCtrl"
assign_names <- list("Pooled", "Vagina", "Endometrium")
smpl_subT <- list(NULL, "Sample_grp", "Sample_grp")
smpl_sub <- list(NULL, "Vagina", "Endometrium")

confounder_list <- list(c(
  "bmi3",
  "ageyear_strat",
  "p01"
))
table_comb <- c()
viobox_comb <- c()


for (cdl in 1:length(confounder_list)) {
  print(confounder_list[[cdl]])
  for (fct in 1:length(assign_names)) {
    print(assign_names[[fct]])
    print(smpl_subT[[fct]])
    print(smpl_sub[[fct]])

    ## Grouptest
    gt_res <- try({
      GroupTest_mod(
        species.table = Toive_phy_16S,
        meta = Toive_phy_16S,
        group = var,
        group_name = "Case (1) vs Control (2)",
        compare.to = "2",
        dir_for_res = here::here(
          "Results",
          "16S",
          "mare",
          "Case_vs_Control"
        ),
        confounders = confounder_list[[cdl]],
        min.prevalence = 0.10,
        min.abundance = 0.05,
        select.by = smpl_subT[[fct]],
        select = smpl_sub[[fct]],
        p.cutoff = 0.05,
        keep.result = T,
        nonzero = F,
        pdf = T,
        show_quartz = F
      )
    })

    if (class(gt_res) == "try-error") {
      gt_res <- c()
    }

    pavg <- plot_data %>%
      dplyr::group_by(taxa, !!tidyvar) %>%
      dplyr::summarise(Abundance = mean(Abundance))

    gt_summ_cmp <- gt_res %>%
      dplyr::select(taxon, model, contains("FDR")) %>%
      tidyr::pivot_longer(
        cols = contains("FDR"),
        names_to = "comp",
        values_to = "pvals"
      ) %>%
      dplyr::group_by(taxon) %>%
      dplyr::filter(pvals <= 0.1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(variable = str_extract(
        comp,
        pattern = paste(levels(factor(plot_data[, var])),
          collapse = "|"
        )
      )) %>%
      dplyr::group_by(taxon, variable) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sig = case_when(
        (pvals > 0 & pvals <= 0.001) ~ "***",
        (pvals > 0.001 & pvals <= 0.01) ~ "**",
        (pvals > 0.01 & pvals <= 0.05) ~ "*",
        (pvals > 0.05 & pvals <= 0.1) ~ "."
      )) %>%
      dplyr::rename(!!tidyvar := variable) %>%
      dplyr::rename("taxa" = taxon) %>%
      dplyr::mutate(comp = assign_names[[fct]])

    table_comb <- table_comb %>%
      dplyr::bind_rows(gt_summ_cmp) %>%
      dplyr::filter(!is.na(sig))
  }
}


## Clean and rearrange data
p2 <- plot_data %>%
  dplyr::select(SampleN, taxa, Abundance, !!tidyvar) %>%
  dplyr::filter(taxa %in% table_comb$taxa) %>%
  dplyr::mutate(split = ifelse(str_ends(SampleN, "A"),
    "Vagina",
    "Endometrium"
  )) %>%
  {
    dplyr::bind_rows(., dplyr::mutate(., split = "Pooled"))
  } %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(CsCtrl = if_else(CsCtrl == 1, "RPL", "Control")) %>%
  dplyr::mutate(CsCtrl = factor(CsCtrl, levels = c("RPL", "Control"))) %>%
  dplyr::mutate(split = factor(split,
    levels = c(
      "Pooled",
      "Vagina",
      "Endometrium"
    )
  ))

## assign colours for plot
pal <- PNWColors::pnw_palette("Sunset", 7) %>%
  .[c(2, 5)] %>%
  purrr::set_names(levels(p2$CsCtrl))

## significance labels
p2_lab <- table_comb %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::rename("split" = comp) %>%
  dplyr::group_by(taxa, split) %>%
  dplyr::distinct(split, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(CsCtrl = "RPL")


## construct viobox
gt_viobox <- ggplot2::ggplot(data = p2, aes(
  x = CsCtrl,
  y = Abundance,
  fill = CsCtrl
)) +

  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = CsCtrl),
    show.legend = F,
    size = 3,
    alpha = 0.8
  ) +

  ## half violin on the left
  gghalves::geom_half_violin(
    mapping = aes(colour = CsCtrl),
    trim = F,
    side = "l",
    draw_quantiles = c(0.25, 0.5, 0.75),
    scale = "area",
    colour = "grey70",
    size = .3,
    alpha = 0.8
  ) +

  ## half box on the right
  gghalves::geom_half_boxplot(
    mapping = aes(colour = CsCtrl),
    side = "r",
    nudge = .05,
    size = .3,
    colour = "grey70",
    outlier.shape = NA,
    alpha = 0.8,
    show.legend = F
  ) +

  ## facet by taxa
  ggplot2::facet_grid(factor(split,
    levels = c(
      "Endometrium",
      "Vagina",
      "Pooled"
    )
  )
  ~ factor(taxa, levels = c(
      "Lactobacillus crispatus",
      "Lactobacillus jensenii",
      "Gardnerella vaginalis",
      "Gardnerella leopoldii"
    )),
  scales = "free"
  ) +

  ## add asterisks labels for significant taxa)
  ggplot2::geom_text(
    data = p2_lab,
    aes(
      x = CsCtrl,
      y = 30,
      label = sig
    ),
    colour = "skyblue2",
    size = 8,
    nudge_x = 0.15,
    inherit.aes = F,
    show.legend = F
  ) +

  ## add colours
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal) +

  ## start new scale
  ggnewscale::new_scale_colour() +

  ## log10 transform y-axis
  ggplot2::scale_y_log10() +

  ## add point for median
  ggplot2::stat_summary(
    fun = "mean",
    geom = "point",
    colour = "red3",
    size = 3.5,
    show.legend = F
  ) +

  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    legend.key.size = unit(1, "cm")
  ) +
  ggplot2::guides(fill = guide_legend(ncol = 4))


gg_objs <- append(gg_objs, list(gt_viobox))


layout_des <- rbind(
  c(1),
  c(1),
  c(2),
  c(2),
  c(2),
  c(2)
)


g <- gridExtra::arrangeGrob(
  grobs = gg_objs,
  layout_matrix = layout_des,
  # top = grid::textGrob(paste("Case vs Controls"),
  #                     gp=grid::gpar(fontsize=20,
  #                                   col = plot_colour)),
  padding = unit(1.5, "cm"),

  # add margins around the plot
  bottom = "\n",
  right = "\n",
  left = "\n"
)
p <- ggpubr::as_ggplot(g) +
  cowplot::draw_plot_label(
    label = c(
      "A)",
      "B)",
      "C)",
      "D)",
      "E)"
    ),
    x = c(
      0.03,
      0.52,
      0.03,
      0.03,
      0.03
    ),
    y = c(
      0.99,
      0.99,
      0.60,
      0.42,
      0.27
    ),
    size = 15,
    colour = "grey50"
  )


ggplot2::ggsave(
  path = path_to_fig,
  plot = p,
  filename = paste("Fig3_Toive",
    ".pdf",
    sep = ""
  ),
  device = "pdf",
  dpi = 320,
  height = 34,
  width = 38,
  units = "cm",
  limitsize = F
)
