pkgs <- list(here::here("Scripts", 
                        "load_packages.R"),
             here::here("Scripts", 
                        "plot_theme.R"))
pkgs <- purrr::map(pkgs, source)



## Read in phyloseq object
## bacteria
Toive_phy_16S <- readRDS(here::here("Data", 
                                    "Toive_phy_16S.rds")) %>%
  phyloseq::subset_samples(Sample_grp == "Blanks")


# Convert raw counts to relative abundance
Toive_phy_16S_rel <- phyloseq::transform_sample_counts(
  Toive_phy_16S, 
  function(OTU) OTU/sum(OTU)
  ) %>%
  phyloseq::transform_sample_counts(
    function(OTU) 
      ifelse(is.na(OTU), 0, OTU)
    )



##### Rearrange data for plotting ##-----
## Extract ASV table and elongate
ASV_table <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(cols = !SampleN, 
                      names_to = "taxa", 
                      values_to = "Abundance") %>%
  dplyr::mutate(reads = "Relative") %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.05, "Others", taxa))

ASV_raw <- data.frame(Toive_phy_16S@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(cols = !SampleN,
                      names_to = "taxa", 
                      values_to = "Abundance_raw") %>%
  dplyr::mutate(taxa = ifelse(taxa %in% ASV_table$taxa, taxa, "Others"), 
                taxa = str_remove(taxa, "X."),
                taxa = str_replace(taxa, "\\.", " ")) 

## Extract meta table
meta_plot <- data.frame(Toive_phy_16S_rel@sam_data) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

## Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot) %>%
  as.data.frame()



##### define plot variables ##-----
path_to_fig <- here::here("Results", 
                          "16S")

## Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}



##### Figure 2: stacked bar plot ##-----

p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl, ends_with("_prc")) %>%
  dplyr::arrange(SampleN) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN)),
                taxa = str_remove(taxa, "X."),
                taxa = str_replace(taxa, "\\.", " "),
                across(starts_with(c("id")) | ends_with("_prc"), 
                       .fns = ~as.numeric(as.character(.x)))) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)))


## define species colours
mycols <- viridis::viridis(3, 
                           begin = 0.30,
                           end = 0.90,
                           option = "mako") %>%
  base::append(viridis::viridis(3, 
                                option = "rocket", 
                                begin = 0.3, 
                                end = 0.90, 
                                direction = -1)) %>%
  
  purrr::set_names(levels(p1$taxa)) %>%
  base::replace("Others", scico::scico(1, 
                                       palette = "grayC", 
                                       begin = 0.65)) %>%
  as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))


p2 <- ASV_raw %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(mycols) %>%
  dplyr::inner_join(p1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = factor(taxa, levels = c(unique(p1$taxa)))) %>%
  dplyr::arrange(taxa) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  as.data.frame()



## Define plot theme
plot_theme(text_size = 12)
plot_colour <- "grey50"

# initiate and construct the plot
stack_bar <- ggplot2::ggplot(p1, aes(x = SampleN, 
                                     y = Abundance, 
                                     fill = taxa)) +
  
  
  # add bar for species
  ggplot2::geom_bar(stat = "identity", 
                    position = position_stack(reverse = T)) +
  
  
  ggplot2::scale_fill_manual(values = unique(p2$colours), name = "Species") +
  
  # Add reads within the bars
  ggplot2::geom_text(data = p2 %>%
                       dplyr::filter(Abundance > 0.05), 
                     aes(x = SampleN, 
                         y = Abundance, 
                         label = Abundance_raw,
                         size = Abundance_raw, 
                         group = taxa,
                         angle = ifelse(Abundance_raw > 500, 0, 270),
                         colour = shades),
                     position = position_stack(0.5, reverse = T),
                     inherit.aes = F,
                     show.legend = F) +
  
  # add discrete text colour based on background shades
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), name = "Species") +
  
  ggplot2::theme(legend.position = "bottom", 
                 strip.text.x = element_text(colour = plot_colour, 
                                             size = 14, 
                                             face = "italic", hjust = 0),
                 axis.text.x = element_text(colour = plot_colour, size = 6, angle = 270),
                 axis.title.x = element_text(hjust = -0.03)) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 2)) 

ggplot2::ggsave(stack_bar, 
                path = path_to_fig,
                filename = "Suppl_Blanks.pdf",
                device = "pdf", 
                height = 12, 
                width = 25,
                units = "cm",
                dpi = 320)

