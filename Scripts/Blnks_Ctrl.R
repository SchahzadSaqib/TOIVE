pkgs <- list.files(here::here("Scripts", 
                              "to_load"), 
                   full.names = TRUE)
pkgs <- purrr::map(pkgs, source)



## Read in phyloseq object
## bacteria
Toive_phy_16S <- readRDS(
  here::here("Data", 
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
mycols <- pick_colrs(p1$taxa) %>%
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
                 axis.title.x = element_text(hjust = -0.05)) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 2)) 

ggplot2::ggsave(stack_bar, 
                path = path_to_fig,
                filename = "Suppl_Blanks.pdf",
                device = "pdf", 
                height = 8, 
                width = 16,
                units = "cm",
                dpi = 320)


##### positive controls -----
EMV_cont_phy <- readRDS(
  here::here(
    "Data", 
    "2_EMV_cont_phy.rds"))

# Convert raw counts to relative abundance
EMV_cont_phy.trf <- phyloseq::transform_sample_counts(
  EMV_cont_phy, 
  function(OTU) OTU/sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) 
      ifelse(is.na(OTU), 
             0, 
             OTU))


##### Rearrange data for plotting ##-----
meta <- data.frame(EMV_cont_phy.trf@sam_data) %>%
  tibble::rownames_to_column(var = "ID")


## Extract ASV table and elongate
plot_data <- data.frame(EMV_cont_phy.trf@otu_table) %>%
  tibble::rownames_to_column(var = "ID") %>%
  dplyr::inner_join(meta) %>% 
  tidyr::pivot_longer(cols = !c(ID, ReadCount), 
                      names_to = "taxa",
                      values_to = "Abundance") %>% 
  dplyr::mutate(group = ifelse(
    str_ends(ID, "_2"),
    "Earlier (as comparison)",
    "Positive controls"), 
    ID = str_remove(ID, "_2"))


##### define plot variables ##-----
## Define path to figures
path_to_fig <- here::here(
  "Results",
  "16S")

## Path to figures
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}


##### Define plot variables ##-----

# Plot aesthetics
plot_colour <- "grey50"
text_size <- 14
plot_theme(text_size = text_size)

## Clean and rearrange data
p1 <- plot_data %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.02, "Others", taxa)) %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  as.data.frame() 

## Define species colours
mycols <- pick_colrs(p1$taxa) %>%
  as.data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  ## appropriate contrast and colour for text
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))

p1 <- p1 %>%
  dplyr::inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  as.data.frame() 

read_labels <- p1 %>%
  dplyr::group_by(group) %>% 
  dplyr::distinct(ID, .keep_all = T)

bar_plot <- ggplot2::ggplot(p1, 
                            aes(x = ID, 
                                y = Abundance, 
                                fill = taxa)) +
  
  ## facet 
  ggplot2::facet_grid(. ~ group) +
  
  ## add bar plot
  ggplot2::geom_bar(stat = "identity", 
                    position = position_stack(reverse = T), 
                    width = 1) + 
  
  ggplot2::coord_flip() +
  
  ## Add label for sample tally  
  ggplot2::geom_text(data = read_labels, 
                     aes(x = ID,
                         y = 1,
                         label = ReadCount, 
                         colour = shades), 
                     angle = 0, 
                     size = text_size*0.42,
                     position = position_stack(vjust = 0.5),
                     inherit.aes = FALSE, 
                     show.legend = FALSE) +
  
  
  ## assign manual set of colours
  ggplot2::scale_fill_manual(values = mycols$colours, 
                             name = "Species") +
  
  ## scale colour for labels
  ggplot2::scale_colour_manual(values = c("grey85", "grey35"), 
                               name = "Species") +
  
  ## convert y-axis scale to percentage
  ggplot2::scale_y_continuous(labels = scales::percent, 
                              limits = c(0, 1.1), 
                              expand = c(0,0)) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::theme(legend.position = "bottom", 
                 axis.text.x = element_text(angle = 90), 
                 axis.title.x = element_text(hjust = -0.05)) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::guides(fill = guide_legend(ncol = 3))

ggplot2::ggsave(bar_plot,
                filename = "Suppl_PstvCtrls.pdf", 
                device = "pdf", 
                path = path_to_fig, 
                width = 10, 
                height = 10,
                units = "in",
                dpi = 320)


grid <- cowplot::plot_grid(stack_bar, 
                   bar_plot, 
                   ncol = 1, rel_heights = c(0.5, 1), 
                   labels = c("A)", "B)"), 
                   label_size = 16,
                   label_colour = "grey50")

ggplot2::ggsave(
  grid, 
  filename = "FigS3_Blnk_Ctrl.pdf",
  path = here::here(
    "Results", 
    "16S"
  ),
  device = "pdf", 
  width = 10, 
  height = 10,
  units = "in",
  dpi = 320)


##### positive controls correlation -----
cor_dat_p <- plot_data %>% 
  dplyr::filter(group == "Positive controls") %>% 
  tidyr::pivot_wider(id_cols = ID, 
                     names_from = taxa, 
                     values_from = Abundance) %>% 
  tibble::column_to_rownames(var = "ID") %>% 
  t()

cor_dat_e <- plot_data %>% 
  dplyr::filter(group == "Earlier (as comparison)") %>% 
  tidyr::pivot_wider(id_cols = ID, 
                     names_from = taxa, 
                     values_from = Abundance) %>% 
  tibble::column_to_rownames(var = "ID") %>% 
  t()


cor(cor_dat_e, 
    cor_dat_p, 
    method = "pearson")

cor.test(cor_dat_p, 
         cor_dat_e, 
         method = "pearson")

cor.test(cor_dat_p, 
         cor_dat_e, 
         method = "spearman")
