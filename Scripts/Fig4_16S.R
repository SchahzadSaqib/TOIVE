#' load packages and scripts
pkgs <- list(here::here("Scripts", 
                        "load_packages.R"),
             here::here("Scripts", 
                        "plot_theme.R"))
pkgs <- purrr::map(pkgs, source)



#' Read in phyloseq object
#' bacteria
Toive_phy_16S <- readRDS(here::here("Data", 
                                    "Toive_phy_16S.rds")) %>%
  phyloseq::subset_samples(!Sample_grp == "Blanks")


#' Convert raw counts to relative abundance
Toive_phy_16S_rel <- phyloseq::transform_sample_counts(
  Toive_phy_16S, 
  function(OTU) OTU/sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) 
      ifelse(is.na(OTU), 0, OTU)
  )


#' Extract L.crispatus abundance in each sample 
L.crispatus_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.crispatus")

#' Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel)$L.crispatus_prc <- L.crispatus_prc$Lactobacillus.crispatus


#' Extract L.iners abundance in each sample 
L.iners_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Lactobacillus.iners")

#' Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel)$L.iners_prc <- L.iners_prc$Lactobacillus.iners


#' Extract L.crispatus abundance in each sample 
G.vag_prc <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  dplyr::select("Gardnerella.vaginalis")

#' Add update L.crispatus abundance to metadata
phyloseq::sample_data(
  Toive_phy_16S_rel)$G.vag_prc <- G.vag_prc$Gardnerella.vaginalis



##### Rearrange data for plotting ##-----
#' Extract ASV table and elongate
ASV_table <- data.frame(Toive_phy_16S_rel@otu_table) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  tidyr::pivot_longer(cols = !SampleN, 
                      names_to = "taxa", 
                      values_to = "Abundance")

#' Extract meta table
meta_plot <- data.frame(Toive_phy_16S_rel@sam_data) %>%
  tibble::rownames_to_column(var = "SampleN") %>%
  dplyr::mutate(across(.cols = everything(), .fns = ~as.factor(.x)))

#' Combine
plot_data <- ASV_table %>%
  dplyr::inner_join(meta_plot) %>%
  dplyr::filter(!CsCtrl == 0) %>%
  as.data.frame()



##### define plot variables ##-----
path_to_fig <- here::here("Results", 
                          "16S")

#' Check or create directory
if (!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}






##### Figure 4: Polar plot + correlations ##-----

sample_tally <- plot_data %>%
  dplyr::distinct(SampleN, .keep_all = T) %>%
  dplyr::group_by(id) %>%
  dplyr::filter(n()>1)


p1 <- plot_data %>%
  dplyr::select(SampleN:CsCtrl, ends_with("_prc")) %>%
  dplyr::filter(!CsCtrl == 0) %>%
  dplyr::filter(id %in% sample_tally$id) %>%
  dplyr::mutate(CsCtrl = case_when((CsCtrl == 1) ~ "RPL",
                                   (CsCtrl == 2) ~ "Control")) %>%
  dplyr::mutate(taxa = if_else(Abundance < 0.05, "Others", taxa)) %>%
  dplyr::arrange(SampleN) %>%
  dplyr::mutate(SampleN = factor(SampleN, levels = unique(SampleN))) %>%
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::mutate(across(starts_with(c("id")) | ends_with("_prc"), 
                       .fns = ~as.numeric(as.character(.x)))) %>%
  dplyr::group_by(taxa) %>%
  dplyr::add_tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = ifelse(taxa == "Others", 0, n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::mutate(value = 1.07)


#' define species colours
mycols <- viridis::viridis(5, 
                           begin = 0.30,
                           end = 0.90,
                           option = "mako") %>%
  base::append(viridis::viridis(5, 
                                option = "rocket", 
                                begin = 0.3, 
                                end = 0.90, 
                                direction = -1)) %>%
  base::append(PNWColors::pnw_palette(n = 5, 
                                      name = "Sailboat")) %>%
  base::append(PNWColors::pnw_palette(n = 5, 
                                      name = "Winter")) %>%
  base::append(PNWColors::pnw_palette(n = 5, 
                                      name = "Mushroom")) %>%
  base::append(PNWColors::pnw_palette(n = nlevels(p1$taxa) - 25, 
                                      name = "Cascades")) %>%
  
  purrr::set_names(levels((p1$taxa))) %>%
  base::replace("Others", scico::scico(1, 
                                       palette = "grayC", 
                                       begin = 0.65)) 


pal <- PNWColors::pnw_palette(name = "Moth", n = 5) %>%
  .[c(1,5)] %>%
  purrr::set_names(unique(p1$Sample_grp))


p1 <- p1 %>%
  dplyr::arrange(Sample_grp, SampleN) 

plot_theme(text_size = 12)

#' Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nbreaks <- nlevels(p1$Sample_grp)
nObs <- nlevels(factor(plot_data$taxa))
to_add <- data.frame(matrix(NA, empty_bar*nbreaks*nObs, ncol(p1)) )
colnames(to_add) <- colnames(p1)
to_add$Sample_grp <- rep(levels(p1$Sample_grp), each=empty_bar*nObs)

endo_ord <- p1 %>%
  dplyr::filter(Sample_grp == "Endometrium") %>%
  dplyr::mutate(taxa_abd = case_when((L.crispatus_prc >= 0.50) ~ 1,
                                     (L.iners_prc >= 0.50) ~ 2,
                                     (G.vag_prc >= 0.50) ~ 3,
                                     TRUE ~ 4)) %>%
  dplyr::arrange(taxa_abd, SampleN) %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::filter(Sample_grp == "Endometrium")

vg_ord <- p1 %>%
  dplyr::filter(Sample_grp == "Vagina") %>%
  dplyr::arrange(match(id, rev(endo_ord$id))) %>%
  dplyr::bind_rows(to_add) %>%
  dplyr::filter(Sample_grp == "Vagina")

pcmpl <- endo_ord %>%
  dplyr::bind_rows(vg_ord) %>%
  dplyr::mutate(labels = rep(seq(1, 
                                 nlevels(SampleN) + empty_bar*nbreaks), 
                             each=nObs)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SampleN = as.character(SampleN)) %>%
  dplyr::mutate(labels = factor(labels, levels = unique(labels))) %>%
  dplyr::mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  dplyr::mutate(Sample_grp = factor(Sample_grp, levels = unique(Sample_grp)))

#' Get the name and the y position of each label
label_data <- pcmpl %>% 
  dplyr::group_by(labels, SampleN, Sample_grp) %>% 
  dplyr::summarise(tot=sum(Abundance)) %>%
  dplyr::mutate(labels = as.numeric(as.character(labels)))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$labels-0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


#' Make the plot
polar_stack <- ggplot2::ggplot(pcmpl) +      
  
  #' add bar for gestational age
  ggplot2::geom_bar(aes(x = labels, 
               y = value, 
               fill = Sample_grp),
           stat = "identity", position = "dodge") +
  
  ggplot2::scale_fill_manual(values = pal, name = "") +
  
  #' start new scale
  ggnewscale::new_scale_fill() +
  
  #' add stacked bar for taxa
  ggplot2::geom_bar(aes(x = labels, 
               y = Abundance, 
               fill = taxa), 
           stat="identity", 
           position = position_stack(reverse = T)) +
  ggplot2::scale_fill_manual(values = mycols, name = "Species") +
  
  
  ggplot2::theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  
  #' Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = max(as.numeric(as.character(pcmpl$labels)))+3, 
                    y = c(0, .25, .50, .75, 1.0), 
                    label = c("0", "25", "50", "75", "100") , 
                    color="grey60", 
                    size=2 , 
                    angle=0, 
                    fontface="bold", 
                    hjust=1) +
  
  ggplot2::guides(fill = guide_legend(ncol = 2, order = 1)) +
  
  ggplot2::ylim(-.30,1.33) + #' control inner circle radius
  ggplot2::coord_polar() +
  ggplot2::scale_x_discrete() +
  
  #' Add labels on top of each bar
  ggplot2::geom_text(data=label_data, 
            aes(x = labels, 
                y = 1.1, 
                label = SampleN, 
                hjust = hjust), 
            color = "grey60", 
            size = 1.5, 
            angle = label_data$angle, 
            inherit.aes = FALSE)



##### Correlations ##-----

sample_tally <- plot_data %>%
  dplyr::distinct(SampleN, .keep_all = T) %>%
  dplyr::group_by(id) %>%
  dplyr::filter(n()>1)

#' vaginal samples
vg_sample <- plot_data %>%
  dplyr::filter(id %in% sample_tally$id) %>%
  dplyr::select(id, taxa, Abundance, Sample_grp) %>%
  dplyr::filter(Sample_grp == "Vagina") %>%
  dplyr::rename("Vagina" = Abundance) %>%
  dplyr::select(1:3) %>%
  tidyr::pivot_wider(names_from = id, values_from = Vagina) %>%
  tibble::column_to_rownames(var = "taxa") %>%
  dplyr::rename_all(~paste(.x, "A", sep = "")) %>%
  as.matrix()

#' endometrial samples
ed_sample <- plot_data %>%
  dplyr::filter(id %in% sample_tally$id) %>%
  dplyr::select(id, taxa, Abundance, Sample_grp) %>%
  dplyr::filter(Sample_grp == "Endometrium") %>%
  dplyr::rename("Endometrium" = Abundance) %>%
  dplyr::select(1:3) %>%
  tidyr::pivot_wider(names_from = id, values_from = Endometrium) %>%
  tibble::column_to_rownames(var = "taxa") %>%
  dplyr::rename_all(~paste(.x, "D", sep = "")) %>%
  as.matrix()


full_corr <- stats::cor.test(vg_sample, 
                             ed_sample, 
                             method = "pearson")

vg_endo_corr <- Hmisc::rcorr(vg_sample, 
                             ed_sample, 
                             type = "pearson")

#' correlations plot table
vg_ed_corr <- vg_endo_corr$r %>%
  as.data.frame() %>%
  dplyr::select(ends_with("A")) %>%
  dplyr::filter(str_ends(rownames(.), "D")) %>%
  dplyr::mutate(Corr = diag(as.matrix(.))) %>%
  dplyr::select(Corr, everything()) %>%
  dplyr::mutate(vg = colnames(.)[2:length(.)]) %>%
  tibble::rownames_to_column(var = "endo") %>%
  dplyr::select(vg, endo, Corr) %>%
  dplyr::mutate(cs_ctrl = if_else(str_starts(endo, "1"), "RPL", "Control"))


gg_objs <- c()

#' lollipop plot
vg_ed_corr_ord <- vg_ed_corr %>%
  dplyr::mutate(vg = str_replace(vg, "A", "")) %>%
  dplyr::mutate(vg = factor(vg, levels = unique(pcmpl$id))) %>%
  dplyr::mutate(cs_ctrl = factor(cs_ctrl, levels = c("RPL", "Control")))

p <- ggplot2::ggplot(vg_ed_corr_ord, 
                     aes(x = vg, 
                         y = Corr, 
                         colour = cs_ctrl)) +
  ggplot2::geom_segment(aes(x = vg, 
                            xend = vg, 
                            y = -Inf, 
                            yend = Corr, 
                            colour = cs_ctrl)) +
  ggplot2::geom_point(size = 3) +
  viridis::scale_colour_viridis(option = "rocket", 
                                discrete = T, 
                                begin = 0.25, 
                                end = 0.75, 
                                alpha = 0.8,
                                name = "") +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Correlation (Pearson)") +
  ggplot2::theme(legend.justification = c(1,-0.7))


gg_objs <- base::append(gg_objs, list(p))

gg_objs <- base::append(gg_objs, list(polar_stack))



#' create grid plot
layout_des <- rbind(c(2),
                    c(2),
                    c(2),
                    c(1))


g <- gridExtra::arrangeGrob(grobs = gg_objs,
                            layout_matrix = layout_des,
                            padding = unit(4, "cm"))


#' save pdf
ggsave(path = here::here(path_to_fig),
       plot = g,
       filename = "Fig4_16S.pdf", 
       device = "pdf", 
       dpi = 320, 
       height = 20, 
       width = 30,
       units = "cm", 
       limitsize = F)

