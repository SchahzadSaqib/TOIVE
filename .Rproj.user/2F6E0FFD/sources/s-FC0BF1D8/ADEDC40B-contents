## Load packages
pkgs <- list(here::here("Scripts", 
                        "load_packages.R"),
             here::here("Scripts", 
                        "plot_theme.R"))
pkgs <- purrr::map(pkgs, source)


## Read in phyloseq object
## Bacteria
TOIVE_phy_bac <- readRDS(here::here("Data", 
                                    "TOIVE_phy_16S.rds"))

# Convert raw counts to relative abundance
TOIVE_phy_bac_rel <- phyloseq::transform_sample_counts(
  TOIVE_phy_bac, 
  function(OTU) OTU/sum(OTU)) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
    )



##### adonis helpers ##-----
#' permutations
permt_to_use <- 99999

#' cleaner
adonis_clean <- function(x) {
  x <- x %>%
    tibble::rownames_to_column(var = "Var") %>%
    dplyr::rename_with(~paste("Pr_F_stat"), starts_with("Pr")) %>%
    dplyr::rename_with(~paste("F_statistic"), starts_with("F")) %>%
    dplyr::mutate(
      Sigif = case_when((Pr_F_stat > 0 & Pr_F_stat <= 0.001) ~ "***",
                        (Pr_F_stat > 0.001 & Pr_F_stat <= 0.01) ~ "**",
                        (Pr_F_stat > 0.01 & Pr_F_stat <= 0.05) ~ "*",
                        (Pr_F_stat > 0.05 & Pr_F_stat <= 0.1) ~ "."),
      Sigif = tidyr::replace_na(.data$Sigif, replace = ""),
      Pr_F_stat = as.numeric(.data$Pr_F_stat)) %>%
    dplyr::filter(Pr_F_stat > 0) %>%
    dplyr::arrange(Pr_F_stat)
}


##### Summary of metadata variables ##----------------------------------------------------------------------------
# view(summarytools::dfSummary(metadata_bac_grpt, 
# plain.ascii = FALSE, 
# style = "grid"))

clin_tab <- readxl::read_xlsx(
  path = here::here("Data", 
                    "TOIVE_clinical_and microbiome_data.xlsx"))


## Define path to figures
path_to_fig <- here::here("Results", 
                          "16S", 
                          "adonis")

## Check or create directory
if(!dir.exists(here::here(path_to_fig))) {
  dir.create(here::here(path_to_fig), recursive = T)
}

## Define plot theme
plot_theme(text_size = 8)


##### Bacteria ##--------------------------------------------------------------------------------------------------
groups <- c("A|D", "A", "D")
csctrl_sub <- c("V1|V2", "V1", "V2")
label_grp <- c("all", "vg", "ed")
label_sub <- c("pool", "cs", "ctrl")
plot_title <- c("All samples", "Vaginal samples", "Endometrial samples")

grp <- 1
step <- 1
for (grp in 1:length(groups)) {
  print(groups[grp])
  print(label_grp[grp])
  grp_table <- c()
  
  for (step in 1:length(csctrl_sub)) {
    print(csctrl_sub[step])
    print(label_sub[step])
    
    ##### Vaginal sample subset bacteria group 1 ##-------------------------------------------------------------------
    data_step <- data.frame(TOIVE_phy_bac_rel@otu_table) %>%
      dplyr::filter(str_ends(rownames(.), groups[grp]) & str_detect(rownames(.), csctrl_sub[step])) %>%
      dplyr::filter(sum(.) != 0) 
    print(rownames(data_step))
    
    # subset metadata
    meta_step <- data.frame(TOIVE_phy_bac_rel@sam_data) %>%
      dplyr::filter(str_ends(rownames(.), groups[grp]) & str_detect(rownames(.), csctrl_sub[step])) %>%
      droplevels() %>%
      dplyr::select(where(
        ~all(nlevels(.x) > 1) | is.numeric(.x)
      ))
    
    # Run adonis
    adonis_out <- data.frame()
    
    for (i in 1:length(meta_step)) {
      set.seed(26)
      
      # display the variable being tested 
      print(colnames(meta_step[i]))
      
      # remove NAs from input data
      data_in <- data_step %>%
        filter(rownames(.) %in% rownames(na.omit(meta_step[i])))
      
      # run adonis
      adonis_hold <- adonis2(formula = data_in ~ na.omit(meta_step[,i]), 
                             data = meta_step, 
                             method = "bray", permutations = permt_to_use)
      rownames(adonis_hold)[1] <- colnames(meta_step[i])
      adonis_out <- rbind(adonis_out, adonis_hold)
    }
    
    # clean output
    adonis_scrub <- adonis_clean(adonis_out) %>%
      dplyr::mutate(comp = label_sub[step])
    
    # write excel file 
    writexl::write_xlsx(adonis_scrub, 
                        path = here::here("Results", 
                                          "16S", 
                                          paste("adonis_toive_bac_", 
                                                label_grp[grp], "_",
                                                label_sub[step],
                                                ".xlsx", 
                                                sep = "")))
    
    grp_table <- grp_table %>%
      bind_rows(adonis_scrub)
  }
  
  clev_plotdata <- grp_table %>%
    dplyr::arrange(Var) %>%
    dplyr::select(Var, Pr_F_stat, comp) %>%
    tidyr::pivot_wider(names_from = comp, values_from = Pr_F_stat) %>%
    dplyr::mutate(across(pool:ctrl, ~ifelse(is.na(.x), 0, .x))) %>%
    dplyr::mutate(Max = pmax(!!!rlang::syms(names(.)[2:4]))) %>%
    dplyr::mutate(across(pool:ctrl, ~ifelse(.x == 0, 10, .x))) %>%
    dplyr::mutate(Min = pmin(!!!rlang::syms(names(.)[2:4]))) %>%
    dplyr::mutate(across(pool:ctrl, ~ifelse(.x == 10, NA, .x))) %>%
    tidyr::pivot_longer(cols = pool:ctrl, names_to = "comp", values_to = "Pr_F_stat") %>%
    dplyr::arrange(desc(Min)) %>%
    dplyr::inner_join(clin_tab) %>%
    dplyr::mutate(Group = factor(Group, levels = c("Outcome", 
                                                   "Lifestyle", 
                                                   "Clinical/medical", 
                                                   "Pregnancies", 
                                                   "Technical"))) %>%
    dplyr::mutate(comp = case_when((comp == "pool") ~ "Pooled",
                                   (comp == "cs") ~ "RPL",
                                   (comp == "ctrl") ~ "Control")) %>%
    dplyr::mutate(comp = factor(comp, levels = c("RPL", "Control", "Pooled"))) %>%
    dplyr::mutate(Name = factor(Name, levels = unique(Name)))
  
  annot_text <- clev_plotdata %>%
    ungroup() %>%
    dplyr::filter(Group == "Outcome") %>%
    dplyr::slice_tail() %>%
    droplevels()
  
  p <- ggplot2::ggplot(clev_plotdata, 
                       aes(x = Pr_F_stat, 
                           y = Name)) +
    ggplot2::geom_segment(aes(x = Min, 
                              xend = Max, 
                              y = Name, 
                              yend = Name), 
                          colour = "grey80") +
    ggplot2::geom_jitter(aes(colour = comp), 
                         size = 4, 
                         height = 0.1, 
                         width = 0.01) +
    ggplot2::geom_vline(xintercept = 0.05, 
                        colour = viridis::rocket(1, begin = 0.4), 
                        linetype = "dashed", 
                        size = 0.2) +
    ggplot2::geom_label(data = annot_text,
                        aes(y = Name, x = 0.05), 
                        label = paste("cutoff: 0.05", 
                                      sep = " "),
                        size = 2, 
                        vjust = 0.5,
                        hjust = 0.52,
                        alpha = 0.6, 
                        colour = plot_colour) +
    ggplot2::facet_grid(Group~., scales = "free_y", space = "free") +
    ggplot2::scale_x_log10(limits = c(min(na.omit(clev_plotdata$Pr_F_stat))*0.5, 
                                      max(clev_plotdata$Pr_F_stat)*1.5)) +
    viridis::scale_colour_viridis(option = "rocket", 
                                  begin = 0.2, 
                                  end = 0.8, 
                                  discrete = T, 
                                  alpha = 0.7, 
                                  name = "") +
    ggplot2::ggtitle(plot_title[grp]) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, 
                                             colour = plot_colour), 
                   strip.text.y = element_text(size = 7, 
                                               colour = plot_colour,
                                               face = "italic")) +
    ggplot2::xlab("p-values (log10)") +
    ggplot2::ylab("Variables")
  
  ggsave(p, 
         path = path_to_fig,
         filename = paste("pvals_bac_", label_grp[grp], ".pdf", sep = ""),
         device = "pdf", 
         height = 30, 
         width = 20, 
         units = "cm",
         dpi = 320,
         limitsize = F)  
  
  assign(paste("adonis_bac_", label_grp[grp], sep = ""), grp_table, .GlobalEnv)
}


adonis_bac_vg_filt <- adonis_bac_vg %>%
  dplyr::inner_join(clin_tab) %>%
  dplyr::select(Name, everything(), -Var) %>%
  dplyr::group_by(comp, Group) %>%
  dplyr::arrange(comp, Group) %>%
  dplyr::group_by(comp) %>%
  tidyr::nest() %>%
  dplyr::pull(name = comp) %>%
  purrr::reduce(.f = function(x, y) dplyr::full_join(x, y, 
                                                     by = "Name", 
                                                     suffix = c("_case", "_control")))

writexl::write_xlsx(adonis_bac_vg_filt, 
                    path = here::here("Results",
                                      "16S",
                                      "adonis", 
                                      "adonis_vg_bac_filt.xlsx"))

adonis_bac_ed_filt <- adonis_bac_ed %>%
  dplyr::inner_join(clin_tab) %>%
  dplyr::select(Name, everything(), -Var) %>%
  dplyr::group_by(comp, Group) %>%
  dplyr::arrange(comp, Group) %>%
  dplyr::group_by(comp) %>%
  tidyr::nest() %>%
  dplyr::pull(name = comp) %>%
  purrr::reduce(.f = function(x, y) dplyr::full_join(x, y, 
                                                     by = "Name", 
                                                     suffix = c("_case", "_control")))

writexl::write_xlsx(adonis_bac_ed_filt, 
                    path = here::here("Results", 
                                      "16S", 
                                      "adonis",
                                      "adonis_ed_bac_filt.xlsx"))
