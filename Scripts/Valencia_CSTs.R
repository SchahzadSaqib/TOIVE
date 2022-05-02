source(
  here::here("Scripts", 
             "load_packages.R")
)

TOIVE_phy_16S <- readRDS(
  here::here("Data", 
             "TOIVE_phy_16S.rds")
)


# Raw data
Valencia_input <- TOIVE_phy_16S@otu_table %>% 
  t() %>%
  data.frame() %>%
  dplyr::mutate(sp = TOIVE_phy_16S@tax_table@.Data[,7]) %>%
  dplyr::mutate(g = TOIVE_phy_16S@tax_table@.Data[,6]) %>%
  dplyr::mutate(across(where(is.factor), 
                       as.character)) %>%
  dplyr::select(c(sp, g), everything()) %>%
  dplyr::mutate(g = ifelse(
    stringr::str_detect(
      sp, 
      "uncultured"), 
    "uncultured", 
    g)
  ) %>%
  dplyr::mutate(g = ifelse(
    is.na(g), 
    sp, 
    g)) %>%
  dplyr::mutate(sp = stringr::str_replace(
    sp, 
    " ", 
    "_")) %>%
  dplyr::mutate(g = paste("g_", g, sep = "")) %>%
  dplyr::mutate(sp = ifelse(
    stringr::str_detect(
      g, 
      "Lactobacillus|Gardnerella|Prevotella|Fannyhessea|Sneathia"), 
    sp, 
    g)) %>%
  dplyr::select(-g) %>%  
  dplyr::group_by(sp) %>%
  dplyr::summarise_all(~ sum(.)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "sp") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sampleID") %>%
  dplyr::mutate(read_count = rowSums(
    dplyr::select(., -sampleID))) %>%
  dplyr::select(sampleID, 
                read_count, 
                everything())

write_csv(Valencia_input,
          file = here::here("Data",
                            "Valencia_input_HC.csv"))

#python3 ~/Documents/Valencia/Valencia.py -r ~/Documents/Valencia/CST_centroids_012920.csv -i Data/Valencia_input_HC.csv -o Data/TOIVE_Val_CSTs -p Data/plot

Valencia <- read.csv(
  here::here("Data",
             "TOIVE_Val_CSTs.csv")) %>%
  dplyr::select(sampleID, CST, subCST) %>%
  dplyr::filter(!str_detect(sampleID, "Blank")) %>%
  dplyr::arrange(CST) 


#write_excel_csv2(path = "Valencia_CSTs.csv")
writexl::write_xlsx(Valencia, 
                    path = here::here("Results", 
                                      "16S",
                                      "Valencia_CSTs.xlsx"))

Valencia <- Valencia %>%
  tibble::column_to_rownames(var = "sampleID")

# Specify color
CST_colours <- viridis::rocket(
  length(unique(Valencia$CST)), 
  begin = 0.3, 
  end = 0.9)
names(CST_colours) <- unique(Valencia$CST)

subCST_colours <- PNWColors::pnw_palette(n = length(unique(Valencia$subCST)), 
                                         name = "Cascades")
#subCST_colours <- sample(subCST_colours)
names(subCST_colours) <- unique(Valencia$subCST)

annot_colours <- list(CST = CST_colours, subCST = subCST_colours)




TOIVE_phy_16S_rel <- phyloseq::transform_sample_counts(
  TOIVE_phy_16S, function(OTU) OTU/sum(OTU)) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  )

hist_plot <- TOIVE_phy_16S_rel@otu_table %>% 
  t() %>%
  data.frame() %>%
  dplyr::select(-c(1:2)) %>%
  dplyr::mutate(sp = TOIVE_phy_16S_rel@tax_table@.Data[,7]) %>%
  dplyr::mutate(across(where(is.factor), 
                       as.character)) %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise_all(~ sum(.)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "sp") %>%
  as.matrix()


h <- pheatmap::pheatmap(hist_plot, 
                        color = viridis::mako(100, 
                                              begin = 0.3,
                                              end = 0.95, 
                                              direction = -1),
                        annotation_col = Valencia,
                        annotation_colors = annot_colours,
                        clustering_method = "ward.D2", 
                        cutree_cols = 5,
                        fontsize = 18, 
                        treeheight_row = 300,
                        treeheight_col = 300,
                        angle_col = 90,
                        border_color = NA,
                        filename = here::here("Results", 
                                              "16S", 
                                              "TOIVE_heatmap_CSTs.png"), 
                        width = 40, 
                        height = 18)
dev.off()
h$gtable$grobs[[5]]$gp = grid::gpar(fontface = "italic", col = "grey40")
h$gtable$grobs[[4]]$gp = grid::gpar(col = "grey40")

h        

tab <- data.frame(TOIVE_phy_16S@sam_data) %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::mutate(CsCtrl = as.numeric(as.character(CsCtrl))) %>%
  dplyr::filter(CsCtrl > 0) %>%
  dplyr::mutate(CsCtrl = ifelse(CsCtrl == 1, "RPL", "Control"))

pair_match <- Valencia %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::left_join(tab) %>%
  dplyr::arrange(id) %>%
  dplyr::group_by(id) %>%
  dplyr::add_tally() %>%
  dplyr::filter(n > 1) %>%
  dplyr::select(-n) %>%
  dplyr::mutate(
    CST_match = ifelse(length(unique(CST)) == 1, "match", "diff"),
    subCST_match = ifelse(length(unique(subCST)) == 1, "match", "diff")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(id, .keep_all = T)


pair_match1 <- pair_match %>%
  dplyr::group_by(CST_match) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

writexl::write_xlsx(pair_match1, 
                    path = here::here("Results", 
                                      "16S", 
                                      "pair_CSTs.xlsx"))

pair_match2 <- pair_match %>%
  dplyr::group_by(subCST_match) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

writexl::write_xlsx(pair_match2, 
                    path = here::here("Results", 
                                      "16S", 
                                      "pair_subCSTs.xlsx"))

full_vg <- Valencia %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::left_join(tab) %>%
  dplyr::filter(Sample_grp == "Vagina")

full_vg1 <- full_vg %>%
  dplyr::group_by(CST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

full_vg2 <- full_vg %>%
  dplyr::group_by(subCST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

CsCtrl_vg1 <- full_vg %>%
  dplyr::group_by(CsCtrl, CST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100) %>%
  tidyr::pivot_wider(id_cols = CST, 
                     names_from = CsCtrl, 
                     values_from = c(n, Prc))

writexl::write_xlsx(CsCtrl_vg1, 
                    path = here::here("Results", 
                                      "16S", 
                                      "VG_CSTs.xlsx"))

CsCtrl_vg2 <- full_vg %>%
  dplyr::group_by(CsCtrl, subCST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100) %>%
  tidyr::pivot_wider(id_cols = subCST, 
                     names_from = CsCtrl, 
                     values_from = c(n, Prc))

writexl::write_xlsx(CsCtrl_ed2, 
                    path = here::here("Results", 
                                      "16S", 
                                      "VG_subCSTs.xlsx"))

full_ed <- Valencia %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::left_join(tab) %>%
  dplyr::filter(Sample_grp == "Endometrium")

full_ed1 <- full_ed %>%
  dplyr::group_by(CST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

full_ed2 <- full_ed %>%
  dplyr::group_by(subCST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100)

CsCtrl_ed1 <- full_ed %>%
  dplyr::group_by(CsCtrl, CST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100) %>%
  tidyr::pivot_wider(id_cols = CST, 
                     names_from = CsCtrl, 
                     values_from = c(n, Prc))

writexl::write_xlsx(CsCtrl_ed1, 
                    path = here::here("Results", 
                                      "16S", 
                                      "Endo_CSTs.xlsx"))

CsCtrl_ed2 <- full_ed %>%
  dplyr::group_by(CsCtrl, subCST) %>%
  dplyr::tally() %>%
  dplyr::mutate(Prc = n/sum(n)*100) %>%
  tidyr::pivot_wider(id_cols = subCST, 
                     names_from = CsCtrl, 
                     values_from = c(n, Prc))

writexl::write_xlsx(CsCtrl_ed2, 
                    path = here::here("Results", 
                                      "16S", 
                                      "Endo_subCSTs.xlsx"))
