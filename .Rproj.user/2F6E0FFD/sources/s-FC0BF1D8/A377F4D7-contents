#' load packages
source(here::here("Scripts", "load_packages.R"))

#' load data
base::load(file = here::here("Data", 
                             "TOIVE_annot_ITS_files.RData"))

#' Set NCBI key 
set_entrez_key("cd9e1966da8433f65cac528b9de11447bd08")
Sys.getenv("ENTREZ_KEY")

#' Align against nt database, restricted by accession IDs 
#' (no environmental samples)

#' BLAST based alignment
Blast_hits_noUC <- taxminer::txm_align(
  input = seqtab.nochim, 
  database_path = here::here("~", 
                             "Documents", 
                             "NCBI_databases", 
                             "nt"), 
  database_name = "nt", 
  output_path = here::here("Results", 
                           "ITS"), 
  output_name = "TOIVE_align_noUC_ITS", 
  accession_list = "Fungi_noENV.seq", 
  accession_path = here::here("~",
                              "Documents",
                              "Accession_lists"), 
  do_acc_check = F,
  Run_Blast = F,
  threads = 5, 
  qcvg = 98, 
  pctidt = 98, 
  max_out = 5000)

#' Align against nt database, restricted by accession IDs 
#' (only environmental samples)

#' BLAST based alignment
Blast_hits_UC <- taxminer::txm_align(
  input = seqtab.nochim, 
  database_path = here::here("~", 
                             "Documents", 
                             "NCBI_databases", 
                             "nt"), 
  database_name = "nt", 
  output_path = here::here("Results", 
                           "ITS"), 
  output_name = "TOIVE_align_UC_ITS", 
  accession_list = "Fungi_ENV.seq", 
  accession_path = here::here("~",
                              "Documents",
                              "Accession_lists"), 
  do_acc_check = F,
  Run_Blast = F,
  threads = 5, 
  qcvg = 98, 
  pctidt = 98, 
  max_out = 5000)

#' Text-mining based filtration of hits
Annot_TOIVE_noUC <- Blast_hits_noUC %>%  
  taxminer::txm_ecosrc(
    filter_host = "human", 
    filter_site = c("vagina+FRS", "gut+oral+clinical"), 
    filter_negate = "non_human", 
    Precomp_tbl  = here::here("~", 
                              "OneDrive - University of Helsinki", 
                              "Eco_lists", 
                              "Fungi_ecosys.rds")) %>%
  
  #' Clean data
  dplyr::group_by(ID) %>%
  dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::distinct(Species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::ungroup() %>%
  base::as.data.frame() 



#' Text-mining based filtration of hits
Annot_TOIVE_UC <- Blast_hits_UC %>%
  taxminer::txm_ecosrc(
    filter_host = "human", 
    filter_site = c("vagina+FRS", "gut+oral+clinical"), 
    filter_negate = "non_human", 
    Precomp_tbl = here::here("~", 
                             "OneDrive - University of Helsinki", 
                             "Eco_lists", 
                             "Fungi_ecosys.rds")
  ) %>%
  
  #' Clean data
  dplyr::group_by(ID) %>%
  dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::distinct(Species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::ungroup() %>%
  base::as.data.frame() %>%
  dplyr::filter(!ID %in% Annot_TOIVE_noUC$ID)


#' Combining annotations and adding lineage
Annot_TOIVE <- Annot_TOIVE_noUC %>%
  dplyr::bind_rows(Annot_TOIVE_UC) %>%
  dplyr::arrange(ID) %>%
  taxminer::txm_lineage(Precomp_tbl = here::here("~", 
                                                 "OneDrive - University of Helsinki", 
                                                 "Eco_lists", 
                                                 "Fungi_lineage.rds"))



#' Create taxonomy table
taxonomy <- Annot_TOIVE %>%
  dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
  dplyr::rename_with(tolower)

#' Filter out sequences that were removed
TOIVE_seqtab_filt <- seqtab.nochim %>%
  base::as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(ids = seq(1:nrow(.))) %>%
  dplyr::filter(ids %in% Annot_TOIVE$ID) %>%
  dplyr::select(-ids) %>%
  t()


#' Assign sequences to rownames of taxonomy  
rownames(taxonomy) <- colnames(TOIVE_seqtab_filt)
taxonomy <- as.matrix(taxonomy)
saveRDS(taxonomy, file = here::here("Data", 
                                    "TOIVE_ITS_taxonomy.rds"))

saveRDS(TOIVE_seqtab_filt, file = here::here("Data", 
                                             "TOIVE_ITS_seqtab_filt.rds"))


#' Add annotated read percentage to tracking
track_cmpl <- data.frame(TOIVE_seqtab_filt) %>%
  dplyr::mutate(reads_annotated = rowSums(.)) %>%
  tibble::rownames_to_column("samples") %>%
  dplyr::select(samples, reads_annotated) %>%
  dplyr::inner_join(track) %>%
  dplyr::relocate(reads_annotated, .after = last_col()) #%>%
  dplyr::mutate(Perc_reads_annotated = reads_annotated/rowSums.seqtab.nonchim.*100) %>%
  dplyr::arrange(str_detect(samples, "\\+"), desc(Perc_retained_afterfilt))


writexl::write_xlsx(track_cmpl, 
                    path = here::here("Results",
                                      "ITS",
                                      "TOIVE_readtrack_annot.xlsx"))


