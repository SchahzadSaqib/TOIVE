source(here::here("Scripts", "load_packages.R"))


##### Load and clean data - bacteria ##------
#' ASV table
TOIVE_16S_ASVs <- readr::read_rds(
  here::here("Data", 
             "TOIVE_16S_seqtab_filt.rds") 
) 

#' taxonomy
TOIVE_16S_taxonomy <- readr::read_rds(
  here::here("Data", 
             "TOIVE_16S_taxonomy.rds")) %>%
  base::as.data.frame() %>%
  dplyr::mutate(species = paste(species, " ", sep = "")) %>%
  dplyr::mutate(species = ifelse(str_ends(species, ".1"), 
                                 "uncultured bacteria", 
                                 species)) %>%
  dplyr::mutate(species = str_extract(species, "^(?:[^ ]+ ){2}")) %>%
  dplyr::mutate(species = str_replace_all(species, " ", ".")) %>%
  dplyr::mutate(species = str_replace(species, "\\.", " ")) %>%
  dplyr::mutate(species = str_remove_all(species, "\\.")) %>%
  base::as.matrix()

#' clean meta file
metadata <- readxl::read_xlsx(
  here::here("Data", 
             "TOIVE_clinical_and microbiome_data.xlsx")
) %>%
  dplyr::select(where(
    ~!all(is.na(.x))
  )
  ) %>%
  dplyr::mutate(
    across(where(is.numeric), 
           .fns = ~tidyr::replace_na(
             data = .x, replace = 0))
  ) %>%
  dplyr::mutate(
    across(where(is.character),
           .fns = ~tidyr::replace_na(
             data = .x, replace = "unknown"
           ))
  ) %>%
  dplyr::select(1:education) %>%
  dplyr::select(where(
    ~ is.numeric(.x) && sum(.x) !=0 || is.character(.x))
  ) %>%
  dplyr::mutate(across(where(is.numeric), ~ifelse(.x == 99, NA, .x)))


#' Combined meta file
meta_cmb <- data.frame(rownames(TOIVE_16S_ASVs)) %>%
  purrr::set_names("Sample_ID") %>%
  dplyr::filter(!is.na(Sample_ID)) %>%
  dplyr::mutate(
    id = str_extract(Sample_ID, pattern = "\\d.\\d*"),
    id = ifelse(is.na(id), paste(1:sum(is.na(id))), id),
    id = as.numeric(id)
  ) %>%
  dplyr::full_join(metadata) %>%
  dplyr::mutate(
    Sample_grp = case_when((stringr::str_starts(Sample_ID, "Bl")) ~ "Blanks",
                           (stringr::str_ends(Sample_ID, "A")) ~ "Vagina",
                           (stringr::str_ends(Sample_ID, "D")) ~ "Endometrium"),
    CsCtrl = case_when((stringr::str_starts(Sample_ID, "Bl")) ~ 0,
                       (stringr::str_detect(Sample_ID, "V1")) ~ 1,
                       (stringr::str_detect(Sample_ID, "V2")) ~ 2),
    dg = if_else(dg1 == 0, 0, 1),
    cause0 = if_else(cause == 2, 0, cause),
    ageyear_strat = case_when((ageyear < 30) ~ 1,
                              (ageyear >= 30 & ageyear < 35) ~ 2,
                              (ageyear >= 35) ~ 3),
    bmi = round(bmi, digits = 1),
    bmi3 = case_when((bmi < 24.9) ~ 1, 
                     (bmi >= 24.9 & bmi < 29.9) ~ 2,
                     (bmi >= 29.9) ~ 3),
    p01 = if_else(p == 0, 0, 1),
    g_strat = if_else(g > 6, 6, g),
    p_strat = if_else(p > 2, 2, p),
    alcdose_01 = if_else(alcdose == 0, 0, 1),
    probiotics_01 = dplyr::case_when((is.na(probiotics)) ~ probiotics,
                                     (probiotics %in% c(1, 2)) ~ 1,
                                     (probiotics %in% c(3,4)) ~ 0),
    Richness = vegan::specnumber(TOIVE_16S_ASVs),
    Diversity = vegan::diversity(TOIVE_16S_ASVs, index = "simpson"),
    ReadCount = rowSums(TOIVE_16S_ASVs)) %>%
  dplyr::select(id, Sample_grp, CsCtrl, sort(names(.))) %>%
  dplyr::mutate(across(-c(height, weight, bmi, Richness, 
                          Diversity, ReadCount, id, ageyear,
                          agemonth, lastpdate, lastptimem, lastptimed), 
                       .fns = ~as.factor(.x))) %>%
  tibble::column_to_rownames(var = "Sample_ID")



##### Construct phyloseq object ##-----
#' Phyloseq object 
TOIVE_phy_16S <- phyloseq::phyloseq(otu_table(TOIVE_16S_ASVs, 
                                              taxa_are_rows = F), 
                                    sample_data(meta_cmb), 
                                    tax_table(TOIVE_16S_taxonomy)) %>%
  phyloseq::tax_glom(taxrank = "species") 


#' Constructing a phylogenetic tree
#seqs <- phyloseq::taxa_names(TOIVE_phy_16S)
#names(seqs) <- phyloseq::taxa_names(TOIVE_phy_16S)
#alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = T)
#phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
#dm <- phangorn::dist.ml(phangAlign)
#treeNJ <- phangorn::NJ(dm)
#fit = phangorn::pml(treeNJ, data = phangAlign)
#fitGTR <- stats::update(fit, k = 4, inv = 0.2)
#fitGTR <- phangorn::optim.pml(fitGTR, model = "GTR", optInv = T, optGamma = T, rearrangement = "stochastic",
#                              control = phangorn::pml.control(trace = 0))
#saveRDS(fitGTR, file = here::here("Data", "TOIVE_phy_16S_tree.rds"))
fitGTR <- base::readRDS(here::here("Data", 
                                   "TOIVE_phy_16S_tree.rds"))


#' Merge with phyloseq object
TOIVE_phy_16S <- phyloseq::merge_phyloseq(TOIVE_phy_16S, 
                                          phy_tree(fitGTR$tree)) %>%
  phyloseq::filter_taxa(function(x){sum(x > 0) > (0.01*length(x))}, 
                        prune = TRUE) %>%
  phyloseq::filter_taxa(function(x){sum(x) > 100}, 
                        prune = TRUE)

#' Convert ASVs sequences to corresponding species names
phyloseq::taxa_names(TOIVE_phy_16S) <- as.character(
  TOIVE_phy_16S@tax_table@.Data[,7]
)

#' adjust readcounts
reads <- data.frame(TOIVE_phy_16S@otu_table) %>%
  dplyr::summarise(rowSums(.)) %>%
  purrr::set_names("ReadCount")
phyloseq::sample_data(TOIVE_phy_16S)$ReadCount <- reads$ReadCount

TOIVE_phy_16S <- phyloseq::subset_samples(TOIVE_phy_16S, ReadCount > 500)
TOIVE_phy_16S


if (identical(as.vector(rowSums(data.frame(TOIVE_phy_16S@otu_table))), TOIVE_phy_16S@sam_data$ReadCount)) {
  print("reads match")
} else {
  stop("ReadCounts are not identical")
}

#' Save phyloseq object
base::saveRDS(TOIVE_phy_16S, file = here::here("Data", 
                                               "TOIVE_phy_16S.rds"))

#' Save table 
TOIVE_16S_table <- data.frame(TOIVE_phy_16S@otu_table) %>%
  tibble::rownames_to_column(var = "Sample_IDs")

writexl::write_xlsx(TOIVE_16S_table, 
                    path = here::here("Data", "TOIVE_asv_table_Raw.xlsx"))


# Convert raw counts to relative abundance
TOIVE_phy_16S_rel <- phyloseq::transform_sample_counts(
  TOIVE_phy_16S, function(OTU) OTU/sum(OTU)) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  )

#' Save table 
TOIVE_16S_table <- data.frame(TOIVE_phy_16S_rel@otu_table) %>%
  tibble::rownames_to_column(var = "Sample_IDs")

writexl::write_xlsx(TOIVE_16S_table, 
                    path = here::here("Data", "TOIVE_asv_table_Rel.xlsx"))





##### Load and clean data - Fungi ##-----
#' seqtab
TOIVE_fungi_seqtab <- readr::read_rds(
  here::here("Data", "TOIVE_ITS_seqtab_filt.rds") 
) %>%
  as.data.frame() %>%
  dplyr::filter(!str_detect(rownames(.), "Blank")) %>%
  as.matrix()


#' taxonomy
TOIVE_fungi_taxonomy <- readr::read_rds(
  here::here("Data", "TOIVE_ITS_taxonomy.rds")
)

#' clean meta file
metadata <- readxl::read_xlsx(
  here::here("Data", 
             "TOIVE_clinical_and microbiome_data.xlsx")
) %>%
  dplyr::select(where(
    ~!all(is.na(.x))
  )
  ) %>%
  dplyr::mutate(
    across(where(is.numeric), 
           .fns = ~tidyr::replace_na(
             data = .x, replace = 0))
  ) %>%
  dplyr::mutate(
    across(where(is.character),
           .fns = ~tidyr::replace_na(
             data = .x, replace = "unknown"
           ))
  ) %>%
  dplyr::select(1:education) %>%
  dplyr::select(where(
    ~ is.numeric(.x) && sum(.x) !=0 || is.character(.x))
  ) %>%
  dplyr::mutate(across(where(is.numeric), ~ifelse(.x == 99, NA, .x)))


#' metadata
metadata_fungi_grpt <- data.frame(rownames(TOIVE_fungi_seqtab)) %>%
  purrr::set_names("Sample_ID") %>%
  dplyr::filter(!is.na(Sample_ID)) %>%
  dplyr::mutate(
    id = str_extract(Sample_ID, pattern = "\\d.\\d*"),
    id = as.numeric(id)) %>%
  dplyr::inner_join(metadata) %>%
  dplyr::mutate(
    Sample_grp = if_else(str_ends(Sample_ID, "A"), "Vagina", "Endometrium"),
    CsCtrl = if_else(str_starts(id, "1"), 1, 2),
    dg = if_else(dg1 == 0, 0, 1),
    cause0 = if_else(cause == 2, 0, cause),
    ageyear_strat = case_when((ageyear < 30) ~ 1,
                              (ageyear >= 30 & ageyear < 35) ~ 2,
                              (ageyear >= 35) ~ 3),
    bmi = round(bmi, digits = 1),
    bmi3 = case_when((bmi < 24.9) ~ 1, 
                     (bmi >= 24.9 & bmi < 29.9) ~ 2,
                     (bmi >= 29.9) ~ 3),
    p01 = if_else(p == 0, 0, 1),
    g_strat = if_else(g > 6, 6, g),
    p_strat = if_else(p > 2, 2, p),
    alcdose_01 = if_else(alcdose == 0, 0, 1),
    probiotics_01 = dplyr::case_when((is.na(probiotics)) ~ probiotics,
                                     (probiotics %in% c(1, 2)) ~ 1,
                                     (probiotics %in% c(3,4)) ~ 0),
    Richness = vegan::specnumber(TOIVE_fungi_seqtab),
    Diversity = vegan::diversity(TOIVE_fungi_seqtab, index = "simpson"),
    ReadCount = rowSums(TOIVE_fungi_seqtab)) %>%
  dplyr::select(id, Sample_grp, CsCtrl, sort(names(.))) %>%
  dplyr::mutate(across(-c(height, weight, bmi, Richness, 
                          Diversity, ReadCount, id, ageyear,
                          agemonth, lastpdate, lastptimem, lastptimed), 
                       .fns = ~as.factor(.x))) %>%
  tibble::column_to_rownames(var = "Sample_ID")



##### Construct phyloseq object ##-----
#' Phyloseq object 
TOIVE_phy_fungi <- phyloseq::phyloseq(
  otu_table(TOIVE_fungi_seqtab, taxa_are_rows = F), 
  sample_data(metadata_fungi_grpt), 
  tax_table(TOIVE_fungi_taxonomy)) %>%
  phyloseq::tax_glom(taxrank="species") 


#' Constructing a phylogenetic tree
#seqs <- phyloseq::taxa_names(TOIVE_phy_fungi)
#names(seqs) <- phyloseq::taxa_names(TOIVE_phy_fungi)
#alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose = T)
#phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
#dm <- phangorn::dist.ml(phangAlign)
#treeNJ <- phangorn::NJ(dm)
#fit = phangorn::pml(treeNJ, data = phangAlign)
#fitGTR <- stats::update(fit, k=4, inv=0.2)
#fitGTR <- phangorn::optim.pml(fitGTR, model = "GTR", optInv = T, optGamma = T, rearrangement = "stochastic",
#                              control = phangorn::pml.control(trace = 0))

#saveRDS(fitGTR, file = here::here("Data", "TOIVE_phy_fungi_tree.rds"))
fitGTR <- readRDS(here::here("Data", 
                             "TOIVE_phy_fungi_tree.rds"))


TOIVE_phy_fungi <- phyloseq::merge_phyloseq(TOIVE_phy_fungi, 
                                            phy_tree(fitGTR$tree)) %>%
  phyloseq::filter_taxa(function(x){sum(x > 0) > (0.01*length(x))}, 
                        prune = TRUE) %>%
  phyloseq::filter_taxa(function(x){sum(x) > 100}, 
                        prune = TRUE)

phyloseq::taxa_names(TOIVE_phy_fungi) <- as.character(
  TOIVE_phy_fungi@tax_table@.Data[,7]
  )

#' adjust readcounts
reads <- data.frame(TOIVE_phy_fungi@otu_table) %>%
  dplyr::summarise(rowSums(.)) %>%
  purrr::set_names("ReadCount")
phyloseq::sample_data(TOIVE_phy_fungi)$ReadCount <- reads$ReadCount

TOIVE_phy_fungi <- phyloseq::subset_samples(TOIVE_phy_fungi, ReadCount > 500)
TOIVE_phy_fungi


if (identical(as.vector(rowSums(data.frame(TOIVE_phy_fungi@otu_table))), TOIVE_phy_fungi@sam_data$ReadCount)) {
  print("reads match")
} else {
  stop("ReadCounts are not identical")
}

#' Save phyloseq object
saveRDS(TOIVE_phy_fungi, file = here::here("Data", 
                                           "TOIVE_phy_fungi.rds"))
