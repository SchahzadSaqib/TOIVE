pkgs <- list(here::here("Scripts", 
                  "load_packages.R"),
             here::here("Scripts",
                        "Split_Libraries_cutadapt.R"))
pkgs <- purrr::map(pkgs, source)


##### Remove Ns from the sequences ##-----

#' list files in the directory
fnFs <- sort(list.files(here::here("Data", 
                                   "Raw"),
                        pattern = "_R1_001.fastq", 
                        full.names = T))
fnRs <- sort(list.files(here::here("Data", 
                                   "Raw"),
                        pattern = "_R2_001.fastq", 
                        full.names = T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), 
                       function(x) paste(x[1], collapse = '_'))


#' define separate directory for cleaned files
Nrmv_Fs <- file.path(here::here("Data", 
                                "Remove_Ns"), 
                                paste0(sample.names, 
                                       "_F_Nrmv.fastq.gz"))
Nrmv_Rs <- file.path(here::here("Data", 
                                "Remove_Ns"), 
                     paste0(sample.names, 
                            "_R_Nrmv.fastq.gz"))
names(Nrmv_Fs) <- sample.names
names(Nrmv_Rs) <- sample.names

#' dada2  based N removal              
outN <- dada2::filterAndTrim(fwd = fnFs, 
                             filt = Nrmv_Fs, 
                             rev = fnRs, 
                             filt.rev = Nrmv_Rs, 
                             maxN = 0, 
                             compress = T, 
                             multithread = T, 
                             verbose = T) 
head(outN)
rownames(outN) <- sample.names

#' list new files without Ns
Nrmv_Fs <- sort(list.files(here::here("Data", 
                                      "Remove_Ns"), 
                           pattern = "_F_", 
                           full.names = T))
Nrmv_Rs <- sort(list.files(here::here("Data", 
                                      "Remove_Ns"),
                           pattern = "_R_", 
                           full.names = T))
sample.names <- c(sapply(strsplit(basename(Nrmv_Fs), "_"), 
                         function(x) paste(x[1], collapse = '_')))
sample.names



##### quality profiles ##-----
Qualp_Nrmv_F <- dada2::plotQualityProfile(Nrmv_Fs, 
                                   aggregate = T) # Quality of forward reads
Qualp_Nrmv_F
Qualp_Nrmv_R <- dada2::plotQualityProfile(Nrmv_Rs, 
                                   aggregate = T) # Quality of reverse reads
Qualp_Nrmv_R
# Grey: Heatmap of frequency of quality scores at each base position
# Green: Median Quality scores
# Orange: Quartiles of the quality score
# Red: Scaled proportion of reads, useful for technologies other than illumina


##### run cutadapt to split libraries into 16S and ITS reads ##-----

# Run library split
split.res = splitLibrariesByPrimers(
  fnFs = list.files(path = here::here("Data", 
                                      "Remove_Ns"), 
                    pattern = "_F_", 
                    full.names = T), # Path for raw data
  fnRs = list.files(path = here::here("Data", 
                                      "Remove_Ns"), 
                    pattern = "_R_", 
                    full.names = T),
  path.16S = here::here("Data", 
                        "Cutadapt_clean", 
                        "16S"), # Path for 16S library
  path.ITS = here::here("Data", 
                        "Cutadapt_clean", 
                        "ITS"), # Path for ITS library
  FWD.16S = "CCTACGGGNGGCWGCAG",
  RWD.16S = "GACTACHVGGGTATCTAATCC",
  FWD.ITS = "GGTCATTTAGAGGAAGTAA",
  RWD.ITS = "GCTGCGTTCTTCATCGATGC",
  cutadapt = "python -m cutadapt")




##### 16S pipeline for pre-processing ##-----

#' list 16S files
Cutadpt_Fs <- sort(list.files(here::here("Data",
                                   "Cutadapt_clean",
                                   "16S"),
                        pattern = "_F_Nrmv", 
                        full.names = T))
Cutadpt_Rs <- sort(list.files(here::here("Data",
                                   "Cutadapt_clean",
                                   "16S"),
                        pattern = "R_Nrmv", 
                        full.names = T))
sample.names <- sapply(strsplit(basename(Cutadpt_Fs), "_"), "[", 2)


#' define separate directory for trimmed/truncated files
QC_trim_Fs <- file.path(here::here("Data",
                                   "dada2_QC_trim", 
                                   "16S"), 
                        paste0(sample.names, "_F_16S_QC_trim.fastq.gz"))
QC_trim_Rs <- file.path(here::here("Data",
                                   "dada2_QC_trim", 
                                   "16S"), 
                        paste0(sample.names, "_R_16S_QC_trim.fastq.gz"))
names(QC_trim_Fs) <- sample.names
names(QC_trim_Rs) <- sample.names


#' dada2 filter and trim
out_QC_trim <- dada2::filterAndTrim(fwd = Cutadpt_Fs, 
                                    filt = QC_trim_Fs, 
                                    rev = Cutadpt_Rs, 
                                    filt.rev = QC_trim_Rs, 
                                    maxN = 0, 
                                    truncQ = 2, 
                                    maxEE = 2, 
                                    truncLen = c(260, 220),
                                    minLen = 150, 
                                    rm.phix = T, 
                                    compress = T, 
                                    multithread = T, 
                                    matchIDs = T, 
                                    verbose = T) 

c(mean(out_QC_trim[,2]), sd(out_QC_trim[,2]))
head(out_QC_trim)


#' list new trimmed/truncated files
QC_trim_Fs <- sort(list.files(here::here("Data", 
                                         "dada2_QC_trim", 
                                         "16S"),
                              pattern = "_F_", full.names = T))

QC_trim_Rs <- sort(list.files(here::here("Data", 
                                         "dada2_QC_trim", 
                                         "16S"),
                              pattern = "_R_", full.names = T))
sample.names <- sapply(strsplit(basename(QC_trim_Fs), "_"), "[", 1)
rownames(out_QC_trim) <- sample.names


#' check quality profiles
Qualp_QC_trim_F <- dada2::plotQualityProfile(QC_trim_Fs, aggregate = T)
Qualp_QC_trim_F
Qualp_QC_trim_R <- dada2::plotQualityProfile(QC_trim_Rs, aggregate = T)
Qualp_QC_trim_R

g <- gridExtra::arrangeGrob(Qualp_Nrmv_F, 
                            Qualp_Nrmv_R, 
                            Qualp_QC_trim_F, 
                            Qualp_QC_trim_R, 
                            nrow = 2, 
                            ncol = 2)
ggplot2::ggsave(g, 
                path = here::here("Results"),
                filename = paste("TOIVE_QC_aggregate", 
                                 ".pdf", sep = ""), 
                device = "pdf", 
                width = 40, 
                height = 30, 
                units = "cm", 
                dpi = 320)


#' learn the error rates
errors_F <- dada2::learnErrors(QC_trim_Fs, 
                               nbases = 1e+09, 
                               randomize = T, 
                               multithread = T, 
                               verbose = T)

errors_R <- dada2::learnErrors(QC_trim_Rs, 
                               nbases = 1e+09, 
                               randomize = T, 
                               multithread = T, 
                               verbose = T)


#' plot errors for visual inspection
plot_err_F <- dada2::plotErrors(errors_F, nominalQ = T)
plot_err_F
ggplot2::ggsave(path = here::here("Results"), 
                filename = "TOIVE_error_rates_F.pdf", 
                plot_err_F, 
                device = "pdf", 
                dpi = 320)

plot_err_R <- dada2::plotErrors(errors_R, nominalQ = T)
plot_err_R
ggplot2::ggsave(path = here::here("Results"), 
                filename = "TOIVE_error_rates_R.pdf", 
                plot_err_R, 
                device = "pdf", 
                dpi = 320)
# Black: Observed errors
# Red: Expected error rate



#' sample inference dada2
SInf_Fs <- dada2::dada(QC_trim_Fs, 
                       err = errors_F, 
                       multithread = T, 
                       verbose = T)

SInf_Rs <- dada2::dada(QC_trim_Rs, 
                       err = errors_R, 
                       multithread = T, 
                       verbose = T)



#' merge paired reads
mergers <- dada2::mergePairs(SInf_Fs, 
                             QC_trim_Fs, 
                             SInf_Rs, 
                             QC_trim_Rs, 
                             verbose = T)
head(mergers[[1]])



#' Construct sequence table 
seqtab <- dada2::makeSequenceTable(mergers)
dim(seqtab)
table(nchar(dada2::getSequences((seqtab))))



#' Remove Chimeras
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, 
                                           method = "consensus", 
                                           multithread = T, 
                                           verbose = T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
rownames(seqtab.nochim) <- sample.names



#' Tracking reads through the pipeline
getN <- function(x) sum(dada2::getUniques(x))
track <- outN %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column(var = "samples") %>%
  dplyr::inner_join(out_QC_trim <- data.frame(out_QC_trim) %>%
                      tibble::rownames_to_column(var = "samples"),
                    by = "samples") %>%
  dplyr::inner_join(SInf_Fs <- data.frame(sapply(SInf_Fs, getN)) %>%
                      dplyr::mutate(samples = sample.names),
                    by = "samples") %>%
  dplyr::inner_join(SInf_Rs <- data.frame(sapply(SInf_Rs, getN)) %>%
                      dplyr::mutate(samples = sample.names),
                    by = "samples") %>%
  dplyr::inner_join(mergers <- data.frame(sapply(mergers, getN)) %>%
                      dplyr::mutate(samples = sample.names),
                    by = "samples") %>%
  dplyr::inner_join(nochim <- data.frame(rowSums(seqtab.nochim)) %>%
                      dplyr::mutate(samples = sample.names),
                    by = "samples") %>%
  purrr::set_names(c("samples", 
                     "inputN", 
                     "filteredN", 
                     "inputF", 
                     "filteredF", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", 
                     "nonchim")) %>%
  dplyr::arrange(desc(nonchim)) %>%
  dplyr::mutate(Perc_retained_afterfilt = nonchim/filteredF*100)


#' save image for annotations
base::save(list = c("seqtab.nochim", "track"), 
           file = here::here("Data", 
                             "TOIVE_annot_files.RData"))

#' save pre-processing data
base::save.image(file = here::here("Data", 
                                   "TOIVE_pre_processing.RData"))
