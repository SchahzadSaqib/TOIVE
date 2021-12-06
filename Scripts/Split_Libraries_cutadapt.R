## cutadapt installation:
# pip3 install --user --upgrade cutadapt

## Multicore zip installation
# LINUX: sudo apt-get install pigz
# MAC: brew install pigz

## DADA2 installation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")


# RUN
# system.time(split.res11 <- splitLibrariesByPrimers())



splitLibrariesByPrimers <- function(fnFs = list.files(path = "~/Documents/16S-ITS_trial", 
                                                      pattern = "R1_001.fastq.gz", 
                                                      full.names = TRUE), # Path for raw data
                                    fnRs = list.files(path = "~/Documents/16S-ITS_trial", 
                                                      pattern = "R2_001.fastq.gz", 
                                                      full.names = TRUE),
                                    path.16S = "16S", # Path for 16S library
                                    path.ITS = "ITS", # Path for ITS library
                                    FWD.16S = "CCTACGGGNGGCWGCAG",
                                    RWD.16S = "GACTACHVGGGTATCTAATCC",
                                    FWD.ITS = "GGTCATTTAGAGGAAGTAA",
                                    RWD.ITS = "GCTGCGTTCTTCATCGATGC",
                                    cutadapt = "/Users/schahzad/miniconda3/bin/cutadapt") # cutadapt path
{
  
  # Full paths for output folders if subdirectory
  if (dirname(path.16S) == ".") {
    path.16S <- file.path(dirname(fnFs[1]), path.16S)
  } 
  if (dirname(path.ITS) == ".") {
    path.ITS <- file.path(dirname(fnFs[1]), path.ITS)
  }
  
  path.ITS.withPrimers <- file.path(path.ITS, "withPrimers") # Remove to get rid of the untrimmed ITS reads
  
  # Create output folders if necessary
  if (!dir.exists(path.16S)) dir.create(path.16S, recursive = T)
  if (!dir.exists(path.ITS)) dir.create(path.ITS, recursive = T)
  if (!dir.exists(path.ITS.withPrimers)) dir.create(path.ITS.withPrimers, recursive = T)
  
  # 16S output files
  fnFs.16S <- file.path(path.16S, paste0("16S_", basename(fnFs)))
  fnRs.16S <- file.path(path.16S, paste0("16S_", basename(fnRs)))
  
  # ITS output files
  ## ITS with primers
  fnFs.ITS.wp <- file.path(path.ITS.withPrimers, basename(fnFs))
  fnRs.ITS.wp <- file.path(path.ITS.withPrimers, basename(fnRs))
  ## Final ITS
  fnFs.ITS <- file.path(path.ITS, paste0("ITS_", basename(fnFs)))
  fnRs.ITS <- file.path(path.ITS, paste0("ITS_", basename(fnRs)))
  
  
  # Create results data.frame
  split.result <- data.frame(row.names = fnFs,
                             "Sample" = gsub("_.*", "", basename(fnFs)), #Sample name = filename to first "_"
                             "ReadCountOrig" = rep(0,length(fnFs)),
                             "ReadCount16S" = rep(0,length(fnFs)),
                             "ReadCountITS" = rep(0,length(fnFs)),
                             "Perc.16S" = rep(0,length(fnFs)),
                             "Perc.ITS" = rep(0,length(fnFs)),
                             "Perc.Lost" = rep(0,length(fnFs)),
                             "RC.Lost" = rep(0,length(fnFs)),
                             "File16S" = fnFs.16S,
                             "FileITS" = fnFs.ITS,
                             "Primer16SFWD" = rep(FWD.16S,length(fnFs)),
                             "Primer16SRWD" = rep(RWD.16S,length(fnFs)),
                             "PrimerITSFWD" = rep(FWD.ITS,length(fnFs)),
                             "PrimerITSRWD" = rep(RWD.ITS,length(fnFs))
  )
  
  
  
  #### 16S ####
  FWD <- FWD.16S
  REV <- RWD.16S
  FWD.RC <- paste0(dada2:::rc(FWD)) # Reverse complement
  REV.RC <- paste0(dada2:::rc(REV)) # Reverse complement
  R1.flags <- paste("-g", FWD, "-a", REV.RC)
  R2.flags <- paste("-G", REV, "-A", FWD.RC)

  for (i in seq_along(fnFs)) {
    
    cutadapt_split <- rstudioapi::terminalExecute(paste(cutadapt, R1.flags, R2.flags, "-n", 2, "--cores=0 --discard-untrimmed --pair-filter=both",
                               "-o", fnFs.16S[i], "-p", fnRs.16S[i], # output files
                               fnFs[i], fnRs[i])) # input files
    
    while (is.null(rstudioapi::terminalExitCode(cutadapt_split))) {
      Sys.sleep(0.1)
    }
    rstudioapi::terminalKill(cutadapt_split)
  }
  
  #### ITS ####
  FWD <- FWD.ITS
  REV <- RWD.ITS
  FWD.RC <- paste0(dada2:::rc(FWD)) # Reverse complement
  REV.RC <- paste0(dada2:::rc(REV)) # Reverse complement
  R1.flags <- paste("-g", FWD, "-a", REV.RC)
  R2.flags <- paste("-G", REV, "-A", FWD.RC)

  for (i in seq_along(fnFs)) {
    
    cutadapt_split <- rstudioapi::terminalExecute(paste(cutadapt, 
                                                        R1.flags, 
                                                        R2.flags, 
                                                        "-n", 
                                                        2, 
                                                        "--cores=0 --discard-untrimmed --pair-filter=both", # Read length cutoff removed for testing "--minimum-length 10" 
                                                        "-o", 
                                                        fnFs.ITS[i], 
                                                        "-p", 
                                                        fnRs.ITS[i], # output files
                                                        fnFs[i], 
                                                        fnRs[i])) # input files # You can put fn(F/R)s.ITS.wp[i] here, if you want to exclude 16S reads and get rid of possible overlap in libraries.
    
    while (is.null(rstudioapi::terminalExitCode(cutadapt_split))) {
      Sys.sleep(0.1)
    }
    rstudioapi::terminalKill(cutadapt_split)
    
    # Get read counts
    RC.orig <- as.numeric(system(paste0("echo $(zcat < ", fnFs[i], "|wc -l)/4|bc"), intern = TRUE))
    RC.16S <- as.numeric(system(paste0("echo $(zcat < ", fnFs.16S[i], "|wc -l)/4|bc"), intern = TRUE))
    RC.ITS <- as.numeric(system(paste0("echo $(zcat < ", fnFs.ITS[i], "|wc -l)/4|bc"), intern = TRUE))
    PERC.16S <- RC.16S / RC.orig * 100
    PERC.ITS <- RC.ITS / RC.orig * 100
    RC.Lost <- RC.orig - RC.16S - RC.ITS
    PERC.Lost <- RC.Lost / RC.orig * 100
    
    # Save read counts to a data.frame
    split.result[fnFs[[i]], c("ReadCountOrig", "ReadCount16S", "ReadCountITS", "Perc.16S", "Perc.ITS", "Perc.Lost", "RC.Lost")] <- c(RC.orig, RC.16S, RC.ITS, PERC.16S, PERC.ITS, PERC.Lost, RC.Lost)
  }
  
  # Write stats to Excel file (date and hour included in filename)
  xlsx::write.xlsx2(split.result, file = file.path(dirname(fnFs.16S[1]), paste0("Split_Stats_ITS_16S_", gsub(" |-|:.*", "", Sys.time()), ".xlsx")), asTable = T, row.names = T)

  print(paste0("Input folder: ", dirname(fnFs[1])))
  print(paste0("16S output folder: ", dirname(fnFs.16S[1])))
  print(paste0("ITS output folder: ", dirname(fnFs.ITS[1])))
  print(paste0("Split stats: ",file.path(dirname(fnFs.16S[1]), paste0("Split_Stats_ITS_16S_", gsub(" |-|:.*", "", Sys.time()), ".xlsx"))))
    
  # Return stats
  return(split.result)
}
