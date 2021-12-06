GroupTest_mod <- function(species.table=NULL, 
                          genus.table=NULL, 
                          family.table=NULL, 
                          order.table=NULL,
                          class.table=NULL, 
                          phylum.table=NULL, 
                          meta, 
                          group,
                          group_name = NULL, 
                          compare.to = NULL, 
                          add_rwns = NULL, 
                          dir_for_res = NULL, 
                          readcount.cutoff = 0, 
                          confounders = NULL, 
                          subject.ID = NULL,
                          outlier.cutoff = 3, 
                          p.cutoff = 0.05, 
                          select.by = NULL, 
                          select = NULL, 
                          pdf = T,
                          min.prevalence = 0.05, 
                          min.abundance = 0, 
                          label.direction = 1, 
                          keep.result = F,
                          nonzero = T, 
                          relative = T, 
                          show_quartz = F) {
  
  ##### Create directory for results ##-------------------------------------------------------------------------
  if(!dir.exists(here::here(dir_for_res, "plots"))) {
    dir.create(here::here(dir_for_res, "plots"), recursive = T)
  }
  
  if(!dir.exists(here::here(dir_for_res, "tables"))) {
    dir.create(here::here(dir_for_res, "tables"), recursive = T)
  }
  
  ##### Check system for quartz ##------------------------------------------------------------------------------
  if(Sys.info()[['sysname']] == "Linux") {
    quartz <- function() {X11()}
  }
  if(Sys.info()[['sysname']] == "Windows") {
    quartz <- function() {X11()}
  }
  
  ##### Load and prepare data ##--------------------------------------------------------------------------------
  # metadata
  
  #meta <- meta
  if (class(meta) == "phyloseq") {
    meta <- data.frame(meta@sam_data)
  } else {
    meta <- read.delim(meta)
  }
  
  if (!is.null(add_rwns)) {
    meta <- meta %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(var = add_rwns)
  }
  
  if(!relative) meta[,"ReadCount"]<-1
  
  # taxonomic tables
  if (class(species.table) == "phyloseq") {
    species <- data.frame(species.table@otu_table)
  } else {
    species <- read.delim(species.table)
  }
  
  #if(length(species.table)>0) species <- read.delim(species.table) else species <- NA
  if(length(genus.table)>0) genus <- read.delim(genus.table) else genus <- NA
  if(length(family.table)>0) family <- read.delim(family.table) else family <- NA
  if(length(order.table)>0) order <- read.delim(order.table) else order <- NA
  if(length(class.table)>0) class <- read.delim(class.table) else class <- NA
  if(length(phylum.table)>0) phylum <- read.delim(phylum.table) else phylum <- NA
  
  #if(min.prevalence<0.3) min.prevalence <- 0.3
  taxa <- data.frame(cbind(species,genus,family,order,class,phylum))
  
  # Subset by selection criteria
  if (length(select.by) != 0) {
    meta <- meta %>%
      dplyr::mutate(selection = .[,select.by])
    
    if (select %in% meta$selection) {
      # subset meta by selection
      meta <- meta %>%
        dplyr::filter(selection %in% select) %>%
        droplevels()
      # subset taxa by selection
      taxa <- taxa %>%
        dplyr::filter(rownames(.) %in% rownames(meta))
    } else {
      stop("Selection variable not found in selected column")
    }
  }
  
  # Subset taxa that are not NAs
  taxa <- taxa[,names(colSums(taxa)[!is.na(colSums(taxa))])]
  
  # Clean column names for taxa
  colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
  colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))])
  #colnames(taxa)[grepl(pattern="Erysipelotrichi_",x=colnames(taxa))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(taxa)[grepl(pattern="Erysipelotrichi_",x=colnames(taxa))])
  
  for(i in grep(pattern="_IncertaeSedis",x=colnames(taxa))){
    colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_IncertaeSedis",fixed=T,
                              replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
  }
  
  for(i in grep(pattern="_incertaesedis",x=colnames(taxa))){
    colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_incertaesedis",fixed=T,
                              replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
  }
  
  for(i in grep(pattern="_uncultured",x=colnames(taxa))){
    colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_uncultured",fixed=T,
                              replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
  }
  
  # Subset taxa and meta with read cutoff threshold
  taxa <- taxa[meta$ReadCount > readcount.cutoff, ]
  meta <- meta[meta$ReadCount > readcount.cutoff, ]
  
  # 
  reltaxa <- (1 + taxa)/meta$ReadCount
  for (i in names(reltaxa)) {
    for (j in 1:nrow(taxa)) {
      reltaxa[j, i][reltaxa[j, i] > (mean(reltaxa[, i]) + outlier.cutoff * 
                                       sd(reltaxa[, i]))] <- mean(reltaxa[, i]) + outlier.cutoff * 
        sd(reltaxa[, i])
    }
  }
  
  
  taxa <- round(reltaxa * meta$ReadCount - 1)
  taxa[taxa<0]<-0
  #taxa <- taxa[, colSums(taxa/ meta$ReadCount > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)] %>%
  #  rownames_to_column(var = "Sample") %>%
  #  filter(Sample %in% rownames(meta)) %>%
  #  #filter(Sample %in% meta$Sample) %>%
  #  column_to_rownames(var = "Sample")
  
  tidyvar <- sym(group)
  
  taxa <- taxa %>%
    
    #' create anchor variables
    tibble::rownames_to_column(var = "SampleN") %>%
    dplyr::mutate(!!tidyvar := meta[[tidyvar]],
                  ReadCounts = meta$ReadCount) %>%
    
    #' reshape to longer format
    tidyr::pivot_longer(cols = -c(!!tidyvar, 
                                  SampleN, 
                                  ReadCounts), 
                        names_to = "taxa",
                        values_to = "Abundance") %>%
    
    #' group by grouping variable and get tally of samples
    #' within each group
    dplyr::group_by(!!tidyvar) %>%
    dplyr::mutate(group_tally = n_distinct(SampleN)) %>%
    
    #' group by samples and convert to relative abundance
    dplyr::group_by(SampleN) %>%
    dplyr::mutate(test_abd = Abundance/sum(Abundance)) %>%
    
    #' group by grouping variable and taxa to determine the 
    #' frequency of occurrence of each taxon above the specified 
    #' minimum abundance threshold
    dplyr::group_by(!!tidyvar, taxa) %>%
    dplyr::mutate(freq = sum(test_abd >= min.abundance),
                  freq = freq/group_tally) %>%
    
    #' filter out taxa that do not meet the minimum prevalence 
    #' threshold in either group
    dplyr::group_by(taxa) %>%
    dplyr::filter(any(freq >= min.prevalence)) %>%
    dplyr::ungroup() %>%
    
    #' clean and return to original dimensions 
    dplyr::select(taxa, Abundance, SampleN) %>%
    tidyr::pivot_wider(id_cols = SampleN, 
                       names_from = "taxa", 
                       values_from = "Abundance") %>%
    dplyr::filter(SampleN %in% rownames(meta)) %>%
    tibble::column_to_rownames(var = "SampleN")
  
  
  
  ##### Prepare dataset ##-----------------------------------------------------------------------------------------
  
  if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
  if(ncol(taxa)>0) {   
    
    dataset <- data.frame(meta, taxa)
    dataset[,group] <- as.character(dataset[,group])
    dataset <- dataset[!is.na(dataset[,group]),]
    dataset[,group][dataset[,group]==""] <- "nogroup"
    for(i in unique(dataset[,group])) if(table(dataset[,group])[i] < 3) dataset[,group][dataset[,group]==i] <- "toofewcases"
    dataset <- dataset[dataset[,group]!="toofewcases",]
    dataset[, group] <- dataset[, group][drop = T]
    # stop if there are less than 2 levels within the grouping variable
    if (nlevels(as.factor(dataset[,group])) < 2) stop("Not enough cases")
    dataset$G <- as.factor(dataset[,group])
    #if (compare.to != "0") dataset[, group][dataset[, group] == "0" & !is.na(dataset[, group])] <- "group0"
    #dataset[, group][dataset[, group] == compare.to & !is.na(dataset[, group])] <- "0"
    dataset[, group] <- as.factor(dataset[, group])
    if (!is.null(compare.to)) {
      dataset <- dataset %>%
        dplyr::mutate(across(.cols = all_of(group), .fns = ~forcats::fct_relevel(.x, compare.to)))
    }
    if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
    print(levels(dataset[,group]), sep = "")
    
    for(i in c(1:ncol(dataset))[-c(1:ncol(meta),ncol(dataset))]){
      if(min(dataset[,i]) == max(dataset[,i]))dataset[,i] <- NA
    }
    
    taxa <- taxa[,names(colSums(dataset[,-c(1:ncol(meta),ncol(dataset))])[!is.na(colSums(dataset[,-c(1:ncol(meta),ncol(dataset))]))])]
    dataset <- data.frame(dataset[,c(1:ncol(meta),ncol(dataset))],
                          dataset[,names(colSums(dataset[,-c(1:ncol(meta),ncol(dataset))])[!is.na(colSums(dataset[,-c(1:ncol(meta),ncol(dataset))]))])])
    
    
    confounders_formula <- paste(confounders, collapse = " + ")
    assign("confounders_formula", paste(confounders, collapse = " + "), envir = .GlobalEnv)
    assign("group", group, envir = .GlobalEnv)
    group_test <- data.frame(array(dim = c(length(names(taxa)),length(c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"),
                                                                        paste("pvals", levels(dataset[,group])[-1], sep = "_"))))))
    rownames(group_test) <- names(taxa)
    names(group_test) <- c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"), 
                           paste("pvals", levels(dataset[,group])[-1], sep = "_"))
    group_test$taxon <- rownames(group_test)
    group_test_nonzero <- group_test
    model <- list()
    updatedmodel <- list()
    
    if (length(subject.ID) != 0) {
      dataset$ID <- as.factor(dataset[, subject.ID])
      modeldata <- na.omit(dataset[, c(names(taxa), group,"G", paste(confounders, sep = ", "), "ID", 
                                       "ReadCount")[c(names(taxa), group, "G",paste(confounders, sep = ", "), "ID", "ReadCount") != ""]])
      modeldata[, group] <- modeldata[, group][drop=T]
      
      resids <- modeldata
      resids[,names(taxa)]<-NA
      modeldata2 <- modeldata
      modeldata2[,names(taxa)]<- modeldata2[,names(taxa)]+1
      print(levels(modeldata[,group]), sep = "")
      print(levels(modeldata2[,group]), sep = "")
      
      for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
        model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                   i, ")~", confounders_formula, "+ ", group, 
                                                                   "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                                                  data = modeldata, 
                                                  admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                               error = function(e) NULL)
        
        if(length(model[[i]])==0|max(model[[i]]$b)==0){
          model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                     i, "+1)~", confounders_formula, "+ ", group, 
                                                                     "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                                                    data = modeldata, 
                                                    admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                                 error = function(e) NULL)
        }
        
        if(length(model[[i]])==0|max(model[[i]]$b)==0){
          model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                     i, ")~", confounders_formula, "+ ", group, 
                                                                     "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "poisson", 
                                                    data = modeldata, 
                                                    admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                                 error = function(e) NULL)
        }
        
        if(length(model[[i]])==0|max(model[[i]]$b)==0){
          if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.255){
            model[[i]] <- tryCatch(nlme::lme(as.formula(paste("((", i, ")/ReadCount)~", confounders_formula, "+ ",group)),random = ~1 | ID,  data = modeldata2),
                                   error = function(e) NULL)
            
          } else{
            model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders_formula, "+ ",group)),random = ~1 | ID,  data = modeldata2),
                                   error = function(e) NULL)
          }}  
        
        tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
        tmp$resdev <- abs(tmp$res) 
        tmp$group <- modeldata[,group]
        
        if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
           anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          
          if(length(summary(model[[i]])$tTable)==0){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$coef),pattern=group,replacement=""))){
              group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j,sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)] }   
          } else {
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
              group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)] }
          }
        } else {
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
            if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                               data = modeldata2, weights = nlme::varExp(),random = ~1 | ID,  
                                               control = nlme::glsControl(maxIter=1000)),
                                     error = function(e) NULL)
              
            } else {
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                               data = modeldata2, weights = nlme::varExp(),random = ~1 | ID,  
                                               control = nlme::glsControl(maxIter=1000)),
                                     error = function(e) NULL)
            }
          } else {
            if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                               data = modeldata2, weights = nlme::varIdent(form=as.formula(paste("~1|",group))),
                                               random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                                     error = function(e) NULL)
            } else {
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                               data = modeldata2, weights = nlme::varIdent(form=as.formula(paste("~1|",group))),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                                     error = function(e) NULL)  
            }
          }    
          if(length(model[[i]])>0){
            tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
            tmp$resdev <- abs(tmp$res) 
            tmp$group <- modeldata[,group]
            if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
               anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
              for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
                group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)]  }
            } else {
              if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
                model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, random = ~1 | ID, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                                                 control = nlme::glsControl(maxIter=1000)),
                                       error = function(e) NULL)
              } else {
                model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, random = ~1 | ID, 
                                                 weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                                                 control = nlme::glsControl(maxIter=1000)),
                                       error = function(e) NULL)  
              }
              if(length(model[[i]])>0){
                tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
                tmp$resdev <- abs(tmp$res) 
                tmp$group <- modeldata[,group]
                if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
                   anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                  for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
                    group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)]  }  
                }
              }
            }
          }
        } 
        
        if(length(model[[i]])>0){
          group_test[i, "model"] <- strsplit(as.character(summary(model[[i]])$call),split="(",fixed=T)[[1]][1] 
          
          if(group_test[i, "model"]== "lme.formula"){
            if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)>0.24999){
              group_test[i, "model"] <- paste("log",group_test[i, "model"])
            } 
          }
          
          for(j in levels(modeldata[,group])[-1]){
            
            if(group_test[i, "model"]== "glmmADMB::glmmadmb"){     
              group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
            } 
            
            if(group_test[i, "model"] == "log lme.formula"){     
              group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
            }
            
            if(group_test[i, "model"] == "lme.formula"){     
              group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange((summary(model[[i]])$tTable[1,1]+summary(model[[i]])$tTable[paste(group,j,sep=""),1]), 
                                                                                         summary(model[[i]])$tTable[1,1])
            }
          }
          
          if(length(confounders[confounders!=""])>0){
            assign("i", paste(i), envir = .GlobalEnv)
            
            resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), 
                                          type="pearson"),
                                    error = function(e) NA)
          }
          
        }
      }
      if(nonzero){       
        residsNonzero <- modeldata
        confounders2 <- confounders
        for(h in seq_along(confounders[confounders!=""])) if(min(table(modeldata[modeldata[,i]>0,confounders[h]]))==0) confounders2[h] <- ""
        
        for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
          if(min(table(modeldata[modeldata[,i]>0,group]))>2){        
            model[[paste(i,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                                        i, ")~", confounders_formula, "+ ", group, 
                                                                                        "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                                                                       data = modeldata[modeldata[,i]>0,], 
                                                                       admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                                                    error = function(e) NULL)
            
            if(length(model[[paste(i,"nonzero")]])==0){
              model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                         i, "+1)~", confounders_formula, "+ ", group, 
                                                                         "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                                                        data = modeldata[modeldata[,i]>0,], 
                                                        admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                                     error = function(e) NULL)
            }
            
            if(length(model[[paste(i,"nonzero")]])==0){
              model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                                                                         i, ")~", confounders_formula, "+ ", group, 
                                                                         "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "poisson", 
                                                        data = modeldata[modeldata[,i]>0,], 
                                                        admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                                     error = function(e) NULL)
            }
            
            if(length(model[[paste(i,"nonzero")]])==0){
              if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){
                model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("((", i, ")/ReadCount)~", confounders_formula, "+ ",group)),random = ~1 | ID,  data = modeldata[modeldata[,i]>0,]),
                                                        error = function(e) NULL)
              } else {
                model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders_formula, "+ ",group)),random = ~1 | ID,  data = modeldata[modeldata[,i]>0,]),
                                                        error = function(e) NULL)
              }}  
            
            tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
            tmp$resdev <- abs(tmp$res) 
            tmp$group <- modeldata[modeldata[,i]>0,group]
            
            if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
               anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
              
              if(length(summary(model[[paste(i,"nonzero")]])$tTable)==0){
                for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$coef),pattern=group,replacement=""))){
                  group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)] }   
              } else {
                for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
                  group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)] }
              }
            } else {
              if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
                if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                } else {
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                }
              } else {
                if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){  
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~",confounders_formula, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], 
                                                                    weights = nlme::varIdent(form=as.formula(paste("~1|",group))),random = ~1 | ID,  
                                                                    control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                } else {
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form=as.formula(paste("~1|",group))),
                                                                    random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)  
                }
              }    
              if(length(model[[paste(i,"nonzero")]])>0){
                tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
                tmp$resdev <- abs(tmp$res) 
                tmp$group <- modeldata[modeldata[,i]>0,group]
                if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
                   anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                  for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
                    group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)]  }
                } else {
                  if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){ 
                    model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                                      data = modeldata[modeldata[,i]>0,], random = ~1 | ID, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                                                                      control = nlme::glsControl(maxIter=1000)),
                                                            error = function(e) NULL)
                  } else {
                    model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                                      data = modeldata[modeldata[,i]>0,], random = ~1 | ID, 
                                                                      weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                                                                      control = nlme::glsControl(maxIter=1000)),
                                                            error = function(e) NULL)  
                  }
                  if(length(model[[paste(i,"nonzero")]])>0){
                    tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
                    tmp$resdev <- abs(tmp$res) 
                    tmp$group <- modeldata[modeldata[,i]>0,group]
                    if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
                       anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                      for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
                        group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)]  }  
                    }
                  }
                }
              }
            }
            
            if(length(model[[paste(i,"nonzero")]])>0){
              group_test_nonzero[i, "model"] <- strsplit(as.character(summary(model[[paste(i,"nonzero")]])$call),split="(",fixed=T)[[1]][1] 
              
              if(group_test_nonzero[i, "model"]== "lme.formula"){
                if(abs(mean(modeldata[,i]/modeldata$ReadCount)-median(modeldata[,i]/modeldata$ReadCount))/median(modeldata[,i]/modeldata$ReadCount)>0.24999){
                  group_test_nonzero[i, "model"] <- paste("log",group_test_nonzero[i, "model"])
                } 
              }
              
              for(j in levels(modeldata[,group])[-1]){
                if(group_test_nonzero[i, "model"]== "glmmADMB::glmmadmb"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
                } 
                
                if(group_test_nonzero[i, "model"] == "log lme.formula"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
                }
                
                if(group_test_nonzero[i, "model"] == "lme.formula"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange((summary(model[[paste(i,"nonzero")]])$tTable[1,1]+summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1]), 
                                                                                                     summary(model[[paste(i,"nonzero")]])$tTable[1,1])
                }
              }
            }  
            
            if(length(confounders[confounders!=""])>0){
              residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), 
                                                                  type="pearson"),
                                                            error = function(e) NA)
            }
            
          }  
          
        }
      }
    } else {
      modeldata <- na.omit(dataset[, c(names(taxa), group, "G", paste(confounders, sep = ", "),  
                                       "ReadCount")[c(names(taxa), group,"G", paste(confounders, sep = ", "), "ReadCount") != ""]])
      modeldata[, group] <- modeldata[, group][drop=T]
      resids <- modeldata
      modeldata2 <- modeldata
      modeldata2[,names(taxa)]<- modeldata2[,names(taxa)]+1
      print(levels(modeldata[,group]), sep = "")
      print(levels(modeldata2[,group]), sep = "")
      
      for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
        model[[i]] <- tryCatch(MASS::glm.nb(formula(paste("round(", i, ")~", confounders_formula, "+ ",group,"+", "offset(log(ReadCount))")), 
                                            data = modeldata, control = glm.control(maxit = 1000,epsilon=1e-12)),
                               error = function(e) NULL)
        
        ## check whether removing the grouping variable breaks the model (crucial for downstream resid-update)
        model_test <- tryCatch(MASS::glm.nb(formula(paste("round(", i, ")~", confounders_formula,"+", "offset(log(ReadCount))")), 
                                            data = modeldata, control = glm.control(maxit = 1000,epsilon=1e-12)),
                               error = function(e) NULL)
        
        if(length(model[[i]])==0 | is.null(model_test)){
          model[[i]] <- tryCatch(MASS::glm.nb(formula(paste("round(", i, "+1)~", confounders_formula, "+ ",group,"+", "offset(log(ReadCount))")), 
                                              data = modeldata, control = glm.control(maxit = 1000,epsilon=1e-12)),
                                 error = function(e) NULL)
        }
        if(length(model[[i]])==0){
          model[[i]] <- tryCatch(stats::glm(formula(paste("round(", i, ")~", confounders_formula, "+ ",group,"+", "offset(log(ReadCount))")), family=poisson,
                                            data = modeldata,control = glm.control(maxit = 1000,epsilon=1e-12)),
                                 error = function(e) NULL)
        }
        if(length(model[[i]])==0){
          if(abs(mean(modeldata[,i]/modeldata$ReadCount)-median(modeldata[,i]/modeldata$ReadCount))/median(modeldata[,i]/modeldata$ReadCount)<0.25){
            model[[i]] <- tryCatch(stats::lm(formula(paste("(", i, "+1)/ReadCount~", confounders_formula, "+ ",group)), 
                                             data = modeldata),
                                   error = function(e) NULL)   
          } else {
            model[[i]] <- tryCatch(stats::lm(formula(paste("log(", i, "+1)/ReadCount~", confounders_formula, "+ ",group)), 
                                             data = modeldata),
                                   error = function(e) NULL)   
          }
        }  
        if(length(model[[i]])>0){
          tmp <- data.frame(res = stats::resid(model[[i]],type="pearson"), fitted = stats::predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
          
          if(summary(stats::lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(stats::lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             stats::anova(stats::lm(res ~ group, data=tmp))$Pr[1]>0.1 & stats::anova(stats::lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$coef),pattern=group,replacement=""))){
              group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)] }   
            
          } else {
            
            if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
              if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){   
                model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                                       error = function(e) NULL)
              } else {
                model[[i]] <- tryCatch(nlme::gls(formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, weights = nlme::varExp(), control = nlme::glsControl(maxIter=5000)),
                                       error = function(e) NULL) 
              }
              
            } else {
              if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){    
                model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, weights = nlme::varIdent(form = formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                                       error = function(e) NULL)
              } else(
                model[[i]] <- tryCatch(nlme::gls(formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                 data = modeldata2, weights = nlme::varIdent(form = formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                                       error = function(e) NULL)
              )
            }    
            if(length(model[[i]])>0){
              tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
              tmp$resdev <- abs(tmp$res) 
              tmp$group <- modeldata[,group]
              if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
                 anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
                  group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,4)]}
              } else {
                if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){    
                  model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                   data = modeldata2, 
                                                   weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = formula(paste(" ~1|",group)))), 
                                                   control = nlme::glsControl(maxIter=1000)),
                                         error = function(e) NULL)
                } else{
                  model[[i]] <- tryCatch(nlme::gls(formula(paste("log(", i, "/ReadCount)~", confounders_formula, "+ ",group)), 
                                                   data = modeldata2, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = formula(paste(" ~1|",group)))), 
                                                   control = nlme::glsControl(maxIter=1000)),
                                         error = function(e) NULL)
                }
                
                if(length( model[[i]])>0){
                  tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
                  tmp$resdev <- abs(tmp$res) 
                  tmp$group <- modeldata[,group]
                  if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
                     anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                    for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
                      group_test[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,4)]
                    } 
                  }  
                }
              }
            }
          }
          if(length(model[[i]])>0){
            group_test[i, "model"] <- strsplit(as.character(summary(model[[i]])$call),split="(",fixed=T)[[1]][1] 
            if(group_test[i, "model"]== "lm" | group_test[i, "model"]== "nlme::gls"){
              if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)>0.24999){
                group_test[i, "model"] <- paste("log",group_test[i, "model"])
              } 
            }
            
            for(j in levels(modeldata[,group])[-1]){
              
              if(group_test[i, "model"]== "MASS::glm.nb" | group_test[i, "model"]== "glm"){     
                group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
                if(length(confounders[confounders!=""])>0){
                  assign("i", paste(i), envir = .GlobalEnv)
                  
                  resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                                       init.theta = 5, control = glm.control(epsilon = 1e-12, maxit = 2500, trace = T)), type="pearson"),
                                          error = function(e) NA)
                  
                  if(is.na(sum(resids[,i]))){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                                         control = glm.control(epsilon = 1e-20, maxit = 1000, trace = T)), type="pearson"),
                                            error = function(e) NA)
                    
                  }
                  if(is.na(sum(resids[,i]))){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                                         init.theta = 10,control = glm.control(epsilon = 1e-20, maxit = 2500, trace = T)), type="pearson"),
                                            error = function(e) NA)
                    
                  }
                }
              }
              
              if(group_test[i, "model"] == "nlme::gls"){     
                group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(fitted(model[[i]])[modeldata$G==j]), 
                                                                                           mean(fitted(model[[i]])[modeldata$G==compare.to]))#(summary(model[[i]])$tTable[1,1]+summary(model[[i]])$tTable[paste(group,j,sep=""),1])/summary(model[[i]])$tTable[1,1]
                if(length(confounders[confounders!=""])>0){
                  assign("i", paste(i), envir = .GlobalEnv)
                  
                  resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                          error = function(e) NA)
                }   
              } 
              
              if(group_test[i, "model"] == "log nlme::gls") {     
                group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(exp(fitted(model[[i]])[modeldata$G==j])), 
                                                                                           mean(exp(fitted(model[[i]])[modeldata$G==compare.to])))#exp(summary(model[[i]])$tTable[1,1]+(summary(model[[i]])$tTable[paste(group,j,sep=""),1]))/exp(summary(model[[i]])$tTable[1,1])
                if(length(confounders)>0){
                  assign("i", paste(i), envir = .GlobalEnv)
                  
                  resids[,i] <- tryCatch(resid(update(model[[i]], as.formula(paste(".~. -", group, sep=""))), type="pearson"),
                                         error = function(e) NA)
                }   
              } 
              
              if(group_test[i, "model"] == "lm"){
                
                group_test[i, paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(fitted(model[[i]])[modeldata$G==j]), 
                                                                                            mean(fitted(model[[i]])[modeldata$G==compare.to]))#(summary(model[[i]])$coef[1,1]+summary(model[[i]])$coef[paste(group,j,sep=""),1])/summary(model[[i]])$coef[1,1]
                if(length(confounders[confounders!=""])>0){
                  assign("i", paste(i), envir = .GlobalEnv)
                  
                  resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                          error = function(e) NA)
                }   
              } 
              
              if(group_test[i, "model"] == "log lm"){
                assign("i", paste(i), envir = .GlobalEnv)
                
                group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(exp(fitted(model[[i]])[modeldata$G==j])), 
                                                                                           mean(exp(fitted(model[[i]])[modeldata$G==compare.to])))#(exp(summary(model[[i]])$coef[1,1]+summary(model[[i]])$coef[paste(group,j,sep=""),1]))/exp(summary(model[[i]])$coef[1,1])
                if(length(confounders[confounders!=""])>0){
                  assign("i", paste(i), envir = .GlobalEnv)
                  
                  resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                          error = function(e) NA)
                  
                  
                }  
              } 
            }
          }
        }
      }
      
      if(nonzero){ 
        
        residsNonzero <- modeldata
        
        for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
          if(min(table(modeldata[modeldata[,i]>0,group]))>3){    
            confounders2 <- confounders
            for(h in seq_along(confounders[confounders!=""])) if(min(table(modeldata[modeldata[,i]>0,confounders[h]]))==0) confounders2[h] <- ""
            confounders_formula2 <- paste(confounders2[which(!confounders2 == "")], collapse = " + ")
            assign("confounders_formula2", confounders_formula2, envir = .GlobalEnv)
            model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, ")~", confounders_formula2, "+ ",group,"+", "offset(log(ReadCount))")), 
                                                                 data = modeldata[modeldata[,i]>0,], control = glm.control(maxit = 1000)),
                                                    error = function(e) NULL)
            if(length(model[[paste(i,"nonzero")]])==0){
              model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, "+1)~", confounders_formula2, "+ ",group,"+", "offset(log(ReadCount))")), 
                                                                   data = modeldata[modeldata[,i]>0,], control = glm.control(maxit = 1000)),
                                                      error = function(e) NULL)
            }
            if(length(model[[paste(i,"nonzero")]])==0){
              model[[i]] <- tryCatch(glm(formula(paste("round(", i, ")~", confounders_formula, "+ ",group,"+", "offset(log(ReadCount))")), family=poisson,
                                         data = modeldata[modeldata[,i]>0,]),
                                     error = function(e) NULL)
            }
            if(length(model[[paste(i,"nonzero")]])==0){
              if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){
                model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("(", i, "+1)/ReadCount~", confounders_formula2, "+ ",group)), 
                                                           data = modeldata[modeldata[,i]>0,]),
                                                        error = function(e) NULL)   
              } else {
                model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("log(", i, "+1)/ReadCount~",confounders_formula2, "+ ",group)), 
                                                           data = modeldata[modeldata[,i]>0,]),
                                                        error = function(e) NULL)   
              }
            }  
            
            tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
            tmp$resdev <- abs(tmp$res) 
            tmp$group <- modeldata[modeldata[,i]>0,group]
            
            if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]!="NaN" & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]!="NaN" &
               summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
               anova(lm(res ~ group, data=tmp))$Pr[1]!="NaN" & anova(lm(resdev ~ group, data=tmp))$Pr[1]!="NaN" &
               anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
              
              for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$coef),pattern=group,replacement=""))){
                group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)] }   
              
            } else {
              
              if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
                if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~", confounders_formula2, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                } else {
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula2, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL) 
                }
                
              } else if(anova(lm(res ~ group, data=tmp))$Pr[1]<0.01 | anova(lm(resdev ~ group, data=tmp))$Pr[1]<0.01){
                if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){  
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~", confounders_formula2, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form = as.formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                } else(
                  model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula2, "+ ",group)), 
                                                                    data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form = as.formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                                                          error = function(e) NULL)
                )
              }    
              if(length( model[[paste(i,"nonzero")]])>0){
                tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
                tmp$resdev <- abs(tmp$res) 
                tmp$group <- modeldata[modeldata[,i]>0,group]
                if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
                   anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                  for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
                    group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,4)]}
                } else {
                  if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){   
                    model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~",confounders_formula2, "+ ",group)), 
                                                                      data = modeldata[modeldata[,i]>0,], 
                                                                      weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = as.formula(paste(" ~1|",group)))), 
                                                                      control = nlme::glsControl(maxIter=1000)),
                                                            error = function(e) NULL)
                  } else{
                    model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log(", i, "/ReadCount)~", confounders_formula2, "+ ",group)), 
                                                                      data = modeldata[modeldata[,i]>0,], weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = as.formula(paste(" ~1|",group)))), 
                                                                      control = nlme::glsControl(maxIter=1000)),
                                                            error = function(e) NULL)
                  }
                  
                  if(length( model[[paste(i,"nonzero")]])>0){
                    tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
                    tmp$resdev <- abs(tmp$res) 
                    tmp$group <- modeldata[modeldata[,i]>0,group]
                    if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
                       anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
                      for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
                        group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste("pvals", j, sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,4)]
                      } 
                    }  
                  }
                }
              }
            }    
            
            if(length(model[[paste(i,"nonzero")]])>0){
              group_test_nonzero[i, "model"] <- strsplit(as.character(summary(model[[paste(i,"nonzero")]])$call),split="(",fixed=T)[[1]][1] 
              if(group_test_nonzero[i, "model"]== "lm" | group_test_nonzero[i, "model"]== "nlme::gls"){
                if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])>0.24999){
                  group_test_nonzero[i, "model"] <- paste("log",group_test_nonzero[i, "model"])
                } 
              }    
              
              for(j in levels(modeldata[,group])[-1]){
                if(group_test_nonzero[i, "model"]== "MASS::glm.nb" |group_test_nonzero[i, "model"]== "glm"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
                  
                  if(length(confounders[confounders!=""])>0){
                    assign("i", paste(i), envir = .GlobalEnv)
                    try( residsNonzero[modeldata[,i]>0,i] <-  resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                                                           init.theta = 5, control = glm.control(epsilon = 1e-12, maxit = 2500, trace = T)), type="pearson"))#,
                    #error = function(e) NA)
                    if(is.na(sum(residsNonzero[,i]))){
                      assign("i", paste(i), envir = .GlobalEnv)
                      
                      residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                                                                 control = glm.control(epsilon = 1e-20, maxit = 1000, trace = T)), type="pearson"),
                                                                    error = function(e) NA)
                    }
                    if(is.na(sum(residsNonzero[,i]))){
                      assign("i", paste(i), envir = .GlobalEnv)
                      
                      residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                                                                 init.theta = 10,control = glm.control(epsilon = 1e-20, maxit = 2500, trace = T)), type="pearson"),
                                                                    error = function(e) NA)
                      
                    }  }  } 
                
                if(group_test_nonzero[i, "model"] == "nlme::gls"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0]), 
                                                                                                     mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0]))#(summary(model[[paste(i,"nonzero")]])$tTable[1,1]+summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1])/summary(model[[paste(i,"nonzero")]])$tTable[1,1]
                  if(length(confounders[confounders!=""])>0){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                                                  error = function(e) NA)
                  }   } 
                
                if(group_test_nonzero[i, "model"] == "log nlme::gls"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0])), 
                                                                                                     mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0])))#exp(summary(model[[paste(i,"nonzero")]])$tTable[1,1]+(summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1]))/exp(summary(model[[paste(i,"nonzero")]])$tTable[1,1])
                  if(length(confounders[confounders!=""])>0){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                                                  error = function(e) NA)
                  }   } 
                
                if(group_test_nonzero[i, "model"] == "lm"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <-  gtools::foldchange(mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0]), 
                                                                                                      mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0]))#(summary(model[[paste(i,"nonzero")]])$coef[1,1]+summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),1])/summary(model[[paste(i,"nonzero")]])$coef[1,1]
                  if(length(confounders[confounders!=""])>0){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                                                  error = function(e) NA)
                  }   } 
                
                if(group_test_nonzero[i, "model"] == "log lm"){     
                  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- gtools::foldchange(mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0])), 
                                                                                                     mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0])))#(exp(summary(model[[paste(i,"nonzero")]])$coef[1,1]+summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),1]))/exp(summary(model[[paste(i,"nonzero")]])$coef[1,1])
                  if(length(confounders[confounders!=""])>0){
                    assign("i", paste(i), envir = .GlobalEnv)
                    
                    residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                                                                  error = function(e) NA)
                  }  } 
              }
              
              
              
              
            }
            
            if(length(confounders2[confounders2!=""])>0){
              assign("i", paste(i), envir = .GlobalEnv)
              
              residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update( model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), 
                                                                  type="pearson"),
                                                            error = function(e) NA)
            }      
          }     
          
        }
      }
      
    }
    
    if(nrow(na.omit(group_test))>0){
      
      for (k in names(group_test)[grepl(pattern="^pvals",x=names(group_test))]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,k], "fdr")
      tmp <- group_test
      tmp[, -1][is.na(tmp[, -1])]<-1 
      sig <- as.character(rownames(tmp)[sapply(data.frame(t(tmp[, grepl(pattern="_FDR$",x=names(group_test))])), min,na.rm=T) < p.cutoff])
      
      if(nonzero){
        for (k in names(group_test_nonzero)[grepl(pattern="^pvals",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test_nonzero[,k], "fdr")
        tmp2 <- group_test_nonzero
        tmp2[, -1][is.na(tmp2[, -1])]<-1 
        sigNonzero <- as.character(rownames(tmp2)[sapply(data.frame(t(tmp2[, grepl(pattern="_FDR$",x=names(group_test_nonzero))])), min,na.rm=T) < p.cutoff & sapply(data.frame(t(tmp2[, grepl(pattern="_FDR$",x=names(group_test_nonzero))])), min,na.rm=T) < sapply(data.frame(t(tmp[, grepl(pattern="_FDR$",x=names(group_test))])), min,na.rm=T)])
      }  
      
      if(length(confounders[confounders!=""])>0){
        
        if(nonzero){ 
          
          names(group_test_nonzero)[-1] <-  paste("nonzero",names(group_test_nonzero)[-1],sep="_")
          group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)
          rownames(group_test)<-group_test$taxon
          
          resids2 <- data.frame(cbind(resids[,setdiff(colnames(resids),colnames(taxa))],
                                      resids[,setdiff(sig,sigNonzero)],
                                      residsNonzero[,sigNonzero]))
          colnames(resids2)<-c(setdiff(colnames(resids),colnames(taxa)),setdiff(sig,sigNonzero),sigNonzero)
        } else {
          sigNonzero <- NULL  
          resids2 <- data.frame(cbind(resids[,setdiff(colnames(resids),colnames(taxa))],resids[,sig]))
          colnames(resids2)<-c(setdiff(colnames(resids),colnames(taxa)),sig)
        }
        
        ##### function to remove outliers ##-----
        normalizeRange <- function(x){
          Q1 <- quantile(x, 0.25)
          Q3 <- quantile(x, 0.75)
          iqR <- IQR(x)
          upper.limit <- Q3 + (1.5*iqR)
          lower.limit <- Q1 - (1.5*iqR)
          x <- x %>% 
            as.data.frame() %>%
            mutate(across(.cols = everything(), ~if_else(.x < lower.limit, lower.limit, 
                                                         if_else(.x > upper.limit, upper.limit, .x)))) %>%
            as.matrix()
        }
        
        dataset2 <- data.frame(resids2) %>%
          dplyr::select(where(~all(!is.na(.x)))) %>%
          dplyr::mutate(across(.cols = any_of(colnames(taxa)), .fns = ~normalizeRange(.x)))
        dataset2[,group] <- dataset2$G
        
      } else {    
        for(i in levels(modeldata$G)){
          for(j in group_test$taxon){
            group_test[j,paste("Mean",i,sep="_")] <-   mean(modeldata[modeldata$G==i,j]/modeldata$ReadCount[modeldata$G==i],na.rm=T)
          }}
        for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
          group_test[,paste("FoldChange",i,sep="_")] <- gtools::foldchange(group_test[,paste("Mean",i,sep="_")],
                                                                           group_test[,paste("Mean",compare.to,sep="_")])
        }
        
        if(nonzero){
          for(i in levels(modeldata$G)){
            for(j in group_test_nonzero$taxon){
              group_test_nonzero[j,paste("Mean",i,sep="_")] <-   mean(modeldata[modeldata$G==i&modeldata[,j]>0,j]/modeldata$ReadCount[modeldata$G==i&modeldata[,j]>0],na.rm=T)
            }}
          
          for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
            group_test_nonzero[,paste("FoldChange",i,sep="_")] <- gtools::foldchange(group_test_nonzero[,paste("Mean",i,sep="_")],
                                                                                     group_test_nonzero[,paste("Mean",compare.to,sep="_")])
          }
          names(group_test_nonzero)[-1] <- paste("nonzero",names(group_test_nonzero)[-1], sep="_") 
          group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)
        } else { 
          sigNonzero <- NULL
        }
        names(group_test)[grepl(pattern="group0",x=names(group_test))] <-gsub(pattern="group0",replacement="0",x=names(group_test)[grepl(pattern="group0",x=names(group_test))])
        
        
        dataset2 <- modeldata
        if(relative) dataset2[,names(taxa)]<-(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
        if(nonzero) for(i in setdiff(sigNonzero,sig)) dataset2[modeldata[,i]==0,i] <- NA
      }
      
      group_test <- group_test %>%
        dplyr::arrange(taxon) %>%
        dplyr::select(taxon, 
                      starts_with("Mean"),
                      starts_with("Fold"),
                      ends_with("estimate"),
                      contains("FoldChange"),
                      starts_with("pvals"),
                      ends_with("FDR"),
                      model)
      
      
      
      writexl::write_xlsx(group_test, 
                          path = 
                            here::here(dir_for_res, "tables",
                                       paste("GroupTest_", group, "_", compare.to, "_", 
                                             select.by, "_", select, 
                                             "_prev_", min.prevalence*100, "_abd_", 
                                             min.abundance*100, ".xlsx", sep = "")))
      
      if(length(confounders[confounders!=""])>0){  
        writexl::write_xlsx(resids, 
                            path = 
                              here::here(dir_for_res, "tables",
                                         paste("GroupTest_residuals_", group, "_", 
                                               compare.to, "_", select.by, "_", select,
                                               "_prev_", min.prevalence*100, 
                                               "_abd_", min.abundance*100,
                                               ".xlsx", sep = "")))
      }
      #------------
      sig <- unique(c(sig,sigNonzero))
      sig <- sig[sig %in% colnames(dataset2)]
      if(length(sig)>0){    
        
        #------------
        if (length(sig) > 1) {  
          sig <- na.omit(sig[apply(dataset2[,sig],MARGIN=2,FUN=sum,na.rm=T)!=0])
        }
        
        lphy <- length(sig)
        
        if(lphy < 5) {
          ncols <- lphy
          nrows <- 1
        } else {
          ncols <-  round(sqrt(lphy)) + 1
          nrows <-  floor(sqrt(lphy))
        }
        if(relative) yaxis <- round((100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))) else yaxis <- round((seq(0,max(dataset[,sig]),max(dataset[,sig])/10)))
        
        if (pdf) {
          pdf(width=ncols*4,height=nrows*4, 
              file = here::here(dir_for_res, "plots",
                                paste("GroupTest_", group, "_", compare.to, "_", select.by, 
                                      "_", select, 
                                      "_prev_", min.prevalence*100, "_abd_", 
                                      min.abundance*100, "_Boxplot.pdf", sep = "")))
          par(mfrow = c(nrows, ncols))
          for(i in sig[order(sig)]){
            if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
              if(length(lab)==0) lab <- "Unassigned taxa"
            } else{
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
              if(lab==" deviance") lab <- "Unassigned taxa deviance"
            }
            boxplot(dataset2[, i] ~ dataset2[, "G"], 
                    ylab="", 
                    main =  lab, 
                    xlab=group_name, 
                    las = label.direction, #yaxt="n",
                    outline = F,
                    col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
                            "#A65628", "#F781BF", "#999999","blue","firebrick4",
                            'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
                            'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
                    outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
                                           "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
                                           'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
          }
          dev.off()
          
          pdf(width=ncols*4,height=nrows*4, 
              file = here::here(dir_for_res, "plots",
                                paste("GroupTest_", group, "_", compare.to, "_", select.by, "_", 
                                      select,
                                      "_prev_", min.prevalence*100, "_abd_", 
                                      min.abundance*100, "_Beanplot.pdf", sep = "")))
          par(mfrow = c(nrows, ncols))
          for(i in sig[order(sig)]){
            if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
              if(length(lab)==0) lab <- "Unassigned taxa"
            } else {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
              if(lab==" deviance") lab <- "Unassigned taxa deviance"
            }
            tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], 
                                        xlab=group_name, 
                                        las = label.direction,#yaxt="n",
                                        ll = 0.1, 
                                        ylab = "", 
                                        main=lab, 
                                        beanlines="median",
                                        col=list(c('#E41A1C','black','black','black'),
                                                 c('orange','black','black','black'),
                                                 c('#377EB8','black','black','black'),
                                                 c('skyblue','black','black','black'),
                                                 c("#4DAF4A",'black','black','black'),
                                                 c('#984EA3','black','black','black'),
                                                 c('#FFFF33','black','black','black'),
                                                 c('#A65628','black','black','black'),
                                                 c('#F781BF','black','black','black'),
                                                 c('#999999','black','black','black'),
                                                 c('blue','black','black','black'),
                                                 c('firebrick4','black','black','black'),
                                                 c('yellowgreen','black','black','black'),
                                                 c('pink','black','black','black'),
                                                 c('turquoise2','black','black','black'),
                                                 c('plum','black','black','black'),
                                                 c('darkorange','black','black','black'),
                                                 c('lightyellow','black','black','black'),
                                                 c('gray','black','black','black')),
                                        border = "black"), error = function(e) NULL)
          } 
          dev.off()
        }
        
        if(show_quartz) {
          quartz()
          par(mfrow = c(nrows, ncols))
          for(i in sig[order(sig)]){
            if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
              if(length(lab)==0) lab <- "Unassigned taxa"
            } else{
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
              if(lab==" deviance") lab <- "Unassigned taxa deviance"
            }
            boxplot(dataset2[, i] ~ dataset2[, "G"], 
                    ylab="", 
                    main =  lab, 
                    xlab=group_name, 
                    las = label.direction, #yaxt="n",
                    outline = F,
                    col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
                            "#A65628", "#F781BF", "#999999","blue","firebrick4",
                            'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
                            'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
                    outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
                                           "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
                                           'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
            #  if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=exp(yaxis)) } else  axis(side=2)
          }
          quartz()
          par(mfrow = c(nrows, ncols))
          for(i in sig[order(sig)]){
            if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
              if(length(lab)==0) lab <- "Unassigned taxa"
            } else {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
              if(lab==" deviance") lab <- "Unassigned taxa deviance"
            }
            tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], 
                                        xlab=group_name, 
                                        las = label.direction,#yaxt="n",
                                        ll = 0.1, ylab = "", 
                                        main=lab, beanlines="median",
                                        col=list(c('#E41A1C','black','black','black'),
                                                 c('orange','black','black','black'),
                                                 c('#377EB8','black','black','black'),
                                                 c('skyblue','black','black','black'),
                                                 c("#4DAF4A",'black','black','black'),
                                                 c('#984EA3','black','black','black'),
                                                 c('#FFFF33','black','black','black'),
                                                 c('#A65628','black','black','black'),
                                                 c('#F781BF','black','black','black'),
                                                 c('#999999','black','black','black'),
                                                 c('blue','black','black','black'),
                                                 c('firebrick4','black','black','black'),
                                                 c('yellowgreen','black','black','black'),
                                                 c('pink','black','black','black'),
                                                 c('turquoise2','black','black','black'),
                                                 c('plum','black','black','black'),
                                                 c('darkorange','black','black','black'),
                                                 c('lightyellow','black','black','black'),
                                                 c('gray','black','black','black')),
                                        border = "black"), error = function(e) NULL)
            #   if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=signif(exp(yaxis),digits=1)) } else  axis(side=2)
          }
        }
      }
    
      
      assign("model", model, .GlobalEnv)
      rm(list = c("i", "group", "confounders_formula", "confounders_formula2"), envir = .GlobalEnv)
      
      if(keep.result) return(group_test)
    }
  }
}

