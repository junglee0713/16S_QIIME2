library(tidyverse)
library(ggsci)

###=============
### permutation based p-value
### originally from Kyle
###=============

permutation_pvalue <- function (observed_stat, permutation_stats) {
  (sum(permutation_stats >= observed_stat) + 1) / (length(permutation_stats) + 1)
}

###=============
### detect outlier
### outlier definition consistent with ggplot
###=============

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

###=============
### given a taxa count matrix, identify most dominant taxa based on mean abundance
### least_prevalence: taxa should be present in 100 x least_prevalence % of samples
### top_taxa_n (<20): how many taxa to show? all other taxa are aggregated as "Other" 
### returns a pie chart, a list of top taxa and aggregated count matrix
###=============

getPrevalentTaxa <- function(summed_counts, least_prevalence = 0.1, top_taxa_n = 9) {
  
  taxa_prevalence <- apply(summed_counts, 1, function(x) sum(x > 0) / length(x)) 
  prevalent_taxa <- names(taxa_prevalence)[taxa_prevalence > least_prevalence]
  
  mean_prop <- sweep(summed_counts, 2, colSums(summed_counts), "/") %>%
    rowMeans() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    setNames(c("Taxa", "MeanProp")) %>%
    filter(Taxa %in% prevalent_taxa) %>%
    arrange(desc(MeanProp))
  
  top_taxa <- mean_prop$Taxa[1:top_taxa_n]
  
  mean_prop <- mean_prop %>%
    mutate(newTaxa = ifelse(Taxa %in% top_taxa, Taxa, "Other"))
  
  aggregated_mean_prop <- aggregate(MeanProp ~ newTaxa, data = mean_prop, sum) %>%
    arrange(desc(MeanProp)) %>%
    mutate(newTaxa = factor(newTaxa, levels = newTaxa)) %>%
    mutate(newTaxa = fct_relevel(newTaxa, "Other", after = Inf)) %>%
    arrange(newTaxa)
  
  pie <- ggplot(aggregated_mean_prop, aes(x = "", y = MeanProp, fill = newTaxa)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_d3(palette = "category20") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          aspect.ratio = 1) + 
    labs(fill = "Taxa", 
         x = NULL, 
         y = NULL,
         title = paste0("Mean relative abundance"))
  
  pie <- pie + 
    coord_polar(theta = "y", start = 0, direction = -1)
  
  newRowname <- ifelse(rownames(summed_counts) %in% top_taxa,
                       rownames(summed_counts),
                       "Other")
  aggregated_counts <- rowsum(summed_counts, newRowname)
  
  out <- list(pie = pie, top_taxa = top_taxa, aggregated_counts = aggregated_counts) 
  
  return(out)  
}

###=============
### Thanks to Ceylan
### Obtain tidy output from permanova test result
###=============

tidy_permanova <- function(anov){
  data.frame(Term = rownames(anov$aov.tab), anov$aov.tab, row.names = NULL) %>%
    rename(p.value = Pr..F.)
}

###=============
### Thanks to Ceylan
### Take care of PERMANOVA with repeated measure
###=============

permanova_with_shuffle_1_group <- function(dist_matrix, s_toTest, group_label, rep_mes_label, perm, is_within=F){
  dist_toTest <- dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", group_label)
  a_ixn <- adonis(as.formula(form1), data=s_toTest, permutations=perm)
  f_ixn <- a_ixn$aov.tab[1, 4]
  set.seed(1)
  
  fs_permuted <- replicate(perm, {
    s_permuted <- s_toTest
    if (is_within){
      s_permuted[,group_label] <- shuffle_within_groups(s_permuted[,group_label], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label] <- shuffle_between_groups(s_permuted[,group_label], s_permuted[,rep_mes_label])
    }
    a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
    a_permuted$aov.tab[1, 4]
  })
  p_ixn <- sum(c(f_ixn, fs_permuted) >= f_ixn) / (length(fs_permuted) + 1)
  a_ixn$aov.tab[1,6] <- p_ixn
  tidy_permanova(a_ixn)
}

###=============
### Thanks to Ceylan
### Take care of PERMANOVA when interaction is considered
###=============

permanova_with_shuffle_2_groups <- function(dist_matrix, s_toTest, group_label1, group_label2, rep_mes_label, covariates, perm, first_within=F, second_within=F){
  set.seed(1)
  dist_toTest <- dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", group_label1, " * ", group_label2)
  if (!is.na(covariates)) {
    form1 <- paste(form1, " + ", covariates)
  }
  a_ixn_orj <- adonis(as.formula(form1), data=s_toTest, permutations=perm)
  
  terms_perm <- c(group_label1, group_label2, paste0(group_label1, ":", group_label2))
  f_ixn_all <- tidy_permanova(a_ixn_orj) %>%
    filter(Term %in% terms_perm) %>%
    pull(F.Model)
  #select(Term, F.Model)
  
  fs_permuted <- replicate(perm, {
    s_permuted <- s_toTest
    
    if (first_within) {
      s_permuted[,group_label1] <- shuffle_within_groups(s_permuted[,group_label1], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label1] <- shuffle_between_groups(s_permuted[,group_label1], s_permuted[,rep_mes_label])
    }
    
    if (second_within) {
      s_permuted[,group_label2] <- shuffle_within_groups(s_permuted[,group_label2], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label2] <- shuffle_between_groups(s_permuted[,group_label2], s_permuted[,rep_mes_label])
    }
    
    a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
    
    tidy_permanova(a_permuted) %>%
      filter(Term %in% terms_perm) %>%
      pull(F.Model)
    #c(a_permuted_g1$aov.tab[1, 4], a_permuted_g2$aov.tab[1, 4], a_permuted$aov.tab[3, 4])
  })
  
  p_ixn <- rowSums(cbind(f_ixn_all, fs_permuted) >= f_ixn_all, na.rm = T) / (dim(fs_permuted)[2] + 1)
  
  tidy_output <- tidy_permanova(a_ixn_orj)
  
  tidy_output[match(terms_perm, tidy_output$Term),"p.value"] <- p_ixn
  tidy_output  
}

###=============
### Thanks to Chunyu
### Read ko assignment table
###=============

read_ko_assign <- function(filePath){
  ko <- read.table(filePath, sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                   quote = "", comment = "",row.names = 1)
  ko <- sweep(ko, 2, colSums(ko), "/")
  colnames(ko) <- colnames(ko) %>%
    sub("^PCMP_", "", ., perl=T)
  return(ko)
}

###=============
### Get pairwise adonis plots
###=============

getPairBeta <- function(pair_adonis_result, pcoa_df, grp = "study_group") {
  all_df <- NULL
  for (i in 1:nrow(pair_adonis_result)) {
    curr_pair <- as.character(pair_adonis_result$Pairs[i])
    member <- unlist(strsplit(as.character(pair_adonis_result$Pairs[i]), " vs ", fixed = T))
    df <- data.frame(curr_pair, member)
    all_df <- rbind(all_df, df)
  }
  
  out_df <- merge(pcoa_df, all_df, by.x = grp, by.y = "member") %>%
    mutate(curr_pair = gsub(" vs ", "\nvs\n", curr_pair))
  
  g <- ggplot(out_df, aes(Axis.1, Axis.2, color = get(grp))) + 
    theme_bw() +
    ggtitle(dist_name) +
    geom_point() +
    theme(aspect.ratio = 1, axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(strip.text.x = element_text(size = 7)) +
    facet_wrap(~curr_pair) +
    labs(color = grp)
}

###=============
### differential abundance test comparing two groups
###=============

diff_abundance_taxa_t_test <- function(props_toTest, row_to_test, s_toTest, group_var, nominal_p_cut = 0.05) {
  combine <- props_toTest %>%
    as.data.frame() %>%
    rownames_to_column(var = "row") %>%
    filter(row %in% row_to_test) %>%
    gather(key = SampleID, value = prop, -row) %>%
    mutate(prop = prop + 1e-10) %>%
    mutate(log_prop = log(prop)) %>%
    left_join(s_toTest, by = "SampleID")
  
  model <- combine %>%
    group_by(row) %>%
    do(out = t.test(log_prop ~ get(group_var), data = .))
  
  summaries <- lapply(1:length(model$out), 
                      function(x) {
                        out <- model$out[[x]]
                        data.frame(row = model$row[[x]],
                                   estimate = t(out$estimate),
                                   pvalue = out$p.value,
                                   stringsAsFactors = F)
                      })
  
  summaries_df <- bind_rows(summaries) 
  summaries_df <- summaries_df %>%
    setNames(gsub("estimate.mean.in.group.", "mean_", colnames(summaries_df))) %>%
    mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    mutate(Sig.Label = case_when(FDR < 0.0001 ~ "***",
                                 FDR < 0.001 ~ "**",
                                 FDR < 0.01 ~ "*", 
                                 FDR < 0.1 ~ ".",
                                 T ~ "")) %>%
    filter(pvalue < nominal_p_cut) %>%
    arrange(pvalue) %>%
    rename(`$p$-value` = pvalue)
  
  return_list <- list(summaries_df = summaries_df, 
                      combine = combine, 
                      group_var = group_var)
  return(return_list)
}

###=============
### draw box plots after differential abundance test
### landscape mode
### ```{r, fig.height = 7, fig.width = 9}
###=============

draw_6_plots_per_page <- function(res, convert_strip_df = NULL, jitter = T) {
  row_to_plot <- unique(res$summaries_df$row)
  num_plot_groups <- ceiling(length(row_to_plot)/6)
  prop_df <- res$combine %>%
    mutate(row = factor(row, levels = row_to_plot))
  
  if (!is.null(convert_strip_df)) {
    prop_df <- left_join(prop_df, convert_strip_df, by = "row") %>%
      mutate(rowStrip = paste0(row, "\n", Taxa))
  } else {
    prop_df <- prop_df %>%
      mutate(rowStrip = row)
  }
  
  if (num_plot_groups > 0) {
    for (i in 1:num_plot_groups) {
      row_group_to_plot <- row_to_plot[(6*(i-1)+1):min(6*i, length(row_to_plot))]
      g <- prop_df %>%
        filter(row %in% row_group_to_plot) %>%
        ggplot(aes(x = get(res$group_var), y = prop)) 
        
      if (jitter == T) {
        g <- g + geom_boxplot(outlier.alpha = 0) + geom_jitter(width = 0.2, height = 0)
      } else {
        g <- g + geom_boxplot()
      }
      
      g <- g +  
        theme(aspect.ratio = 1) +
        scale_y_log10() +
        labs(x = res$group_var) +
        facet_wrap(~rowStrip, scales = "free_y", ncol = 3) +
        theme(strip.text.x = element_text(size = 6)) 
      
      ## Determine x.axis.text angle
      LabelTooMany <- length(unique(prop_df[, res$group_var])) > 3 # more than 3 labels
      LabelTooLong <- max(nchar(as.character(unique(prop_df[, res$group_var])))) > 9 # one label has more than 9 characters
      
      if (LabelTooMany | LabelTooLong) {
        g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
    
      print(g)
    }
  }
}

###=============
### draw box plots after differential abundance test
### landscape mode
### ```{r, fig.height = 7, fig.width = 9}
###=============

draw_6 <- function(summaries_df, props_toTest, group_var = "genotype", jitter = T) {
  row_to_plot <- unique(summaries_df$Taxa)
  num_plot_groups <- ceiling(length(row_to_plot)/6)
  prop_df <- props_toTest
  
  if (num_plot_groups > 0) {
    for (i in 1:num_plot_groups) {
      row_group_to_plot <- row_to_plot[(6*(i-1)+1):min(6*i, length(row_to_plot))]
      g <- prop_df %>%
        filter(Taxa %in% row_group_to_plot) %>%
        mutate(Taxa = factor(Taxa, levels = unique(Taxa))) %>%
        ggplot(aes(x = get(group_var), y = prop)) 
      
      if (jitter == T) {
        g <- g + geom_boxplot(outlier.alpha = 0) + geom_jitter(width = 0.2, height = 0)
      } else {
        g <- g + geom_boxplot()
      }
      
      g <- g +  
        theme(aspect.ratio = 1) +
        scale_y_log10() +
        labs(x = group_var) +
        facet_wrap(~Taxa, scales = "free_y", ncol = 3) +
        theme(strip.text.x = element_text(size = 6)) 
      
      ## Determine x.axis.text angle
      LabelTooMany <- length(unique(prop_df[, group_var])) > 3 # more than 3 labels
      LabelTooLong <- max(nchar(as.character(unique(prop_df[, group_var])))) > 9 # one label has more than 9 characters
      
      if (LabelTooMany | LabelTooLong) {
        g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
      
      print(g)
    }
  }
}

###=============
### differential abundance test comparing multiple groups
###=============

diff_abundance_taxa_pairwise_t_test <- function(props_toTest, row_to_test, s_toTest, group_var, nominal_p_cut = 0.05) {
  combine <- props_toTest %>%
    as.data.frame() %>%
    rownames_to_column(var = "row") %>%
    filter(row %in% row_to_test) %>%
    gather(key = SampleID, value = prop, -row) %>%
    mutate(prop = prop + 1e-10) %>%
    mutate(log_prop = log(prop)) %>%
    left_join(s_toTest, by = "SampleID")
  
  model <- combine %>%
    mutate(GroupVarPairwise = get(group_var)) %>%
    group_by(row) %>%
    do(out = pairwise.t.test(.$log_prop, .$GroupVarPairwise, p.adjust.method = "none"))
  
  summaries <- lapply(1:length(model$out), 
                      function(x) 
                        data.frame(
                          row = model$row[[x]],
                          model$out[[x]]$p.value,
                          Group1 = rownames(model$out[[x]]$p.value),
                          check.names = F,
                          stringsAsFactors = F
                        )) 
  
  summaries_df <- bind_rows(summaries) 
  summaries_df <- summaries_df %>%
    gather(key = "Group2", value = "$p$-value", -row, -Group1) %>%
    filter(!is.na(`$p$-value`)) %>%
    mutate(FDR = p.adjust(`$p$-value`, method = "BH")) %>%
    mutate(Sig.Label = case_when(FDR < 0.0001 ~ "***",
                                 FDR < 0.001 ~ "**",
                                 FDR < 0.01 ~ "*", 
                                 FDR < 0.1 ~ ".",
                                 T ~ "")) %>%
    filter(`$p$-value` < nominal_p_cut) %>%
    arrange(`$p$-value`) 
  
  return_list <- list(summaries_df = summaries_df, 
                      combine = combine, 
                      group_var = group_var)
  return(return_list)
}

###=============
### get pretty heatmap
###=============

get_pheatmap <- function(res, anno_color = NULL, grps = NULL, satu_limit = 0.4) {
  row_to_plot <- unique(res$summaries_df$row)
  if (is.null(grps)) {
    grps <- c("Sampling_Site", res$group_var)
  }
  s_Heat <- res$combine %>%
    select(SampleID, grps) %>%
    unique() %>%
    arrange_(.dots = grps)
  
  anno <- s_Heat %>% select(-SampleID)
  rownames(anno) <- s_Heat$SampleID
  colnames(anno) <- grps
  
  heatmap_prop_df <- res$combine %>%
    filter(row %in% row_to_plot) %>%
    select(row, SampleID, prop) %>%
    spread(key = SampleID, value = prop)
  
  heatmap_prop <- heatmap_prop_df[,rownames(anno)]
  rownames(heatmap_prop) <- heatmap_prop_df$row
  
  color = saturated_rainbow(101, saturation_limit = satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  g <- pheatmap(heatmap_prop, annotation = anno, cluster_cols = F, cluster_rows = F, 
                color = color, breaks = breaks,  
                annotation_colors = anno_color, 
                cellwidth = 8, cellheight = 8, fontsize_col = 8, fontsize_row = 8)
  
  return(g)
}

###=============
### get pretty heatmap for co-occurrence analysis
###=============

heatmap_co_occur <- function(summed_COUNT_16S, summed_COUNT_ITS, s_Heat, grps, anno_color = NULL,
                             prop_cut = 0.01, prev_cut = 0.4, satu_limit = 0.4) {
  
  summed_PROP_16S <- sweep(summed_COUNT_16S, 2, colSums(summed_COUNT_16S), "/")
  summed_PROP_ITS <- sweep(summed_COUNT_ITS, 2, colSums(summed_COUNT_ITS), "/")
  
  s_Heat <- s_Heat %>%
    select(SampleID, grps) %>%
    arrange_(.dots = grps)
  
  anno <- s_Heat %>% select(-SampleID)
  rownames(anno) <- s_Heat$SampleID
  colnames(anno) <- grps
  
  # which 16S taxa to show?
  subset_PROP_16S <- summed_PROP_16S[, rownames(anno)]
  taxa_prop_cut_pass <- rownames(subset_PROP_16S)[apply(subset_PROP_16S, 1, max) >= prop_cut]
  taxa_prev_cut_pass <- rownames(subset_PROP_16S)[apply(subset_PROP_16S, 1, function(x) {sum(x>0)/length(x)}) >= prev_cut]
  taxa_pass <- intersect(taxa_prop_cut_pass, taxa_prev_cut_pass)
  heat_16S_df <- subset_PROP_16S[taxa_pass, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Taxa") %>%
    mutate(Label = "bacteria")
  count_16S_df <- summed_COUNT_16S[taxa_pass, rownames(anno)]
  
  # which ITS taxa to show?
  subset_PROP_ITS <- summed_PROP_ITS[, rownames(anno)]
  taxa_prop_cut_pass <- rownames(subset_PROP_ITS)[apply(subset_PROP_ITS, 1, max) >= prop_cut]
  taxa_prev_cut_pass <- rownames(subset_PROP_ITS)[apply(subset_PROP_ITS, 1, function(x) {sum(x>0)/length(x)}) >= prev_cut]
  taxa_pass <- intersect(taxa_prop_cut_pass, taxa_prev_cut_pass)
  heat_ITS_df <- subset_PROP_ITS[taxa_pass, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Taxa") %>%
    mutate(Label = "fungi")
  count_ITS_df <- summed_COUNT_ITS[taxa_pass, rownames(anno)]
  
  heat_comb_df <- rbind(heat_16S_df, heat_ITS_df)
  heat_comb <- heat_comb_df %>%
    select(-Taxa, -Label) %>%
    as.matrix()
  rownames(heat_comb) <- heat_comb_df$Taxa
  anno_row = heat_comb_df %>%
    select(Label)
  rownames(anno_row) <- heat_comb_df$Taxa
  
  color = saturated_rainbow(101, saturation_limit = satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  g <- pheatmap(heat_comb, annotation_col = anno, annotation_row = anno_row,
                color = color, breaks = breaks, 
                annotation_colors = anno_color,
                fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,
                cellheight = 8, cellwidth = 8, gaps_row = nrow(heat_16S_df))
  
  return(list(g = g, count_16S_df = count_16S_df, count_ITS_df = count_ITS_df))
}

###=============
### get heatmap for co-occurence indices
###=============

co_occur_index <- function(res, dist_name = "Dice") {
  require(ade4)
  tab_16S <- res$count_16S_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "Taxa") %>% 
    gather(key = SampleID, value = Count, -Taxa) %>%
    mutate(Label = "bacteria")
  
  tab_ITS <- res$count_ITS_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "Taxa") %>% 
    gather(key = SampleID, value = Count, -Taxa) %>%
    mutate(Label = "fungi")
  
  tab_comb <- rbind(tab_16S, tab_ITS) %>%
    select(-Label) %>%
    spread(key = Taxa, value = Count)
  rownames(tab_comb) <- tab_comb$SampleID
  tab_comb <- tab_comb %>% 
    select(-SampleID)  
  
  col.df <- rbind(tab_16S, tab_ITS) %>%
    select(Taxa, Label) %>%
    unique() 
  
  if (dist_name == "Dice") {
    dist.mat <- dist.binary(t(tab_comb), method = 5)
    dist.mat <- as.matrix(dist.mat)
  } else if (dist_name == "Bray-Curtis") {
    dist.mat <- vegdist(t(tab_comb))
    dist.mat <- as.matrix(dist.mat)
  } else {
    stop("dist_name should be either 'Dice' or 'Bray-Curtis'")
  }
  
  annc <- col.df[order(match(col.df$Taxa, colnames(dist.mat))),"Label"] %>% as.data.frame()
  rownames(annc) <- colnames(dist.mat)
  colnames(annc) <- "Label"
  
  ann_colors <- list(Label = c(bacteria=brewer.pal(n = 8, name = "Accent")[5], fungi=brewer.pal(n = 8, name = "Accent")[6]))
  
  g <- pheatmap(dist.mat, 
           color = colorRampPalette(brewer.pal(n = 7, name = "YlGnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), 
           show_rownames = T, show_colnames = F, 
           annotation = annc, annotation_colors = ann_colors, 
           cellheight = 8, cellwidth = 8,
           main = dist_name)
  
  return(g)
}

###=============
### a simple heatmap
###=============

heat <- function(s_Heat, props_Heat, grps, satu_limit = satu_limit, prop_cut = prop_cut, 
                 annotation_legend = T,
                 must_add_row = NULL,
                 gaps_col = NULL,
                 gaps_row = NULL,
                 show_colnames = T,
                 use_annotation = T,
                 annotation_colors = NA) {
  
  s_Heat <- s_Heat %>%
    dplyr::select(SampleID, grps) %>%
    unique() %>%
    arrange(SampleID) %>%
    arrange_(.dots = grps)
  
  anno <- s_Heat %>% select(-SampleID)
  rownames(anno) <- s_Heat$SampleID
  colnames(anno) <- grps
  
  props_Heat <- props_Heat[, s_Heat$SampleID]
  row_to_show <- rownames(props_Heat)[apply(props_Heat, 1, max) >= prop_cut]
  row_to_show <- c(row_to_show, must_add_row)
  
  props_Heat <- props_Heat[row_to_show, rownames(anno)]
  
  color = saturated_rainbow(101, saturation_limit = satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  if (use_annotation) {
    g <- pheatmap(props_Heat, annotation = anno, cluster_cols = F, cluster_rows = F, 
                  color = color, breaks = breaks, annotation_colors = annotation_colors,
                  gaps_col = gaps_col, gaps_row = gaps_row,
                  show_colnames = show_colnames, annotation_legend = annotation_legend, 
                  cellwidth = 8, cellheight = 8, fontsize_col = 8, fontsize_row = 8)
  } else {
    g <- pheatmap(props_Heat, cluster_cols = F, cluster_rows = F, 
                  color = color, breaks = breaks, annotation_colors = annotation_colors,
                  gaps_col = gaps_col, gaps_row = gaps_row,
                  show_colnames = show_colnames, annotation_legend = annotation_legend, 
                  cellwidth = 8, cellheight = 8, fontsize_col = 8, fontsize_row = 8)
  }
  
  return(g)
}

###=============
### lmer diff. abund. test
###=============

diff_abundance_taxa_lmer <- function(props_toTest, row_to_test, s_toTest, group_var, nominal_p_cut) {
  combine <- props_toTest %>%
    as.data.frame() %>%
    rownames_to_column(var = "row") %>%
    filter(row %in% row_to_test) %>%
    gather(key = SampleID, value = prop, -row) %>%
    mutate(prop = prop + 1e-6) %>%
    mutate(log_prop = log(prop)) %>%
    left_join(s_toTest, by = "SampleID")
  
  model <- combine %>%
    group_by(row) %>%
    do(fit = summary(lmer(log_prop ~ get(group_var) + (1|SubjectID), data = .)))
  
  summaries <- lapply(1:length(model$fit),
                      function(x) {
                        data.frame(row = model$row[[x]],
                                   model$fit[[x]]$coefficients,
                                   Term = rownames(model$fit[[x]]$coefficients),
                                   stringsAsFactors = F)
                      })
  
  summaries_df <- bind_rows(summaries) %>%
    mutate(Term = as.character(Term)) %>%
    filter(Term != "(Intercept)") %>%
    mutate(Term = gsub("get(group_var)", "Reference $\\rightarrow$ ", Term, fixed = T)) %>%
    rename(`$p$-value` = Pr...t..) %>%
    mutate(FDR = p.adjust(`$p$-value`, method = "BH")) %>%  
    mutate(Sig.Label = case_when(FDR < 0.0001 ~ "***",
                                 FDR < 0.001 ~ "**",
                                 FDR < 0.01 ~ "*", 
                                 FDR < 0.1 ~ ".",
                                 T ~ "")) %>%
    select(row, Term, Estimate, `$p$-value`, FDR, Sig.Label) %>%
    arrange(`$p$-value`) %>%
    filter(`$p$-value` < nominal_p_cut) %>%
    as.data.frame() 
  
  return_list <- list(summaries_df = summaries_df, 
                      combine = combine, 
                      group_var = group_var)
  return(return_list)
}



