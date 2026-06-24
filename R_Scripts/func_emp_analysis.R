
#Cleaning up raw data for downstream processing
process_data <- function(df) {
  # List of known reference sample names
  ref_samples <- c("Dmel_OregonR_non-inbred", "Dmel_Helsinki_inbred_F11", "Dmel_Samarkand_inbred")
  
  #Filter reference rows
  if ("sample" %in% names(df)) {
    df <- df %>% filter(!sample %in% ref_samples)
  } 
  
  #Cleaning up and separating columns
  if ("sample" %in% names(df)) {
    df <- df %>% separate_wider_delim(sample, "_", names = c("Sample", "Temperature", "generation", "replicate"))
  }
  df <- df%>%dplyr::select(!Temperature) #removing temp columns
  if ("replicate" %in% names(df)) {
    df <- df %>% mutate(replicate = gsub("", " ", replicate))
  }  
  if ("generation" %in% names(df)) {
    df <- df %>% mutate(generation = gsub("", " ", generation))
  }
  if ("generation" %in% names(df)) {
    df <- df %>% separate_wider_delim(generation, " ", names = c("F", "F2", "gen"), too_many = "merge")
  }
  if ("replicate" %in% names(df)) {
    df <- df %>% separate_wider_delim(replicate, " ", names = c("r", "r2", "rep"), too_many= "merge")
  }
  drop_cols2 <- intersect(c(6,7,9,10), seq_along(df))
  df <- df[-drop_cols2]
  if ("gen" %in% names(df)) {
    df <- df %>% mutate(gen = gsub(" ", "", gen))
  }
  if ("rep" %in% names(df)) {
    df <- df %>% mutate(rep = gsub(" ", "", rep))
  }
  
  #Renaming columns
  rename_cols <- c("Dmel_OregonR_non-inbred" = "ore_af",
                   "Dmel_Samarkand_inbred" = "sam_af",
                   "Dmel_Helsinki_inbred_F11" = "hel_af") 
  existing <- intersect(names(rename_cols), names(df))
  rename_list <- setNames(existing, rename_cols[existing])
  df <- df %>% dplyr::rename(!!!rename_list)
  
  # Filter specific generation 10
  if ("gen" %in% names(df)) {
    df$gen <- as.numeric(df$gen)
    df <- df[df$gen %in% c(10), ]
  }
  
  #Generating a mpos column as the median position of the window
  df<- df%>%dplyr::select(!c(window, pos))
  df$window<- paste(df$chr, df$win, sep = ":")
  df<- df%>% separate_wider_delim(win, "-", names=c("start", "end"))
  df$start<- as.numeric(df$start)
  df$end<- as.numeric(df$end)
  df$mpos<- (df$start+df$end)/2
  df <- df %>% dplyr::select(-all_of(c("start", "end")))
  return(df)
}

#Calculating logit-transformed AF
logit_freq<- function(samp_df){
  line_a<- colnames(samp_df)[ncol(samp_df)-1]
  line_b<- colnames(samp_df)[ncol(samp_df)]
  col_1<- paste0("logit_", line_a)
  col_2<- paste0("logit_", line_b)
  samp_df<- samp_df
  temp_df <- samp_df %>% dplyr::mutate(
    !!col_1 := log(samp_df[[line_a]] / samp_df[[line_b]]),
    !!col_2 := -log(samp_df[[line_a]] / samp_df[[line_b]])
  )
  temp_df<- temp_df%>%pivot_longer(cols =c(all_of(col_1), all_of(col_2)), 
                                   names_to = "line", values_to = "s")%>%
    separate_wider_delim(line, delim ="_", names = c("logit", "line", "drop"), )%>%
    dplyr::select(c(chr, window, mpos, Sample, rep, s, line, gen))%>%dplyr::rename("group"=Sample)
  return(temp_df)
}

#Hypothesis testing and computing FDR to classifying a region as significantly deviating (1) or not (0)
compute_p<- function(data_merge,  data_sample, ref){
  data_merge<- data_merge%>% dplyr::select(c("window", "group", "rep", "s", "mpos","chr"))
  data_merge$window<- as.factor(data_merge$window)
  data_merge$rep<- as.factor(data_merge$rep)
  data_merge<-data_merge%>%mutate(group = fct_relevel(group, ref)) 
  
  #Estimating as the focal sample's selection response as mean logit-AF
  data_sample<- data_sample%>% dplyr::select(c("window", "group", "rep", "s", "mpos","chr"))
  data_sample<- data_sample %>% 
    group_by(window) %>%
    summarise_at(vars(s),list(val=mean) )
  
  #Hypothesis testing for each window 
  model_offset<- data_merge%>% group_by(window, mpos,chr)%>%
    do({
      model <- lm(s ~ group, data = .)
      s<- data_sample$val[data_sample$window == unique(.$window)]
      window_val<- data_sample$window[data_sample$window == unique(.$window)]
      em <- emmeans(model, pairwise ~ group)
      a<-test(em,null=s)
      contrast_df <- as.data.frame(em$contrasts)
      p<- as.numeric(a[["contrasts"]][7])
      contrast_value <- contrast_df$estimate[1]
      tibble(p_value = p, contrast = contrast_value, s_val=s)
    })
  
  #Multiple testing coreection using Benjamini-Hochberg
  model_offset <- model_offset %>%
    ungroup() %>%
    mutate(adj_P = p.adjust(p_value, method = "BH"))
  merge_df<- model_offset
  merge_df$group<- NA
  
  #Annotate with class: 1 significant, 0 non significant
  merge_df <- merge_df %>%
    mutate(
      group = case_when(
        adj_P <= 0.05 ~ 1,
        adj_P >=0.05 ~ 0,
        TRUE ~ 2
      )
    )
  return(merge_df)
}

#Combining significance information for all the three tests
find_pattern<- function( merge_df){
  
  merge_df<- merge_df %>%
    mutate(
      pattern = paste(group_HS, group_HO, group_SO, sep = "-")
    )
  
  merge_df<- merge_df%>%mutate(group=case_when(
    pattern=="1-1-1" ~ 3,
    pattern=="0-0-0" ~ 0,
    pattern %in% c("1-0-0","0-1-0", "0-0-1") ~ 1,
    pattern %in% c("1-1-0","0-1-1", "1-0-1") ~ 2,
    TRUE ~ 4
  ))
}
