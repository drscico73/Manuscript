### Loading the libraries ###
library("tidyverse")
library("ggplot2")
library("dplyr")
library("tidyr")
library("stringr")
library("patchwork")
library("scales")
library("poolSeq")
library("here")
library("emmeans")
library("ggh4x")

### Load functions ###
source("Desktop/manuscript_files/functions_empirical_analysis.R")
source("Desktop/manuscript_files/functions_visualisation.R")

### Read RDS files having window AF estimates ###
rds_hel_ho<- readRDS("~/Desktop/manuscript_files/rds_files/hel_ho_out_500.rds")
rds_hel_hs<- readRDS("~/Desktop/manuscript_files/rds_files/hel_hs_out_500.rds")
rds_ore_ho<- readRDS("~/Desktop/manuscript_files/rds_files/ore_ho_out_500.rds")
rds_ore_so<- readRDS("~/Desktop/manuscript_files/rds_files/ore_so_out_500.rds")
rds_sam_so<- readRDS("~/Desktop/manuscript_files/rds_files/sam_so_out_500.rds")
rds_sam_hs<- readRDS("~/Desktop/manuscript_files/rds_files/sam_hs_out_500.rds")

### Data Cleaning and AF resestimation ###
#Cleaning raw data files
data_hel_ho<- process_data(rds_hel_ho)
data_hel_hs<- process_data(rds_hel_hs)
data_ore_ho<- process_data(rds_ore_ho)
data_ore_so<- process_data(rds_ore_so)
data_sam_so<- process_data(rds_sam_so)
data_sam_hs<- process_data(rds_sam_hs)

#Merging cleaned files
#Note: AF for each sample was estimated twice; once from markers of line A and once from markers of line B.
data_ho<- merge(data_hel_ho, data_ore_ho,  by=c("chr","mpos","Sample",
                                                "gen","rep","window"), suffixes = c("_hel", "_ore"))
data_so<- merge(data_sam_so,data_ore_so,by=c("chr","mpos","Sample",
                                             "gen","rep","window"), suffixes = c("_sam", "_ore"))
data_hs<- merge(data_sam_hs,data_hel_hs,by=c("chr","mpos","Sample",
                                             "gen","rep","window"), suffixes = c("_sam", "_hel"))

#Restimating AF as the average of the window AF estimated from 2 different marker sets
data_ho$ore_af<- (data_ho$ore_af_hel+data_ho$ore_af_ore)/2
data_ho$hel_af<- (data_ho$hel_af_hel+data_ho$hel_af_ore)/2
data_so$ore_af<- (data_so$ore_af_sam+data_so$ore_af_ore)/2
data_so$sam_af<- (data_so$sam_af_sam+data_so$sam_af_ore)/2
data_hs$sam_af<- (data_hs$sam_af_hel+data_hs$sam_af_sam)/2
data_hs$hel_af<- (data_hs$hel_af_hel+data_hs$hel_af_sam)/2

### Selection Response ###
#Estimating logit transformed AF as response to selection
sel_ho<- logit_freq(data_ho)%>%group_by(chr,window, mpos, group, rep, gen)%>%pivot_wider(names_from = line, values_from = s)%>%ungroup()%>%
  dplyr::rename("ore_s_ho"=ore, "hel_s_ho"=hel)%>%filter(gen=="10")%>%dplyr::select(!gen)
sel_so<- logit_freq(data_so)%>%group_by(chr, mpos, group, rep, gen)%>%pivot_wider(names_from = line, values_from = s)%>%ungroup()%>%
  dplyr::rename("ore_s_so"=ore, "sam_s_so"=sam)%>%filter(gen=="10")%>%dplyr::select(!gen)
sel_hs<- logit_freq(data_hs)%>%group_by(chr, mpos, group, rep, gen)%>%pivot_wider(names_from = line, values_from = s)%>%ungroup()%>%
  dplyr::rename("sam_s_hs"=sam, "hel_s_hs"=hel)%>%filter(gen=="10")%>%dplyr::select(!gen)

### Testing for deviations from additive assumption ###
#Estimating deviation as the difference between the predicted logit AF and logit AF calculated from the sample

#S(SO-OS)-H(HO-OH)=S(HS-SH)
sel_sam_so<- sel_so%>%dplyr::select(-ore_s_so)%>%dplyr::rename("s"=sam_s_so)
sel_hel_ho<- sel_ho%>%dplyr::select(-ore_s_ho)%>%dplyr::rename("s"=hel_s_ho)
model_df_hs_10<- rbind(sel_sam_so, sel_hel_ho)
sel_sam_hs<- sel_hs%>%dplyr::select(-hel_s_hs)%>%dplyr::rename("s"=sam_s_hs)%>%
  mutate(group="HS-SH")
merge_hs_10<- compute_p(model_df_hs_10, sel_sam_hs, "SO-OS")
merge_hs_10$dev<- merge_hs_10$s_val-merge_hs_10$contrast

#H(HS-SH)-O(SO-OS)=H(HO-OH)
sel_ore_so<- sel_so%>%dplyr::select(-sam_s_so)%>%dplyr::rename("s"=ore_s_so)
sel_hel_hs<- sel_hs%>%dplyr::select(-sam_s_hs)%>%dplyr::rename("s"=hel_s_hs)
model_df_ho_10<- rbind(sel_ore_so, sel_hel_hs)
sel_hel_ho<- sel_ho%>%dplyr::select(-ore_s_ho)%>%dplyr::rename("s"=hel_s_ho)
merge_ho_10<- compute_p(model_df_ho_10, sel_hel_ho, "HS-SH")
merge_ho_10$dev<- merge_ho_10$s_val-merge_ho_10$contrast

#O(HO-OH)-S(HS-SH)=O(SO-OS)
sel_ore_so<- sel_ho%>%dplyr::select(-hel_s_ho)%>%dplyr::rename("s"=ore_s_ho)
sel_sam_hs<- sel_hs%>%dplyr::select(-hel_s_hs)%>%dplyr::rename("s"=sam_s_hs)
model_df_so_10<- rbind(sel_ore_so, sel_sam_hs)
sel_ore_so<- sel_so%>%dplyr::select(-sam_s_so)%>%dplyr::rename("s"=ore_s_so)
merge_so_10<- compute_p(model_df_so_10, sel_ore_so, "HO-OH")
merge_so_10$dev<- merge_so_10$s_val-merge_so_10$contrast

#Merging the results from all the three comparisons
merge_dev_10<- merge(merge_hs_10%>%dplyr::select(window, chr, mpos, dev), merge_ho_10%>%dplyr::select(window, chr, mpos, dev),
                     by=c("window", "mpos", "chr"), suffixes = c("_HS", "_HO"))
merge_dev_10<- merge(merge_dev_10, merge_so_10%>%dplyr::select(window, chr, mpos, dev), by=c("window", "mpos", "chr"))%>%
  dplyr::rename("dev_SO"=dev)
merge_dev_10<- merge_dev_10%>% pivot_longer(cols=c(dev_HS, dev_SO, dev_HO), names_to="sample", values_to="dev")
merge_comp_10<- merge(merge_ho_10%>%dplyr::select(window, group),
                      merge_hs_10%>%dplyr::select(window, group), by="window",
                      suffixes = c("_HO", "_HS"))

merge_comp_10<- merge(merge_comp_10,  merge_so_10%>%dplyr::select(window, group),
                      by="window")%>%dplyr::rename("group_SO"=group)
merge_df_10<- find_pattern(merge_comp_10)
sup_4<-plot_class(merge_dev_10, merge_df_10, "1", "10")
print(sup_4)
