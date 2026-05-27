### Loading the libraries ###
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("stats")
library("stringr")
library("patchwork")
library("scales")
library("ggpubr")
library("cowplot")
library("foreach")
library("poolSeq")
library("boot")

### Allele Frequency Plots ###
# data loading and cleaning
df <- readRDS("/Path/to/file/AlleleFreq.rds")
df2 <-df %>%
  separate_wider_delim(sample, "_", names = c("Cross", "O", "Temperature", "generation", "replicate")) %>%
  mutate(across(c(replicate, generation), ~ gsub(" ", "", .))) %>%
  mutate(replicate = gsub("", " ", replicate),
         generation = gsub("", " ", generation)) %>%
  separate_wider_delim(replicate, " ", names = c("r", "r2", "rep", "x"), too_many = "merge") %>%
  select(-c(2,5,6,8,9,11)) %>%
  mutate(generation = gsub(" ", "", generation))
df2$chr <- factor(df2$chr, levels = c("X", "2L", "2R", "3L", "3R"))

# Logit transform frequencies
logit_freq <- df2 %>%
  group_by(chr, pos, Cross, generation, rep)%>%
  summarise_at(vars(`Dmel_OregonR_non-inbred`),
               list(logit_Freq = logit))

data_F1 <- filter(logit_freq, generation == "F1")
data_F10 <- filter(logit_freq, generation == "F10")

# Calculate the selection strengths
FreqDiff <- data.frame(
  chr = data_F1$chr,
  pos = data_F1$pos,
  Cross = data_F1$Cross,
  rep = data_F1$rep,
  generation <- "F10"
  FreqDiff10 = data_F10$logit_Freq - data_F1$logit_Freq,
  s = (data_F10$logit_Freq - data_F1$logit_Freq)/10
)

# Performing t-test to find windows differing in Selection response
t_test <- compare_means(s ~ Cross, data = FreqDiff, method = "t.test",
                            group.by = c("chr", "pos"), p.adjust.method = "fdr") %>% mutate(generation = "F10")

t_test_results <- bind_rows(t_test) %>%
  mutate(
    significance = ifelse(p <= 0.05, "Significant", "Not significant"),
    log_pval = -log10(p),
    chr = factor(chr, levels = c("X", "2L", "2R", "3L", "3R"))
  )

# Allele Frequency plot
p1 <- ggplot(df2) + 
  stat_summary(aes(x = pos, y = `Dmel_OregonR_non-inbred`, color = Cross, group = Cross),
               geom = "line", fun = "median",  linewidth = 1, alpha = 1) + 
  stat_summary(aes(x = pos, y = `Dmel_OregonR_non-inbred`, fill = Cross, group = Cross),
               geom = "ribbon", fun.data = "mean_sd",  linewidth = 0.5, alpha = 0.3, show.legend = FALSE) +
  facet_grid(generation~chr, scales = "free_x") +
  scale_color_manual(values = c("darkgreen", "orange1")) +
  scale_fill_manual(values = c("darkgreen", "orange1")) +
  labs(x = "Position (Mb)", y = "Allele Frequency", color = "Population") + 
  theme_bw()+
  theme(
      axis.text.x = element_text(vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1, "cm"),
      legend.position = "bottom",
      strip.text.x = element_text(size = 10)
      ) +
  scale_x_continuous(labels = label_number(scale = 1/1e6), sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Generation", breaks = NULL, labels = NULL), limits = c(0,1)) 

# Selection strength plot
p2 <- ggplot(FreqDiff) +
  stat_summary(aes(pos, s, color = Cross), fun = "mean", geom = "line", linewidth = 1) +
  stat_summary(aes(pos, s, fill = Cross), fun.data = "mean_sd", geom = "ribbon", linewidth = 0.5,
               alpha = 0.3, show.legend = FALSE) +
  labs(x = "Position (Mb)", y = "Selection strength", color = "Population") +
  scale_color_manual(values = c("darkgreen", "orange1")) +
  scale_fill_manual(values = c("darkgreen", "orange1")) +
  facet_grid(generation~chr, scales = "free_x") +
  theme_bw()+
  theme(
      axis.text.x = element_text(vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1, "cm"),
      legend.position = "bottom",
      strip.text.x = element_text(size = 10)
      ) +
  scale_x_continuous(labels = label_number(scale = 1/1e6),
                     sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Generation", breaks = NULL, labels = NULL)) 

# t-test results
p3 <- ggplot(subset(t_test_results, .y.=="s"), aes(x = pos, y = log_pval, color = significance)) +
    geom_point(alpha = 0.7, size=2.5) +
    facet_grid(generation~chr, scales = "free_x") +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "cyan4")) +
    labs(x = "Position (Mb)", y = "-log10(p)", color = "Significance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.7) +
    theme_bw()+
    theme(
      axis.text.x = element_text(vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(1, "cm"),
      legend.position = "bottom",
      strip.text.x = element_text(size = 8)
      ) +
    scale_x_continuous(labels = label_number(scale = 1/1e6),
                       sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = "Generation", breaks = NULL, labels = NULL)) +
    guides(color = guide_legend(override.aes = list(alpha = 0.5)))
