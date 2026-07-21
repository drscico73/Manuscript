theme_common <- theme_bw() +
  theme(
    axis.text.x  = element_text(vjust = 0.5, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    plot.title   = element_text(size = 20, face = "bold"),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14)
  )

theme_set(theme_common)

### Plotting AF ###

plot_af_gen<- function(df_wide, line){
  
  col_1 <- switch(line,
                  "OregonR"   = "ore_af",
                  "Samarkand" = "sam_af",
                  "Helsinki"  = "hel_af",
                  stop("Wrong input for 'line'. Use: OregonR, Samarkand, or Helsinki"))
  sample<- unique(df_wide$Sample)
  gen_samp<- unique(df_wide$gen)
  df<- df_wide%>%pivot_longer(cols=all_of(col_1), 
                              names_to = "line", values_to = "AF")
  custom_colors <- c("HO-OH" = "purple3", "SO-OS" = "green4", "HS-SH" = "#619CFF")
  df$color_group <- with(df, ifelse(sample == "HO-OH", "HO-OH",
                                    ifelse(sample == "HS-SH", "HS-SH", "SO-OS")))
  df$chr <- factor(df$chr, levels = c("X", sort(unique(df$chr[df$chr != "X"]))))
  df$gen<- as.factor(df$gen)
  p <- ggplot(df, aes(x = mpos, y = AF, color = Sample, fill = Sample)) +
    stat_summary(aes(fill = sample), fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "ribbon", alpha = 0.2, color = NA, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    labs(x = "Position (Mb)",
      y = paste0(line, " AF"),
      color = "Sample") +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors)+
    facet_grid( ~ chr, scales = "free_x") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18, face="bold"),
      legend.key.size = unit(1.5, "cm"))+
    ylim(min = 0.0, max = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +
    scale_x_continuous(labels = scales::label_number(scale = 1/1e6)) 
  return(p)
}

### Plotting the predicted architecture along the genome ###

plot_class<- function(df1, df2, gen_beg, gen_end){
  custom_colors <- c("add" = "navyblue", "dev" = "magenta")
  group_colours<- c("1/3"="#FFBFAA", "2/3"="#FF7B5A", "0/3"="#f9f1f1", "3/3"="#FF0000")
  df1<- df1%>%dplyr::select(chr, mpos, dev, sample)%>% 
    mutate(add= 0)%>%
    pivot_longer(cols = c("add", "dev"), names_to = "group", values_to = "Deviation")%>%
    mutate(color_sample =ifelse(group == "add", "add","dev"))
  df1$chr <- factor(df1$chr, levels = c("X", sort(unique(df1$chr[df1$chr != "X"]))))
  
  d<- df2%>%dplyr::select("window", "group")%>%
    separate_wider_delim("window", ":", names=c("chr", "pos"))%>%
    separate_wider_delim("pos", "-", names=c("low", "high"))%>%
    mutate(color_group = case_when(
      group == 0 ~ "0/3",
      group == 1 ~ "1/3",
      group == 2 ~ "2/3",
      TRUE ~ "3/3"
    ),
    low = as.numeric(low),
    high = as.numeric(high)
    )
  d$chr <- factor(d$chr, levels = c("X", sort(unique(d$chr[d$chr != "X"]))))
  
  plot_dev<- ggplot() +
    geom_rect(
      data = d,
      inherit.aes = FALSE,
      aes(xmin = low, xmax = high, ymin = -Inf, ymax = Inf, fill = color_group),
      alpha = 0.5
    ) +
    geom_line(data = df1[df1$sample=="dev_HS",], aes(x = mpos, y = Deviation, color = color_sample)) +
    geom_line(data = df1[df1$sample=="dev_SO",], aes(x = mpos, y = Deviation, color = color_sample)) +
    geom_line(data = df1[df1$sample=="dev_HO",], aes(x = mpos, y = Deviation, color = color_sample)) +
    stat_summary(data = df1, aes(x = mpos, y = Deviation, color = color_sample),
                 geom = "line"  ,fun = mean, linewidth=1.5)+
    ylim(-1.75, 1.75)+
    facet_wrap( ~chr, scales = "free_x", ncol = 5) +
    labs(#title = paste0("Additive vs Non-additive regions ", gen_beg, "-", gen_end),
      x = "Position (Mb)", y = "Predicted - Calculated (logit AF)", color = "", 
      fill="Significance" 
    ) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 18, face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(1.3, "cm")
    ) +
    scale_x_continuous(labels = scales::label_number(scale = 1 / 1e6)) +
    scale_color_manual(
      values = custom_colors,
      labels = c("add" = "Expected Difference", "dev" = "Calculated Difference")
    )+
    scale_fill_manual(values=group_colours)
  return(plot_dev)
}
