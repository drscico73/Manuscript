theme_common<- theme_bw()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 17, face = "bold"),
        axis.title.y = element_text(size = 17, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.position = "right",
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

### Plotting predicted vs true architecture ###

plot_sig<- function(df, arch_type){
  df$sign<- as.factor(df$sign)
  result_color<- c("0"="green2", "1"="red2")
  p_dist_plot<- ggplot(df, aes(x = diff_pos,
                               y = -log10(p_value), colour=sign)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(
      x = "Position 1 - Position 2 (Mb)",
      y = "-log10(P-Value)",
      colour = "Predicted Architecture",
      title = paste0("Model Result vs Distance between targets:", arch_type, " architecture")
    ) +
    theme_common+
    scale_color_manual(values=result_color, labels=c("0"="Additive", "1"="Non-additive"))+
    scale_x_continuous(labels = label_number(scale = 1/1e6)) 
  
  return(p_dist_plot)
}

### Plotting the difference between the predicted and calculated selection response as a function of distance between the two selected targets ###

plot_diff_dist<- function(df, suffix){
  arch_color<- c("add"="yellow", "non-add"="pink")
  s_dist_plot<- ggplot(df, aes(x = diff_pos,
                               y = diff_s, color=arch)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(
      x = "Position 1 - Position 2 (Mb)",
      y = paste0("Predicted - ", suffix),
      title = paste("Predicted - ", suffix, 
                    " logit frequency vs Distance between targets"),
      color="Architecture Type"
    ) +
    scale_color_manual(values=arch_color, labels=c("add"="Additive", "non-add"="Non-Additive"))+
    theme_common+ 
    scale_x_continuous(labels = label_number(scale = 1/1e6)) 
  
  return(s_dist_plot)
}

### Plotting roc curve ###

plot_roc<- function(df, thresh){
  plot_color = '#233f7d'  
  text_color = 'black' 
  highlight_df <- df %>%
    filter(round(threshold, 2) == thresh)
  
  plot_roc<- ggplot(data = df) +
    geom_line(aes(x = fpp, y = tpp),
              size = 1, color = plot_color) + 
    geom_ribbon(aes(x=fpp,ymin=0,ymax=tpp),alpha=0.2, fill=plot_color) +
    annotate("text", alpha = 0.6, x = 0.7, y = 0.1,
             label = paste("AUC:", auc),
             size = 8, color = text_color) + 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1.01), expand = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "False Positive Rate",
         y = "True Positive Rate",
         title = "ROC curve of the model") +
    theme_common
  print(plot_roc)
}
