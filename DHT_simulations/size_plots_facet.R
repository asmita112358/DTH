
library(ggplot2)
library(data.table)
library(ggpubr)
library(Haplin)
setwd("/Users/asmitaroy/Downloads/Dispersion/pvalue_combinations_v2")
##All Bar plots for size, faceted
p = list()
p[[1]] <- read.csv("norm_size_2_euc_bal.csv")[,-c(1, 6,7)]
p[[2]] <- read.csv("norm_size_3_euc_bal.csv")[,-c(1, 6,7)]
p[[3]] <- read.csv("norm_size_5_euc_bal.csv")[,-c(1, 6,7)]
p[[4]] <- read.csv("norm_size_2_euc_unbal.csv")[,-c(1, 6,7)]
p[[5]] <- read.csv("norm_size_3_euc_unbal.csv")[,-c(1, 6,7)] 
p[[6]] <- read.csv("norm_size_5_euc_unbal.csv")[,-c(1, 6,7)] 

method_names = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT")
s = matrix(nrow = 6, ncol = 5)
for(i in 1:6)
{
  colnames(p[[i]]) <- method_names
  rownames(p[[i]]) <- NULL
  s[i,] = colMeans(p[[i]] <= 0.05)
}
colnames(s) <- method_names
design_type <- rep(c("balanced", "unbalanced"), each = 3)
group = rep(c(2,3,5), times = 2)
group = paste("G =", group)
size = data.frame(group, design_type, s)
data_long <- size %>%
  pivot_longer(
    cols = c(GO, GO.perm, PERMANOVA, BETADISPER, DHT),
    names_to = "Method",
    values_to = "Size"
  )%>%
  mutate(
    Method = factor(
      Method,
      levels = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT") # Desired order
    )
  )
data_long$se_size = sqrt(data_long$Size*(1-data_long$Size)/500)

all_size = ggplot(data_long, aes(x = Method, y = Size, fill = Method)) + 
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = Size - 1.96 * se_size, ymax = Size + 1.96 * se_size), 
                width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(design_type), cols = vars(group)) +
  geom_hline(yintercept = 0.05) +
  labs(title = "Empirical size for Normal data", x = NULL, y = "Empirical size") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),   # Increase axis title size
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("barp_norm_all.png", all_size, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)



##All QQplots, faceted
qqp = list()
title = rep(c("G = 2", "G = 3", "G = 5"), times = 2)
title = paste0(title,", ", design_type)
for(case in seq_along(p)){
  x_loc <- -log10(0.05)
  y_loc <- -log10(0.05)
  mc = p[[case]]
  N_rep = nrow(mc)
  log.pvals = -log10(apply(mc, 2, sort))
  #cat(max(log.pvals))
  log.theo = -log10((1:N_rep)/(N_rep+1))
  data = data.frame(log.theo,log.pvals)
  colnames(data) <- c("theo_q", method_names)
  data.plot <- data.table(data)
  #.lim <- c(0, max(log.pvals)) * 1.05
  qqp[[case]] <- data.plot %>%
    melt(id.vars = 1,
         variable.name = "Method") %>%
    ggplot(aes(x = theo_q, y = value,
               group = Method, color = Method, shape = Method)) +
    scale_shape_manual(values = 0:5) +
    geom_point() +
    geom_abline() +
    geom_segment(aes(x = x_loc, y = 0, xend = x_loc, yend = x_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, y = y_loc, xend = y_loc, yend = y_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    labs(title = title[case],
         x = "", y = "") +
    scale_x_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
    scale_y_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
    theme_minimal()+
    theme(legend.text = element_text(size = 12),   # Adjust legend text size
          legend.title = element_text(size = 14, face = "bold"),  # Adjust legend title size
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(face = "bold",hjust = 0.5),
          plot.title.position = "panel")
  
  qqp[[case]]
}
qqp_all = ggarrange(qqp[[1]], qqp[[2]], qqp[[3]],qqp[[4]], qqp[[5]], qqp[[6]], nrow = 2,ncol = 3, align = "h", legend = "bottom", common.legend = T)
final_qqp <- annotate_figure(qqp_all, 
                             top = text_grob("QQ Plots for Normal outcome", face = "bold", size = 14),
                             bottom = text_grob("Expected P-value (-log10 scale)", size = 14),left = text_grob("Observed P-value (-log10 scale)",size = 14,rot = 90 ))
final_qqp


ggsave("qqp_norm_all.png", final_qqp, device = "png",
       width = 7, height = 5, units = "in", bg = "white",
       dpi = 1200)

##############################
##############################
##############################
##############################
## NB data ##################
##############################
##############################
##############################
##############################

##All Bar plots for size, faceted
p = list()
p[[1]] <- read.csv("NB_size_2_bray_bal.csv")[,-c(1, 6,7)]
p[[2]] <- read.csv("NB_size_3_bray_bal.csv")[,-c(1, 6,7)]
p[[3]] <- read.csv("NB_size_5_bray_bal.csv")[,-c(1, 6,7)]
p[[4]] <- read.csv("NB_size_2_bray_unbal.csv")[,-c(1, 6,7)]
p[[5]] <- read.csv("NB_size_3_bray_unbal.csv")[,-c(1, 6,7)] 
p[[6]] <- read.csv("NB_size_5_bray_unbal.csv")[,-c(1, 6,7)] 

method_names = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT")
s = matrix(nrow = 6, ncol = 5)
for(i in 1:6)
{
  colnames(p[[i]]) <- method_names
  rownames(p[[i]]) <- NULL
  s[i,] = colMeans(p[[i]] <= 0.05)
}
colnames(s) <- method_names
design_type <- rep(c("balanced", "unbalanced"), each = 3)
group = rep(c(2,3,5), times = 2)
group = paste("G =", group)
size = data.frame(group, design_type, s)
data_long <- size %>%
  pivot_longer(
    cols = c(GO, GO.perm, PERMANOVA, BETADISPER, DHT),
    names_to = "Method",
    values_to = "Size"
  )%>%
  mutate(
    Method = factor(
      Method,
      levels = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT") # Desired order
    )
  )
data_long$se_size = sqrt(data_long$Size*(1-data_long$Size)/500)

all_size = ggplot(data_long, aes(x = Method, y = Size, fill = Method)) + 
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = Size - 1.96 * se_size, ymax = Size + 1.96 * se_size), 
                width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(design_type), cols = vars(group)) +
  geom_hline(yintercept = 0.05) +
  labs(title = "Empirical size for Negative Binomial data", x = NULL, y = "Empirical size") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),   # Increase axis title size
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("barp_NB_all.png", all_size, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)



##All QQplots, faceted
qqp = list()
title = rep(c("G = 2", "G = 3", "G = 5"), times = 2)
title = paste0(title,", ", design_type)
for(case in seq_along(p)){
  x_loc <- -log10(0.05)
  y_loc <- -log10(0.05)
  mc = p[[case]]
  N_rep = nrow(mc)
  log.pvals = -log10(apply(mc, 2, sort))
  #cat(max(log.pvals))
  log.theo = -log10((1:N_rep)/(N_rep+1))
  data = data.frame(log.theo,log.pvals)
  colnames(data) <- c("theo_q", method_names)
  data.plot <- data.table(data)
  #.lim <- c(0, max(log.pvals)) * 1.05
  qqp[[case]] <- data.plot %>%
    melt(id.vars = 1,
         variable.name = "Method") %>%
    ggplot(aes(x = theo_q, y = value,
               group = Method, color = Method, shape = Method)) +
    scale_shape_manual(values = 0:5) +
    geom_point() +
    geom_abline() +
    geom_segment(aes(x = x_loc, y = 0, xend = x_loc, yend = x_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, y = y_loc, xend = y_loc, yend = y_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    labs(title = title[case],
         x = "", y = "") +
    scale_x_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
    scale_y_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
    theme_minimal()+
    theme(legend.text = element_text(size = 12),   # Adjust legend text size
          legend.title = element_text(size = 14, face = "bold"),  # Adjust legend title size
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(face = "bold",hjust = 0.5),
          plot.title.position = "panel")
  
  qqp[[case]]
}
qqp_all = ggarrange(qqp[[1]], qqp[[2]], qqp[[3]],qqp[[4]], qqp[[5]], qqp[[6]], nrow = 2,ncol = 3, align = "h", legend = "bottom", common.legend = T)
final_qqp <- annotate_figure(qqp_all, 
                             top = text_grob("QQ Plots for Negative Binomial outcome", face = "bold", size = 14),
                             bottom = text_grob("Expected P-value (-log10 scale)", size = 14),left = text_grob("Observed P-value (-log10 scale)",size = 14,rot = 90 ))
final_qqp


ggsave("qqp_NB_all.png", final_qqp, device = "png",
       width = 7, height = 5, units = "in", bg = "white",
       dpi = 1200)
