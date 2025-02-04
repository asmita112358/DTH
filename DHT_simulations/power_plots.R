##Power plots for parametric sims

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(data.table)
library(ggpattern)
setwd("/Users/asmitaroy/Downloads/Dispersion/pvalue_combinations_v2")
p1 <- read.csv("NB_pow_2_bray_bal.csv")[,-c(1, 7, 8)]
p2 <- read.csv("NB_pow_3_bray_bal.csv")[,-c(1, 7, 8)]
p3 <- read.csv("NB_pow_5_bray_bal.csv")[,-c(1, 7, 8)]
p4 <- read.csv("NB_pow_2_bray_unbal.csv")[,-c(1, 7, 8)]
p5 <- read.csv("NB_pow_3_bray_unbal.csv")[,-c(1, 7, 8)]
p6 <- read.csv("NB_pow_5_bray_unbal.csv")[,-c(1, 7, 8)]



##All power plots together
design_type = rep(c("balanced", "unbalanced"), each = nrow(p1)*3)
G = rep(rep(c(2,3,5), each = nrow(p1)), times = 2)
G = paste("G =",G)
pow_data = data.frame(design_type, G, rbind(p1, p2, p3, p4, p5, p6))
data_long <- pow_data %>%
  pivot_longer(
    cols = c(GO, GO.perm, PERMANOVA, BETADISPER, DHT),
    names_to = "Method",
    values_to = "Power"
  )%>%
  mutate(
    Method = factor(
      Method,
      levels = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT") # Desired order
    )
  )

# Create the plot
all_pow <- ggplot(data_long, aes(x = theta, y = Power, color = Method, shape = Method)) +
  scale_shape_manual(values = 0:6) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  facet_grid(rows = vars(design_type), cols = vars(G)) +
  labs(
    title = "Power plots for Negative Binomial data",
    x = "Theta",
    y = "Power",
    color = "Method"
  ) +
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),   # Increase legend text size
        legend.title = element_text(size = 14, face = "bold"),  # Increase legend title size
        axis.title = element_text(size = 14),   # Increase axis title size
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))
all_pow
ggsave("pow_NB_all.png", all_pow, device = "png",
       width = 8, height = 5, units = "in", bg = "white",
       dpi = 1200)



###Power plots for Normal data
p1 <- read.csv("norm_pow_2_euc_bal.csv")[,-c(1, 7, 8)]
p2 <- read.csv("norm_pow_3_euc_bal.csv")[,-c(1, 7, 8)]
p3 <- read.csv("norm_pow_5_euc_bal.csv")[,-c(1, 7, 8)]
p4 <- read.csv("norm_pow_2_euc_unbal.csv")[,-c(1, 7, 8)]
p5 <- read.csv("norm_pow_3_euc_unbal.csv")[,-c(1, 7, 8)]
p6 <- read.csv("norm_pow_5_euc_unbal.csv")[,-c(1, 7, 8)]



##All power plots together
design_type = rep(c("balanced", "unbalanced"), each = nrow(p1)*3)
G = rep(rep(c(2,3,5), each = nrow(p1)), times = 2)
G = paste("G =",G)
pow_data = data.frame(design_type, G, rbind(p1, p2, p3, p4, p5, p6))
data_long <- pow_data %>%
  pivot_longer(
    cols = c(GO, GO.perm, PERMANOVA, BETADISPER, DHT),
    names_to = "Method",
    values_to = "Power"
  )%>%
  mutate(
    Method = factor(
      Method,
      levels = c("GO", "GO.perm", "PERMANOVA", "BETADISPER", "DHT") # Desired order
    )
  )

# Create the plot
all_pow <- ggplot(data_long, aes(x = theta, y = Power, color = Method, shape = Method)) +
  scale_shape_manual(values = 0:6) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  facet_grid(rows = vars(design_type), cols = vars(G)) +
  labs(
    title = "Power plots for Normal data",
    x = "Theta",
    y = "Power",
    color = "Method"
  ) +
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),   # Increase legend text size
        legend.title = element_text(size = 14, face = "bold"),  # Increase legend title size
        axis.title = element_text(size = 14),   # Increase axis title size
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))
all_pow
ggsave("pow_norm_all.png", all_pow, device = "png",
       width = 8, height = 5, units = "in", bg = "white",
       dpi = 1200)

