##plot codes for IBD simulations and continuous covariates simulation
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpattern)
library(ggpubr)


pval = read.csv("cont_size.csv")[,c(-1, -6,-7)]
len = nrow(pval)
size = apply(pval <= 0.05, 2, mean)
se_size = sqrt(size*(1-size)/len)
x_loc <- -log10(0.05)
y_loc <- -log10(0.05)
theo_q = (1:len)/(len + 1)
theo_q = -log10(theo_q)
p.sort = -log10(apply(pval, 2, sort))
data = data.frame(theo_q, p.sort)
colnames(data) = c("theo_q", "GO","GO.perm","PERMANOVA", "Betadisper", "DHT")
data.plot <- data.table(data)
qqp <- data.plot %>%
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
  labs(title = "QQplot",
       x = "Expected P-value (-log10 scale)", y = "Observed P-value (-log10 scale)") +
  scale_x_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), breaks = 0:3, expand = c(0,0)) +
  theme_minimal()+
  theme(legend.text = element_text(size = 12),   # Adjust legend text size
        legend.title = element_text(size = 14, face = "bold"),  # Adjust legend title size
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(face = "bold",hjust = 0.5),
        plot.title.position = "panel")


p1 = read.csv("cont_pow.csv")[,c(-1, -7, -8)]

colnames(p1)[1] <- "theta"
data_long <- reshape2::melt(p1, id.vars = "theta", variable.name = "Method", value.name = "Power")
pow <- ggplot(data_long, aes(x = theta, y = Power, color = Method, shape = Method)) +
  scale_shape_manual(values = 0:6) +
  #scale_x_reverse() + 
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  labs(title = "Empirical Power",
       x = "theta",
       y = "Power") +
  theme_minimal()+
  theme(legend.text = element_text(size = 12),   # Adjust legend text size
        legend.title = element_text(size = 14, face = "bold"),  # Adjust legend title size
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(face = "bold",hjust = 0.5),
        plot.title.position = "panel")

pow
allp <- ggarrange(qqp,pow,ncol = 2, legend = "right", common.legend = T)
allp
ggsave("cont_all.png", allp, device = "png",
       width = 8, height = 4, units = "in", bg = "white",
       dpi = 1200)
