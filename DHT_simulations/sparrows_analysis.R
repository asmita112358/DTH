setwd("~/Downloads/Dispersion/GO-data&code")
source("~/Downloads/Dispersion/Dispersion_codes/all_funs.R")
library(vegan)
#Sparrows and coral data analysis for all competing methods

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # The Bumpuses’s sparrows data  # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
DATA <- read.table("sparrows.csv", sep = ";", dec = ".", header=TRUE);


# Standardized data
DATA[,1:5] <- scale(DATA[,1:5])
ident.sparrows <- DATA$status
group1 = DATA %>% filter(status == "Survived") %>% select(-status)
group2 = DATA %>% filter(status == "Non-survived")%>% select(-status)
n.sample = c(21, 28)
pval = disp.pval(list(group1, group2), n.perm = 999, distance = "euclidean", binary = FALSE)
names(pval) = c("GO","GO.perm","PERMANOVA", "PERMDISP", "KS", "Wasserstein", "DHT")


phi = list()
rel.abund = list(group1, group2)
n.g = 2
for(i in 1:n.g)
{
  phi[[i]] = as.vector(unlist(vegdist(rel.abund[[i]], method = "euclidean", binary = FALSE)))
}


#n.sample = c(nrow(rel1), nrow(rel2))
Survival_status = factor(rep(c("Survived", "Non-survived"), times = n.sample*(n.sample - 1)/2))
data_vi1 = data.frame(disp = c(phi[[1]], phi[[2]]), Survival_status)
p_smoke <- ggplot(data_vi1, aes(x=Survival_status, y=disp, fill = Survival_status)) + 
  geom_violin(trim=FALSE)   + scale_fill_brewer(palette=4)+ 
  geom_boxplot(width=0.1)+
  labs(
    title = "Bumpus' Sparrows",
    x = "Survival Status of Sparrows",
    y = "Within group distance"
  )+
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 12,face = "bold", hjust = 0.5),
        panel.border = element_rect(color = "black", size = 1, fill = NA),legend.position = "none")
p_smoke
ggsave("sparrows_violin.png", p_smoke, device = "png",
       width = 4, height = 4, units = "in", bg = "white",
       dpi = 1200)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # The Tikus Island corals data  # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # 



load("tikus.RData")
group = list()
group[[1]] = tikus$abund[tikus$x$time == 81,]
group[[2]] = tikus$abund[tikus$x$time == 83,]
group[[3]] = tikus$abund[tikus$x$time == 84,]
group[[3]] = tikus$abund[tikus$x$time == 85,]
group[[4]] = tikus$abund[tikus$x$time == 87,]
group[[5]] = tikus$abund[tikus$x$time == 84,]
pval = disp.pval(group, n.perm = 999, distance = "bray", binary = FALSE)
names(pval) = c("GO","PERMANOVA", "PERMDISP", "KS", "Wasserstein")