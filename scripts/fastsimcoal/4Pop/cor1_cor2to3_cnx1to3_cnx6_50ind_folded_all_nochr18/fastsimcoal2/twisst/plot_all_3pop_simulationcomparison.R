install.packages("ggtern")
library(ape, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")
library(ggtern)
library(ggplot2)
library(dplyr)
library(scales)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst")
source("plot_twisst.R")
weights_file <- "all_3pop.weights.edited.csv.gz"
window_data_file <- "all_3pop.phyml_bionj.w50.data.tsv"
w1 <- read.csv(weights_file, skip = 3, header=TRUE, sep = "\t")
t1 <- read.csv(window_data_file, header=TRUE, sep = "\t")
dat <- cbind(t1,w1)
dattern <- mutate(dat, type = case_when (scaffold == "scaffold_78" | scaffold == "scaffold_60" ~ "barrier",
                                         scaffold != "scaffold_78" | scaffold != "sacffold_60" ~ "neutral"))

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst")
w3popsim7 <- read.csv("./3PopModel7xtree_run2/3PopModel7x_run2.weights.edited.csv.gz", skip = 3, header=TRUE, sep = "\t")
t3popsim7 <- read.csv("./3PopModel7xtree_run2/3PopModel7x_run2.win.data.tsv", header=TRUE, sep = "\t")
dat3popsim7 <- cbind(t3popsim7 ,w3popsim7)
w3popsim8 <- read.csv("./3PopModel8xtree_run2/3PopModel8x_run2.weights.edited.csv.gz", skip = 3, header=TRUE, sep = "\t")
t3popsim8 <- read.csv("./3PopModel8xtree_run2/3PopModel8x_run2.win.data.tsv", header=TRUE, sep = "\t")
dat3popsim8 <- cbind(t3popsim8 ,w3popsim8)

p1 <- ggtern(dattern %>% filter(type == "neutral"),aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) +
  stat_tri_tern(bins=20)
p2 <- ggtern(dat3popsim7,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) +
  stat_tri_tern(bins=20)
p3 <- ggtern(dat3popsim8,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) +
  stat_tri_tern(bins=20)

#grid.arrange(p1,p2,p3,nrow=1)

emp <- ggplot_build(p1)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))
sim7 <- ggplot_build(p2)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))
sim8 <- ggplot_build(p3)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))

comp <- cbind(emp, "sim7"=sim7$count, "sim8"=sim8$count) %>%
  mutate(sim7diff=(sim7-count), sim8diff=(sim8-count))


sim7diff <- data.frame(comp$x,comp$y,comp$z,comp$sim7diff) %>%
  #uncount(abs(comp$sim7diff)) %>%
  rename("x"="comp.x","y"="comp.y","z"="comp.z","count"="comp.sim7diff")
sim8diff <- data.frame(comp$x,comp$y,comp$z,comp$sim8diff) %>%
  rename("x"="comp.x","y"="comp.y","z"="comp.z","count"="comp.sim8diff")

p7 <- ggtern(data = sim7diff, aes(x=x,y=y,z=z)) +
  geom_point(size=4, aes(color=count)) +
  theme_bw() +
  scale_colour_gradient2(low = muted("darkred"),mid = "white",high = muted("darkblue"),midpoint = 0, limits = c(-10000, 10000)) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(color  = "Count",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")  +
  ggtitle("Exp(GWS)-Obs") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

p8 <- ggtern(data = sim8diff, aes(x=x,y=y,z=z)) +
  geom_point(size=4, aes(color=count)) +
  theme_bw() +
  scale_colour_gradient2(low = muted("darkred"),mid = "white",high = muted("darkblue"),midpoint = 0, limits = c(-10000, 10000)) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(fill   = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")  +
  ggtitle("Exp(LSI)-Obs") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

grid.arrange(p7,p8,nrow=1)
