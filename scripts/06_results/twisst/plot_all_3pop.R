install.packages("ggtern")

library(ape, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")
library(ggtern)
library(ggplot2)
library(dplyr)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst")
source("plot_twisst.R")

hash = paste("#topo1 (O,((EURn,SPA),EURc));",
             "#topo2 (O,((EURn,EURc),SPA));",
             "#topo3 (O,(EURn,(EURc,SPA)));", sep="\n")

weights_file <- "all_3pop.weights.edited.csv.gz"
window_data_file <- "all_3pop.phyml_bionj.w50.data.tsv"
w1 <- read.csv(weights_file, skip = 3, header=TRUE, sep = "\t")
t1 <- read.csv(window_data_file, header=TRUE, sep = "\t")
dat <- cbind(t1,w1)


# writeLines(hash, con = "all_3pop.weights.edited.name.csv")
# write.table(dat[7:length(dat)],file="all_3pop.weights.edited.name.csv", sep= "\t", quote=FALSE, append=TRUE)
# system("gzip all_3pop.weights.edited.name.csv")

# twisst_data <- import.twisst(weights_files="all_3pop.weights.edited.name.csv.gz", window_data_files=window_data_file)
# pdf("all_3pop_plotsummary_edited.pdf") #6x14
# plot.twisst.summary(twisst_data, lwd=3, cex=1)
# dev.off()

######################### identifying scaffold78 and 60 region ###############################

dat1 <- dat %>% group_by(scaffold) %>%
  filter(scaffold != "scaffold_78" & scaffold != "scaffold_60")

w1edit <- dat1[7:length(dat)]
t1edit <- dat1[1:6]

# writeLines(hash, con = "all_3pop.weights.edited.nochr18.csv")
# write.table(w1edit,file="all_3pop.weights.edited.nochr18.csv", sep= "\t", quote=FALSE, append=TRUE)
# system("gzip all_3pop.weights.edited.nochr18.csv")
# write.table(t1edit,file="all_3pop.phyml_bionj.w50.data.nochr18.tsv", sep= "\t", quote=FALSE)

weights_file_nochr18 <- "all_3pop.weights.edited.nochr18.csv.gz"
window_data_file_nochr18 <- "all_3pop.phyml_bionj.w50.data.nochr18.tsv"
# twisst_data_nochr18 <- import.twisst(weights_files=weights_file_nochr18, window_data_files=window_data_file_nochr18)
# plot.twisst.summary(twisst_data_nochr18, lwd=3, cex=1)

dat2 <- dat %>% group_by(scaffold) %>%
  filter(scaffold == "scaffold_78" | scaffold == "scaffold_60")

w2edit <- dat2[7:length(dat)]
t2edit <- dat2[1:6]

# writeLines(hash, con = "all_3pop.weights.edited.chr18.csv")
# write.table(w2edit,file="all_3pop.weights.edited.chr18.csv", sep= "\t", quote=FALSE, append=TRUE)
# system("gzip all_3pop.weights.edited.chr18.csv")
# write.table(t2edit,file="all_3pop.phyml_bionj.w50.data.chr18.tsv", sep= "\t", quote=FALSE)

weights_file_chr18 <- "all_3pop.weights.edited.chr18.csv.gz"
window_data_file_chr18 <- "all_3pop.phyml_bionj.w50.data.chr18.tsv"
# twisst_data_chr18 <- import.twisst(weights_files=weights_file_chr18, window_data_files=window_data_file_chr18)
# plot.twisst.summary(twisst_data_chr18, lwd=3, cex=1)


########################### ternary plot ###########################
dattern <- mutate(dat, type = case_when (scaffold == "scaffold_78" | scaffold == "scaffold_60" ~ "barrier",
                                         scaffold != "scaffold_78" | scaffold != "sacffold_60" ~ "neutral"))

#### MAIN ####
# pdf("ternaryplot_3pop_all.pdf")
ggtern(data=dattern,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density of sites",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")
# dev.off()
##only neutral sites
emp7 <- ggtern(dattern %>% filter(type == "neutral"),aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density of sites",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  ggtitle("Observed neutral sites") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
# pdf("ternaryplot_3pop_neutral.pdf")
emp7
# dev.off()

#### ternary contour ####
## conducting density estimation on barrier (scaffold78+60) and neutral sites separately
c1 <- ggtern(dattern, aes(x = topo1,y = topo2, z = topo3, color = type)) +
  stat_density_tern(aes(alpha = ..level.., fill = type), 
                    geom = 'polygon', 
                    bins = 50,
                    color = "grey",bdl=0.005) +
  #geom_point(alpha = 0.4, size=1) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  theme_hidegrid_major()
## conducting density estimation on all loci collectively
c2 <- ggtern(dattern, aes(x = topo1,y = topo2, z = topo3)) +
  stat_density_tern(aes(alpha = ..level..), 
                    geom = 'polygon', 
                    bins = 50,
                    color = "grey",bdl=0.005) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(fill = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  theme_hidegrid_major()

# pdf("ternarycontour_3pop.pdf",width = 14, height = 8.5)
grid.arrange(c2,c1,ncol=2)
# dev.off()

############ simulated twisst ############
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst")
weights_file_3popsim7 <- "./3PopModel7xtree_run2/3PopModel7x_run2.weights.csv.gz"
weights_file_3popsim8 <- "./3PopModel8xtree_run2/3PopModel8x_run2.weights.csv.gz"
window_data_file_3popsim7 <- "./3PopModel7xtree_run2/3PopModel7x_run2.win.data.tsv"
window_data_file_3popsim8 <- "./3PopModel8xtree_run2/3PopModel8x_run2.win.data.tsv"

#edit files due to name changes
w3popsim7 <- read.csv(weights_file_3popsim7, skip = 3, header=TRUE, sep = "\t")
colnames(w3popsim7) <- c("topo3","topo1","topo2")
w3popsim7 <- relocate(w3popsim7,topo3, .after = topo2)
t3popsim7 <- read.csv(window_data_file_3popsim7, header=TRUE, sep = "\t")
dat3popsim7 <- cbind(t3popsim7 ,w3popsim7)
w3popsim8 <- read.csv(weights_file_3popsim8, skip = 3, header=TRUE, sep = "\t")
colnames(w3popsim8) <- c("topo3","topo1","topo2")
w3popsim8 <- relocate(w3popsim8,topo3, .after = topo2)
t3popsim8 <- read.csv(window_data_file_3popsim8, header=TRUE, sep = "\t")
dat3popsim8 <- cbind(t3popsim8 ,w3popsim8)

#writeLines(hash, con = "./3PopModel7jaathatree_run1/3PopModel7jaatha_run1.weights.edited.csv")
#write.table(w3popsim7,file="./3PopModel7jaathatree_run1/3PopModel7jaatha_run1.weights.edited.csv", sep= "\t", quote=FALSE, append=TRUE)
#system("gzip ./3PopModel7jaathatree_run1/3PopModel7jaatha_run1.weights.edited.csv")
#writeLines(hash, con = "./3PopModel8jaathatree_run1/3PopModel8jaatha_run1.weights.edited.csv")
#write.table(w3popsim8,file="./3PopModel8jaathatree_run1/3PopModel8jaatha_run1.weights.edited.csv", sep= "\t", quote=FALSE, append=TRUE)
#system("gzip ./3PopModel8jaathatree_run1/3PopModel8jaatha_run1.weights.edited.csv")

twisst_dat3popsim7 <- import.twisst(weights_files=weights_file_3popsim7, window_data_files=window_data_file_3popsim7)
#pdf("./3PopModel7xtree_run2/3PopModel7x_plotsummary2.pdf")
#plot.twisst.summary(twisst_dat3popsim7, lwd=3, cex=1)
#dev.off()

twisst_dat3popsim7 <- import.twisst(weights_files=weights_file_3popsim8, window_data_files=window_data_file_3popsim8)
#pdf("./3PopModel8xtree_run2/3PopModel8x_plotsummary2.pdf")
#plot.twisst.summary(twisst_dat3popsim8, lwd=3, cex=1)
#dev.off()

sim7 <- ggtern(dat3popsim7,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill = "Density",
    Tarrow = "ILS + Introgression",
    Larrow = "ILS",
    Rarrow = "Barrier")+
  ggtitle(paste("Simulated neutral sites from", "genome-wide introgression model", sep = "\n"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15))

sim8 <- ggtern(dat3popsim8,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")+
  ggtitle(paste("Simulated neutral sites from", "locus-specific introgression model",sep="\n")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

pdf("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst/ternaryplot_emp_sim7_sim8_jaatha.pdf",width=16,height=6)
grid.arrange(emp7,sim7,sim8, nrow=1)
dev.off()

#### ternary count for neutral sites ####
p1 <- ggtern(dattern %>% filter(type == "neutral"),aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) 
  #stat_tri_tern(bins=4,geom='text',aes(label=sprintf("%.0f",after_stat(count))),size=2, color='white',centroid=TRUE)
p2 <- ggtern(dat3popsim7,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count)))
  #stat_tri_tern(bins=20,geom='text',aes(label=sprintf("%.0f",after_stat(count))),size=2, color='white',centroid=TRUE)
p3 <- ggtern(dat3popsim8,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) 
  #stat_tri_tern(bins=20,geom='text',aes(label=sprintf("%.0f",after_stat(count))),size=2, color='white',centroid=TRUE)


## using the count to compare against simulated results
emp <- ggplot_build(p1)$data[[1]] %>% group_by(group,count)  %>% 
  summarize() 
plot(ggplot_build(p1))
sim7 <- ggplot_build(p2)$data[[1]] %>% group_by(group,count)  %>% 
  summarize() 

sim8 <- ggplot_build(p3)$data[[1]] %>% group_by(group,count)  %>% 
  summarize() 

comp <- cbind(emp, sim7$count, sim8$count)
colnames(comp) <- c("group","emp","sim7","sim8")
compdiff <- comp %>% mutate(sim7diff=abs(sim7-emp), sim8diff=abs(sim8-emp)) 
sum(compdiff$sim7diff,na.rm=TRUE)
sum(compdiff$sim8diff,na.rm=TRUE)

######################## subset to only specific regions #########################
pdf(file="all_3pop_plotweightings_smooth_windows_indscaff.pdf")
par(mfrow = c(2,2))
for (i in  c(0:15, 17:20, 22:31, 33:38, 40:45, 47:53, 55:63, 66:70, 72:73, 75:100)) {
  regions <- paste("scaffold_",i,sep="")
  twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions)
  twisst_data_subset_smooth <- smooth.twisst(twisst_data_subset, span_bp = 20000, spacing = 1000)
  plot.weights(weights_dataframe=twisst_data_subset_smooth$weights[[1]], positions=twisst_data_subset_smooth$pos[[1]],
               line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)
  mtext(side=3,text=paste("scaffold_",i, sep=""))
}
dev.off()

### chr18 ####
pdf("all_3pop_plotweightings_smooth_windows_chr18.pdf")
par(mfrow = c(1,1))
for (i in c("scaffold_60", "scaffold_70", "scaffold_78")) {
  regions <- i
  twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions)
  twisst_data_subset_smooth <- smooth.twisst(twisst_data_subset, span_bp = 80000, spacing = 1000)
  plot.weights(weights_dataframe=twisst_data_subset_smooth$weights[[1]], positions=twisst_data_subset_smooth$pos[[1]],
               line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)
  mtext(side=3,text=paste(i))
}
dev.off()

#### chr1A ####
pdf(file="all_3pop_plotweightings_smooth_windows_chr1A.pdf")
par(mfrow = c(1,1))
for (i in c("scaffold_11", "scaffold_22", "scaffold_48", "scaffold_62", "scaffold_77", "scaffold_8")) {
  regions <- i
  twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions)
  twisst_data_subset_smooth <- smooth.twisst(twisst_data_subset, span_bp = 80000, spacing = 1000)
  plot.weights(weights_dataframe=twisst_data_subset_smooth$weights[[1]], positions=twisst_data_subset_smooth$pos[[1]],
               line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)
  mtext(side=3,text=paste(i))
}
dev.off()

#### chr1 ####
pdf(file="all_3pop_plotweightings_smooth_windows_chr1.pdf")
par(mfrow = c(1,1))
for (i in c("scaffold_2", "scaffold_23", "scaffold_30", "scaffold_34", "scaffold_35", "scaffold_37", "scaffold_52", "scaffold_53", "scaffold_75", "scaffold_84", "scaffold_97")) {
  regions <- i 
  twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions)
  twisst_data_subset_smooth <- smooth.twisst(twisst_data_subset, span_bp = 80000, spacing = 1000)
  plot.weights(weights_dataframe=twisst_data_subset_smooth$weights[[1]], positions=twisst_data_subset_smooth$pos[[1]],
               line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)
  mtext(side=3,text=paste(i))
}
dev.off()

#### chr9 ####
pdf(file="all_3pop_plotweightings_smooth_windows_chr9.pdf")
par(mfrow = c(1,1))
for (i in c("scaffold_27", "scaffold_51", "scaffold_56", "scaffold_83")) {
  regions <- i 
  twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions)
  twisst_data_subset_smooth <- smooth.twisst(twisst_data_subset, span_bp = 80000, spacing = 1000)
  plot.weights(weights_dataframe=twisst_data_subset_smooth$weights[[1]], positions=twisst_data_subset_smooth$pos[[1]],
               line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)
  mtext(side=3,text=paste(i))
}
dev.off()