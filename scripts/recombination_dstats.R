library(ggplot2)
library(dplyr)
library(purrr)
library(ggpubr, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")  
library(cowplot, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/dsuite")

# Get list of all chromosome files (assuming they are in a folder named "chrom_files/")
file_list <- list.files(path = "/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr/", pattern = "*.scaffolds", full.names = TRUE)

# Create a lookup table: A dataframe mapping scaffolds to chromosome numbers
chromosome_mapping <- map_df(file_list, function(file) {
  chrom_number <- tools::file_path_sans_ext(basename(file))  # Extract chromosome number from filename
  scaffolds <- readLines(file)  # Read scaffold names from file
  data.frame(scaffold = scaffolds, chromosome = chrom_number)})
# using Vijay's final_rho_table from 2016 paper
recom <- read.csv("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst/final_rho_table", 
                  header=TRUE, sep = "\t")
recom_mean <- recom %>%
  group_by(Scaffold) %>%
  summarise(mean_recomb_rate = mean(cor1, na.rm = TRUE))

recom_mean_chr <- merge(chromosome_mapping, recom_mean, by.x="scaffold", by.y="Scaffold") 

#########################################################################################################################
########### dstats vs chrlength ###########
data1 <- read.table(file="cor1_cor2_cnx3P_localFstats_trios3name_50_25.txt",header=TRUE)
data2 <- read.table(file="cnx6_cnx3P_cor2_localFstats_trios3name_50_25.txt",header=TRUE)

merged_df2 <- merge(data1,recom_mean_chr, by.x="chr",by.y = "scaffold") %>%
  mutate(length=windowEnd-windowStart+1) %>%
  mutate(f_d = if_else(f_d < 0, NA_real_, f_d)) %>%
  mutate(d_f = if_else(d_f < 0, NA_real_, d_f)) %>%
  mutate(f_dM = if_else(f_dM < 0, NA_real_, f_dM)) %>%
  group_by(chromosome) %>%
  summarise(length=sum(length, na.rm=TRUE),
            mean_df = mean(d_f, na.rm = TRUE),
            mean_D = mean(D, na.rm = TRUE),
            mean_fd = mean(f_d, na.rm = TRUE),
            mean_fdm = mean(f_dM, na.rm=TRUE),
            mean_recom=mean(mean_recomb_rate, na.rm=TRUE)) %>%
  filter(!str_detect(chromosome, "chrUN"))

merged_df3 <- merge(data2,recom_mean_chr, by.x="chr",by.y = "scaffold") %>%
  mutate(length=windowEnd-windowStart+1) %>%
  mutate(f_d = if_else(f_d < 0, NA_real_, f_d)) %>%
  mutate(d_f = if_else(d_f < 0, NA_real_, d_f)) %>%
  mutate(f_dM = if_else(f_dM < 0, NA_real_, f_dM)) %>%
  group_by(chromosome) %>%
  summarise(length=sum(length, na.rm=TRUE),
            mean_df = mean(d_f, na.rm = TRUE),
            mean_D = mean(D, na.rm = TRUE),
            mean_fd = mean(f_d, na.rm = TRUE),
            mean_fdm = mean(f_dM, na.rm=TRUE),
            mean_recom=mean(mean_recomb_rate, na.rm=TRUE)) %>%
  filter(!str_detect(chromosome, "chrUN"))

# check correlation between mean recombination and chr length
a1 <- ggplot(merged_df2, aes(x=length, y=mean_recom)) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2,size=0.5) +
  geom_point(aes(x = length, y = mean_recom), color = "red", alpha=0.6, size = 2) +
  labs(y = "Mean recombination \n rate (rho)", x = "Chromosome length (bp)") +
  scale_x_log10(labels = scales::scientific) + 
  scale_y_log10() + 
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.8) +
  stat_cor(aes(label = paste(format.pval(..p.value.., digits = 3))), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.7) +
  geom_text_repel(data=merged_df2,aes(x = length, y = mean_recom, label=chromosome), size=3,show.legend = FALSE) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# check correlation between dfm and chr length
b3 <- ggplot(merged_df2, aes(x = length, y = mean_fdm)) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2,size=0.5) +
  geom_point(aes(x = length, y = mean_fdm), color = "red", alpha=0.6, size = 2) +  
  scale_x_log10(labels = scales::scientific) + 
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.8) +
  stat_cor(aes(label = paste(format.pval(..p.value.., digits = 3))), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.75) +
  labs(y = "Allele sharing between EURw and EURn \n [f_dm of (((SPA,EURw),EURn),O)]", x = "Chromosome length (bp)") +
  geom_text_repel(data=merged_df2,aes(x = length, y = mean_fdm, label=chromosome), size=3,show.legend = FALSE) +
  #geom_text_repel(data=merged_df2%>%filter(chromosome=="chr18"),aes(x = length, y = mean_fdm, label=chromosome), size=3,show.legend = FALSE) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

c3 <- ggplot(merged_df3, aes(x = length, y = mean_fdm)) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2,size=0.5) +
  geom_point(aes(x = length, y = mean_fdm), color = "red", alpha=0.6, size = 2) +  
  scale_x_log10(labels = scales::scientific) + 
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.8) +
  stat_cor(aes(label = paste(format.pval(..p.value.., digits = 3))), method = "pearson", color = "black", size = 3, label.x.npc = 0.8, label.y.npc = 0.75) +
  labs(y = "Allele sharing between EURw and EURn \n [f_dm of (((IRQ,EURn),EURw),O)]", x = "Chromosome length (bp)") +
  geom_text_repel(data=merged_df3,aes(x = length, y = mean_fdm, label=chromosome), size=3,show.legend = FALSE) +
  #geom_text_repel(data=merged_df2%>%filter(chromosome=="chr18"),aes(x = length, y = mean_fdm, label=chromosome), size=3,show.legend = FALSE) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

plot1 <- plot_grid(a1, plot_grid(b3, c3,labels = c("(b)", "(c)")), labels = "(a)", ncol = 1, rel_heights =c(1,2))
pdf(file="recombination_lengths_combined.pdf",width=8,height=8)
plot1
dev.off()
