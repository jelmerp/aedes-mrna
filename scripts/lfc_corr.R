## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",        # Misc. data manipulation and plotting
              "ggpubr",           # Stats for plots
              "here",             # Managing file paths
              "readxl")           # Read Excel files             
pacman::p_load(char = packages)

## Input and output files
qpcr_file <- "results/qpcr/lfc_correlation.xls"
plotfile <- "results/qpcr/qpcr_plot.png"

## Read input file
qpcr <- read_xls(qpcr_file) %>%
  mutate(geneID = sub(" .*", "", geneID))

## Test for normality
shapiro.test(qpcr$lfc_RNA)
shapiro.test(qpcr$lfc_qPCR)

## QQ plots
ggqqplot(qpcr$lfc_RNA)
ggqqplot(qpcr$lfc_RNA)

## Create the plot
p <- ggplot(qpcr) +
  aes(x = lfc_RNA, y = lfc_qPCR) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.8) +
  geom_point(shape = 21, fill = "red", size = 2.5, stroke = 1) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  labs(x = "Log-fold change in RNAseq",
       y = "Log-fold change in qPCR") +
  ggpubr::stat_cor(method = "pearson", size = 4, color = "grey40",
                   label.x.npc = 0.38, label.y.npc = 0.99) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
p

## Save the plot
ggsave(plotfile, p, width = 6, height = 6, bg = "white")
