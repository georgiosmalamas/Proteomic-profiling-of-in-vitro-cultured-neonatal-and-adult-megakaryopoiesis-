######################################################
## Loading neccessary libraries and seed  ##
######################################################
set.seed(42)
## libraries
library(limma)
library(viridis)
library(tidyverse)
library(readxl)
library(ggrepel)
library(ComplexHeatmap)
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(WGCNA)
library(clusterProfiler)
library(ReactomePA)
library(openxlsx)
library(ggtext)
library(RcisTarget)
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp) 
library(DT)
library(tibble)
library(writexl)
library(showtext)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(ggpubr)
library(viridis)
library(eulerr)
library(gridExtra)
library(grid)
library(umap)
library(patchwork)
library(svglite)
library(reshape2)
library(ggalluvial)
library(msigdbr)
######################################################
## QC,  Imputation  ##
######################################################

pr_matrix_path <- "report.pr_matrix.tsv"
pg_matrix_path <- "report.pg_matrix.tsv"

dia_pq <- read.delim(pr_matrix_path) %>% 
  dplyr::select(Protein.Group, Genes, Proteotypic, Precursor.Id, matches("\\.d$")) %>%
  tidyr::gather(key = "Sample", value = "Intensity", matches("\\.d$")) %>%
  filter(!is.na(Intensity)) %>%
  filter(Proteotypic ==1) %>%
  group_by(Sample, Protein.Group) %>%
  distinct(Precursor.Id, .keep_all = TRUE) %>%
  dplyr::summarise(count = n()) %>%
  filter(count > 1) %>%
  dplyr::select(-count)

#  proteins
protein_summary <- read.delim(pg_matrix_path) %>%
  dplyr::select(Protein.Group, Protein.Names, Genes, First.Protein.Description, matches("\\.d$")) %>%
  tidyr::gather(key = "Sample", value = "Intensity", matches("\\.d$")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  mutate(Sample = gsub(pattern = ".*CB_MPB.", replacement = "", x = Sample)) %>%
  dplyr::summarise(count_protein = n())


#  peptides
peptide_summary <- read.delim(pr_matrix_path) %>%
  dplyr::select(Protein.Group, Genes, Proteotypic, Stripped.Sequence, matches("\\.d$")) %>%
  gather(key = "Sample", value = "Intensity", matches("\\.d$")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  distinct(Stripped.Sequence, .keep_all = T) %>%
  mutate(Sample = gsub(Sample, pattern = ".*CB_MPB.", replacement = "")) %>%
  dplyr::summarise(count_peptide = n())

#  precursors
precursors_summary2 <- read.delim(pr_matrix_path) %>%
  dplyr::select(Protein.Group, Genes, Proteotypic, Precursor.Id, matches("\\.d$")) %>%
  gather(key= "Sample", value= "Intensity", matches("\\.d$")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  distinct(Precursor.Id, .keep_all = T) %>%
  mutate(Sample = gsub(Sample, pattern = ".*CB_MPB.", replacement = "")) %>%
  dplyr::summarise(count_precursor = n())

# visualization
protein_summary2 %>% 
  full_join(peptide_summary2, "Sample") %>%
  full_join(precursors_summary2, "Sample") %>%
  mutate(Sample = gsub(Sample, pattern = "\\.d$", replacement = "")) %>%
  mutate(Sample = gsub(Sample, pattern = ".*CB_MPB.", replacement = "")) %>%
  gather("Type", "Count", -Sample) %>% # 
  ggplot(aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "Darkblue") +
  facet_wrap(~Type, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_blank())

# Processing

df_prot_raw <- read.delim(pg_matrix_path) %>%
  # Filter out all values based on 1 or less precursors
  pivot_longer(cols = matches("\\.d$"), names_to = "Sample", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  inner_join(dia_pq, by = c("Sample", "Protein.Group")) %>%
  group_by(Sample) %>%
  filter(!is.na(Intensity)) %>%
  pivot_wider(names_from = Sample, values_from = Intensity) %>%
  # Add UID
  dplyr::mutate(UID = paste0('UID', seq_along(1:n())))%>%   # Clean up column names
  rename_with(~str_remove(., '\\.d$')) %>%
  rename_with(~gsub(".*CB_MPB.", "", .)) %>% 
  
  # log2 transform LFQ data
 dplyr:: mutate(across(where(is.numeric), log2))


# Distirbution check
df_prot_raw2 %>% 
  dplyr::select(-Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description, -UID) %>%
  gather() %>%
  # filter(!grepl(key, pattern = "blank")) %>%
  ggplot(aes(x = key, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "bottom")+
  labs(fill ="", x = "Sample", y = "LFQ intensity")

















