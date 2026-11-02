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

      ### Processing

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
df_prot_raw %>% 
  dplyr::select(-Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description, -UID) %>%
  gather() %>%
  # filter(!grepl(key, pattern = "blank")) %>%
  ggplot(aes(x = key, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "bottom")+
  labs(fill ="", x = "Sample", y = "LFQ intensity")


df_prot_clean <- df_prot_raw %>%
  dplyr::select(-Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>%
  dplyr::select(UID, everything()) 

df_uid_information <- df_prot_raw %>%
  dplyr::select(UID, Protein.Group, Protein.Ids, Protein.Names, Genes, First.Protein.Description)


# AANPASSEN
pheno <- data.frame(samplesnames = colnames(df_prot_clean2 %>%  
                          dplyr::select(-UID))) %>% 
  mutate(samples = gsub(samplesnames, pattern = "GM.*", replacement = "")) %>%
  mutate(samples = gsub(samplesnames, pattern = "S4.*", replacement = "")) %>%
  mutate(ID= samples) %>%
  mutate(ct = factor(ID,levels=unique(ID))) %>%
  ungroup() 

## For cb,mpb and plts 

pheno_parsed <- parse_sample_name(pheno$samplesnames)
pheno <- pheno %>%
  dplyr::mutate(ID = dplyr::case_when(
    pheno_parsed$cell_type %in% params$cell_types_with_day ~ paste0(pheno_parsed$cell_type, "_day", pheno_parsed$day),
    TRUE ~ "PLTs"
  ))


# Process the data
vvs = df_prot_clean %>% 
      column_to_rownames('UID') %>% 
      dplyr::mutate(across(.cols = everything(), ~ !is.na(.))) %>% # is a value present?
      t() %>% 
      as.data.frame() %>% 
      add_column(ct = pheno$ct) %>% 
      group_by(ct) %>% 
      dplyr::summarise_all(sum) %>% 
      column_to_rownames('ct') %>%
      dplyr::mutate(across(.cols = everything(), ~ . / as.numeric(table(list(pheno$ct))))*100) %>% 
      t() %>%
      as.data.frame() %>% 
      dplyr::mutate(No_pass = rowSums(across(.cols = everything(), ~ . >= valid_value_threshold ))) %>% # valid value threshold
      rownames_to_column('UID')                              

##Select valid_values per percentage (>65%) 
valid_value_threshold <- 65
df_prot_valid = df_prot_clean %>% 
              subset(UID %in% vvs$UID[vvs$No_pass >=1])
# plot of the amount of quantified proteins per type
vvs %>% 
    dplyr::select(pheno$ct) %>% 
    gather() %>% 
    filter(value >= valid_value_threshold) %>% # Threshold
    group_by(key=factor(key, levels = levels(pheno$ct))) %>% 
    dplyr::summarise(Freq=n()) %>% 
  ungroup() %>%
    # Visualization
    ggplot(aes(x=key, y=Freq)) + #, fill=as.factor(rev(round(value)))))+
          geom_bar(stat='identity', fill = "darkblue")+
          scale_y_continuous(limits=c(0,8500), breaks=seq(0,7500,1000), expand=c(0,0))+
          scale_fill_brewer()+
          xlab('')+
          ylab('#')+
          ggtitle('Valid values (%)')+
          theme_classic()+
          theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),legend.title = element_blank())
      #300x250



# impute
# define parameters
mann_downshift = imputation_downshift
mann_width = imputation_width

# create imputed values
imputed_values = df_prot_valid %>% 
  dplyr::select(pheno$samplesnames[order(pheno$ct)]) %>% 
  dplyr::mutate(across(.cols = everything(), ~ rnorm(length(.), (mean(na.omit(.) -  mann_downshift*sd(na.omit(.)))), (sd(na.omit(.))*mann_width))))

# replace NA's in DF with imputed values using coalesce (no need for a mask in tidyverse)  
df_prot_imp = df_prot_valid %>% 
  dplyr::select(colnames(imputed_values)) %>%
  # select(pheno$samplesnames[order(pheno$ct)]) %>% 
  map2_df(.,imputed_values2,coalesce) %>% 
  bind_cols(UID = df_prot_valid2$UID)
# dit laatste kan ws beter, bind met een charateristics df?, maar coalsce mag geen rownames hebben... 
                
# check
df_prot_imp %>% 
  dplyr::select(pheno$samplesnames[order(pheno$ct)]) %>% 
  boxplot()








