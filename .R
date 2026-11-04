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
  dplyr::mutate(UID = paste0('UID', seq_along(1:n())))%>%   
  rename_with(~str_remove(., '\\.d$')) %>%
  rename_with(~gsub(".*CB_MPB.", "", .)) %>% 
    dplyr:: mutate(across(where(is.numeric), log2))  # log2 transform LFQ data
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

# PHENO METADATA
pheno <- data.frame(samplesnames = colnames(df_prot_clean %>%  
                          dplyr::select(-UID))) %>% 
  mutate(samples = gsub(samplesnames, pattern = "GM.*", replacement = "")) %>%
  mutate(samples = gsub(samplesnames, pattern = "S4.*", replacement = "")) %>%
  mutate(ID= samples) %>%
  mutate(ct = factor(ID,levels=unique(ID))) %>%
  ungroup() 
## For cb,mpb and plts 
##function for making pheno automated
parse_sample_name <- function(samplesnames,
                               cell_types_with_day = params$cell_types_with_day,
                               cell_types_no_day   = params$cell_types_no_day) {
  samplesnames <- as.character(samplesnames)
  pattern_with_day <- paste0("^(", paste(cell_types_with_day, collapse = "|"), ")([0-9]+)\\.([0-9]+)\\.")
  pattern_no_day   <- paste0("^(", paste(cell_types_no_day, collapse = "|"), ")([0-9]+)")

  m_day    <- stringr::str_match(samplesnames, pattern_with_day)
  m_no_day <- stringr::str_match(samplesnames, pattern_no_day)

  tibble::tibble(
    samplesnames = samplesnames,
    cell_type = dplyr::coalesce(m_day[, 2], m_no_day[, 2]),
    replicate = as.integer(dplyr::coalesce(m_day[, 3], m_no_day[, 3])),
    day       = suppressWarnings(as.integer(m_day[, 4]))
  )
}
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

### Imputation
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
######################################################
## Limma comparison  ##
######################################################
triplicate_parsed <- parse_sample_name(pheno$samplesnames)
pheno <- pheno %>%
  dplyr::mutate(triplicate = dplyr::case_when(
    triplicate_parsed$cell_type %in% params$cell_types_with_day ~
      paste0(triplicate_parsed$cell_type, triplicate_parsed$replicate, "_day", triplicate_parsed$day),
    triplicate_parsed$cell_type %in% params$cell_types_no_day ~
      paste0(triplicate_parsed$cell_type, triplicate_parsed$replicate),
    TRUE ~ "PLTs"
  ))
pheno$time <- factor(pheno$time, levels = c("day5", "day8", "day11", "day14"))
pheno$cell_type <-as.factor(pheno$cell_type)
pheno$combGroup <- paste(pheno$cell_type, pheno$time, sep = ".")
design_matrix <- model.matrix(~ 0 + ct, data = pheno)
colnames(design) <- levels(pheno$ct)
# Blocking factor (e.g. for triplicates)
blocking_factor <- as.factor(pheno$triplicate)
contrast_definitions <- c("MPB_day5 - CB_day5","MPB_day8 - CB_day8","MPB_day11 - CB_day11","MPB_day14 - CB_day14", "MPB_day8 - MPB_day5", "MPB_day11 - MPB_day5", "MPB_day14 - MPB_day5","CB_day8 - CB_day5", "CB_day11 - CB_day5", "CB_day14 - CB_day5")
# Construct contrast matrix
cont.matrix <- makeContrasts(
  # 1. Cell types at each time point
  MPB_vs_CB_day5  = MPB_day5 - CB_day5,
  MPB_vs_CB_day8  = MPB_day8 - CB_day8,
  MPB_vs_CB_day11 = MPB_day11 - CB_day11,
  MPB_vs_CB_day14 = MPB_day14 - CB_day14,
  # 2. Changes over time within each cell type
  CB_day8_vs_day5 = CB_day8 - CB_day5,
  CB_day11_vs_day5 = CB_day11 - CB_day5,
  CB_day14_vs_day5 = CB_day14 - CB_day5,
  MPB_day8_vs_day5 = MPB_day8 - MPB_day5,
  MPB_day11_vs_day5 = MPB_day11 - MPB_day5,
  MPB_day14_vs_day5 = MPB_day14 - MPB_day5,
  # 3. Compare each MPB day14 against platelets (PLTs)
#  MPB_day14_vs_PLTs = MPB_day14 - PLTs,
  levels = design_matrix
)
print(cont.matrix)
# Estimate duplicate correlation
corfit <- duplicateCorrelation(
  df_prot_raw2 %>%
    column_to_rownames("UID") %>%
    dplyr::select(pheno$samplesnames),
  design = design,
  block = blocking_factor)  
fit <- df_prot_imp %>%
  na.omit() %>%
  column_to_rownames("UID") %>%
  dplyr::select(pheno$samplesnames) %>%
  lmFit(
    design = design_matrix,
    block = blocking_factor,
    correlation = corfit$consensus.correlation)
# Fit contrasts, compute empirical Bayes statistics
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrasts)
sig_tests <- decideTests(fit_ebayes, p.value = de_pvalue_cutoff, adjust.method = "BH", lfc = de_lfc_cutoff)
# Extract and annotate significant results
df_sigresults <- sig_tests %>%
   #dplyr::select(-c('CB_day14_vs_day5','MPB_day14_vs_day5','Interaction_d14_vs_d5')) %>%
  as.data.frame() %>%
  filter(if_any(where(is.numeric), ~ .x >= 1 | .x <= -1)) %>%
  rownames_to_column("UID") %>%
  left_join(df_uid_information2, by = "UID")
# Extracte significant results
df_sigresults_day_comparison <- sig_tests %>%
  as.data.frame() %>%
   dplyr::select(-c("MPB_day8 - MPB_day5","MPB_day11 - MPB_day5","MPB_day14 - MPB_day5","CB_day8 - CB_day5","CB_day11 - CB_day5","CB_day14 - CB_day5")) %>%
  as.data.frame() %>%
  filter(if_any(where(is.numeric), ~ .x >= 1 | .x <= -1)) %>%
  rownames_to_column("UID") %>%
  left_join(df_uid_information2, by = "UID")
df_sigresults_mpb_ <- sig_tests %>%
  as.data.frame() %>%
   dplyr::select(-c("MPB_day5 - CB_day5","MPB_day8 - CB_day8","MPB_day11 - CB_day11","MPB_day14 - CB_day14","CB_day8 - CB_day5","CB_day11 - CB_day5","CB_day14 - CB_day5")) %>%
  as.data.frame() %>%
  filter(if_any(where(is.numeric), ~ .x >= 1 | .x <= -1)) %>%
  rownames_to_column("UID") %>%
  left_join(df_uid_information2, by = "UID")
df_sigresults_cb_ <- sig_tests %>%
  as.data.frame() %>%
   dplyr::select(-c("MPB_day5 - CB_day5","MPB_day8 - CB_day8","MPB_day11 - CB_day11","MPB_day14 - CB_day14","MPB_day8 - MPB_day5","MPB_day11 - MPB_day5","MPB_day14 - MPB_day5")) %>%
  as.data.frame() %>%
  filter(if_any(where(is.numeric), ~ .x >= 1 | .x <= -1)) %>%
  rownames_to_column("UID") %>%
  left_join(df_uid_information2, by = "UID")

######################################################
## WGCNA  ##
######################################################







######################################################
## Rcistarget  ##
######################################################





######################################################
## Figure 2  ##
######################################################


##Figure 2B
Mks_annotations <- data.frame(
  Proteins = c("RUNX1", "GATA1","GF1B",       
            "MECOM", "ETV6", "CXCR4", "MPL", "ITGA2B", "GP1BB", "PF4" ),      
  protein_type = c("MKs markers", "MKs markers","MKs markers",
                   "MKs markers", "MKs markers","MKs markers","MKs markers", "MKs markers", "MKs markers", "MKs markers")
)
df_long <- df_prot_valid %>%
  dplyr::select(all_of(pheno$samplesnames), UID) %>%
  gather("samplesnames", "Intensity", -UID) %>%
  filter(Intensity != 0, is.finite(Intensity)) %>%
  left_join(df_uid_information, by = "UID") %>%
  left_join(pheno, by = "samplesnames")

# Calculate mean intensity per UID per ctcol
df_summary <- df_long %>%
  group_by(UID, Genes, Protein.Group, ctcol) %>%
  summarise(Intensity = mean(Intensity), .groups = "drop") %>%
  group_by(ctcol) %>%
  mutate(
    total = sum(Intensity),
    relative = Intensity / total,
    ranks = rank(-relative)
  ) %>%
  left_join(Mks_annotations, by = "Genes") %>%
  mutate(protein_type = ifelse(is.na(protein_type), "Other", protein_type))
# Plot
plot_abundance <- df_summary %>%
  ggplot(aes(x = ranks, y = Intensity)) +
   geom_point(aes(color = protein_type), size = 4, alpha = 0.8) +
    geom_point(
    data = . %>% filter(granule_type != "Other"),
    aes(color = granule_type),
    shape = 21,
    fill = NA,
    stroke = 1.2,
    size = 4
  ) +
   scale_y_continuous(limits = c(10,25)) +
 scale_x_continuous(limits = c(0,8000)) +
  geom_label_repel(
    data = . %>% filter(granule_type != "Other"),
    aes(label = Genes, fill = granule_type),
    color = "black",
    size = 5.5,
    label.size = 0.25,
    label.r = unit(0.2, "lines"),
    box.padding = 0.5,
    max.overlaps = 15,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("MKs markers" = "#B22222", "Other" = "grey75")
  ) +
  scale_fill_manual(
    values = c("MKs markers" = "#F08080", "Other" = "grey80")
  ) +
  labs(
    x = "Relative Rank",
    y = "Mean Protein Abundance",
    color = "Granule Type",
    fill = "Granule Type",
    title = "Abundance Ranking of Proteins by Cell Type"
  ) +
  facet_wrap(~ ctcol) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )
# Save
ggsave(file.path( "abundance_ranking_by_celltype.pdf"), plot_abundance, width = 12, height = 8, dpi = 300)







######################################################
## Figure 3  ##
######################################################
##Figure 3A
#run UMAP
umap_result <- df_prot_imp %>%
  column_to_rownames("UID") %>%
  na.omit() %>%
  t() %>%
  umap()
# Convert result to data frame
umap_df <- as.data.frame(umap_result$layout)
umap_df$samplesnames <- rownames(umap_df)
# Join phenotype data
umap_df <- left_join(umap_df, pheno, by = "samplesnames")
# Create plot
umap_gplot <- ggplot(umap_df, aes(x = V1, y = V2, color = ct, shape = ctcol)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "UMAP Projection of Proteomic Profiles",
    x = "UMAP-1",
    y = "UMAP-2",
    color = "Cell Type",
    shape = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
#Set the right dimensions
umap_gplot <- umap_gplot +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_x_continuous(breaks = -3:3)
# Print the plot
print(umap_gplot)
# Save the plot to a PDF file
pdf(file.path("plt_umap_plot.pdf"), width = 12, height = 6)
print(umap_gplot, column_title = "")
dev.off()

######################################################
## Figure 4  ##
######################################################





######################################################
## Figure 5  ##
######################################################






