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
sig_tests <- decideTests(fit_ebayes, p.value = 0.05, adjust.method = "BH", lfc = 1)
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
#If you like only the signifcant proteins to appear in the WGCNA
df_prot_imp_sig <- df_prot_imp %>% 
  filter((UID %in% df_sigresults$UID))
df_uid_info_sig <- df_uid_information %>% 
  filter(UID %in% df_prot_imp_sig$UID)
#########WGCNA BEGIN
# Call the network topology analysis function
# and adjacency() below need this raw expression matrix, not a correlation matrix.
df_wgcna_cor_pearson <- df_prot_imp_sig %>%
  left_join(df_uid_info_sig, "UID") %>%
  mutate(ID = paste(UID, Genes, sep = "_")) %>%
  distinct(ID, .keep_all = T) %>%
  filter(!is.na(ID)) %>%
  `rownames<-`(NULL) %>%
  column_to_rownames("ID") %>%
  dplyr::select(unique(pheno_statistics$samplesnames)) %>%
  t()
df_wgcna_cor <- df_wgcna_cor_pearson %>%
  Hmisc::rcorr(type = "pearson")
# df_wgcna_cor$r[df_wgcna_cor$n < 30 ] = 0
sft <- df_wgcna_cor_pearson%>%
  WGCNA::pickSoftThreshold(powerVector = c(c(1:10), seq(12,30,2)), verbose = 1, networkType = 'signed')
# look for the optimal parameters
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9
powers = c(c(1:10), seq(12,30,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# set the soft-power, based on above
# according to wgcna website, above .8 is ok, of does not happen use 20
softPower = 10
# calculate the adjacency, Mann uses a signed network
adjacency = df_wgcna_cor_pearson %>%
  WGCNA::adjacency(power = softPower, type='signed')
# turn the adjacency into a topological map, and make a dissimilarity matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# Clustering of the TOM, using hclust: Mann paper used average
protTree = hclust(as.dist(dissTOM), method='average')
# plot the clustering
# sizeGrWindow(12,9)
plot(protTree, xlab='', sub='', main = 'Protein clustering on TOM-based dissimilarity',labels = FALSE, hang = 0.4)
# now select modules using dynamic tree cut
dynamicMods = dynamicTreeCut::cutreeDynamic(dendro=protTree, distM=dissTOM, deepSplit = 4, pamRespectsDendro = F, minClusterSize = 15)
table(dynamicMods)
# dynamicMods is one module label per protein, in the column order of
# df_wgcna_cor_pearson (i.e. the row order before it was transposed for WGCNA) --
# so transpose back to recover "UID_Genes" IDs as rows before attaching ME.
wgcna_correlation <- t(df_wgcna_cor_pearson) %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  separate(ID, into = c("UID", "Genes"), sep ="_", remove = F) %>%
  mutate(ME = dynamicMods) %>%
  dplyr::select(ID, UID, Genes, ME, contains("UID"))
colnames(wgcna_correlation) <- gsub("UID","Protein.Group", colnames(wgcna_correlation))
rownames(wgcna_correlation) <- wgcna_correlation$Protein.Group
#zscores table creation
zscores_df_prot <- df_prot_imp_sig %>%
  dplyr::mutate(Proteins = paste0(UID)) %>%
  #filter(rowSums(abs(results))>0) %>%
  column_to_rownames('UID') %>%
  dplyr::select(matches('_')) %>%
  t() %>%
  scale() %>%
  t() 
zscores_df_prot <- as.data.frame(zscores_df_prot)
colnames(zscores_df_prot)
tz <- t(zscores_df_prot) 
# Renaming logic for heatmap column labels (replicate-collapsed). Day is zero-padded
rename_logic <- function(sample) {
  p <- parse_sample_name(sample)
  case_when(
    p$cell_type %in% params$cell_types_with_day ~ sprintf("%s_day%02d", p$cell_type, p$day),
    TRUE ~ "PLTs" )}
new_colnames <- make.unique(map_chr(colnames(zscores_df_prot), ~ rename_logic(.)))  # rename_logic()
# Assign the new column names to your dataframe
colnames(zscores_df_prot) <- new_colnames
wgcna_impi<- wgcna_correlation[,1:4]
tz <- t(zscores_df_prot)
zscoresMeans1 <- data.frame(matrix(ncol = 0, nrow = ncol(zscores_df_prot)))
rownames(zscoresMeans1) <- colnames(zscores_df_prot)  # one row per sample
   ##Automatize it!
modules <- unique(wgcna_correlation$ME)
for (module in modules) {
  module_columns <- which(grepl(module, wgcna_correlation$ME))
  zscoresMeans <- rowMeans(tz[, module_columns])
 zscoresMeans1[[paste0("ME", module)]] <- zscoresMeans}
zscoresMeans2 <- t(zscoresMeans1)
colnames(zscoresMeans2) <- new_colnames
######################################################
## Rcistarget  ##
######################################################





######################################################
## Figure 2  ##
######################################################

##Figure 2A

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
##Figure 3B
#per celltype comparison
barplot_sig <- df_sigresults_day_comparison %>% 
  gather("Significance_test", "Significant", -UID, -Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>% 
  filter(Significant != 0) %>%
  group_by(Significance_test) %>%
  distinct(UID, Significant) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::mutate(Significance_test = factor(Significance_test, levels = c("MPB_day5 - CB_day5","MPB_day8 - CB_day8","MPB_day11 - CB_day11","MPB_day14 - CB_day14"))) %>%  # Order bars
  ggplot(aes(x = Significance_test, y = count)) +
  geom_bar(stat = "identity", fill = "gray56", width = 0.7) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  scale_y_continuous(limits = c(0,1500)) +
  labs(
    title = "Significance Plot",
    x = "",
    y = "Protein significant count"
  )
print(barplot_sig)
ggsave( "barplot_significance_tests_per_cell_type.pdf", barplot_sig, width = 14, height = 12, dpi = 300)
#per cb
barplot_sig <- df_sigresults_cb %>% 
  gather("Significance_test", "Significant", -UID, -Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>% 
  filter(Significant != 0) %>%
  group_by(Significance_test) %>%
  distinct(UID, Significant) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::mutate(Significance_test = factor(Significance_test, levels = c("CB_day8 - CB_day5","CB_day11 - CB_day5","CB_day14 - CB_day5"))) %>%  # Order bars
  ggplot(aes(x = Significance_test, y = count)) +
  geom_bar(stat = "identity", fill = "gray56", width = 0.7) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  scale_y_continuous(limits = c(0,1500)) +
  labs(
    title = "Significance Plot",
    x = "",
    y = "Protein significant count"
  )
print(barplot_sig)
ggsave( "barplot_significance_tests_per_cb.pdf", barplot_sig, width = 14, height = 12, dpi = 300)
#per mpb
barplot_sig <- df_sigresults_mpb %>% 
  gather("Significance_test", "Significant", -UID, -Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>% 
  filter(Significant != 0) %>%
  group_by(Significance_test) %>%
  distinct(UID, Significant) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::mutate(Significance_test = factor(Significance_test, levels = c("MPB_day8 - MPB_day5","MPB_day11 - MPB_day5","MPB_day14 - MPB_day5"))) %>%  # Order bars
  ggplot(aes(x = Significance_test, y = count)) +
  geom_bar(stat = "identity", fill = "gray56", width = 0.7) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  scale_y_continuous(limits = c(0,1500)) +
  labs(
    title = "Significance Plot",
    x = "",
    y = "Protein significant count"
  )
print(barplot_sig)
ggsave( "barplot_significance_tests_per_mpb.pdf", barplot_sig, width = 14, height = 12, dpi = 300)

##Figure 3C
contrast_name <- colnames(contrast_matrix)[1]   #day5
p_cut   <- 0.05
lfc_cut <- 1
n_label <- 15
highlight_proteins <- NULL                         
col_sig <- c("Up in MPB" = "#C0392B", "Up in CB" = "#2471A3", NS = "grey75")
vdf <- limma::topTable(fit_ebayes, coef = contrast_name, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("UID") %>%
  dplyr::left_join(df_uid_information2, by = "UID") %>%
  dplyr::mutate(
    Gene = sub(";.*", "", dplyr::coalesce(Genes, UID)),
    sig  = dplyr::case_when(
      adj.P.Val < p_cut & logFC >=  lfc_cut ~ "Up in MPB",
      adj.P.Val < p_cut & logFC <= -lfc_cut ~ "Up in CB",
      TRUE ~ "NS"),
    sig = factor(sig, levels = c("Up in MPB", "Up in CB", "NS")))
n_up   <- sum(vdf$sig == "Up in MPB")
n_down <- sum(vdf$sig == "Up in CB")
labels_df <- dplyr::bind_rows(
  vdf %>% dplyr::filter(sig == "Up in MPB")   %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(sig == "Up in CB") %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(Gene %in% highlight_proteins)
) %>% dplyr::distinct(UID, .keep_all = TRUE)
x_lim <- ceiling(max(abs(vdf$logFC), na.rm = TRUE))
y_max <- max(-log10(vdf$adj.P.Val), na.rm = TRUE)
# Visualization
volcano <- ggplot(vdf, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = ~ dplyr::filter(.x, sig == "NS"),
             colour = col_sig["NS"], size = 1.1, alpha = 0.6, stroke = 0) +
  geom_point(data = ~ dplyr::filter(.x, sig != "NS"),
             aes(fill = sig), shape = 21, colour = "white", stroke = 0.25, size = 2.2) +
  ggrepel::geom_text_repel(
    data = labels_df, aes(label = Gene, colour = sig),
    size = 2.8, fontface = "italic",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "grey40",
    box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf, seed = 1,
    show.legend = FALSE
  ) +
  annotate("text", x = x_lim, y = y_max, label = paste0(n_up, " up in MPB"),
           hjust = 1, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in MPB"]) +
  annotate("text", x = -x_lim, y = y_max, label = paste0(n_down, " up in CB"),
           hjust = 0, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in CB"]) +
  scale_fill_manual(values = col_sig, guide = "none") +
  scale_colour_manual(values = col_sig, guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::breaks_pretty(n = 7)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"(adjusted "*italic(p)*")"),
    title = contrast_name
  ) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9, colour = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    plot.margin  = margin(6, 8, 6, 6)
  )
# Exporting
file_stub <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
ggsave(paste0("Volcano_", file_stub, ".pdf"), volcano,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
##Figure 3D
contrast_name <- colnames(contrast_matrix)[2]   #day8
p_cut   <- 0.05
lfc_cut <- 1
n_label <- 15
highlight_proteins <- NULL                         
col_sig <- c("Up in MPB" = "#C0392B", "Up in CB" = "#2471A3", NS = "grey75")
vdf <- limma::topTable(fit_ebayes, coef = contrast_name, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("UID") %>%
  dplyr::left_join(df_uid_information2, by = "UID") %>%
  dplyr::mutate(
    Gene = sub(";.*", "", dplyr::coalesce(Genes, UID)),
    sig  = dplyr::case_when(
      adj.P.Val < p_cut & logFC >=  lfc_cut ~ "Up in MPB",
      adj.P.Val < p_cut & logFC <= -lfc_cut ~ "Up in CB",
      TRUE ~ "NS"),
    sig = factor(sig, levels = c("Up in MPB", "Up in CB", "NS")))
n_up   <- sum(vdf$sig == "Up in MPB")
n_down <- sum(vdf$sig == "Up in CB")
labels_df <- dplyr::bind_rows(
  vdf %>% dplyr::filter(sig == "Up in MPB")   %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(sig == "Up in CB") %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(Gene %in% highlight_proteins)
) %>% dplyr::distinct(UID, .keep_all = TRUE)
x_lim <- ceiling(max(abs(vdf$logFC), na.rm = TRUE))
y_max <- max(-log10(vdf$adj.P.Val), na.rm = TRUE)
# Visualization
volcano <- ggplot(vdf, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = ~ dplyr::filter(.x, sig == "NS"),
             colour = col_sig["NS"], size = 1.1, alpha = 0.6, stroke = 0) +
  geom_point(data = ~ dplyr::filter(.x, sig != "NS"),
             aes(fill = sig), shape = 21, colour = "white", stroke = 0.25, size = 2.2) +
  ggrepel::geom_text_repel(
    data = labels_df, aes(label = Gene, colour = sig),
    size = 2.8, fontface = "italic",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "grey40",
    box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf, seed = 1,
    show.legend = FALSE
  ) +
  annotate("text", x = x_lim, y = y_max, label = paste0(n_up, " up in MPB"),
           hjust = 1, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in MPB"]) +
  annotate("text", x = -x_lim, y = y_max, label = paste0(n_down, " up in CB"),
           hjust = 0, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in CB"]) +
  scale_fill_manual(values = col_sig, guide = "none") +
  scale_colour_manual(values = col_sig, guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::breaks_pretty(n = 7)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"(adjusted "*italic(p)*")"),
    title = contrast_name
  ) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9, colour = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    plot.margin  = margin(6, 8, 6, 6)
  )
# Exporting
file_stub <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
ggsave(paste0("Volcano_", file_stub, ".pdf"), volcano,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
##Figure 3E
contrast_name <- colnames(contrast_matrix)[3]   #day11
p_cut   <- 0.05
lfc_cut <- 1
n_label <- 15
highlight_proteins <- NULL                         
col_sig <- c("Up in MPB" = "#C0392B", "Up in CB" = "#2471A3", NS = "grey75")
vdf <- limma::topTable(fit_ebayes, coef = contrast_name, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("UID") %>%
  dplyr::left_join(df_uid_information2, by = "UID") %>%
  dplyr::mutate(
    Gene = sub(";.*", "", dplyr::coalesce(Genes, UID)),
    sig  = dplyr::case_when(
      adj.P.Val < p_cut & logFC >=  lfc_cut ~ "Up in MPB",
      adj.P.Val < p_cut & logFC <= -lfc_cut ~ "Up in CB",
      TRUE ~ "NS"),
    sig = factor(sig, levels = c("Up in MPB", "Up in CB", "NS")))
n_up   <- sum(vdf$sig == "Up in MPB")
n_down <- sum(vdf$sig == "Up in CB")
labels_df <- dplyr::bind_rows(
  vdf %>% dplyr::filter(sig == "Up in MPB")   %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(sig == "Up in CB") %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(Gene %in% highlight_proteins)
) %>% dplyr::distinct(UID, .keep_all = TRUE)
x_lim <- ceiling(max(abs(vdf$logFC), na.rm = TRUE))
y_max <- max(-log10(vdf$adj.P.Val), na.rm = TRUE)
# Visualization
volcano <- ggplot(vdf, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = ~ dplyr::filter(.x, sig == "NS"),
             colour = col_sig["NS"], size = 1.1, alpha = 0.6, stroke = 0) +
  geom_point(data = ~ dplyr::filter(.x, sig != "NS"),
             aes(fill = sig), shape = 21, colour = "white", stroke = 0.25, size = 2.2) +
  ggrepel::geom_text_repel(
    data = labels_df, aes(label = Gene, colour = sig),
    size = 2.8, fontface = "italic",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "grey40",
    box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf, seed = 1,
    show.legend = FALSE
  ) +
  annotate("text", x = x_lim, y = y_max, label = paste0(n_up, " up in MPB"),
           hjust = 1, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in MPB"]) +
  annotate("text", x = -x_lim, y = y_max, label = paste0(n_down, " up in CB"),
           hjust = 0, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in CB"]) +
  scale_fill_manual(values = col_sig, guide = "none") +
  scale_colour_manual(values = col_sig, guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::breaks_pretty(n = 7)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"(adjusted "*italic(p)*")"),
    title = contrast_name
  ) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9, colour = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    plot.margin  = margin(6, 8, 6, 6)
  )
# Exporting
file_stub <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
ggsave(paste0("Volcano_", file_stub, ".pdf"), volcano,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
##Figure 3F 
contrast_name <- colnames(contrast_matrix)[4]   #day14
p_cut   <- 0.05
lfc_cut <- 1
n_label <- 15
highlight_proteins <- NULL                         
col_sig <- c("Up in MPB" = "#C0392B", "Up in CB" = "#2471A3", NS = "grey75")
vdf <- limma::topTable(fit_ebayes, coef = contrast_name, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("UID") %>%
  dplyr::left_join(df_uid_information2, by = "UID") %>%
  dplyr::mutate(
    Gene = sub(";.*", "", dplyr::coalesce(Genes, UID)),
    sig  = dplyr::case_when(
      adj.P.Val < p_cut & logFC >=  lfc_cut ~ "Up in MPB",
      adj.P.Val < p_cut & logFC <= -lfc_cut ~ "Up in CB",
      TRUE ~ "NS"),
    sig = factor(sig, levels = c("Up in MPB", "Up in CB", "NS")))
n_up   <- sum(vdf$sig == "Up in MPB")
n_down <- sum(vdf$sig == "Up in CB")
labels_df <- dplyr::bind_rows(
  vdf %>% dplyr::filter(sig == "Up in MPB")   %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(sig == "Up in CB") %>% dplyr::slice_min(adj.P.Val, n = n_label, with_ties = FALSE),
  vdf %>% dplyr::filter(Gene %in% highlight_proteins)
) %>% dplyr::distinct(UID, .keep_all = TRUE)
x_lim <- ceiling(max(abs(vdf$logFC), na.rm = TRUE))
y_max <- max(-log10(vdf$adj.P.Val), na.rm = TRUE)
# Visualization
volcano <- ggplot(vdf, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = ~ dplyr::filter(.x, sig == "NS"),
             colour = col_sig["NS"], size = 1.1, alpha = 0.6, stroke = 0) +
  geom_point(data = ~ dplyr::filter(.x, sig != "NS"),
             aes(fill = sig), shape = 21, colour = "white", stroke = 0.25, size = 2.2) +
  ggrepel::geom_text_repel(
    data = labels_df, aes(label = Gene, colour = sig),
    size = 2.8, fontface = "italic",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "grey40",
    box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf, seed = 1,
    show.legend = FALSE
  ) +
  annotate("text", x = x_lim, y = y_max, label = paste0(n_up, " up in MPB"),
           hjust = 1, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in MPB"]) +
  annotate("text", x = -x_lim, y = y_max, label = paste0(n_down, " up in CB"),
           hjust = 0, vjust = 1, size = 3.5, fontface = "bold", colour = col_sig["Up in CB"]) +
  scale_fill_manual(values = col_sig, guide = "none") +
  scale_colour_manual(values = col_sig, guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::breaks_pretty(n = 7)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"(adjusted "*italic(p)*")"),
    title = contrast_name
  ) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9, colour = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    plot.margin  = margin(6, 8, 6, 6)
  )
# Exporting
file_stub <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
ggsave(paste0("Volcano_", file_stub, ".pdf"), volcano,
       width = 89, height = 89, units = "mm", device = cairo_pdf)


######################################################
## Figure 4  ##
######################################################
##Figure 4A
column_order <- order(colnames(zscoresMeans2))
column_split_df <- data.frame(
  names= colnames(zscoresMeans2)
)

df_combined <- column_split_df %>%
  dplyr::mutate(Split = gsub("^([A-Za-z]+)[0-9]+(_day[0-9]+)$", "\\1\\2", names))
colnames(zscoresMeans2)
# Create the heatmap
# order
desired_split_order <- c("CB_day05", "CB_day08", "CB_day11", "CB_day14","MPB_day05", "MPB_day08", "MPB_day11", "MPB_day14")
column_split <- factor(df_combined$Split, levels = desired_split_order)
col_fun_zscoresMeans2 <- colorRamp2(
  c(min(zscoresMeans2, na.rm = TRUE), 0, max(zscoresMeans2, na.rm = TRUE)),
  c("blue", "white", "red")
)

ht1 <- Heatmap(
  (zscoresMeans2),
  name = "z-score",
 col = col_fun_zscoresMeans2,
  column_order = column_order,
 cluster_column_slices = TRUE,
  show_row_names = TRUE,
  show_column_names = F,
  show_column_dend = FALSE,
  show_row_dend = TRUE,
 row_dend_width = unit(3,"cm"),
 # guards against fewer rows (WGCNA modules) than requested clusters -- the module count
 # from dynamicTreeCut isn't fully deterministic run-to-run, so this isn't safe to hardcode
 row_km = min(5, nrow(zscoresMeans2)),
 column_split = column_split,
  column_names_side = 'bottom',
  heatmap_legend_param = list(
    title = "z-score",
    legend_direction = "vertical",
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    legend_width = unit(2, "cm")
  )
)
# Save the heatmap to a PDF file
pdf(file.path(output_dir, "heatmap_zscores_no_pearson.pdf"), width = 8, height = 5)
draw(ht1, column_title = "")
dev.off()

##Figure 4B
df_uid_information <- df_uid_information
enztrezid_gene <- clusterProfiler::bitr(geneID =df_uid_information$Protein.Group, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>%
  magrittr::set_colnames(c("Protein.Group", "ENTREZID"))
#clusters
wgcna_correlation$UID <- wgcna_correlation$Protein.Group
cluster_prepare <- wgcna_correlation %>%
  dplyr::select(UID, ME) %>%
  left_join(df_uid_information) %>%
  left_join(enztrezid_gene, "Protein.Group")
# Convert to data.table
dt_cluster_prepare <- as.data.table(cluster_prepare)
dt_cluster_prepare <- dt_cluster_prepare[!is.na(ENTREZID), .(ME, ENTREZID)]
dt_cluster_prepare[, ME := paste("ME", ME, sep = "")]
dt_cluster_prepare[, row_id := 1:.N, by = ME]
dt_spread <- data.table::dcast(dt_cluster_prepare, row_id ~ ME, value.var = "ENTREZID") 
cluster_annotation <- as.list(dt_spread[, -"row_id", with = FALSE])
cluster_annotation <- lapply(cluster_annotation, function(x) x[!is.na(x)])
#BP
clusters <- clusterProfiler::compareCluster(geneClusters = cluster_annotation, 
                             fun = "enrichGO", ont = "BP",
                              universe = unique(cluster_prepare$ENTREZID), 
                              OrgDb = org.Hs.eg.db,
                              pvalueCutoff = 0.05, 
                               readable = T
                                )
#CC
clusters2 <- clusterProfiler::compareCluster(geneClusters = cluster_annotation, 
                             fun = "enrichGO", ont = "CC",
                              universe = unique(cluster_prepare$ENTREZID), 
                              OrgDb = org.Hs.eg.db,
                              pvalueCutoff = 0.05, 
                               readable = T
                                )
#MF
clusters3 <- clusterProfiler::compareCluster(geneClusters = cluster_annotation,
                             fun = "enrichGO", ont = "MF",
                              universe = unique(cluster_prepare$ENTREZID),
                              OrgDb = org.Hs.eg.db,
                              pvalueCutoff = 0.05,
                               readable = T
                                )

df_bp <- as.data.frame(clusters@compareClusterResult)
df_cc <- as.data.frame(clusters2@compareClusterResult)
df_mf <- as.data.frame(clusters3@compareClusterResult)
df_bp$GOtype <- "BP"
df_cc$GOtype <- "CC"
df_mf$GOtype <- "MF"
df_all <- bind_rows(df_bp, df_cc, df_mf)
 write.csv(df_all, file="GO_enrichemnt_modules_for_cytoscape.csv")     # Export to table for cytoscape
##Figure 4C
##Figure 4D

# Data
df_bar <- data.frame(
  Group = c("a", "b", "c", "d", "e"),
  Numerator = c(7, 17, 11, 14, 12),
  Denominator = c(20, 45, 31, 35, 47)
) %>%
  mutate(Percentage = Numerator / Denominator * 100)
#  colors
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00")
# Plot
TF_barplot <- ggplot(df_bar, aes(x = Group, y = Percentage, fill = Group)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(
    aes(label = paste0(Numerator, "/", Denominator, 
                       " (", round(Percentage, 1), "%)")),
    vjust = -0.3, angle = 90, size = 4
  ) +
  scale_y_continuous(
  limits = c(0, 60),
  breaks = seq(0, 60, 15),  
  expand = c(0, 0)
) +
  scale_fill_manual(values = my_cols) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14)
  ) +
  labs(y = "Percentage (%)")
# Export
ggsave( filename = "TFS_barplot.pdf",  plot = TF_barplot, device = cairo_pdf, width = 8,height = 6,  units = "in")

######################################################
## Thom & Stokhuijzen et al dataset inclusion ##
######################################################


######################################################
## Figure 5  ##
######################################################
##Figure 5A
##Figure 5B
##Figure 5C





