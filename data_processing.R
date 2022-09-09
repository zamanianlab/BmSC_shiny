install.packages("feather")

library(feather)
library(Seurat)


#-----------------------------------------------------------------------------
# new_combined <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BmSinglecell-ms/2_10XGenomics/scmulti_integrated.RDS")
# data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell
# 
# md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information
# 
# counts <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
#   pivot_longer(!gene_id, names_to = "index", values_to = "counts") 
# 
# 
# sc_data <- data %>% 
#   left_join(md) %>% 
#   left_join(counts) %>%
#   filter(counts > 0)
# 
# 
# sc_data <- sc_data %>% 
#   mutate(annotation = case_when(
#     integrated_snn_res.0.5 == "0" ~ "Unannotated",
#     integrated_snn_res.0.5 == "1" ~ "Body Wall Muscle",
#     integrated_snn_res.0.5 == "2" ~ "Unannotated",
#     integrated_snn_res.0.5 == "3" ~ "Unannotated",
#     integrated_snn_res.0.5 == "4" ~ "Unannotated",
#     integrated_snn_res.0.5 == "5" ~ "Coelomocyte",
#     integrated_snn_res.0.5 == "6" ~ "Unannotated",
#     integrated_snn_res.0.5 == "7" ~ "Unannotated",
#     integrated_snn_res.0.5 == "8" ~ "Mesoderm",
#     integrated_snn_res.0.5 == "9" ~ "Unannotated",
#     integrated_snn_res.0.5 == "10" ~ "Motor Neuron",
#     integrated_snn_res.0.5 == "11" ~ "Neuron",
#     integrated_snn_res.0.5 == "12" ~ "Neuron",
#     integrated_snn_res.0.5 == "13" ~ "Canal-associated",
#     integrated_snn_res.0.5 == "14" ~ "Secretory",
#     integrated_snn_res.0.5 == "15" ~ "Unannotated",
#     integrated_snn_res.0.5 == "16" ~ "Mesoderm",
#     integrated_snn_res.0.5 == "17" ~ "Neuron",
#     integrated_snn_res.0.5 == "18" ~ "Body Wall Muscle",
#     integrated_snn_res.0.5 == "19" ~ "Unannotated",
#     integrated_snn_res.0.5 == "20" ~ "Unannotated",
#     integrated_snn_res.0.5 == "21" ~ "Inner Body",
#     integrated_snn_res.0.5 == "22" ~ "Neuron",
#     integrated_snn_res.0.5 == "23" ~ "Neuron",
#     integrated_snn_res.0.5 == "24" ~ "Neuron",
#     integrated_snn_res.0.5 == "25" ~ "Neuron",
#     integrated_snn_res.0.5 == "26" ~ "Neuron"
#   )) %>% 
#   mutate(cluster = case_when(
#     integrated_snn_res.0.5 == "0" ~ 1,
#     integrated_snn_res.0.5 == "1" ~ 2,
#     integrated_snn_res.0.5 == "2" ~ 3,
#     integrated_snn_res.0.5 == "3" ~ 4,
#     integrated_snn_res.0.5 == "4" ~ 5,
#     integrated_snn_res.0.5 == "5" ~ 6,
#     integrated_snn_res.0.5 == "6" ~ 7,
#     integrated_snn_res.0.5 == "7" ~ 8,
#     integrated_snn_res.0.5 == "8" ~ 9,
#     integrated_snn_res.0.5 == "9" ~ 10,
#     integrated_snn_res.0.5 == "10" ~ 11,
#     integrated_snn_res.0.5 == "11" ~ 12,
#     integrated_snn_res.0.5 == "12" ~ 13,
#     integrated_snn_res.0.5 == "13" ~ 14,
#     integrated_snn_res.0.5 == "14" ~ 15,
#     integrated_snn_res.0.5 == "15" ~ 16,
#     integrated_snn_res.0.5 == "16" ~ 17,
#     integrated_snn_res.0.5 == "17" ~ 18,
#     integrated_snn_res.0.5 == "18" ~ 19,
#     integrated_snn_res.0.5 == "19" ~ 20,
#     integrated_snn_res.0.5 == "20" ~ 22,
#     integrated_snn_res.0.5 == "21" ~ 22,
#     integrated_snn_res.0.5 == "22" ~ 23,
#     integrated_snn_res.0.5 == "23" ~ 24,
#     integrated_snn_res.0.5 == "24" ~ 25,
#     integrated_snn_res.0.5 == "25" ~ 26,
#     integrated_snn_res.0.5 == "26" ~ 27
#   ))
# 
# sc_data <- sc_data %>% 
#   select(-"percent.mt", -"Percent.Largest.Gene", -"integrated_snn_res.0.5")
# col_order <- c("gene_id", "counts", "cluster", "annotation", "orig.ident", "index", "nCount_RNA", "nFeature_RNA", "UMAP_1", "UMAP_2")
#dataset <- dataset[,col_order]

#------------------------------------



gene_list <- read.csv("~/Desktop/gene_list.csv")

gene_list <- gene_list %>% 
  select(-"Genome.project", -"Transcript.stable.ID", -"Gene.description")

gene_list <- gene_list[!duplicated(gene_list),]

colnames(gene_list) <- c("gene_id", "gene_name")

combined <- gene_list %>% 
  left_join(dataset, by = "gene_id")

combined <- combined %>% 
  mutate(Treatment = case_when(
    orig.ident == "utBM" ~ "Untreated",
    orig.ident == "tBM" ~ "Ivermectin (1µM)"
  )) 




reduced <- combined %>% 
  select("gene_name", "gene_id", "counts", "cluster", "annotation", "Treatment")
colnames(reduced) <- c("Gene Name", "Gene ID", "Counts", "Cluster", "Annotation", "Treatment")
reduced$Counts <- round(reduced$Counts, digits = 3)


mapping <- combined %>% 
  select(-"nCount_RNA", -"nFeature_RNA", -"orig.ident")
colnames(mapping) <- c("Gene ID", "Gene Name", "Counts", "Cluster", "Annotation", "Index", "UMAP_1", "UMAP_2", "Treatment")
mapping$Counts <- round(mapping$Counts, digits = 3)


col_order <- c("Gene Name", "Gene ID", "Counts", "Cluster", "Annotation", "Treatment", "Index", "UMAP_1", "UMAP_2")
mapping<- mapping[,col_order]

#-----------------------------------------------------------
# Generate dataframe for dot plots
#calculate average gene expression per cluster using seurat's DotPlot function
list <- gene_list$gene_id
dot <- DotPlot(new_combined, features = list,  assay = "RNA", scale = FALSE)
dot <- dot$data
dot <- rownames_to_column(dot, "genes")
dot <- dot %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")

#rename clusters
dot$id<- factor(dot$id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

dot <- dot %>% 
  mutate(ID = case_when(
    id == "1" ~ "Unannotated",
    id == "2" ~ "MS",
    id == "3" ~ "Unannotated",
    id == "4" ~ "Unannotated",
    id == "5" ~ "Unannotated",
    id == "6" ~ "C",
    id == "7" ~ "Unannotated",
    id == "8" ~ "Unannotated",
    id == "9" ~ "MD",
    id == "10" ~ "Unannotated",
    id == "11" ~ "Neuron",
    id == "12" ~ "Neuron",
    id == "13" ~ "Neuron",
    id == "15" ~ "S",
    id == "14" ~ "CA",
    id == "16" ~ "Unannotated",
    id == "17" ~ "MD",
    id == "18" ~ "Neuron",
    id == "19" ~ "MS",
    id == "20" ~ "Unannotated",
    id == "21" ~ "Unannotated",
    id == "22" ~ "IB",
    id == "23" ~ "Neuron",
    id == "24" ~ "Neuron",
    id == "25" ~ "Neuron",
    id == "26" ~ "Neuron",
    id == "27" ~ "Neuron"))


dot$ID <- factor(dot$ID, levels = c("MS","MD", "C", "S", "CA", "IB", "Neuron", "Unannotated"))
dot <- dot %>% select(-"features.plot")

dot <- dot %>% 
  left_join(gene_list)





# export dataframes
write_feather(reduced, "~/Desktop/reduced.feather")
write_feather(mapping, "~/Desktop/mapping.feather")
write_feather(combined, "~/Desktop/complete.feather")
write_feather(dot, "~/Desktop/dot.feather")



