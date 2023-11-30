library(ggplot2)
library(monocle3)
library(Seurat)
library(Matrix)
library(anndata)

cds <- readRDS('/home/levinsj/Fetal_dir/AminData/iPT.new.cds.5.3.2023.rds')

print(cds)

pdf(file = "/home/levinsj/Fetal_dir/Analysis/rScripts/iPT_SPP_PDGFR.pdf", width = 8, heigh = 4)
plot_cells(cds, genes=c("SPP1","PDGFB"), show_trajectory_graph= FALSE, cell_size = 1, label_branch_points = FALSE, label_leaves = FALSE, label_groups_by_cluster = FALSE, label_roots = FALSE, graph_label_size = 0, labels_per_group = 0, label_principal_points = FALSE, )
dev.off()


pdf(file = "/home/levinsj/Fetal_dir/Analysis/rScripts/iPT_TGFB1_TGRBR1.pdf", width = 8, height = 4)
plot_cells(cds, genes=c("TGFB1","TGFBR1"), show_trajectory_graph= FALSE, cell_size = 1, label_branch_points = FALSE, label_leaves = FALSE,  label_groups_by_cluster = FALSE, label_roots = FALSE, graph_label_size = 0, labels_per_group = 0, label_principal_points = FALSE)
dev.off()

pdf(file = "/home/levinsj/Fetal_dir/Analysis/rScripts/iPT_CXCX13_13.pdf", width = 8, height = 4)
plot_cells(cds, genes=c("CXCL12","CXCL13"), show_trajectory_graph= FALSE, cell_size = 1, label_branch_points = FALSE, label_leaves = FALSE,  label_groups_by_cluster = FALSE, label_roots = FALSE, graph_label_size= 0, label_principal_points = FALSE, labels_per_group = 0)
dev.off()
