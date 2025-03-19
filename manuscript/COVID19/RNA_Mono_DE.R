rna_data = readRDS(file = paste0(out_path, "/RDS/rna_data_", run_name, ".rds")) 
rm_ind = readRDS(file = paste0(out_path, "/RDS/rm_ind1_ind2_", run_name, ".rds"))

sce = readH5AD("./scripts/COVID19/data/covid_portal_210320_with_raw.h5ad")

rownames(rna_data) = colnames(sce@assays@data@listData$X[, -rm_ind])
colnames(rna_data) = rownames(sce@assays@data@listData$X)[1:24737]

> which(colnames(rna_data) == "CD38")
[1] 5476
> which(colnames(rna_data) == "FCGR1A")
[1] 1327
> which(colnames(rna_data) == "SIGLEC1")
[1] 21457


rna_data[cell_hd, which(colnames(rna_data) == "CD38")] %>% summary
rna_data[cell_pt, which(colnames(rna_data) == "CD38")] %>% summary

rna_data[cell_hd, which(colnames(rna_data) == "FCGR1A")] %>% summary
rna_data[cell_pt, which(colnames(rna_data) == "FCGR1A")] %>% summary

rna_data[cell_hd, which(colnames(rna_data) == "SIGLEC1")] %>% summary
rna_data[cell_pt, which(colnames(rna_data) == "SIGLEC1")] %>% summary


rna_norm = sce@assays@data@listData$X[1:24737, -rm_ind] %>% t()

## for each mono cell type
rna_mono_summary = c()
for(cell_type_each in c("CD16_mono", "CD14_mono", "CD83_CD14_mono")){
    cell_hd = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary == "Healthy")
    cell_pt = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary %in% c("Asymptomatic", "Mild", "Moderate", "Severe", "Critical"))

    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_hd, which(colnames(rna_norm) == "CD38")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "CD38") %>% rbind(rna_mono_summary)

    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_pt, which(colnames(rna_norm) == "CD38")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "CD38") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_hd, which(colnames(rna_norm) == "FCGR1A")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "FCGR1A") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_pt, which(colnames(rna_norm) == "FCGR1A")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "FCGR1A") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_hd, which(colnames(rna_norm) == "SIGLEC1")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "SIGLEC1") %>% rbind(rna_mono_summary)

    rna_mono_summary = data.frame(
        gene_exp = rna_norm[cell_pt, which(colnames(rna_norm) == "SIGLEC1")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "SIGLEC1") %>% rbind(rna_mono_summary)

}
rna_mono_summary$status = factor(rna_mono_summary$status, levels = c("Patient", "Healthy"))

rna_mono_summary %>% ggplot(aes(x = gene, y = gene_exp, fill = status)) + 
geom_boxplot() + 
facet_wrap(~cell_type, scales = "free_y") + 
theme_bw(base_size = 25) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_brewer(palette = "Set1") +
xlab("") +
ylab("Normalized Gene Expression")

pdf("./results/COVID19/Figures/RNA_Mono_DE.pdf", width = 17, height = 8)
rna_mono_summary %>% dplyr::group_by(cell_type, status, gene) %>% dplyr::summarize(mean_gene_exp = mean(gene_exp), prop_pos_cell = sum(gene_exp > 0)/length(gene_exp) * 100) %>% ggplot(aes(x = gene, y = status, fill = mean_gene_exp, size = prop_pos_cell)) +
geom_point(shape = 21) +
facet_wrap(~cell_type) +
theme_bw(base_size = 25) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_gradient(low = "white", high = "red") +
xlab("") +
ylab("") +
scale_size_continuous(range = c(1, 40)) +
theme(legend.position = "top")
dev.off()



## for each mono cell type
rna_mono_summary = c()
for(cell_type_each in c("CD16_mono", "CD14_mono", "CD83_CD14_mono")){
    cell_hd = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary == "Healthy")
    cell_pt = which(adt_feature$full_clustering == cell_type_each & adt_feature$Status_on_day_collection_summary %in% c("Asymptomatic", "Mild", "Moderate", "Severe", "Critical"))

    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_hd, which(colnames(rna_data) == "CD38")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "CD38") %>% rbind(rna_mono_summary)

    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_pt, which(colnames(rna_data) == "CD38")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "CD38") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_hd, which(colnames(rna_data) == "FCGR1A")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "FCGR1A") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_pt, which(colnames(rna_data) == "FCGR1A")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "FCGR1A") %>% rbind(rna_mono_summary)
    
    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_hd, which(colnames(rna_data) == "SIGLEC1")],
        cell_type = cell_type_each,
        status = "Healthy",
        gene = "SIGLEC1") %>% rbind(rna_mono_summary)

    rna_mono_summary = data.frame(
        gene_exp = rna_data[cell_pt, which(colnames(rna_data) == "SIGLEC1")],
        cell_type = cell_type_each,
        status = "Patient",
        gene = "SIGLEC1") %>% rbind(rna_mono_summary)

}
rna_mono_summary$status = factor(rna_mono_summary$status, levels = c("Patient", "Healthy"))
rna_mono_summary %>% ggplot(aes(x = gene, y = gene_exp, fill = status)) + 
geom_boxplot() + 
facet_wrap(~cell_type, scales = "free_y") + 
theme_bw(base_size = 25) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_brewer(palette = "Set1") +
xlab("") +
ylab("Raw Gene Expression")

wilcox.test(rna_data[cell_hd, which(colnames(rna_data) == "CD38")],
rna_data[cell_pt, which(colnames(rna_data) == "CD38")])

rna_data[cell_hd, which(colnames(rna_data) == "FCGR1A")] %>% summary
rna_data[cell_pt, which(colnames(rna_data) == "FCGR1A")] %>% summary

rna_data[cell_hd, which(colnames(rna_data) == "SIGLEC1")] %>% summary
rna_data[cell_pt, which(colnames(rna_data) == "SIGLEC1")] %>% summary
