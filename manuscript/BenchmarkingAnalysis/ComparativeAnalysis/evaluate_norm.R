## evaluation using ARI, silhouette score on cell types and batches, LISI on batches

library(data.table)
library(gridExtra)
library(mclust)
library(cluster)
library(stringr)
library(harmony)
library(Seurat)
library(pdfCluster)
library(dplyr)
library(umap)
# library(lisi)

mlgFindNeighbors = function(latent){
  rownames(latent) = paste0("cell", 1:nrow(latent))
  NN_graph = FindNeighbors(latent, verbose = FALSE, nn.method = "rann") ## "annoy"
  G = NN_graph$snn
  return(G)
}
mlgARI = function(embed, cluster_label){ ## louvain clustering using resolution 0.8
    cluster = cluster_label
    G = mlgFindNeighbors(embed)
    mlgc = FindClusters(G, resolution = 0.8, random.seed = 20220525, verbose = FALSE)[[1]] #0.8
    mlgARI = adjustedRandIndex(mlgc, cluster)
    return(mlgARI)
}
# evaluationLISI = function(embed, summary){ ## run it using conda environment
#     if ("batch" %in% colnames(summary)) {
#         return(compute_lisi(embed, summary, 'batch')$batch)
#     } else {
#         return(rep(NA, nrow(summary)))
#     }
# }


## latent can be PCA
## latent_type = "PCA" or "UMAP"
## evaluation_method = "ARI", "Si"
## label_type = "cell_type" or "batch"
## cluster_label = c(......) cell types of batches of each cell
evaluation_norm = function(latent = NULL, latent_type = NULL, evaluation_method = NULL, cluster_label = NULL){
    if(evaluation_method == "ARI"){
        return(mlgARI(latent, cluster_label))
    }

    if(evaluation_method == "Si"){
        if(latent_type %in% c("UMAP", "TSNE")){
            clusterNum = factor(cluster_label, levels = unique(cluster_label), label = 1:length(unique(cluster_label))) %>% as.numeric
            score = silhouette(clusterNum, dist(as.matrix(latent), method = "manhattan"))[,3] #
            return(score)
        }else{
            if(latent_type == "PCA"){
                embed = umap(latent)$layout
                clusterNum = factor(cluster_label, levels = unique(cluster_label), label = 1:length(unique(cluster_label))) %>% as.numeric
                score = silhouette(clusterNum, dist(as.matrix(embed)))[,3]
                return(score)
            }
        }
        
    }

    # if(evaluation_method == "LISI"){
    #     return(evaluationLISI(latent, data.frame(batch = cluster_label)))
    # }
}


args = commandArgs(trailingOnly=TRUE)

method = args[1]
evaluation_type = args[2]
study_name = args[3] #"public13Dataset_CITEseq"
out_path = "./manuscript/results/public13Dataset_CITEseq/" 

# adt_data = readRDS(file = paste0(out_path, "/RDS/adt_data_RawCount_", study_name, ".rds"))
if(method %in% c("totalVI_sample_GPU", "totalVI_study_GPU", "sciPENN_sample_GPU", "sciPENN_study_GPU")){
    adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", method, "_", study_name, ".rds"))

}else{
    adt_feature = readRDS(file = paste0(out_path, "/RDS/adt_feature_", study_name, ".rds"))
}

tmp = readRDS( file = paste0(out_path, "/RDS/adt_", method, "_", study_name, "_pca_umap.rds"))
ans_pca = tmp[[1]]
ans_umap = tmp[[2]]
dim(ans_pca)
dim(ans_umap)

seed_list = c(
    20220524, 20230817, 20200502, 20230801, 20190826, 20220525, 20230818, 20200503, 20230829, 20190827, 
    20137910, 20349283, 24918503, 38298493, 150, 2940, 10, 302934245, 3428, 9123784
)

# subsample_index = sample(1:nrow(ans_umap), 50000, replace = FALSE)
generate_index = function(seed){
    set.seed(seed)
    subsample_index = c()
    for(study_each in unique(adt_feature$study_name)){
        study_index = which(adt_feature$study_name == study_each)
        if(length(study_index) < 5000){
            subsample_index = c(subsample_index, study_index)
        }else(
            subsample_index = c(subsample_index, sample(study_index, 5000, replace = FALSE))
        )
    }
    print(length(subsample_index))

    return(subsample_index)
}


## Silhouette score
if(evaluation_type == "Si"){
    print("Silhouette score...")
    si_score_sample = c()
    si_score_study = c()
    si_score_broadCT = c()
    si_score_refineCT = c()

    for (seed in seed_list){
        subsample_index = generate_index(seed)
        # latent = ans_pca[subsample_index, 1:min(ncol(ans_pca), 50)] #ans_umap[subsample_index, ]
        latent = ans_umap[subsample_index, ] #ans_umap[subsample_index, ]
    
        cluster_label = adt_feature$sample[subsample_index] %>% as.character
        si_score_sample = evaluation_norm(latent, "UMAP", "Si", cluster_label) %>% rbind(si_score_sample, .)
        
        cluster_label = adt_feature$study_name[subsample_index] %>% as.character
        si_score_study = evaluation_norm(latent, "UMAP", "Si", cluster_label) %>% rbind(si_score_study, .)

        cluster_label = adt_feature$cell_type_l2[subsample_index] %>% as.character
        cell_select = which(cluster_label != "undefined")
        cluster_label = cluster_label[cell_select]
        si_score_refineCT = evaluation_norm(latent[cell_select, ], "UMAP", "Si", cluster_label) %>% rbind(si_score_refineCT, .)

        cluster_label = adt_feature$cell_type_l1[subsample_index] %>% as.character
        cell_select = which(cluster_label != "undefined")
        cluster_label = cluster_label[cell_select]
        si_score_broadCT = evaluation_norm(latent[cell_select, ], "UMAP", "Si", cluster_label) %>% rbind(si_score_broadCT, .)
        saveRDS(list(si_score_sample, si_score_study, si_score_refineCT, si_score_broadCT), paste0(out_path, "/RDS/adt_", method, "_", study_name, "_umap_si.rds"))

        
    }
    print(rowMeans(si_score_sample))
    print(rowMeans(si_score_study))
    print(rowMeans(si_score_broadCT))
    print(rowMeans(si_score_refineCT))

    # saveRDS(list(si_score_sample, si_score_study, si_score_refineCT, si_score_broadCT), paste0(out_path, "/RDS/adt_", method, "_", study_name, "_umap_si.rds"))

}

## ARI
if(evaluation_type == "ARI"){
    print("ARI score")
    
    ari_sample = c()
    ari_study = c()
    ari_refineCT = c()
    ari_broadCT = c()

    for(seed in seed_list){
        print(seed)
        subsample_index = generate_index(seed)
        latent = ans_pca[subsample_index, 1:min(ncol(ans_pca), 50)]

        cluster_label = adt_feature$sample[subsample_index] %>% as.character
        ari_sample = c(ari_sample, evaluation_norm(latent, "UMAP", "ARI", cluster_label))
        
        cluster_label = adt_feature$study_name[subsample_index] %>% as.character
        ari_study = c(ari_study, evaluation_norm(latent, "UMAP", "ARI", cluster_label))
        
        cluster_label = adt_feature$cell_type_l2[subsample_index] %>% as.character
        ari_refineCT = c(ari_refineCT, evaluation_norm(latent, "UMAP", "ARI", cluster_label))
        
        cluster_label = adt_feature$cell_type_l1[subsample_index] %>% as.character
        ari_broadCT = c(ari_broadCT, evaluation_norm(latent, "UMAP", "ARI", cluster_label))
        saveRDS(list(ari_sample, ari_study, ari_refineCT, ari_broadCT), paste0(out_path, "/RDS/adt_", method, "_", study_name, "_ari.rds"))
        
    }
    print(ari_sample)
    print(ari_study)
    print(ari_refineCT)
    print(ari_broadCT)

    saveRDS(list(ari_sample, ari_study, ari_refineCT, ari_broadCT), paste0(out_path, "/RDS/adt_", method, "_", study_name, "_ari.rds"))
}
