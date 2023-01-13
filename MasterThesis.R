


###################################

# Normalization 

###################################

bea.normalization <- function(dataset.name, annot, patient ,tumor,paper,stage,normal,dataset, path_save, cancer){
    library(dplyr) 
    library(dplyr) 
    library(SeuratObject)
    library(Seurat)
    library(patchwork)
    library(sctransform)

    if(dataset.name == "kim"){dataset[["Cluster"]] <- Idents(dataset)}
    dataset[["Patient"]] <- annot[,patient]
    if(tumor != 0){dataset[["Pri_Met"]] <- annot[, tumor]}
    dataset[["Paper"]] <- annot[,paper]
    dataset[["Stage"]] <- annot[,stage]


    # Remove cancer and normal cells 
    cancer.cell <- rownames(cancer[which(cancer$Epi_Cancer == "Cancer"),])
    #print(length(cancer.cell))

    normal.cell <- rownames(annot[which(annot[,normal] == "Normal"),])
    #print(length(normal.cell))
    dataset <- dataset[,!colnames(dataset) %in% normal.cell]
    dataset <- dataset[,!colnames(dataset) %in% cancer.cell]
    

    # Percentage mitochondrial genes and normalization
    dataset <- PercentageFeatureSet(dataset, pattern = "^MT-", col.name = "percent.mt")
    dataset <- SCTransform(dataset, vars.to.regress = "percent.mt", verbose = FALSE)

    # PCA and UMAP
    dataset <- RunPCA(dataset)
    dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)

    #Clusters
    dataset <- FindNeighbors(dataset, dims = 1:30)
    dataset <- FindClusters(dataset, resolution = 1.62)
    if(dataset.name != "kim"){dataset[["Cluster"]] <- Idents(dataset)}

    pdf(width=10,height=10,file= paste(path_save,
        paste0(dataset.name,"_original_UMAP.pdf"), sep = ""))

    print(DimPlot(dataset, reduction = "umap", label = TRUE, group.by = "Paper"))

    print(DimPlot(dataset, reduction = "umap", label = TRUE, group.by = "Cluster"))

    dev.off()


    return(dataset)

}

###################################

# Seurat Integration 

###################################

bea.integration <- function(dataset.name, dataset, path_save, path_ref, levels){
    library(Seurat)
    library(patchwork)

    # Reference
    hlca2 <- readRDS(path_ref)

    # Anchors 
    dataset.anchors <- FindTransferAnchors(reference = hlca2, query = dataset,
       dims = 1:30, reference.reduction = "pca")
    saveRDS(dataset.anchors, file = paste0(path_save,dataset.name, "_anchors.rds"))

    #dataset.anchors <- readRDS(paste0(path_save,paste0(dataset.name, "_anchors.rds")))  
    
    table <- matrix(,nrow = length(colnames(dataset)), ncol = 0)

    for(lev in levels){
        print(dim(table))
        table <- cbind(bea.annotation.level(lev, dataset, hlca2) , table)
    }

    # Table with all the annotation
    rownames(table) <- colnames(dataset)
    colnames(table) <- c("ann_finest_level", "ann_level_5","ann_level_4", "ann_level_3",
        "ann_level_2" , "ann_level_1") #columns are reversed 

    saveRDS(table , file= paste(path_save,
        paste(dataset.name, "HLCA_annot_table.rds", sep = "_"), sep = ""))

    #UMAP
    dataset[["HLCA_level2"]] <- table["ann_level_2"]
    
    pdf(width=10,height=10,file= paste(path_save,
        paste0(dataset.name,"_IntHLCA_UMAP.pdf"), sep = ""))

    print(DimPlot(dataset, reduction = "umap", label = TRUE, group.by = "HLCA_level2"))

    dev.off()

    pdf(width=8,height=8,file= paste(path_save,
        paste0(dataset.name,"_IntHLCA_UMAP_small.pdf"), sep = ""))

    print(DimPlot(dataset, reduction = "umap", label = TRUE, group.by = "HLCA_level2", repel = TRUE))

    dev.off()

    return(table)
}

bea.annotation.level <- function(level, dataset, hlca2){
    vec = as.character(hlca2@meta.data[[level]])
    names(vec) = rownames(hlca2@meta.data)

    predictions <- TransferData(anchorset = dataset.anchors, refdata = vec,dims = 1:30)
    dataset <- AddMetaData(dataset, metadata = predictions)
    return(as.vector(dataset[["predicted.id"]]))
    #saveRDS(dataset, file = paste(dataset.name, level, "annot.rds", sep = "_")) 

}


###################################

# Jaccard similarity 

###################################

bea.jaccard_matrix_1 <- function(paper_annot, table, paper, level){
	library(jaccard)

	######List of cell type
	paper_type <- unlist(unique(paper_annot[[paper]])) ######## Depends on the dataset 
	met_type <- unique(table[[level]])


	######Matrix 
	jaccard_mat <- matrix(,nrow = length(paper_type), ncol = length(met_type))
	colnames(jaccard_mat) <- met_type
	rownames(jaccard_mat) <- paper_type

	
	######Jaccard index 
	for(i in seq_along(paper_type)){
    	for (j in seq_along(met_type)) {
        	jaccard_mat[i,j] <- jaccard(c(paper_annot[[paper]] == paper_type[i]), 
        		c(table[[level]] == met_type[j]))
        }
    }


	######Ordering 
    or <- apply(jaccard_mat, 1, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[or,]

	or <- apply(jaccard_mat, 2, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[,or]

    return(jaccard_mat)

}

bea.jacc.htmap <- function(path_save, dataset.name, jaccard_list, met, title){
	library(gplots)

	pdf(width=10,height=10,file= paste(path_save,
    	paste0(dataset.name,"_jacc_htmap_",title,".pdf"), sep = ""))

	my_palette <- colorRampPalette(c("white","aquamarine2" ,"darkgreen"))

	for (t in seq_along(jaccard_list)){

		heatmap.2(x=jaccard_list[[t]], 
		Colv=FALSE, 
		Rowv=FALSE,
		scale="none",
		col=my_palette,
		trace="none",
		main=paste("Jaccard", dataset.name, met[t]),
		ylab="Original Annotation",
		xlab= paste("Novel annotation",met[t]),
		margins=c(14,14))

	}

	dev.off()
}


###################################

# Consensus annotation

###################################

consensus <- list("Tcell" = c("CD4_MEMORY_EFFECTOR_T_CELL", "CD4_NAIVE_T_CELL", 
		"CD8_MEMORY_EFFECTOR_T_CELL","CD8_NAIVE_T_CELL", "NATURAL_KILLER_T_CELL",
		"PROLIFERATING_NK_T_CELL", "CD8.T.cell.exhausted","CD4.T.cell.exhausted"
		, "CD8.T.cell.naive" ,"CD8.T.cell.C2.CD28" , "CD8.T.cell.effector"
		,"CD8.T.cell.GZMK","CD8.T.cell.tissue.resident.memory","CD8.T.cell.MAIT"
		,"CD4.T.cell.naive", "CD4.T.cell.central.memory.blood","CD4.T.cell.effector"
		, "CD4.T.cell.CD69","CD4.T.cell.EOMES","CD4.T.cell.GZMA"
		,"CD4.T.cell.Treg.FOXP3.resting.blood", "CD4.T.cell.Treg.suppressive.CTLA4",
		"CD4.T.cells","CD8.T.cells","T.cells.proliferating", "Gamma.delta.T.cells",
		"Natural.killer.T.cells","T.cells","T.cytotoxic.cells","T.follicular.helper.cells",
		"T.helper.cells","T.memory.cells","T.regulatory.cells", "Conventional.T.cells",
		"Regulatory.T.cells","CD8..T.cells", "Naive CD8+ T cells","Naive CD4+ T cells",
		"Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells" ,
		"Effector CD4+ T cells" ,"γδ-T cells", "CD8+ NKT-like cells","CD4+ NKT-like cells",
		"CD4+ T-cells", "CD8+ T-cells", "T_cells", "CD4 T cells", "CD8 T cells"),

	"Bcell" = c("B.cells", "B.cells.memory", "B.cells.naive", "B_CELL", "Pro-B cells",
		"Pre-B cells","Naive B cells","Memory B cells", "B-cells", "B_cell", "Pre-B_cell_CD34-",
		"Pro-B_cell_CD34+", "B cells"),

	"NKcell" = c("NK.cells", "NATURAL_KILLER_CELL", "Natural killer  cells", "NK cells",
		"NK_cell"),

	"MASTcell" = c("Mast.cells", "BASOPHIL_MAST_1_CELL", "BASOPHIL_MAST_2_CELL", "Mast cells"),

	"Macrophage" = c("Alveolar.Macrophages", "Monocyte.derived.macrophages", "Alveolar.Mφ.CCL3.",
		"Alveolar.Mφ.MT.positive", "Alveolar.Mφ.proliferating", "Monocyte.derived.Mφ",
		"Interstitial.Mφ.perivascular", "Macrophages", "Red.pulp.macrophages","MACROPHAGE_CELL", 
		"PROLIFERATING_MACROPHAGE_CELL", "Alveolar macrophages", "Macrophage", "Monocyte-derived Mφ",
		"Alveolar.macrophages"),

	"Monocyte" = c("Monocytes","CLASSICAL_MONOCYTE_CELL", "NONCLASSICAL_MONOCYTE_CELL", 
		"OLR1_CLASSICAL_MONOCYTE_CELL", "Classical Monocytes", "Non-classical monocytes",
		"Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Monocyte", 
		"Classical.monocytes", "Non.classical.monocytes", "Classical monocytes"),

	"Dcell" = c("Myeloid.dendritic.cells", "Migratory.DCs", "DC1", "DC2", "Dendritic.cells","EREG_DENDRITIC_CELL", 
		"IGSF21_DENDRITIC_CELL", "MYELOID_DENDRITIC_TYPE_1_CELL", "TREM2_DENDRITIC_CELL",
		"Myeloid Dendritic cells", "Myeloid Dendritic cells", "DC"),

	"PlasmacytoidDcell" =c("Plasmacytoid.dendritic.cells", "Plasmacytoid.DCs","PLASMACYTOID_DENDRITIC_CELL",
		"Plasmacytoid Dendritic cells", "Plasmacytoid DCs"),

	"PLASMAcell" = c("Plasma.cells", "Plasma B cells", "Plasma cell", "Plasma cells"),

	"Megakaryocyte" = c("Megakaryocytes", "PLATELET_MEGAKARYOCYTE_CELL", "Platelets", "Megakaryocyte"),

	"MDScell" = c("Myeloid.derived.suppressor.cells"),

	"Endothelial" = c("Endothelial.cells", "Lymphatic.endothelial.cells", "EC.aerocyte.capillary",
		"EC.venous.systemic", "Lymphatic.EC.mature", "Lymphatic.EC.differentiating",
		"EC.arterial", "EC.general.capillary", "EC.venous.pulmonary", "Lymphatic.EC.proliferating",
		"Endothelial", "Erythroid-like and erythroid precursor cells", "Endothelial cell",
		"Endothelial cells", "Endothelial_cells", "CAPILLARY_CELL", "EC venous systemic", 
		"EC general capillary"),

	"Epithelial"  = c("Epithelial", "Alveolar.type.2", "Alveolar.type.1", "Basal.cells",
		"Club.cells", "Ciliated.cells", "Basal.resting", "Deuterosomal", "Goblet..subsegmental.",
		"Ionocyte", "SMG.mucous", "AT1", "AT2", "AT2.proliferating", "Suprabasal",
		"Multiciliated..non.nasal.", "Goblet..bronchial.", "Transitional.Club.AT2", "Tuft",
		"SMG.serous..bronchial.", "SMG.duct", "Airway.epithelial.cells", "Clara.cells",
		"Ionocytes", "Pulmonary.alveolar.type.I.cells", "Airway.goblet.cells", 
		"Pulmonary.alveolar.type.II.cells", "Airway epithelial cells", "Airway goblet cells",
		"Basal cells (Airway progenitor cells)", "Ciliated cells", "Clara cells", "Epithelial cells",
		"Pulmonary alveolar type I cells", "Pulmonary alveolar type II cells", "Secretory cell",
		"Epithelial cells", "Keratinocytes", "Melanocytes", "Epithelial_cells", "ALVEOLAR_EPITHELIAL_TYPE_1_CELL",
		"SIGNALING_ALVEOLAR_EPITHELIAL_TYPE_2_CELL", "GOBLET_CELL", "BASAL_CELL", "MUCOUS_CELL", "CILIATED_CELL",
		"CLUB_CELL", "Multiciliated (non-nasal)", "Basal resting", "Transitional Club-AT2", "Goblet (nasal)",
		"Club..non.nasal.","Mesothelial.cells", "Mesothelial cells", "PROXIMAL_CILIATED_CELL"),

	"Fibroblast" = c("Stromal", "Myofibroblasts", "Fibroblasts", "Smooth.muscle.cells", "Adventitial.fibroblasts",
		"Pericytes", "Fibromyocytes", "Mesothelium", "Peribronchial.fibroblasts",
		"Alveolar.fibroblasts", "Subpleural.fibroblasts", "Smooth.muscle", "SM.activated.stress.response", 
		"Myocytes", "Pericytes", "Smooth muscle", "Skeletal muscle", "Smooth_muscle_cells",
		"Chondrocytes", "VASCULAR_SMOOTH_MUSCLE_CELL", "ADVENTITIAL_FIBROBLAST_CELL", "Adventitial fibroblasts",
		"Alveolar fibroblasts", "ALVEOLAR_FIBROBLAST_CELL", "FIBROMYOCYTE_CELL"),

	"Neuron_Glia" = c("Neuroendocrine", "Neurons", "Astrocytes", "Astrocyte", "Neuroepithelial_cell", "NEUROENDOCRINE_CELL"),
	"Basophils" = c("Basophils"),
	"Eosinophils" = c("Eosinophils"),
	"Neutrophils" = c("Neutrophils", "Granulocytes", "Myelocyte", "Pro-Myelocyte"),
	"Nuocytes" = c("Nuocytes"),
	"RedBloobCell" = c("Erythroid-like and erythroid precursor cells", "Erythrocytes", "Erythroblast"),
	"HematopoieticStemCell" = c("HSC/MPP cells", "HSC_-G-CSF", "HSC_CD34", "HSC"), 
	"PrecusorCell" = c("Progenitor cells", "MSC", "MEP", "BM & Prog.", "Tissue_stem_cells", 
		"BM", "iPS_cells", "GMP", "CMP"), 
	"VariousCell" = c("Mesangial", "Adipocytes", "Gametocytes", "Embryonic_stem_cells", 
		"Osteoblasts", "Hepatocytes"),
	"Immune" = c("Immune system cells", "ISG expressing immune cells"))


bea.consensus.annot = function(dataset, path_save, pred_full, consensus){
	pred_full <- pred_full[ !pred_full$ref %in% c("level_1", "level_2", "level_3", "level_5", "level_4", "Guo"), ] 
	pred_full$met_ref <- paste0(pred_full$method, "_", pred_full$ref)

	cluster <- as.numeric(levels(Idents(dataset)))
	prediction <- data.frame("cluster" = cluster, "pred" = NA)

	for(cl in cluster){
		con <- c()
		annot = pred_full[pred_full[,"cluster"] == cl,"ct"]
		for (ct in annot) {
			for (i in 1:length(consensus)) {
				if(ct %in% consensus[[i]]){

					con = c(con, names(consensus[i]))
				}
			}
		}
		prediction[prediction[,"cluster"] == cl,"pred"] = names(which.max(table(con)))
	}

	return(prediction)
}


bea.consensus.umap = function(dataset, path_save, prediction){

	new.cluster.ids <- c()

	for (i in 1:dim(dataset)[2]) {
		new.cluster.ids <- c(new.cluster.ids, prediction[prediction$cluster == dataset@meta.data[i,"Cluster"],"pred"])
	}

	dataset[["Consensus"]] <- new.cluster.ids

	pdf(width=7,height=7,file= paste(path_save,
	    paste0(dataset.name,"_consensus_UMAP.pdf"), sep = ""))


	print(DimPlot(dataset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Consensus", repel = TRUE))

	#print(DimPlot(object = dataset, cells.highlight = rownames(tcell), cols.highlight = "red", cols = "gray"))
	
	dev.off()

	return(dataset)

}

###################################

# T cells refinement 

###################################

bea.tcell.id <- function(dataset, prediction){

	tcell.id <- c()
	for (i in 1:dim(dataset)[2]) {
		if(dataset@meta.data[i, "Consensus"] == "Tcell"){
			tcell.id <- c(tcell.id, colnames(dataset)[i])
		}
	}
	return(tcell.id)

}

bea.TCell <- function(dataset, dataset.name, tcell.id, path_save){ #Calls all the other functions 
	library(dplyr) 
	library(SeuratObject)
	library(Seurat)
	library(patchwork)
	library(ggplot2)
	library(sctransform)
	library(jaccard)
	library(gplots)
	library(RColorBrewer)
	library(UCell)
	library(ProjecTILs)
	library(stringr)

	#querydata <- readRDS(file= paste("/mnt/ndata/beatrice/Temp/",str_to_title(dataset.name),"/",
    #paste(dataset.name, "level2.rds", sep = "_"), sep = ""))


    querydata <- dataset[,colnames(dataset) %in% tcell.id]

    ref <- load.reference.map()

    #table <- readRDS(table , file= paste(path_save,
    	#paste(dataset.name, "HLCA_annot_table.rds", sep = "_"), sep = "")) 

	#PROJECTION 
	query.projected <- make.projection(querydata, ref=ref, filter.cells = FALSE) #Filter T cells 

	#PREDICTION 
	query.projected <- cellstate.predict(ref=ref, query=query.projected)
	print(table(query.projected$functional.cluster))

	tcell <- matrix(,nrow = length(colnames(query.projected)), ncol = 0)
	rownames(tcell) <- colnames(query.projected)
	tcell <- cbind(tcell, query.projected$functional.cluster)
	colnames(tcell) <- "TILs"



	pdf(width=7,height=7,file= paste(path_save,
		   paste0(dataset.name, "_TCell_filtered.pdf"), sep = ""))

	#REF PLOT 
	refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", 
		"#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
	print(DimPlot(ref,label = T, cols = refCols))

	#PROJECTION PLOT
	print(plot.projection(ref, query.projected, cols = refCols))

	#HISTOGRAM
	print(plot.statepred.composition(ref, query.projected,metric = "Percent") + theme_bw()+ theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)))


	dev.off()


	######List of cell type
	dataset_type <- unlist(unique(querydata[["Paper"]]))  
	til_type <- unique(tcell[,"TILs"])

	######Matrix 
	jaccard_mat <- matrix(,nrow = length(dataset_type), ncol = length(til_type))
	colnames(jaccard_mat) <- til_type
	rownames(jaccard_mat) <- dataset_type


	######Jaccard index 
	for(i in seq_along(dataset_type)){
    	for (j in seq_along(til_type)) {
        	jaccard_mat[i,j] <- jaccard(c(querydata[["Paper"]] == dataset_type[i]), 
        		tcell[,"TILs"] == til_type[j])
    		}}

    ######Ordering 
    or <- apply(jaccard_mat, 1, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[or,]

	or <- apply(jaccard_mat, 2, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[,or]

	pdf(width=10,height=10,file= paste(path_save,
    	paste0(dataset.name,"_jacc_htmap_tcell.pdf"), sep = ""))

	my_palette <- colorRampPalette(c("white","aquamarine2" ,"darkgreen"))

	heatmap.2(x=jaccard_mat, 
	Colv=FALSE, 
	Rowv=FALSE,
	scale="none",
	col=my_palette,
	trace="none",
	main=paste("Jaccard", dataset.name, "ProjecTILs"),
	ylab="Dataset Annotation",
	xlab= "ProjecTILs",
	margins=c(14,14))


	dev.off()

	saveRDS(tcell , file= paste(path_save,
		paste(dataset.name,"_tcell.rds", sep = "_"), sep = ""))

	return(tcell)

}



bea.TILs.incorporation <- function(dataset, tcell){

	new.cluster.ids <- dataset@meta.data[,"Consensus"]
	names(new.cluster.ids) <- colnames(dataset)
	for (i in seq_along(new.cluster.ids)) {
		if(new.cluster.ids[i] == "Tcell"){
			new.cluster.ids[i] <- tcell[names(new.cluster.ids)[i],"TILs"]
		}
	}
	dataset[["ConsensusTILs"]] <- new.cluster.ids

	pdf(width=7,height=7,file= paste(path_save,
	    paste0(dataset.name,"_consensusTILs_UMAP.pdf"), sep = ""))

	print(DimPlot(dataset, reduction = "umap", label = TRUE, pt.size = 0.5, 
		group.by = "ConsensusTILs", repel = TRUE)
		+ guides(color=guide_legend(override.aes = list(size=4),ncol =1)))

	dev.off()

	return(dataset)
}

###################################

# Jaccard similarity 

###################################

bea.methods_jaccard <- function(pred_full, dataset, prediction, paper, consensus){
	#load(path_pred) # named pred_full
	pred_full <- pred_full[ !(pred_full$ref %in% c("level_1", "level_2", "level_3", "level_5", "level_4", "Guo")), ] 
	pred_full$met_ref <- paste0(pred_full$method, "_", pred_full$ref)
	pred_full$common <- NA

	for(j in seq_along(pred_full$ct)){
		for (i in 1:length(consensus)) {
			if(pred_full$ct[j] %in% consensus[[i]]){
				pred_full$common[j] <- names(consensus[i])
			}
		}
	}

	common_annot <- data.frame(matrix(,nrow = dim(dataset[["seurat_clusters"]])[1], ncol = 3 + length(unique(pred_full$met_ref)), 
		dimnames = list(NULL, c("cluster", "paper", unique(pred_full$met_ref), "consensus"))))
	if(dataset.name != "kim"){common_annot$cluster <- unlist(dataset[["seurat_clusters"]])}
	if(dataset.name == "kim"){common_annot$cluster <- unlist(dataset[["Cluster"]])}
	common_annot$paper <- unlist(dataset[[paper]])
	rownames(common_annot) <- rownames(dataset[[]])

	for(met in unique(pred_full$met_ref)){
		for(k in seq_along(common_annot$cluster)){
			common_annot[k,met] <- pred_full[pred_full$met_ref == met & pred_full$cluster == common_annot$cluster[k],"common"]
		}
	}

	for(k in seq_along(common_annot$cluster)){
		common_annot[k,"consensus"] <- prediction[prediction$cluster == common_annot$cluster[k],"pred"]
	}

	if(dataset.name == "kim") {
		translation <- list("Myeloid cells" = c("Megakaryocyte", "Eosinophils", "Neutrophils",
				"Dcell", "Monocyte", "Macrophage", "Basophils", "PlasmacytoidDcell","MDScell"),
			"MAST cells" = c("MASTcell"),
			"T lymphocytes" = c("Tcell"),
			"NK cells" = c("NKcell"),
			"B lymphocytes" = c("Bcell", "PLASMAcell"),
			"Fibroblasts" = c("Fibroblast"),
			"Epithelial cells" = c("Epithelial"),
			"Endothelial cells" = c("Endothelial"),
			"Neuron_Glia" = c("Neuron_Glia"),
			"Other" = c("Immune", "PrecusorCell","VariousCell", "HematopoieticStemCell", 
				"RedBloobCell", "Nuocytes", "PrecusorCell"))
	}

	if(dataset.name == "qian") {
		translation <- list("Myeloid" = c("Megakaryocyte", "Eosinophils", "Neutrophils",
				"Dcell", "Monocyte", "Macrophage", "Basophils", "PlasmacytoidDcell","MDScell"),
			"MAST-cells" = c("MASTcell"),
			"T-cells" = c("Tcell"),
			"NK-cells" = c("NKcell"),
			"B-cells" = c("Bcell", "PLASMAcell"),
			"Fibroblasts" = c("Fibroblast"),
			"Epithelial" = c("Epithelial"),
			"Endothelial" = c("Endothelial"),
			"Unknown" = c("Immune", "PrecusorCell","VariousCell", "HematopoieticStemCell", 
				"RedBloobCell", "Nuocytes", "PrecusorCell"))
	}
	for(met in colnames(common_annot)[3:length(common_annot)]){
		for(j in seq_along(common_annot[,met])){
			for (i in 1:length(translation)) {
				if(common_annot[j,met] %in% translation[[i]]){
					common_annot[j,met] <- names(translation[i])
				}
			}
		}
	}
	return(common_annot)
}


bea.jaccard_matrix <- function(common_annot, met, Paper){
	

	######List of cell type
	paper_type <- unique(common_annot[, Paper]) ######## Depends on the dataset 
	met_type <- unique(common_annot[, met])
	

	######Matrix 
	jaccard_mat <- matrix(,nrow = length(paper_type), ncol = length(met_type))
	colnames(jaccard_mat) <- met_type
	rownames(jaccard_mat) <- paper_type


	######Jaccard index 
	for(i in seq_along(paper_type)){
    	for (j in seq_along(met_type)) {
        	jaccard_mat[i,j] <- jaccard(c(common_annot[, Paper] == paper_type[i]), 
        		common_annot[, met] == met_type[j])
    	}
   	}

    or <- apply(jaccard_mat, 1, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[or,]


	or <- apply(jaccard_mat, 2, max, na.rm = TRUE)
	or <- order(or, decreasing = TRUE)

	jaccard_mat <- jaccard_mat[,or]

    return(jaccard_mat)

}

###################################

# Comparaison of different methods 

###################################


bea.best_match <- function(jaccard_list, remove_ct, met){
	
	#paper_type <- rownames(jaccard_mat)[!rownames(jaccard_mat) %in% remove_ct]
	paper_type <- rownames(jaccard_list[[1]])[!rownames(jaccard_list[[1]]) %in% remove_ct]


	best_match <- data.frame(matrix(,nrow = length(paper_type), ncol = length(met)))
	colnames(best_match) <- met
	rownames(best_match) <- paper_type

	for (i in seq_along(jaccard_list)) {
		jaccard_mat = jaccard_list[[i]]

		for (ct in paper_type) {
			if(ct == "Oligodendrocytes"){
				if(ct %in% colnames(jaccard_mat)){
					best_match[ct,met[i]] <- jaccard_mat[ct,"Neuron_Glia"]
				} else {
					best_match[ct,met[i]] <- 0
				}
			} else {
				if(ct %in% colnames(jaccard_mat)){
					best_match[ct,met[i]] <- jaccard_mat[ct,ct]
				}
				else {
					best_match[ct,met[i]] <- 0
				}
			}
		}
	}


	return(best_match)

}



bea.best_match_boxplot <- function(best_match, path_save, dataset.name){

	library(tidyverse)
	library(hrbrthemes)
	library(ggrepel)
	
	# Median
	data <- data.frame(
	  method = as.factor( rep(colnames(best_match),each = nrow(best_match)) ),
	  best_match = data.frame(best_match=unlist(best_match))[,"best_match"],
	  ct = as.factor( rep(rownames(best_match),ncol(best_match)) ))

	
	data$method <- with(data, reorder(method, best_match, median))

	pdf(width=20,height=10,file= paste(path_save,
	    paste0(dataset.name,"_best_match_boxplot.pdf"), sep = ""))

	 print(ggplot( data, mapping = aes(x=method, y=best_match, fill=method, label = ct) )
	    + geom_boxplot() 
	    + geom_jitter(color="black", size=0.5, alpha=0.9)
	    + theme_bw()
	    + theme(legend.position="none",
	      axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1))
	    + ggtitle("Cell type equivalent jaccard index") 
	    + xlab("Methods and references")
	    + ylab("Jaccard index"))


  	print(ggplot( data, mapping = aes(x=method, y=best_match, fill=method, label = ct) )
	    + geom_boxplot() 
	    + geom_jitter(color="black", size=0.5, alpha=0.9)
	    + geom_text_repel()
	    + theme_bw()
	    + theme(legend.position="none",
	      axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1))
	    + ggtitle("Cell type equivalent jaccard index") 
	    + xlab("Methods and references")
	    + ylab("Jaccard index"))


  	dev.off()

  	pdf(width=10,height=7,file= paste(path_save,
	    paste0(dataset.name,"_best_match_boxplot_small.pdf"), sep = ""))

	 print(ggplot( data, mapping = aes(x=method, y=best_match, fill=method, label = ct) )
	    + geom_boxplot() 
	    + geom_jitter(color="black", size=0.5, alpha=0.9)
	    + theme_bw()
	    + theme(legend.position="none",
	      axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1))
	    + ggtitle("Cell type equivalent jaccard index") 
	    + xlab("Methods and references")
	    + ylab("Jaccard index"))


  	print(ggplot( data, mapping = aes(x=method, y=best_match, fill=method, label = ct) )
	    + geom_boxplot() 
	    + geom_jitter(color="black", size=0.5, alpha=0.9)
	    + geom_text_repel()
	    + theme_bw()
	    + theme(legend.position="none",
	      axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1))
	    + ggtitle("Cell type equivalent jaccard index") 
	    + xlab("Methods and references")
	    + ylab("Jaccard index"))


  	dev.off()

  	return(data)

}

###################################

# Histologic patterns 

###################################

bea.integration.cancer = function(dataset, annot, normal,patient ,tumor,paper,stage,seu, dataset.name){

	#normal.cell <- rownames(annot[which(annot$SampleType == "Normal"),])
	normal.cell <- rownames(annot[which(annot[,normal] == "Normal"),])
	
	annot <- annot[!rownames(annot) %in% normal.cell,]
  new.cluster.ids <- rep("Cancer", nrow(annot))
  names(new.cluster.ids) <- rownames(annot)

	ConsensusTILS <- dataset@meta.data[,"ConsensusTILs"]
	names(ConsensusTILS) <- colnames(dataset)
	
	dataset <- seu 
	if(dataset.name == "qian"){
		dataset <- CreateSeuratObject(counts = seu)
		not_common <- colnames(dataset)[!colnames(dataset) %in% intersect(colnames(dataset),rownames(annot))]
		dataset <- dataset[,!colnames(dataset) %in% not_common]
	}

  dataset <- dataset[,!colnames(dataset) %in% normal.cell]

	for (i in seq_along(new.cluster.ids)) {
		if(names(new.cluster.ids[i]) %in% names(ConsensusTILS)){
			new.cluster.ids[i] <- ConsensusTILS[names(new.cluster.ids[i])]
		}
	}
	
	dataset[["Patient"]] <- annot[,patient]
	if(tumor != 0){
		dataset[["Pri_Met"]] <- annot[, tumor]
	}
  dataset[["Paper"]] <- annot[,paper]
  dataset[["Stage"]] <- annot[,stage]
	dataset[["ConsensusTILS"]] <- new.cluster.ids

	if(tumor != 0){
		for (i in seq_along(new.cluster.ids)) {
			if(new.cluster.ids[i] == "Cancer"){
				new.cluster.ids[i] <-  dataset@meta.data[names(new.cluster.ids[i]),"Pri_Met"]
			}
		}
		dataset[["ConsensusPM"]] <- new.cluster.ids
	}


  # Percentage mitochondrial genes and normalization
  dataset <- PercentageFeatureSet(dataset, pattern = "^MT-", col.name = "percent.mt")
  dataset <- SCTransform(dataset, vars.to.regress = "percent.mt", verbose = FALSE)

  # PCA and UMAP
  dataset <- RunPCA(dataset)
  dataset <- RunUMAP(dataset, dims = 1:30, verbose = FALSE)

  pdf(width=8,height=8,file= paste(path_save,
        paste0(dataset.name,"_final_annot_UMAP.pdf"), sep = ""))

  print(DimPlot(dataset, reduction = "umap", label = TRUE, group.by = "ConsensusTILS")
  	+ guides(color=guide_legend(override.aes = list(size=4),ncol =1)))

  dev.off()

  return(dataset)
}


bea.remove.noncancer.patient = function( dataset.name, dataset, annot, idents, cellID){
	nb_ct_patient <- bea.nb.ct.patient(dataset.name, dataset, idents)
	noncancer_patient <- rownames(nb_ct_patient[which(nb_ct_patient$Cancer == 0),])
	dataset <- dataset[,!colnames(dataset) %in% annot[annot$Patient %in% noncancer_patient, cellID]]

	return(dataset)
}

bea.nb.ct.patient = function( dataset.name, dataset, idents){

	annot_patient <- data.frame(annotation = dataset[[idents]], 
		patient = dataset[["Patient"]])

	nb_ct_patient <- as.data.frame(matrix(,nrow = length(unlist(unique(dataset[["Patient"]]))),
		ncol = length(unlist(unique(dataset[[idents]])))))
	rownames(nb_ct_patient) <- unlist(unique(dataset[["Patient"]]))
	colnames(nb_ct_patient) <- unlist(unique(dataset[[idents]]))

	for(pat in unlist(unique(dataset[["Patient"]]))){
		value <- c()
		for(ct in unlist(unique(dataset[[idents]]))){
			value <- c(value, 
			nrow(annot_patient[annot_patient[,idents] == ct & annot_patient$Patient == pat,]))
		}
		nb_ct_patient[rownames(nb_ct_patient) == pat,] <- value
	}
	return(nb_ct_patient)
}

bea.score.histo = function(scores, annot, col, normal, tumor, patient){

	#normal.cell <- rownames(annot[which(annot$SampleType == "Normal"),])
	normal.cell <- rownames(annot[which(annot[,normal] == "Normal"),])
	cancer.cell <- names(which(dataset$ConsensusTILS == "Cancer"))

	scores <- scores[!rownames(scores) %in% normal.cell,]
	scores <- scores[cancer.cell,]
	scores <- na.omit(scores)
	scores$diff_score <- scores$solid_up - scores$lepidic_up
	
	# CUT OFF QUANTILES 

	diff_score <- scores$solid_up - scores$lepidic_up
	class <- c()
	color <- c()
	for(sco in diff_score){
		if(sco > quantile(diff_score, probs = c(0.125,0.375,0.625,0.875))[4]){
			class <- c(class, "Solid")
			color <- c(color, "firebrick3")
		}
		else if(sco < quantile(diff_score, probs = c(0.125,0.375,0.625,0.875))[1]){
			class <- c(class, "Lepidic")
			color <- c(color, "dodgerblue4")
		}
		else if(sco < quantile(diff_score, probs = c(0.125,0.375,0.625,0.875))[4] & sco > quantile(diff_score, probs = c(0.125,0.375,0.625,0.875))[1] ){
			class <- c(class, "Transition")
			color <- c(color, "gray69")
		}
	}
	scores$histo <- class
	scores$color <- color
	if(tumor != 0){
		scores$annot_tumor <- annot[rownames(scores), tumor]
	}

	scores$stage <- NA
	for(i in seq_along(rownames(scores))){
		scores$stage[i] <- annot[rownames(scores)[i],col]
	}

	scores$Patient <- NA
	scores$patient <- NA
	for(i in seq_along(rownames(scores))){
		scores$Patient[i] <- annot[rownames(scores)[i],patient]
		scores$patient[i] <- paste0("Patient_",annot[rownames(scores)[i],patient])
	}

	return(scores)
}

bea.tumor.annotation = function(dataset, annot, col_name, annot_name){
	#col_name = name of col with tumor annotation
	#annot_name = name of new annotation col
	dataset[[annot_name]] <- dataset[["ConsensusTILS"]]
	for(coln in colnames(dataset)){
		if(dataset[[annot_name]][coln,1] == "Cancer"){
			dataset[[annot_name]][coln,1] <- annot[coln, col_name]
		}
	}
	return(dataset)
}

bea.HISTO.scatter = function(scores, path_save, dataset.name){

	scores$histo <- factor(scores$histo, levels = c("Solid" , "Transition", "Lepidic"))


	pdf(width=10,height=10,file= paste(path_save,
				paste0(dataset.name,"_HISTO_scatterplot.pdf"), sep = ""))

	print(ggplot(scores, aes(solid_up, lepidic_up)) +
	  geom_point(aes(color = diff_score))+
	  theme_bw() +
	  xlab("Solid score") +
	  ylab("Lepidic score") +
	  scale_color_gradient2(low = "dodgerblue4", mid = "gray69",
	                            high = "firebrick3", space = "Lab" ))

	print(ggplot(scores, aes(solid_up, lepidic_up, color = histo)) +
	  geom_point()+
	  theme_bw() +
	  xlab("Solid score") +
	  ylab("Lepidic score") +
	  scale_color_manual(values = c("firebrick3","gray69","dodgerblue4")))

	#print(ggplot(scores, aes(solid_up, lepidic_up, color = patient)) +
	#  geom_point()+
	#  theme_bw() +
	#  xlab("Solid score") +
	#  ylab("Lepidic score") +
	#  scale_color_manual(values = c("midnightblue","orchid3","mediumvioletred")))

	print(ggplot(scores, aes(x=solid_up, y=lepidic_up)) + 
		geom_point() +
		stat_cor(method="spearman"))	

	dev.off()
}


bea.barplot.MP.LTS = function(scores, dataset.name, path_save){
	stage <- c(rep("Metastasis" , 3), rep("Primary" , 3))
	Condition <- rep(c("Lepidic","Solid" , "Transition") , 2)
	value <- c(table(scores[scores$annot_tumor == "Metastasis","histo"]), table(scores[scores$annot_tumor == "Primary","histo"]))
	data <- data.frame(stage,Condition,value)

	data$Condition <- factor(data$Condition, levels = c("Solid" , "Transition", "Lepidic"))

	pdf(width=3,height=3, file= paste(path_save,
			paste0(dataset.name,"_HISTO_barplot.pdf"), sep = ""))

	#Stackef + Percent
	print(ggplot(data, aes(fill=Condition, y=value, x=stage)) + 
    	geom_bar(position="fill", stat="identity") +
    	theme_bw() +
	    xlab( "Sample Type" ) + ylab( "Proportions" ) +
	    scale_fill_manual(values = c("firebrick3","gray69","dodgerblue4")))

	stage <- c(rep(unique(scores$stage) , each = 3))
	Condition <- rep(c("Lepidic","Solid" , "Transition") , 3)
	value <-c()
	for(sta in unique(scores$stage)){
		value <- c(value, table(scores[scores$stage == sta,"histo"]))
	}
	data <- data.frame(stage,Condition,value)

	data$Condition <- factor(data$Condition, levels = c("Solid" , "Transition", "Lepidic"))

	print(ggplot(data, aes(fill=Condition, y=value, x=stage)) + 
    	geom_bar(position="fill", stat="identity") +
    	theme_bw() +
	    xlab( "Tumor stage" ) + ylab( "Proportions" ) +
	    scale_fill_manual(values = c("firebrick3","gray69","dodgerblue4")))


	dev.off()
}


bea.consensus.aggregate  =function(dataset, common_annot, annot_name_2, col_name_2){

	dataset[[annot_name_2]] <- dataset[["ConsensusHISTO"]]
	for(coln in colnames(dataset)){
		if(!dataset[[annot_name_2]][coln,1] %in% c("Solid", "Transition", "Lepidic")){
			dataset[[annot_name_2]][coln,1] <- common_annot[coln, col_name_2]
		}
	}	
	return(dataset)
}


bea.barplot.LTS.immune.pertype = function(dataset, dataset.name, path_save, idents){
	
	nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, idents)

	patient_list <- rownames(nb_ct_patient)

	ct = colnames(nb_ct_patient)
	Patient <- rep(patient_list, each = length(ct))
	Condition <- 	rep(ct,length(patient_list))
	value <- c()
	for(pat in patient_list){
		value <- c(value, as.numeric(nb_ct_patient[pat,ct]))
	}
	data <- data.frame(Patient,Condition,value)

	percent <- c()
	for (pat in patient_list) {
		for(cell in ct){
			total <- sum(data[data$Patient == pat, "value"])
			percent <- c(percent, data[data$Condition == cell & data$Patient == pat,"value"]*100/total)
		}
	}

	data$percent <- percent

	percent_LTS <- c()
	for (pat in patient_list) {
		for(cell in ct){
			if(cell %in% c("Solid", "Transition", "Lepidic")){
				total <- sum(data[data$Patient == pat & data$Condition %in% c("Solid", "Transition", "Lepidic"), "value"])
				percent_LTS <- c(percent_LTS, data[data$Condition == cell & data$Patient == pat,"value"]*100/total)
			}
			else {
				percent_LTS <- c(percent_LTS, NA)
			}
		}
	}

	data$percent_LTS <- percent_LTS

	order_patient <- unique(data$Patient)[order(data[data$Condition == "Lepidic", "percent_LTS"], decreasing = TRUE)]
	solid_oder <- unique(data$Patient)[order(data[data$Condition == "Solid", "percent_LTS"])]
	lepidic_nonull <- data[data$Condition == "Lepidic" & data$value != 0, "Patient"]
	solid_oder <- solid_oder[!solid_oder %in% lepidic_nonull]
	order_patient <- order_patient[!order_patient %in% solid_oder]
	order_patient <- c(order_patient, solid_oder)

	subset <- data[data$Condition %in% c("Lepidic" , "Solid" , "Transition"),]

	subset$Patient <- factor(subset$Patient, levels = order_patient)
	subset$Condition <- factor(subset$Condition, levels = c("Solid" , "Transition" , "Lepidic"))


	pdf(width=3,height=3.5,file= paste(path_save,
			paste0(dataset.name,"_barplot_LTS_immune_pertype.pdf"), sep = ""))

	par(mfrow = c(2, 1))

	#Stacked 
	print(ggplot(subset, aes(fill=Condition, y=value, x=Patient)) + 
    	geom_bar(position="fill", stat="identity") +
	    xlab( "Patient" ) + ylab( "Histologic pattern proportion" ) +
	    theme_bw() +
	    scale_fill_manual(values = c("firebrick3","gray69","dodgerblue4")) + 
	    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)) +
	    guides(fill=guide_legend(title="Histologic pattern")))

	#subset <- data[!data$Condition %in% c("Lepidic" , "Solid" , "Transition"),]
	#subset$Patient <- factor(subset$Patient, levels = order_patient)

	data$Patient <- factor(data$Patient, levels = order_patient)
	
	#Stacked 
	print(ggplot(data, aes(fill=Condition, y=value, x=Patient)) + 
    	geom_bar(position="fill", stat="identity") +
	    xlab( "Patient" ) + ylab( "Cell type proportion" ) +
	    theme_bw() + 
	    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)) +
	    guides(fill=guide_legend(title="Cell type")))

	dev.off()

	return(data)
}


bea.LTS.immune.corr = function(dataset.name, path_save, LTS.immune){
	library(corrplot)
	library(grid)

	LTS.immune$Condition <- gsub(" ", ".", LTS.immune$Condition)
	
	ct = unique(LTS.immune$Condition)[!unique(LTS.immune$Condition) %in% c("Solid", "Transition", "Lepidic")]
	patient_list <- rownames(LTS.immune)

	corr_mat <- matrix(,nrow = 3, ncol =  length(ct), 
		dimnames = list(c("Lepidic", "Transition", "Solid"), ct))
	corr_mat_p <- matrix(,nrow = 3, ncol =  length(ct), 
		dimnames = list(c("Lepidic", "Transition", "Solid"), ct))

	for(col in colnames(corr_mat)){
		for(row in rownames(corr_mat)){
			# either percent or value ?????
			cor.result <- cor.test(as.vector(LTS.immune[LTS.immune$Condition == col,"percent"]),
				as.vector(LTS.immune[LTS.immune$Condition == row,"percent_LTS"]),
				 method = "spearman")

			corr_mat[row, col] <- cor.result$estimate
			corr_mat_p[row, col] <- cor.result$p.value
		}
	}

	
	pdf(width=8,height=8,file= paste(path_save,
		paste0(dataset.name,"_LTS_corr_mat.pdf"), sep = ""))

	print(corrplot(corr_mat, type="upper", p.mat = corr_mat_p, sig.level = 0.05, 
		insig = "blank", title= "Spearman correlation of cell type proportion",
		mar=c(0,0,10,0)))

	print(corrplot(corr_mat, type="upper", p.mat = corr_mat_p, sig.level = 0.05, 
		insig = "blank", title= "Spearman correlation of cell type proportion",
		mar=c(0,0,10,0), method = "number"))

	myTable <- tableGrob(round(corr_mat, 4), theme=ttheme_minimal(base_size = 8))
	plot.new()
	print(grid.draw(myTable))
	myTable <- tableGrob(round(corr_mat_p, 4), theme=ttheme_minimal(base_size = 8))
	plot.new()
	print(grid.draw(myTable))

	dev.off()

	mat <- list(corr_mat, corr_mat_p)

	return(mat)
}



###############################################

# Cell cell interactions analysis 

################################################


bea.liana.patient = function( dataset.name, dataset, idents){
	library(plyr)
	library(dplyr) 
	library(Seurat)
	library(patchwork)
	library(ggplot2)
	library(sctransform)
	library(gplots)
	library(RColorBrewer)
	library(celldex)

	library(tidyverse)
	library(magrittr)
	library(liana)
	library(circlize)

	# Remove non immune cells 
	non.immune <- c()
	for(cell in colnames(dataset)){
		if(dataset@meta.data[cell, idents] %in% c("Epithelial", "Fibroblast", "Endothelial",
			"Neuron_Glia")){
			non.immune <- c(non.immune, cell)
		}
	}
	#length(non.immune) #5432 cells for Kim

	dataset <- dataset[,!colnames(dataset) %in% non.immune]

	nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, idents)
	Idents(dataset) <- dataset[[idents]]

	liana_list <- list()
	for(i in 1:length(unique(dataset@meta.data$Patient))){
		print(i)
		ten_cells <- colnames(nb_ct_patient)[which(nb_ct_patient[rownames(nb_ct_patient) == unique(dataset@meta.data$Patient)[i],] > 10)]

		liana_test <- liana_wrap(dataset[,Idents(dataset) %in% ten_cells & dataset@meta.data$Patient == unique(dataset@meta.data$Patient)[i]],
			idents_col = idents)
		liana_test <- liana_test %>% liana_aggregate()
		liana_list[[i]] <- liana_test
	}
	return(liana_list)

}

bea.liana.aggregate = function( dataset.name, dataset, idents, path_save, liana_list, sig){


	nb_patient <- length(unique(dataset@meta.data$Patient))

	liana_aggregate <- data.frame(matrix(,nrow = 0, ncol = 4 + nb_patient, 
		dimnames = list(NULL, c("source", "target", "ligand", "receptor", unique(dataset@meta.data$Patient)))))
	
	for(j in 1:nb_patient){
		print(j)
		liana_list[[j]] <- as.data.frame(liana_list[[j]])
		for(k in 1:dim(liana_list[[j]][liana_list[[j]]$aggregate_rank <= sig ,])[1]){
			rowname <- paste(liana_list[[j]]$source[k], liana_list[[j]]$target[k], 
				liana_list[[j]]$ligand.complex[k], liana_list[[j]]$receptor.complex[k], sep = "_")
			
			if(rowname %in% rownames(liana_aggregate)){
				liana_aggregate[rowname,colnames(liana_aggregate) == unique(dataset@meta.data$Patient)[j]] <- 1
			}
			
			else {
				temp <- data.frame(matrix(,nrow = 1, ncol = 4 + nb_patient, 
					dimnames = list(rowname, c("source", "target", "ligand", "receptor", unique(dataset@meta.data$Patient)))))
				temp[1,] <- c(liana_list[[j]][k,1:4], rep(0, nb_patient))
				liana_aggregate <- rbind(liana_aggregate, temp)
				liana_aggregate[rowname,unique(dataset@meta.data$Patient)[j]] <- 1

			}

		}
	}

	sum_patient <- c()
	for(i in 1:nrow(liana_aggregate)){
	sum_patient <- c(sum_patient, sum(liana_aggregate[i,5:dim(liana_aggregate)[2]]))
	}
	liana_aggregate <- cbind(liana_aggregate, sum_patient)


	save(liana_aggregate, file = paste0(path_save ,dataset.name,"_",idents,"_liana_aggregate.RData"))
	return(liana_aggregate)

}


bea.liana.heatmap = function( dataset.name, liana_aggregate, at_least, idents, path_save, patient, patient_list){

	if(patient == "no"){
		limit_pos <- c(which(liana_aggregate$sum_patient >= at_least))
		sub <- liana_aggregate[limit_pos,]
	} 
	if(patient == "yes"){
		limit_pos <- c(which(liana_aggregate$sum_patient >= at_least))
		sub <- liana_aggregate[limit_pos,c("source", "target", "ligand", "receptor",patient_list, "sum_patient")]
		sub$sum_patient <- rowSums(sub[,5:(dim(sub)[2]-1)])
		limit_pos <- c(which(sub$sum_patient >= at_least))
		sub <- sub[limit_pos,]
	}
	if(patient != "no" & patient != "yes") {
		limit_pos <- rownames(liana_aggregate[liana_aggregate[,patient] ==1 & liana_aggregate$sum_patient >=2,])
		sub <- liana_aggregate[limit_pos,]
	}

	int_mat <- matrix(,nrow = length(unique(sub[,"source"]))
		, ncol = length(unique(sub[,"target"])))
	colnames(int_mat) <- unlist(unique(sub[,"target"]))
	rownames(int_mat) <- unlist(unique(sub[,"source"]))

	for(sour in unique(unique(sub[,"source"]))){
		value <- c()
		for(tar in unlist(unique(sub[,"target"]))){
			value <- c(value, dim(sub[sub$source == sour & sub$target %in% tar,])[1])
		}
		int_mat[sour,] <- value
	}


	if(ncol(int_mat) >= 2 & nrow(int_mat) >= 2){
			or <- apply(int_mat, 1, max, na.rm = TRUE)
		or <- order(or, decreasing = TRUE)

		int_mat <- int_mat[or,]


		or <- apply(int_mat, 2, max, na.rm = TRUE)
		or <- order(or, decreasing = TRUE)

		int_mat <- int_mat[,or]

		######PDF file 
		my_palette <- colorRampPalette(c("white","aquamarine2" ,"darkgreen"))

		if(patient == "no"  | patient == "yes"){
			pdf(width=10,height=10,file= paste(path_save,
		    	paste0(dataset.name,"_liana_htmap_atleast_",at_least ,"_",idents,".pdf"), sep = ""))

			print(heatmap.2(x=int_mat, 
				Colv=FALSE, 
				Rowv=FALSE,
				scale="none",
				col=my_palette,
				trace="none",
				main= "Number of interactions",
				ylab="Ligand",
				xlab="Receptor",
				margins=c(14,14)))

			dev.off()
		} else {

			print(heatmap.2(x=int_mat, 
				Colv=FALSE, 
				Rowv=FALSE,
				scale="none",
				col=my_palette,
				trace="none",
				main= paste("Number of interactions - ", patient),
				ylab="Ligand",
				xlab="Receptor",
				margins=c(14,14)))	
		}
	}
	else {
		print("There is only one source or target")
	}


	return(int_mat)
}

bea.pair_ct_interactions <- function(int_mat, dataset.name, path_save, annotation, patient){

	n = max(length(colnames(int_mat)), length(rownames(int_mat)))

	#nb_pairs = (n*(n-1))/2 -> NROW SHOULD BE THIS NUMBER + N


	pair_int <- matrix(,nrow = n,ncol = n)
	colnames(pair_int) <- colnames(int_mat)
	rownames(pair_int) <- colnames(int_mat)

	for(col in colnames(int_mat)){
		for( row in rownames(int_mat)){
			if(col != row){
				pair_int[col, row] <- int_mat[col,row] + int_mat[row,col]
				pair_int[row, col] <- int_mat[col,row] + int_mat[row,col]
			}
			else {
				pair_int[col, row] <- int_mat[col,row]*2 #ADDED TIMES 2
				#pair_int[row, col] <- int_mat[col,row]*2 #ADDED TIMES 2
			}
		}
	}

	my_palette <- colorRampPalette(c("white","aquamarine2" ,"darkgreen"))

	if(patient == "no"){
		pdf(width=15,height=15,pointsize=20,file= paste(path_save,
	    	paste0(dataset.name,"_liana_",annotation,"_pair_ct_htmap.pdf"), sep = ""))

		print(heatmap.2(x=pair_int, 
			Colv=FALSE, 
			Rowv=FALSE,
			scale="none",
			col=my_palette,
			trace="none",
			main= "Cell type pair interactions",
			ylab="CT1",
			xlab="CT2",
			margins=c(14,14)))

	  	dev.off()
	} else {
		print(heatmap.2(x=pair_int, 
			Colv=FALSE, 
			Rowv=FALSE,
			scale="none",
			col=my_palette,
			trace="none",
			main= paste("Cell type pair interactions -", patient),
			ylab="CT1",
			xlab="CT2",
			margins=c(14,14)))
	}

	return(pair_int)
}

bea.patient.cor.htmap <- function(dataset.name, path_save, dataset, pair_int, title){


	mat_cor <- matrix(,nrow = length(unique(dataset@meta.data$Patient)), ncol = length(unique(dataset@meta.data$Patient)), 
			dimnames = list(unique(dataset@meta.data$Patient), unique(dataset@meta.data$Patient)))

	mat_p <- matrix(,nrow = length(unique(dataset@meta.data$Patient)), ncol = length(unique(dataset@meta.data$Patient)), 
			dimnames = list(unique(dataset@meta.data$Patient), unique(dataset@meta.data$Patient)))

	for(i in seq_along(unique(dataset@meta.data$Patient))){
		for(j in seq_along(unique(dataset@meta.data$Patient))){
			co_ct <- intersect(colnames(pair_int[[i]]), colnames(pair_int[[j]]))


			subset1 <- pair_int[[i]][co_ct,co_ct]
			subset2 <- pair_int[[j]][co_ct,co_ct]

			lower1 <- lower.tri(subset1, diag = FALSE)
			lower2 <- lower.tri(subset2, diag = FALSE)

			vec1 <- subset1[lower1]
			vec2 <- subset2[lower2]

			cor.result <- cor.test(vec1,vec2,method = "spearman")

			mat_cor[i,j] <- cor.result$estimate
			mat_p[i,j] <- cor.result$p.value
		}
	}

	distribution <- as.data.frame(mat_cor[lower.tri(mat_cor, diag = FALSE)])
	colnames(distribution) <- "Spearman"

	pdf(width=10,height=10,file= paste(path_save,
		paste0(dataset.name,"_patient_corr_mat_",title,".pdf"), sep = ""))

	print(corrplot(mat_cor, type="upper", p.mat = mat_p, sig.level = 0.05, 
		insig = "blank", title= "Spearman correlation between patient",
		mar=c(0,0,5,0)))

	print(corrplot(mat_cor, type="upper", p.mat = mat_p, sig.level = 0.05, 
		insig = "blank", title= "Spearman correlation between patient",
		mar=c(0,0,5,0), method = "number"))

	print(ggplot(distribution, aes(x=Spearman)) + 
  			geom_histogram() + 
  			xlim(-1,1) + 
  			theme_bw() + 
  			xlab("Spearman coefficient") + 
    		ylab("Density") + 
    		theme(plot.margin = margin(5, 5, 5, 5, "cm"), text = element_text(size = 20)) +
    		guides(fill=guide_legend(title="Distribution of Spearman correlation between patient")))

	dev.off()

	mat <- list(mat_cor, mat_p)

	return(mat)
}

nb_cell_nb_int_patient_scatter <- function(liana_aggregate, dataset.name, path_save, idents, dataset, annotation){
	library(ggpubr)

	non.immune <- c()
	for(cell in colnames(dataset)){
		if(dataset@meta.data[cell, idents] %in% c("Epithelial", "Fibroblast", "Endothelial",
			"Neuron_Glia")){
			non.immune <- c(non.immune, cell)
		}
	}
	#length(non.immune) #5432 cells for Kim

	dataset <- dataset[,!colnames(dataset) %in% non.immune]

	nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, idents)
	nb_ct <- dim(nb_ct_patient)[2]
	nb_patient <- dim(nb_ct_patient)[1]

	nb_cell_nb_int <- data.frame("patient" = rep(rownames(nb_ct_patient), each = nb_ct),
		"ct" = rep(colnames(nb_ct_patient),nb_patient), "nb_cell" = NA, "nb_int" = NA)

	for(pat in nb_cell_nb_int[,"patient"]){
		for (ct in nb_cell_nb_int[,"ct"]) {
			nb_cell_nb_int[nb_cell_nb_int$patient == pat & nb_cell_nb_int$ct == ct,"nb_cell"] = nb_ct_patient[pat,ct]
			nb_cell_nb_int[nb_cell_nb_int$patient == pat & nb_cell_nb_int$ct == ct,"nb_int"] = sum(liana_aggregate[liana_aggregate$source == ct | liana_aggregate$target == ct & liana_aggregate$sum_patient > 1,pat])
		}
	}

	nb_cell_nb_int <- nb_cell_nb_int[nb_cell_nb_int$nb_int > 10,]

	pdf(width=7,height=7, file= paste(path_save,
    	paste0(dataset.name,"_",annotation,"_nb_cell_nb_int_patient_scatter.pdf"), sep = ""))

	print(ggplot(nb_cell_nb_int, aes(x=log(nb_cell), y=log(nb_int), color=ct)) + 
		geom_point()+
		theme_bw()+
		xlab("Log Number of cells" ) + 
		ylab("Log Number of interactions"))

	print(ggplot(nb_cell_nb_int, aes(x=log(nb_cell), y=log(nb_int))) + 
		geom_point() +
		stat_cor(method="spearman"))

  	dev.off()

	return(nb_cell_nb_int)
}


bea.ctpair_histo_inter_boxplot <- function(liana_aggregate, dataset.name, path_save, idents, dataset ){
	
	limit_pos <- c(which(liana_aggregate$sum_patient >= 2))
	liana_aggregate <- liana_aggregate[limit_pos,]

	non.immune <- c()
	for(cell in colnames(dataset)){
		if(dataset@meta.data[cell, idents] %in% c("Epithelial", "Fibroblast", "Endothelial",
			"Neuron_Glia")){
			non.immune <- c(non.immune, cell)
		}
	}
	#length(non.immune) #5432 cells for Kim

	dataset <- dataset[,!colnames(dataset) %in% non.immune]

	ct_pair_histo_nb_int <- data.frame(matrix(,nrow = 0, ncol = 3, 
		dimnames = list(NULL, c("ct_pair", "nb_int", "patient"))))

	patient <- unlist(unique(dataset[["Patient"]]))

	nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, idents)

	histo_pair <- list(c("Transition", "Lepidic"), c("Transition", "Solid"), c("Solid", "Lepidic"))

	for(pat in patient){
		for (pair in histo_pair) {
			if(nb_ct_patient[pat, pair[1]] > 10 & nb_ct_patient[pat, pair[2]] > 10){
				nb_int <- sum(liana_aggregate[liana_aggregate$source == pair[1] & liana_aggregate$target == pair[2]| liana_aggregate$source == pair[2] & liana_aggregate$target == pair[1],pat])
				temp <- data.frame(matrix(,nrow = 1, ncol = 3, dimnames = list(NULL, c("ct_pair", "nb_int", "patient"))))
				temp[1,] <- c(paste0(pair[1],"_",pair[2]),nb_int,pat)
				ct_pair_histo_nb_int <- rbind(ct_pair_histo_nb_int, temp)
			}
			
		}
	}

	ct_pair_histo_nb_int <- ct_pair_histo_nb_int[order(ct_pair_histo_nb_int$ct_pair),]

	ct_pair_histo_nb_int$nb_int <- as.numeric(ct_pair_histo_nb_int$nb_int)

	pdf(width=5,height=5,file= paste(path_save,
    	paste0(dataset.name,"_ctpair_HISTO_inter_boxplot.pdf"), sep = ""))

	print(ggplot(ct_pair_histo_nb_int, aes(x=ct_pair, y=nb_int, fill=ct_pair)) +
	    geom_boxplot() +
	    geom_jitter(color="black", size=0.4, alpha=0.9) +
	    theme_bw() +
	    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)) +
	    xlab("Histologic pattern") + 
	    ylab("Number of interactions per patient"))

	
	dev.off()

	return(ct_pair_histo_nb_int)

}


bea.LTS.interactions <- function(liana_aggregate, dataset.name, path_save, idents, dataset, at_least){

	nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, idents)
	#nb_ct_patient <- bea.nb.ct.patient( dataset.name, dataset, "ConsensusHISTO")
	nb_patient <- dim(nb_ct_patient)[1]
	ct <- colnames(nb_ct_patient)
	ct <- ct[!ct %in% c("Solid", "Lepidic", "Transition", "Endothelial", "Epithelial", "Fibroblast", "Neuron_Glia")]
	nb_ct <- length(ct)
	pat <- rownames(nb_ct_patient)
	if(dataset.name == "qian"){
		pat <- paste0("Patient_", rownames(nb_ct_patient))
		rownames(nb_ct_patient) <- pat
	}



	variety <- rep(ct, each= nb_patient*3)
	treatment <- rep(c(rep(c("Solid","Transition", "Lepidic"),each= nb_patient)),nb_ct)
	note <- rep(NA, length(variety))
	patient <- rep(pat , nb_ct*3)
	LTS_int <- data.frame(variety, treatment ,  note, patient)

	LTS_int$treatment <- factor(LTS_int$treatment, levels = c("Lepidic" , "Transition", "Solid"))

	for(pat in patient){
		for (type in ct) {
			for (LTS in c("Solid", "Transition", "Lepidic")) {
				if(nb_ct_patient[pat, type] > 10 & nb_ct_patient[pat, LTS] > 10){
					nb_int <- sum(liana_aggregate[liana_aggregate$source == type & liana_aggregate$target == LTS| liana_aggregate$source == LTS & liana_aggregate$target == type & liana_aggregate$sum_patient >= at_least,pat])
					LTS_int[LTS_int$variety == type & LTS_int$treatment == LTS & LTS_int$patient == pat,"note"] <- nb_int

				}
			}
		}
	}
 	
	
	pdf(width=10,height=8,file= paste(path_save,
    	paste0(dataset.name,"_LTS_boxplot.pdf"), sep = ""))


	# grouped boxplot
	print(ggplot(LTS_int, aes(x=variety, y=note, fill=treatment)) + 
    	geom_boxplot()
    	+ theme_bw()
    	+ scale_fill_manual(values = c("dodgerblue4","gray69","firebrick3"))
    	+ theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)) 
    	+ xlab("Cell types")
    	+ ylab("Number of interactions")
    	+ guides(fill=guide_legend(title="Histologic pattern")))



	dev.off()

	return(LTS_int)
}




