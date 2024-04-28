library("clusterProfiler")
library("sceasy")
library(Seurat)
library(ggplot2)
library(sceasy)
library(Seurat)
library(reticulate)

# Environment setup and data format conversion.
use_condaenv("my_env", required=TRUE)
sc <- import("scanpy", convert = FALSE)
pathology <- sc$read_h5ad("[pathology_adata].h5ad")
convertFormat("pathology_adata.h5ad", from="anndata", to="seurat", outFile='pathology_adata.rds')
data <- readRDS("pathology_adata.rds")



# heatmap
library(reshape2)
clinical <- read.csv("Clinico-pathology.csv",header = T)
clinical$cogn_global_lv <- -clinical$cogn_global_lv
clinical_norm <- clinical
clinical_norm[,-c(1,9)] <- scale(clinical[,-c(1,9)])
clinical_norm$id <- 1:48
clinical_norm <- clinical_norm[,-c(1,9)]
colnames(clinical_norm)[c(1,3,4)] <- c("Amyloid","NFT","Tangles")
clinical_norm_stack <- melt(clinical_norm,id.vars = "id")

ggplot(clinical_norm_stack, aes(x = id, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#0073B1")+
  geom_rect(aes(xmin=0.5,xmax=24.45,ymin=0.5,ymax=7.5),
            col="#D4A2BE",alpha=0,)+
  geom_rect(aes(xmin=24.55,xmax=39.47,ymin=0.5,ymax=7.5),
            col="#226089",alpha=0)+
  geom_rect(aes(xmin=39.53,xmax=48.46,ymin=0.5,ymax=7.5),
            col="#900C27",alpha=0)+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size=unit(1.8, "cm"))+
  coord_flip()
ggsave("heatmap.pdf",width = 4,height = 24)

# UMAP
data <- ScaleData(data)
data <- RunPCA(data)
dat] <- RunUMAP(data, dims = 1:10)
p21 <- DimPlot(data,reduction = "umap",
                group.by = [feature_for_group], pt.size = 0.4)
p21+theme(axis.text = element_blank(),
           axis.ticks = element_blank(),
           axis.line = element_blank(),
           panel.border = element_rect(colour = "black"),
           legend.position = c(0.01,0.95),
           legend.spacing.x = unit(0.1, 'cm'))+
  labs(title = "data",x="UMAP1",y="UMAP2")



#read biclustering results
read_txt_to_list <- function(file_path) {
    data <- readLines(file_path)
    dict_list <- list()
    current_dict <- list()
    for (line in data) {
        if (startsWith(line, "bc")) {
            if (length(current_dict) > 0) {
                dict_list <- c(dict_list, list(current_dict))
                current_dict <- list()
              }
          } else if (grepl("^\\d+", line)) {
              current_dict$c <- as.integer(unlist(strsplit(line, " ")))
            } else {
                current_dict$r <- unlist(strsplit(line, " "))
              }
      }
    if (length(current_dict) > 0) {
        dict_list <- c(dict_list, list(current_dict))
      }
    return(dict_list)
  }

#venn plot
get_specific_and_common_gene <- function(res,L){
    common <- res[[1]]$r
    for(i in 2:L){
        common <- intersect(common,res[[i]]$r)
      }
    specific <- list()
    for(i in 1:L){
        specific[[paste0("bc",i)]] <- setdiff(res[[i]]$r,common)
      }
    return(list(common=common,specific=specific))
}

getgene <- function(res,L){
    gs=list()
    for(i in 1:L){
        gs[[paste0("bc",i)]] <- res[[i]]$r
      }
    gs
  }
library(VennDetail)
gene_list <- getgene([bicluster_result],6)
venn <- venndetail(gene_list)
getSet(venn,subset = c("bc1_bc2"))
library(RColorBrewer)
mycol <- brewer.pal(6,"Set1")
pdf("upset.pdf",width = 18,height = 6)
plot(venn,type = "upset")
dev.off()

# calculate os and iou
overlap_ratio <- function(ref,my,includeself=T){
    m <- length(ref)
    n <- length(my)
    res <- data.frame(refbc=rep(NA,n),IoU=rep(NA,n))
      for(i in 1:n){
        overlap <- rep(0,m)
        for(j in 1:m){
            overlap[j] <- length(intersect(my[[i]],ref[[j]]))/
                length(union(my[[i]],ref[[j]]))
            if(!includeself&j==i)overlap[j] <- 0
          }
        res[i,1] <- which.max(overlap)
        res[i,2] <- max(overlap)
      }
    res
  }

overlap_exclusive_self <- function(ref){
    library(VennDetail)
    m <- length(ref)
    name <- names(ref)
    venn <- venndetail(ref)
    res <- data.frame(refbc=rep(NA,m),score=rep(NA,m))
    for(i in 1:m){
        overlap <- rep(0,m)
        for(j in 1:m){
            if(i<j){
                p <- i
                q <- j
              }
            else{
                p <- j
                q <- i
              }
            overlap[j] <- nrow(getSet(venn,subset = paste0(name[p],"_",name[q])))/max(c(length(ref[[i]]),length(ref[[j]])))
            if(j==i)overlap[j] <- 0
          }
        res[i,1] <- which.max(overlap)
        res[i,2] <- max(overlap)
      }
    res
  }

overlap_exclusive_score <- function(ref,my){
    library(VennDetail)
    name <- names(ref)
    m <- length(ref)
    n <- length(my)
    res <- data.frame(refbc=rep(NA,n),score=rep(NA,n))
    for(i in 1:n){
        ref$mybc <- my[[i]]
        venn <- venndetail(ref)
        overlap <- rep(0,m)
        for(j in 1:m){
            overlap[j] <- nrow(getSet(venn,subset = paste0(name[j],"_mybc")))/
                max(c(length(my[[i]]),length(ref[[j]])))
          }
        res[i,1] <- which.max(overlap)
        res[i,2] <- max(overlap)
      }
    res
  }

# Extract genes that are shared and unique across different biclusters
gene_analysis <- get_specific_and_common_gene([bicluster_result],6)
gene_ref <- gene_analysis$specific
write.table(overlap_exclusive_self(gene_ref),"bc_align.txt",sep = ",",quote = F)
length(intersect(gene_ref[1],gene_ref[3]))/length(union(gene_ref[1],gene_ref[3]))
write.table(gene_analysis$specific[[1]],"specific_and_common_gene_set.txt",row.names = rep("bc1",length(gene_analysis$specific[[1]])),col.names = F)
for(i in 2:6)
  {
    write.table(gene_analysis$specific[[i]],"specific_and_common_gene_set.txt",row.names = rep(paste0("bc",i),length(gene_analysis$specific[[i]])), append = T,col.names = F)
    }
write.table(gene_analysis$common,"specific_and_common_gene_set.txt",row.names = rep("common",length(gene_analysis$common)), col.names = F,append = T)
test <- read.table("specific_and_common_gene_set.txt",header = F)
write.csv(test,"specific_and_common_gene_set.csv",row.names = F)

# Merge biclusters into FGM(Adjust as needed based on the actual situation.)
combine_gene <- union(gene_analysis1$specific[[1]],gene_analysis1$specific[[2]])
combine_gene <- union(gene_analysis1$specific[[3]],combine_gene)
combine_gene <- union(gene_analysis2$specific[[1]],combine_gene)
combine_gene <- union(gene_analysis2$specific[[2]],combine_gene)
combine_gene <- union(gene_analysis3$specific[[2]],combine_gene)

go <- enrichGO(combine_gene, OrgDb = Db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05,keyType = "SYMBOL")
write.csv(go@result,file="FGM1.csv",row.names=F)
dotplot(go,title="FGM1")
ggsave("FGM1.pdf")

no_combine_cell <- union(no_bicluster_result[[1]]$c, no_bicluster_result[[2]]$c)
no_combine_cell <- union(no_bicluster_result[[3]]$c, no_combine_cell)

early_combine_cell <- union(early_bicluster_result[[1]]$c, early_bicluster_result[[2]]$c)

late_combine_cell <- late_bicluster_result[[2]]$c

# Assign dual labels of cell type and FGM to each cell based on the cell_type
no_cell_type <- unlist(read.table("[no_cell_type].txt"))
len_no <- 0
for (i in 1:6) {
    len_no <- len_no + length(no_bicluster_result[[i]]$c)
  }
FGM_info <- character(len_no)
FGM_info[no_combine_cell] <- 'FGM1'
FGM_info[no_combine_cell2] <- 'FGM2'
FGM_info[no_combine_cell3] <- 'FGM3'
no_cell <- union(no_combine_cell, no_combine_cell2)
no_cell <- union(no_combine_cell3, no_cell)
no_cell_info <- data.frame(
  index = sort(no_cell),
  cell_type = no_cell_type,
  fgm_info = FGM_info
)
write.csv(no_cell_info, "[no_cell_info.csv]", row.names = FALSE)

# FGM perturbation analysis

increase_set <- combine_gene2
decrease_set <- combine_gene
increase <- setdiff(increase_set,decrease_set)
write.csv(increase,file = "comparison_enrich\\Ast\\Ast_early_no_inc_gene.csv",row.names = F)
go <- enrichGO(increase, OrgDb = Db, ont='BP',pAdjustMethod = 'BH', pvalueCutoff = 0.05,keyType = "SYMBOL")
write.csv(go@result,file = "comparison_enrich\\Ast\\Ast_early_no_inc.csv",row.names = T)
increase_table <- -log10(go@result$p.adjust[1:8])
dotplot(go,title="increase")
decrease <- setdiff(decrease_set,increase_set)
write.csv(decrease,file = "comparison_enrich\\Ast\\Ast_early_no_dec_gene.csv",row.names = F)
go2 <- enrichGO(decrease, OrgDb = Db, ont='BP',pAdjustMethod = 'BH', pvalueCutoff = 0.05,keyType = "SYMBOL")
write.csv(go2@result,file = "comparison_enrich\\Ast\\Ast_early_no_dec.csv",row.names = T)
dotplot(go2,title="decrease")

go_dec <- read.csv("[Ex_early_no_dec.csv]", header = TRUE, row.names = 1)
pathways_of_dec <- c("gliogenesis","myelination","ensheathment of neurons")
go_inc <- read.csv("[Ex_early_no_inc.csv]", header = TRUE, row.names = 1)
pathways_of_inc <- c("positive regulation of MAPK cascade","regulation of plasma membrane bounded cell projection assembly","modulation of chemical synaptic transmission")
selected_pathways <- data.frame(pathway = character(), p = numeric(), GeneRatio = numeric(), stringsAsFactors = FALSE)
for (pathway in pathways_of_inc) {
    p_value <- go_inc[go_inc$Description == pathway, "pvalue"]
    GeneRatio <- go_inc[go_inc$Description == pathway, "GeneRatio"]
    pathway_df <- data.frame(pathway = pathway, p = p_value, GeneRatio=GeneRatio)
    selected_pathways <- rbind(selected_pathways, pathway_df)
  }
for (pathway in pathways_of_dec) {
    p_value <- go_inc[go_dec$Description == pathway, "pvalue"]
    GeneRatio <- go_inc[go_dec$Description == pathway, "GeneRatio"]
    pathway_df <- data.frame(pathway = pathway, p = p_value, GeneRatio=GeneRatio)
    selected_pathways <- rbind(selected_pathways, pathway_df)
  }
write.csv(selected_pathways,file = "comparison_enrich\\Ex\\Ex_early_no.csv",row.names = T)

go_dec <- read.csv("[Ex_early_no_dec.csv]", header = TRUE, row.names = 1)
pathways_of_dec <- c("gliogenesis","myelination","ensheathment of neurons")
go_inc <- read.csv("[Ex_early_no_inc.csv]", header = TRUE, row.names = 1)
pathways_of_inc <- c("positive regulation of MAPK cascade","regulation of plasma membrane bounded cell projection assembly","modulation of chemical synaptic transmission")
selected_pathways <- data.frame(pathway = character(), p = numeric(), GeneRatio = numeric(), stringsAsFactors = FALSE)
for (pathway in pathways_of_inc) {
    p_value <- go_inc[go_inc$Description == pathway, "pvalue"]
    GeneRatio <- go_inc[go_inc$Description == pathway, "GeneRatio"]
    pathway_df <- data.frame(pathway = pathway, p = p_value, GeneRatio=GeneRatio)
    selected_pathways <- rbind(selected_pathways, pathway_df)
  }
for (pathway in pathways_of_dec) {
    p_value <- go_inc[go_dec$Description == pathway, "pvalue"]
    GeneRatio <- go_inc[go_dec$Description == pathway, "GeneRatio"]
    pathway_df <- data.frame(pathway = pathway, p = p_value, GeneRatio=GeneRatio)
    selected_pathways <- rbind(selected_pathways, pathway_df)
  }
write.csv(selected_pathways,file = "comparison_enrich\\Ex\\Ex_early_no.csv",row.names = T)



