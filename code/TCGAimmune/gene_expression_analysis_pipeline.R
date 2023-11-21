library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("RColorBrewer")
library("genefilter")
library("immunedeconv")
library("DESeq2")
library("ggplot2")
library("corrplot")
library("ggpubr")
library("ggpubr")
library(rstatix)
library(dplyr)
library(tidyr)
library(stringr)


GDCprojects = getGDCprojects() 
head(GDCprojects[c("project_id", "name")])
immunedeconv::timer_available_cancers
save_dir = "gene_expression_analysis/"

if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

target_genes = c("ARID1A","IFNA1", "IFNB1", "IL6", "CXCL9", "CXCL10")

get_assay <- function(cancertype){
  # Query platform Illumina HiSeq with a list of barcode 
  query_TCGA = TCGAbiolinks::GDCquery(
    project = paste0("TCGA-", cancertype),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
  )
  
  # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
  TCGAbiolinks::GDCdownload(query_TCGA)
  
  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
  # rsem.genes.results as values
  
  tcga <- GDCprepare(query_TCGA)
  assy <- SummarizedExperiment::assay(tcga, "tpm_unstrand")
  
  #Available cancer type 
  assy = as.matrix(assy)

  rownames(assy) = tcga@rowRanges[, "gene_name"]$gene_name
  

  return(assy)
}


rt_info <-function(cancertype, barcodes){
  # Clinical data
  query <- GDCquery(
    project = paste0("TCGA-", cancertype),
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
  )
  GDCdownload(query)
  clinical_data <- GDCprepare(query)
  names(clinical_data)
  clinical_radio_therapy = as.data.frame(clinical_data[[paste0("clinical_radiation_",tolower(cancertype))]])
  clinical_radio_therapy = clinical_radio_therapy[3:nrow(clinical_radio_therapy),]
  
  # return the overlapping patients
  return(intersect(barcodes, clinical_radio_therapy$bcr_patient_barcode))
  
}

calculate_groups <- function(original_assay,gene_name){
  # K-means clusting to make four groups
  arid1a = unlist(original_assay[gene_name,])
  groups = kmeans(arid1a, centers=4, iter.max = 50)$cluster
  
  medians = unlist(lapply(unique(groups), function(x){median(arid1a[groups==x])}))
  
  groups = as.character(groups)
  group_num = 1
  for(g in unique(groups)[order(medians)]){
    groups[groups==g] = paste0("g", as.character(group_num))
    group_num = group_num+1
  }
  
  names(groups) = colnames(original_assay)
  return(groups)
}


deconvolution <- function(original_assay, cancertype, res_dir, gene_name, groups, method="timer"){
  # K-means w.r.t ARID1A expression
  
  assay = original_assay
  arid1a = assay[gene_name,]

  if(length(groups) != dim(assay)[2]){
    print("Dimensional error!")
    return()
  }
  
  # Run immune infiltration analysis
  if(method == "timer"){
    print("RUN timer")
    res_timer = immunedeconv::deconvolute(as.matrix(assay), "timer", rep(cancertype, ncol(assay)), tumor = TRUE)
  }else if(method == "xcell"){
    res_timer = immunedeconv::deconvolute(as.matrix(assay), "xcell", tumor = TRUE)
  }else if(method == "quantiseq"){
    res_timer = immunedeconv::deconvolute(as.matrix(assay), "quantiseq")
  }else if(method == "epic"){
    res_timer = immunedeconv::deconvolute(as.matrix(assay), "epic", tumor = TRUE)
  }else if(method == "mcp_counter"){
    res_timer = immunedeconv::deconvolute(as.matrix(assay), "mcp_counter", tumor = TRUE)
  } 
  
  res_timer = as.data.frame(res_timer)
  res_timer = res_timer[,c("cell_type", colnames(assay))]
  row.names(res_timer) = res_timer[,"cell_type"]
  res_timer = res_timer[,-c(1)]
  res_timer[gene_name,] = assay[gene_name,]
  
  res_timer = t(res_timer)
  res_timer=as.data.frame(res_timer)
  
  # integrate the result
  list_df <- lapply(rownames(res_timer), function(r){
    df_r = res_timer[r, ]
    data.frame(group=groups[r],#rep(groups[r], ncol(df_r)),
               infiltration=unlist(df_r),
               cell_types=colnames(df_r),
               sample=r,
               gene_expression=rep(assay[gene_name, r], ncol(df_r)))
  })
  
  df =  dplyr::bind_rows(list_df)
  print(head(df))
  
  # gene expression box plot
  p = ggplot(data=df[df$cell_types != gene_name,], aes(x=group , y=infiltration, fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter( size=0.1, color="black", alpha=0.2)+
    theme_bw()+
    ggtitle(paste0(cancertype, " (# subjects: ", as.character(nrow(res_timer)), ")"))+
    geom_signif(comparisons = list(c("g1", "g4")), 
                map_signif_level=TRUE, 
                vjust=1.5)+
    facet_wrap(~cell_types, ncol=6, scales = "free")
  ggsave(plot=p, filename =paste0(res_dir, "/", cancertype, "_",method, "_boxplot.pdf"), device = "pdf", dpi = 300, width = 15, height = 2 * (length(unique(df$cell_types))/6))
  
  # ARID1A expression between two groups 
  p = ggplot(data=df[df$cell_types == df$cell_types[1],], aes(x=group, y=gene_expression, fill=group))+
    geom_violin() +
    geom_boxplot(fill="white",width=0.1)+
    ggtitle(paste0(cancertype, " (# subjects: ", as.character(nrow(res_timer)), ")"))+
    #geom_point(position=position_jitterdodge(), size=0.5, color="black", alpha=0.8)+
    theme_bw()+
    geom_signif(comparisons = list(c("g1", "g4")), 
                map_signif_level=TRUE, 
                vjust=1)
  ggsave(plot=p, filename =paste0(res_dir, "/", cancertype, "_geneexp_boxplot.pdf"), device = "pdf", dpi = 300, width=5)
  
  return(df[df$cell_types != gene_name,])
  
}

# List to save all cancer type results
list_df <- list()

for(cancer_name in grep("TCGA", unique(GDCprojects$id), value = T)){
  cancertype = unlist(strsplit(cancer_name, "-"))[2]
  # if the result already exist
  if(cancertype  %in% unique(unlist(list_df["Tumour"]))){
    next
  }
  
  print(paste0("Analyse ", cancertype))
  
  # Collect gene expression data
  assy <- get_assay(cancertype)
  
  # Divide patients into four groups based on ARID1A expression
  groups <- calculate_groups(original_assay = assy, gene_name = "ARID1A")
  
  # Only select high and low expression of ARID1A
  groups = groups[groups %in% c("g1","g4")]
  assy <- assy[,names(groups)]
  
  # Estimate immune cell-type compositions 
  df = deconvolution(original_assay = assy, cancertype = cancertype, res_dir = save_dir, gene_name = "ARID1A", groups = groups, method = "xcell")
  
  df["Tumour"] = cancertype
  
  subset_assy = assy[target_genes,]
  subset_assy = data.frame(t(subset_assy))
  subset_assy = as.data.frame(lapply(subset_assy, as.numeric))
  subset_assy$group = groups
  rownames(subset_assy) = colnames(assy)
  
  df_subset_assy = subset_assy %>%  pivot_longer(values_to = "tpm", names_to = "gene", cols=!group)
  # gene expression box plot
  p=ggplot(df_subset_assy, aes(y=tpm, fill=group, x=group)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter( size=0.1, color="black", alpha=0.2)+
    facet_wrap(.~gene, scales = "free")+
    geom_signif(comparisons = list(c("g1", "g4")), 
                map_signif_level=TRUE, 
                vjust=1.5)+theme_bw()+
    ggtitle(paste0(cancertype, "(# subjects: ", as.character(dim(subset_assy)[1]), ")"))
  ggsave(plot=p, filename =paste0(save_dir, "/", cancertype, "_gex_boxplot.pdf"), device = "pdf", dpi = 300)
  
  # calculate correlation for the comparison of all cancer types
  corr_analysis = lapply(target_genes[2:length(target_genes)], function(x){
                                                            c = cor.test(subset_assy[["ARID1A"]], subset_assy[[x]],  method = "pearson")
                                                            return(data.frame(corr=c$estimate, p_val=c$p.value))
                                                          })
  corr_analysis = bind_rows(corr_analysis)
  corr_analysis$gene = target_genes[2:length(target_genes)]
  corr_analysis$tumour = cancertype
  
  barcodes = unlist(lapply(colnames(assy), function(x){
    x = unlist(strsplit(x,"-"))
    return(paste(x[1:3], collapse = "-"))}))
  
  # Get radio therapy information from TCGA
  rt <- rt_info(cancertype = cancertype, barcodes=barcodes)
  
  if(length(rt) > 0){
    # Analyse patients with radio therapy 
    
    # Find patients existing in both data sets
    rt_subset_assy = subset_assy[barcodes %in% rt,]
    
    df_subset_assy = rt_subset_assy %>%  pivot_longer(values_to = "tpm", names_to = "gene", cols=!group)
    
    # gene expression boxplot for rt patients
    p=ggplot(df_subset_assy, aes(y=tpm, fill=group, x=group)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter( size=0.1, color="black", alpha=0.2)+
      facet_wrap(.~gene, scales = "free")+
      geom_signif(comparisons = list(c("g1", "g4")), 
                  map_signif_level=TRUE, 
                  vjust=1.5)+theme_bw()+
      ggtitle(paste0(cancertype, "(# RT patients: ", as.character(dim(rt_subset_assy)[1]), ")"))
    ggsave(plot=p, filename =paste0(save_dir, "/rt_", cancertype, "_gex_boxplot.pdf"), device = "pdf", dpi = 300)
    
    sub_df = df[df$sample %in% rownames(rt_subset_assy),]
    # immune infiltration results for rt patients
    p = ggplot(data=sub_df[df$cell_types != "ARID1A",], aes(x=group , y=infiltration, fill=group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter( size=0.1, color="black", alpha=0.2)+
      theme_bw()+
      ggtitle(paste0(cancertype, " (# RT patients: ", as.character(length(unique(sub_df$sample))), ")"))+
      geom_signif(comparisons = list(c("g1", "g4")), 
                  map_signif_level=TRUE, 
                  vjust=1.5)+
      facet_wrap(~cell_types, ncol=6, scales = "free")
    ggsave(plot=p, filename =paste0(save_dir, "/rt_", cancertype, "_xcell_boxplot.pdf"), device = "pdf", dpi = 300, width = 15, height = 2 * (length(unique(df$cell_types))/6))
  }
  
  subset_assy$tumour = cancertype
  subset_assy$RT = barcodes %in% rt
  
  if(length(list_df) == 0){
    list_df = df
    list_corr_df = corr_analysis
    list_subset_df = subset_assy
  }else{
    list_df = rbind(list_df, df)
    list_corr_df = rbind(list_corr_df, corr_analysis)
    list_subset_df = rbind(list_subset_df, subset_assy)
  }
}

# Correlation heatmap
write.table(x = list_subset_df, file = paste0(save_dir, "/total_res_gex.csv"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(x = list_df, file = paste0(save_dir, "/total_res_infiltration.csv"), sep = "\t", quote = F, row.names = F, col.names = T)

# load those if you laready have files 
#list_df <- read.table(paste0(save_dir, "/total_res_infiltration.csv"), sep = "\t", header = T, stringsAsFactors = F)#, row.names = F, col.names = T)
#list_subset_df <- read.table( paste0(save_dir, "/total_res_gex.csv"), sep = "\t", header = T, stringsAsFactors = F)#, row.names = F, col.names = T)
#list_corr_df <- read.table(file = paste0(save_dir, "/total_cor_gex.csv"), sep = "\t", header = T, stringsAsFactors = F)

my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)

# Calculate the correlation between infiltration level and gene expresion in each cell type and each cancer type
corr_df = list_df %>% group_by(Tumour,cell_types) %>% summarise(corr = cor.test(infiltration, gene_expression)$estimate, 
                                                                p_val = cor.test(infiltration, gene_expression)$p.value)
corr_mat_pivot = pivot_wider(corr_df[,-c(4)], names_from = Tumour, values_from = corr)
corr_mat = as.matrix(corr_mat_pivot[,-c(1)])
rownames(corr_mat) = corr_mat_pivot$cell_types

p = pheatmap::pheatmap(corr_mat, main = "Correlation between ARID1A expr and infiltration")
ggsave(plot=p, filename =paste0(save_dir, "/total_corr_heatmap_infiltration.pdf"), device = "pdf", dpi = 300, width = 10)


corr_df[corr_df$p_val > 0.05, "corr"]=NA
p = ggplot(corr_df, aes(x=Tumour, y=cell_types, fill=corr))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-.8,.8))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("TCGA Tumour")+
  ylab("Immune cells")+
  ggtitle("Correlation between ARID1A expr and infiltration\n(grey - p.val > 0.05)")

ggsave(plot=p, filename =paste0(save_dir, "/total_corr_heatmap_infiltration_pvalue.pdf"), device = "pdf", dpi = 300, width = 10)

# Calculate correlation between selected genes (expression) and ARID1A in each cancer type
write.table(x = list_corr_df, file = paste0(save_dir, "/total_cor_gex.csv"), sep = "\t", quote = F, row.names = F, col.names = T)
a = pivot_wider(list_corr_df[,-c(2)], names_from = tumour, values_from = corr, )
a_gene_names = a$gene
a = as.matrix(a[,-c(1)])
rownames(a) = a_gene_names
p = pheatmap::pheatmap(a, main="Correlation between ARID1A expr and other gene expr")
ggsave(plot=p, filename =paste0(save_dir, "/total_corr_heatmap_gene_expression.pdf"), device = "pdf", dpi=300, width=10)

list_corr_df[list_corr_df$p_val > 0.05, "corr"]=NA
p = ggplot(list_corr_df, aes(x=tumour, y=gene, fill=corr))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-.5,.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("TCGA Tumour")+
  ylab("Genes")+
  ggtitle("Correlation between ARID1A expr and other gene expr\n(grey - p.val > 0.05)")

ggsave(plot=p, filename =paste0(save_dir, "/total_corr_heatmap_gene_expression_pval.pdf"), device = "pdf", dpi = 300, width = 10)

print(table(list_subset_df[list_subset_df$RT,"tumour"]))

