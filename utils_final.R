
library(here)

#==============================================================================
# function that grabs gene meta data from biomaRt
getGeneMetaData <- function(){
  
  library(biomaRt)
  
  # load tidy expression data
  reads_file = here("data","tidy_reads.csv")
  reads_data = read.csv(reads_file)
  
  # grab gene meta data
  ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
  gene_ids = reads_data$gene_name
  
  attributes2get = c("ensembl_gene_id","entrezgene_id","external_gene_name",
                     "description","chromosome_name","start_position",
                     "end_position","strand","band")
  
  # "transcript_start","transcript_end","transcription_start_site","transcript_length"
  # "hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name",
  # "hsapiens_homolog_ensembl_peptide","hsapiens_homolog_chromosome",
  # "hsapiens_homolog_chrom_start","hsapiens_homolog_chrom_end"
  
  gene_meta_data = getBM(attributes = attributes2get,
                         filters = 'mgi_symbol',
                         values = gene_ids,
                         mart = ensembl)
  write.csv(gene_meta_data, file = here("data","tidy_gene_meta_data.csv"))
  
} # function getGeneMetaData

#==============================================================================
# function that grabs mouse to human homolog gene meta data from biomaRt
getHumanHomolog <- function(){
  
  library(biomaRt)
  
  # load tidy expression data
  reads_file = here("data","tidy_reads.csv")
  reads_data = read.csv(reads_file)
  
  # grab gene meta data
  ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
  gene_ids = reads_data$gene_name
  
  attributes2get = c("ensembl_gene_id","external_gene_name",
                     "hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name",
                     "hsapiens_homolog_ensembl_peptide","hsapiens_homolog_chromosome",
                     "hsapiens_homolog_chrom_start","hsapiens_homolog_chrom_end")
  
  gene_meta_data = getBM(attributes = attributes2get,
                         filters = 'mgi_symbol',
                         values = gene_ids,
                         mart = ensembl)
  write.csv(gene_meta_data, file = here("data","tidy_human_homolog_gene_meta_data.csv"))
  
} # function getHumanHomolog


#==============================================================================
# function to run SVA
runSVA <- function(expr_data, meta_data, form2use, nsv, 
                   svestMethod = "leek"){
  
  library(sva)
  
  # specify models
  tmp_form = form2use
  full_model = model.matrix(tmp_form, data=meta_data)
  null_model = model.matrix(~1,data=meta_data)
  
  # estimate how many surrogate variables there should be
  n.sv = num.sv(as.matrix(expr_data), full_model, method = svestMethod)
  if (n.sv==0) {
    n.sv = nsv
    svobj = sva(as.matrix(expr_data), full_model, null_model, n.sv=n.sv)
  } else {
    svobj = sva(as.matrix(expr_data), full_model, null_model, n.sv=n.sv)
  } # if (n.sv==0) {
  
  sv_names = character()
  sv_form = character()
  for (isv in 1:n.sv){
    sv_names[isv] = sprintf("sv%03d",isv)
    if (isv == 1){
      sv_form = sv_names[isv]
    } else{
      sv_form = sprintf("%s + %s",sv_form,sv_names[isv])
    } # if (isv == 1){
  } # for (isv in 1:n.sv){
  colnames(svobj$sv) = sv_names
  
  meta_data = cbind(meta_data, svobj$sv)
  tmp = as.character(form2use)
  new_form = as.formula(sprintf("~ %s + %s",tmp[2], sv_form))

  res = list(meta_data = meta_data, nsv = n.sv, 
             form2use = new_form, sv_names = sv_names)
  return(res)
  
} # function runSVA


#==============================================================================
# function to run DE modeling
runDEmodel <- function(expr_data, meta_data, form2use, gene_meta_data){
  
  library(limma)
  library(WGCNA)
  
  # Construct contrast matrices 
  full_model = model.matrix(form2use, data = meta_data)
  contrast.matrix = cbind(group = as.numeric(colnames(full_model)=="groupWT"),
                          sex = as.numeric(colnames(full_model)=="sexM"),
                          interaction = as.numeric(colnames(full_model)=="groupWT:sexM"))
  
  # fit DE limma model 
  fit_DE = lmFit(expr_data,full_model)
  fit_coef = fit_DE$coefficients
  
  # fit contrasts
  fitContrasts = contrasts.fit(fit_DE, contrast.matrix)
  
  # empirical Bayes
  eb_res = eBayes(fitContrasts)
  
  # F-stats table
  Fstat_table = topTable(eb_res, number=dim(eb_res)[1])
  # compute Storey FDR as the adj.P.Val
  Fstat_table$adj.P.Val = qvalue(as.numeric(Fstat_table$P.Value))$qvalues
  Fstat_table$ensembl_id = rownames(Fstat_table)
  # add gene meta data info to table
  gene_meta_data_sub = subset(gene_meta_data, is.element(gene_meta_data$ensembl_gene_id, Fstat_table$ensembl_id))
  Fstat_table = merge(Fstat_table, gene_meta_data_sub, by.x = "ensembl_id", by.y = "ensembl_gene_id")
  # # add human homolog gene meta data info to table
  # hsap_gene_meta_data_sub = subset(hsap_homolog_gene_meta_data, 
  #                                  is.element(hsap_homolog_gene_meta_data$ensembl_gene_id, Fstat_table$ensembl_id))
  # Fstat_table = merge(Fstat_table, hsap_gene_meta_data_sub, by.x = "ensembl_id", by.y = "ensembl_gene_id")

  # Group*Sex Interaction Effect table
  interaction_table = topTable(eb_res, coef="interaction", number=dim(eb_res)[1])
  # compute Storey FDR as the adj.P.Val
  interaction_table$adj.P.Val = qvalue(as.numeric(interaction_table$P.Value))$qvalues
  interaction_table$ensembl_id = rownames(interaction_table)
  # add gene meta data info to table
  gene_meta_data_sub = subset(gene_meta_data, is.element(gene_meta_data$ensembl_gene_id, interaction_table$ensembl_id))
  interaction_table = merge(interaction_table, gene_meta_data_sub, by.x = "ensembl_id", by.y = "ensembl_gene_id")
  
  # Group Main Effect table
  group_table = topTable(eb_res, coef="group", number=dim(eb_res)[1])
  # compute Storey FDR as the adj.P.Val
  group_table$adj.P.Val = qvalue(as.numeric(group_table$P.Value))$qvalues
  group_table$ensembl_id = rownames(group_table)
  # add gene meta data info to table
  gene_meta_data_sub = subset(gene_meta_data, is.element(gene_meta_data$ensembl_gene_id, group_table$ensembl_id))
  group_table = merge(group_table, gene_meta_data_sub, by.x = "ensembl_id", by.y = "ensembl_gene_id")

  # Sex Main Effect table
  sex_table = topTable(eb_res, coef="sex", number=dim(eb_res)[1])
  # compute Storey FDR as the adj.P.Val
  sex_table$adj.P.Val = qvalue(as.numeric(sex_table$P.Value))$qvalues
  sex_table$ensembl_id = rownames(sex_table)
  # add gene meta data info to table
  gene_meta_data_sub = subset(gene_meta_data, is.element(gene_meta_data$ensembl_gene_id, sex_table$ensembl_id))
  sex_table = merge(sex_table, gene_meta_data_sub, by.x = "ensembl_id", by.y = "ensembl_gene_id")

  # pack into list to return as output  
  res2use = list(Fstat_table = Fstat_table, 
                 interaction_table = interaction_table, 
                 group_table = group_table, 
                 sex_table = sex_table,
                 fit_coef = fit_coef,
                 full_model = full_model)
  return(res2use)
  
} # function runDEmodel

#==============================================================================
# function to run Chromosome permutation analysis
chromosome_perm <- function(gene_meta_data_sub, de_df, 
                            nperm=10000, 
                            fdr_qThresh=0.05,
                            rand_seed = 999){
  
  set.seed(rand_seed)
  
  gene_info = data.frame(table(gene_meta_data_sub$chromosome_name))
  colnames(gene_info) = c("Chromosome","total")
  
  tmp_df = subset(de_df, de_df$adj.P.Val<fdr_qThresh)
  
  chr_res = data.frame(table(tmp_df$chromosome_name))
  colnames(chr_res) = c("Chromosome","freq")
  chr_res = merge(chr_res, gene_info, by.x = "Chromosome", by.y = "Chromosome")
  chr_res$percent = chr_res$freq/chr_res$total
  
  mask = is.element(gene_meta_data$ensembl_gene_id, tmp_df$ensembl_id)
  n_de_genes = dim(tmp_df)[1]
  
  final_res = data.frame(matrix(nrow = dim(chr_res)[1], ncol = nperm+2))
  rownames(final_res) = chr_res$Chromosome
  perm_labels = character()
  for (iperm in 1:nperm){
    perm_labels[iperm] = sprintf("perm%09d",iperm)
  }
  colnames(final_res) = c("Chromosome","real_pct",perm_labels) 
  final_res[,"Chromosome"] = chr_res$Chromosome
  final_res[,"real_pct"] = chr_res$percent
  
  for (iperm in 1:nperm){
    rand_perm = sample(dim(gene_meta_data_sub)[1])
    idx2use = rand_perm[1:n_de_genes]
    perm_df = gene_meta_data_sub[idx2use,]
    perm_res = data.frame(table(perm_df$chromosome_name))
    colnames(perm_res) = c("Chromosome","freq")
    perm_res = merge(perm_res, gene_info, by.x = "Chromosome", by.y = "Chromosome")
    perm_res$percent = perm_res$freq/perm_res$total
    perm_res = subset(perm_res, perm_res$Chromosome!="Y")
    
    final_res[,iperm+2] = perm_res$percent
  }
  
  res2use = final_res[,c("Chromosome","real_pct")]
  colnames(res2use)[2] = "percent"
  res2use$pval = (rowSums(final_res[,2:ncol(final_res)]>=res2use$percent, na.rm=TRUE))/(nperm+1)
  res2use$fdr = p.adjust(res2use$pval, method = "fdr")
  
  return(res2use)
  
} # function chromosome_perm


#==============================================================================
# function to adjust expression data for covariates
adjExprData <- function(expr_data, nuisance_covnames, coeffs, full_model){
  
  coeffs2use = coeffs[,nuisance_covnames]
  model2use = full_model[,nuisance_covnames]
  expr_data_adj = expr_data - coeffs2use %*% t(model2use)
  return(expr_data_adj)
  
} # function adjExprData


#==============================================================================
# function to run gene set enrichment analysis
#
# genelistOverlap.R
#
# Calculate enrichment odds ratio and p-value from hypergeometric test to
# answer the question of whether genes from one list are enriched in genes
# from another list.
#
# INPUT
#	list1 or list2 = excel file, tab or comma delimited file with gene IDs
#					 assuming each list has a header.
#	backgroundTotal = specify number of genes in background pool
#
# mvlombardo 21.12.2016

genelistOverlap <- function(list1,list2,backgroundTotal, print_result = TRUE, header = FALSE) {
  
  # Read in libraries and set options
  options(stringsAsFactors = FALSE)
  require(readxl)
  require(tools)
  
  if (is.character(list1)){
    # get the file extension of list1
    ext1 = file_ext(list1)
    
    if (is.element(ext1,c("xlsx","xls","txt","csv"))){
      if (ext1=="xlsx" | ext1=="xls") {
        genes1 = read_excel(list1)
      } else if (ext1=="txt") {
        genes1 = read.delim(list1, header = header)
      } else if (ext1=="csv") {
        genes1 = read.csv(list1, header = header)
      }# if (ext1=="xlsx" | ext1=="xls") {
    } else {
      genes1 = data.frame(list1)
    }# if (is.element(ext1,c("xlsx","xls","txt","csv"))){
  } else if (is.data.frame(list1)){
    genes1 = list1
  }# if (is.character(list1)){
  
  if (is.character(list2)){
    # get the file extension of list1
    ext2 = file_ext(list2)
    
    if (is.element(ext2,c("xlsx","xls","txt","csv"))){
      if (ext2=="xlsx" | ext2=="xls") {
        genes2 = read_excel(list2)
      } else if (ext2=="txt") {
        genes2 = read.delim(list2, header = header)
      } else if (ext1=="csv") {
        genes2 = read.csv(list2, header = header)
      }# if (ext1=="xlsx" | ext1=="xls") {
    } else {
      genes2 = data.frame(list2)
    }# if (is.element(ext1,c("xlsx","xls","txt","csv"))){
  } else if (is.data.frame(list2)){
    genes2 = list2
  }# if (is.character(list2)){
  
  # Find overlapping genes
  gene_mask = is.element(genes1[,1],genes2[,1])
  overlapping_genes = genes1[gene_mask,1]
  gene_overlap = sum(gene_mask)
  ngenes1 = length(genes1[,1])
  ngenes2 = length(genes2[,1])
  
  # Calculate odds ratio
  A = gene_overlap;
  B = ngenes1-gene_overlap
  if (ngenes2==gene_overlap){
    # add 0.5 to ngenes2 to avoid OR = Inf
    C = (ngenes2+0.5)-gene_overlap
  } else {
    C = ngenes2-gene_overlap
  }
  D = backgroundTotal-C
  OR = (A*D)/(B*C)
  
  # Calculate p-value from hypergeometric test
  hypergeo_p = sum(dhyper(gene_overlap:ngenes2,ngenes1,backgroundTotal-ngenes1,ngenes2))
  
  # pack into result
  result = vector(mode = "list", length = 1)
  result[[1]]$list1 = list1
  result[[1]]$list2 = list2
  result[[1]]$backgroundTotal = backgroundTotal
  result[[1]]$OR = OR
  result[[1]]$hypergeo_p = hypergeo_p
  result[[1]]$percent_overlap_list1 = gene_overlap/ngenes1
  result[[1]]$percent_overlap_list2 = gene_overlap/ngenes2
  result[[1]]$gene_overlap = gene_overlap
  result[[1]]$ngenes1 = ngenes1
  result[[1]]$ngenes2 = ngenes2
  result[[1]]$overlapping_genes = overlapping_genes
  
  # print result to the screen and then return result
  if (print_result){
    print(sprintf("OR = %f, p = %f",OR,hypergeo_p))
  }
  return(result)
} # function genelistOverlap 


#==============================================================================
# function to Allen scRNA cell type enrichment results
#
# see this paper: https://www.biorxiv.org/content/10.1101/2020.03.30.015214v1
# allen_meta table is from https://www.biorxiv.org/content/biorxiv/early/2020/03/31/2020.03.30.015214/DC3/embed/media-3.csv?download=true
#
annotate_celltypes <- function(celltype_data){
  
  allen_meta = read.csv(here("data","allen_celltype_metadata.csv"))
  mouse_terms = grep("Mouse", celltype_data$Term, value = TRUE)
  
  celltype_data = subset(celltype_data, is.element(celltype_data$Term, mouse_terms))
  
  up_terms = grep(" up", celltype_data$Term, value = TRUE)
  celltype_data_sub = subset(celltype_data, is.element(celltype_data$Term, up_terms))
  celltype_data_sub$cluster_id = NA
  celltype_data_sub$class_label = NA
  celltype_data_sub$subclass_label = NA
  celltype_data_sub$supertype_label = NA
  celltype_data_sub$neighborhood_label = NA
  celltype_data_sub$final_label = NA
  celltype_data_sub$final_label2 = NA
  celltype_data_sub$plotcolor2use = NA
  
  for (i in 1:dim(allen_meta)[1]){
    
    mask = is.element(celltype_data_sub$Term, grep(sprintf("Mouse %d ",allen_meta[i,"cluster_id"]), celltype_data_sub$Term, value=TRUE))
    celltype_data_sub$cluster_id[mask] = allen_meta[i,"cluster_id"]
    celltype_data_sub$class_label[mask] = allen_meta[i,"class_label"]
    celltype_data_sub$subclass_label[mask] = allen_meta[i,"subclass_label"]
    celltype_data_sub$supertype_label[mask] = allen_meta[i,"supertype_label"]
    celltype_data_sub$neighborhood_label[mask] = allen_meta[i,"neighborhood_label"]
    celltype_data_sub$plotcolor2use[mask] = allen_meta[i,"plotcolor2use"]
    celltype_data_sub$final_label[mask] = sprintf("%s - %s - %s - %s",
                                                  allen_meta[i,"class_label"],
                                                  allen_meta[i,"neighborhood_label"],
                                                  allen_meta[i,"subclass_label"],
                                                  allen_meta[i,"supertype_label"])
    celltype_data_sub$final_label2[mask] = sprintf("%s %s %s %d",
                                                  allen_meta[i,"class_label"],
                                                  allen_meta[i,"subclass_label"],
                                                  allen_meta[i,"supertype_label"],
                                                  allen_meta[i,"cluster_id"])
                                                  # allen_meta[i,"supertype_num"])
# celltype_data_sub$final_label[mask] = sprintf("%s - %s - %s",
    #                                               allen_meta[i,"neighborhood_label"],
    #                                               allen_meta[i,"subclass_label"],
    #                                               allen_meta[i,"supertype_label"])
  } # for (i in 1:dim(allen_meta)[1]){

  celltype_data_sub = celltype_data_sub[,c("X","cluster_id","Term","class_label",
                                           "neighborhood_label","subclass_label",
                                           "supertype_label","final_label",
                                           "final_label2","Overlap","P.value",
                                           "Adjusted.P.value","Old.P.value","Old.Adjusted.P.value",
                                           "Odds.Ratio","Combined.Score","Genes","plotcolor2use")]
  return(celltype_data_sub)

} # function annotate_celltypes


#==============================================================================
# function to get colors from ggplot colorwheel
get_ggColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} # function get_ggColorHue

#==============================================================================



#==============================================================================
# function to get colors from ggplot colorwheel
volcanoPlot <- function(data2use, fname2save, fdr_qThresh, plotTopGenes=FALSE, 
                        colors2use=NULL, fontSize=40, dotSize=1, 
                        plotWidth=10, plotHeight=8) {
  
  # dotSize = 1
  # fontSize = 40
  # data2use = data
  
  data2use$sig = "NS"
  data2use$sig[data2use$adj.P.Val<=fdr_qThresh] = "SIG"
  data2use$fc_thresh = "NS"
  data2use$fc_thresh[abs(data2use$logFC)>=1] = "SIG"
  
  data2use$col4plot = "NS"
  data2use$col4plot[data2use$sig=="SIG" & data2use$logFC>0] = "Pos"
  data2use$col4plot[data2use$sig=="SIG" & data2use$logFC<0] = "Neg"
  data2use$col4plot = factor(data2use$col4plot, levels = c("Pos","NS","Neg"))
  
  data2use$delabel = NA
  mask = (-log10(data2use$P.Value)>3 & data2use$sig=="SIG") | (data2use$logFC>1 & data2use$sig=="SIG") 
  data2use$delabel[mask] <- data2use$external_gene_name[mask]
  
  p = ggplot(data = data2use, aes(x = logFC, y = -log10(P.Value), colour = col4plot, label=delabel))
  
  if (plotTopGenes){
    p = p + geom_point(size=dotSize) + geom_text_repel() + xlab("log2FC") + ylab("-log10(p)") 
  }
  
  p = p + geom_point(size=dotSize) + xlab("log2(FC)") + ylab("-log10(p)") + 
    geom_hline(yintercept = -log10(max(data2use$P.Value[data2use$sig=="SIG"])), linetype = "dashed") + 
    # geom_vline(xintercept = 1, linetype = "dashed") + geom_vline(xintercept = -1, linetype = "dashed") + 
    guides(colour=FALSE)  +
    theme(axis.text.x = element_text(size=fontSize),
          axis.text.y = element_text(size=fontSize),
          axis.title.x = element_text(size=fontSize),
          strip.text.x = element_text(size=fontSize),
          axis.title.y = element_text(size=fontSize))
  
  if (!is.null(colors2use)){
    p = p + scale_colour_manual(values = colors2use)
  }
  ggsave(filename = fname2save, width=plotWidth, height=plotHeight) 
  return(p)
} # function volcanoPlot

