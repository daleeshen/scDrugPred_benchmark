library(tidyverse)

mapping_file<-"/media/shenli/Expansion/scDrug_benchmark/scDEAL/scDEAL_dataset_drug.csv"
out_root<-"/media/shenli/Expansion/DRMREF/LLM/scGPT-main/scdata"

rounds<-1:3
cell_counts<-c(500,1000,2000,3000,5000)
ratios<-c("1_3","1_10","1_30","1_100")

mapping<-read.csv(mapping_file,check.names=FALSE,stringsAsFactors=FALSE)
datasets<-unique(mapping$Dataset)

make_scgpt_input<-function(count_file,meta_file,out_file,dataset,scenario){
  count<-read.csv(count_file,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  meta<-read.csv(meta_file,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  
  if(!"condition"%in%colnames(meta)) stop("No condition column: ",meta_file)
  
  # metadata cell names统一 . -> -
  rownames(meta)<-gsub("\\.","-",rownames(meta))
  
  # 分别测试count的行/列哪一个更像cell
  row_ids<-gsub("\\.","-",rownames(count))
  col_ids<-gsub("\\.","-",colnames(count))
  
  row_match<-sum(row_ids%in%rownames(meta))
  col_match<-sum(col_ids%in%rownames(meta))
  
  # 如果cell在列：只修改列名，然后转置
  if(col_match>row_match){
    colnames(count)<-col_ids
    count<-as.data.frame(t(count),check.names=FALSE)
    
    # 如果cell在行：只修改行名
  }else if(row_match>0){
    rownames(count)<-row_ids
    
  }else{
    stop("No matched cells: ",dataset," | ",count_file)
  }
  
  common<-intersect(rownames(meta),rownames(count))
  if(length(common)==0) stop("No matched cells after normalization: ",dataset)
  
  common<-rownames(meta)[rownames(meta)%in%common]
  meta<-meta[common,,drop=FALSE]
  count<-count[common,,drop=FALSE]
  
  out<-data.frame(Cell_barcode=common,Condition=meta$condition,count,check.names=FALSE)
  
  dir.create(dirname(out_file),recursive=TRUE,showWarnings=FALSE)
  write.table(out,out_file,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  
  cat("DONE:",scenario,"|",dataset,"| cells =",nrow(out),"| genes =",ncol(count),"\n")
  
  tibble(dataset=dataset,scenario=scenario,n_cells=nrow(out),n_genes=ncol(count),output=out_file)
}

res<-list()

# ============================================================
# Balanced
# ============================================================

for(n in rounds){
  for(cell_n in cell_counts){
    
    # 1000 cells
    if(cell_n==1000){
      data_dir<-paste0("/media/shenli/Expansion/DRMREF/",n,"_round_raw_1000cells/data_for_SCAD")
      count_prefix<-"balanced_sampled_"
      meta_prefix<-"balanced_sampled_"
      
      # 500 / 2000 / 3000 / 5000 cells
    }else{
      data_dir<-paste0("/media/shenli/Expansion/DRMREF/",n,"_round_raw_",cell_n,"cells/metadata_count")
      count_prefix<-""
      meta_prefix<-""
    }
    
    for(dataset in datasets){
      count_file<-file.path(data_dir,paste0(count_prefix,dataset,"_raw_count_matrix.csv"))
      meta_file<-file.path(data_dir,paste0(meta_prefix,dataset,"_meta_data.csv"))
      out_file<-file.path(out_root,"balanced",cell_n,n,paste0(dataset,".txt"))
      
      if(!file.exists(count_file) || !file.exists(meta_file)){
        cat("SKIP balanced:",cell_n,"| round",n,"|",dataset,"\n")
        next
      }
      
      res[[length(res)+1]]<-make_scgpt_input(
        count_file,meta_file,out_file,dataset,
        paste0("balanced_",cell_n,"_round",n)
      )
    }
  }
}

# ============================================================
# Imbalanced
# ============================================================

for(n in rounds){
  for(ratio in ratios){
    
    data_dir<-paste0(
      "/media/shenli/Expansion/DRMREF/unbalanced_original/",
      n,"_round_unbalanced_sampling/",ratio,"/metadata_count"
    )
    
    for(dataset in datasets){
      count_file<-file.path(data_dir,paste0(dataset,"_raw_count_matrix.csv"))
      meta_file<-file.path(data_dir,paste0(dataset,"_meta_data.csv"))
      out_file<-file.path(out_root,"imbalance",ratio,n,paste0(dataset,".txt"))
      
      if(!file.exists(count_file) || !file.exists(meta_file)){
        cat("SKIP imbalance:",ratio,"| round",n,"|",dataset,"\n")
        next
      }
      
      res[[length(res)+1]]<-make_scgpt_input(
        count_file,meta_file,out_file,dataset,
        paste0("imbalance_",ratio,"_round",n)
      )
    }
  }
}

summary_df<-bind_rows(res)
write.csv(summary_df,file.path(out_root,"scGPT_input_summary.csv"),row.names=FALSE)

print(summary_df)
cat("\nTotal generated files:",nrow(summary_df),"\n")