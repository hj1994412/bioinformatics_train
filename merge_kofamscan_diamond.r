#!/share/home/zhiwenl/miniconda3/envs/meta/bin/Rscript
# auther: huangjing
# use code: Rscript merge_kofam_KEGG.r <work_path> <redun_kofam_file_name> <KEGG_file_name> <out_path> <kofam_strict:TRUE or FALSE> <pident_thres>

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)
if (!requireNamespace("vroom", quietly = TRUE)) {
  install.packages("vroom")
}
library(vroom)
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

# parameters
option_list <- list(
  make_option(c("-d", "--work_dir"), 
              type="character", default=".", dest="work_dir",
              help="Working directory, used as a parameter for the setwd().", metavar="character"),
  make_option(c("-f", "--file_path"), 
              type="character", default="anno_file_path", dest="file_path",
              help="kofam detail-tsv filter file and/or diamond KEGG file path.\n
              The kofam file path is labelled as 'kofam_path', while the diamond KEGG file path is labelled as 'KEGG_path'.\n
              The kofamscan output file must have column names and must include 'score' and 'thrshld'.\n
              The diamond KEGG output file must have column names and must include 'pident'.", 
              metavar="character"),
  make_option(c("-m", "--max_pident_thres"), 
              type="double", default="90", dest="max_pident_thres",
              help="The max pident threshold used for filtering KEGG annotation results.\n
              The default value is 90.", metavar="double"), 
  make_option(c("-e", "--min_pident_thres"), 
              type="double", default="80", dest="min_pident_thres",
              help="The min pident threshold used for filtering KEGG annotation results.\n
              The default value is 80.", metavar="double"),                        
  make_option(c("-o", "--out_path"), 
              type="character", default=".",  dest="out_path",
              help="output directory", metavar="character"),           
  make_option(c("-s", "--suffix"), 
              type="character", default="merge.filter",  dest="suffix",
              help="output file suffix.\n
                    The default value is 'merge.filter'.", 
              metavar="character")                      
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
args

work_dir <- args$work_dir
file_path <- args$file_path
out_path <- args$out_path
max_pident_thres <- args$max_pident_thres
min_pident_thres <- args$min_pident_thres
suffix  <- args$suffix

setwd(work_dir)
dir.create(out_path)

hmm_diamond_merge <- function(diamond_file_path, hmm_file_path,min_pident_thres,max_pident_thres,suffix,out_path... ) {
## 两种注释方法取交集
    cat(paste(basename(diamond_file_path),"\t",basename(hmm_file_path),"\n",sep = ""));
    temp_kofam <- vroom(hmm_file_path, col_names = TRUE, delim = "\t");
    temp_kegg <- vroom(diamond_file_path, col_names = TRUE, delim = "\t");

    temp_kofam %>%
        select(gene_name,KO) %>%
        inner_join((temp_kegg %>% select(gene_name,KO)),by = c("gene_name","KO")) %>%
        mutate(source = "intersection") -> temp_inner
#temp_inner %>% dim
#temp_inner %>% head

    cat(paste("The intersection of",basename(diamond_file_path),"and",basename(hmm_file_path),"contains", (temp_inner %>% nrow), "genes.\n"))

## 剩余kofam与diamond KEGG注释结果的合并过滤
### diamond KEGG中不存在注释结果，但kofamscan中存在注释结果的gene。
    temp_kofam %>%
        select(gene_name) %>%
        anti_join((temp_kegg %>% select(gene_name)),by=c("gene_name")) -> temp_kofam.uniq
    temp_kofam.uniq %>%
        filter(duplicated(gene_name)) %>% nrow -> kofam_duplicated_num
    cat(paste("Have",dim(temp_kofam.uniq)[1],"gene not in diamond KEGG annotation results.\nHave",kofam_duplicated_num,"duplicated gene.\n"))
### kofamscan中不存在注释结果，但diamond KEGG中存在注释结果的gene。
    temp_kegg %>%
        select(gene_name) %>%
        anti_join((temp_kofam %>% select(gene_name)),by=c("gene_name")) -> temp_kegg.uniq    
    temp_kegg.uniq  %>%
        filter(duplicated(gene_name)) %>% nrow -> diamond_duplicated_num
    cat(paste("Have",dim(temp_kegg.uniq)[1],"gene not in Kofamscan KEGG annotation results.\nHave",diamond_duplicated_num,"duplicated gene.\n"))
### 两个注释的余集结果都保留，如果有存在重复注释的gene则删除。
    temp_inner %>%
        bind_rows((temp_kofam %>% 
               select(gene_name,KO) %>% 
               filter(gene_name %in% temp_kofam.uniq$gene_name,
                      !duplicated(gene_name)) %>%
                mutate(source = "only_hmm_kofam"))) %>%
        bind_rows((temp_kegg %>% 
               select(gene_name,KO) %>% 
               filter(gene_name %in% temp_kegg.uniq$gene_name,
                      !duplicated(gene_name)) %>%
                mutate(source = "only_diamond_KEGG"))) %>%
        distinct() -> temp_combine
    #temp_combine %>% head
    #(temp_combine %>% nrow() -> temp_combine.nrows)

### 将kofamscan与diamond KEGG都存在注释结果，且注释结果不同的gene合并。
    temp_kofam %>%
        filter(!gene_name %in% temp_combine$gene_name,
           !gene_name %in% temp_kofam.uniq$gene_name) %>%
        select(-KO_definition) %>%
        rename(hmm_KO = KO) %>%
        full_join((temp_kegg %>%
               filter(!gene_name %in% temp_combine$gene_name,
                      !gene_name %in% temp_kegg.uniq$gene_name) %>%
               select(-KO_definition) %>%
               rename(diamond_KO = KO)),relationship = "many-to-many",by = "gene_name") -> temp_multi
    temp_multi$gene_name %>% unique %>% length -> uncertain_gene_num
    cat(paste("The annotation result of",uncertain_gene_num,"genes is uncertain.\n"))

### 根据注释质量指数以及提供的KO list进行过滤。
#### 1. kofamscan的score >= thrshld且diamond KEGG的pident < 98.65)
### KO list可以是prokaryote.hal或者eukaryote.hal中的Knum。
#max_pident_thres = 90
#min_pident_thres = 50
    temp_multi %>%
        mutate(KO = if_else((thrshld !=0 & score >= thrshld & pident < max_pident_thres),
                        hmm_KO, NA),
           source = if_else((thrshld !=0 & score >= thrshld & pident < max_pident_thres),
                        paste("score >= thrshld & pident <",max_pident_thres), NA)) %>%
        mutate(KO = if_else((pident >= max_pident_thres),
                        diamond_KO, KO),
           source = if_else((pident >= max_pident_thres),
                        paste("pident >=",max_pident_thres), source)) %>%         
        mutate(KO = if_else((thrshld == 0 & (min_pident_thres <= pident & pident < max_pident_thres)),
                        diamond_KO, KO),
           source = if_else((thrshld == 0 & (min_pident_thres <= pident & pident < max_pident_thres)),
                        paste("thrshld == 0 &",min_pident_thres,"<= pident <", max_pident_thres), source)) %>%     
        mutate(KO = if_else((thrshld != 0 & (max_pident_thres > pident & pident >= min_pident_thres)),
                        diamond_KO, KO),
           source = if_else((thrshld != 0 & (max_pident_thres> pident & pident >= min_pident_thres)),
                        paste("thrshld != 0 &", min_pident_thres,"<= pident <", max_pident_thres), source)) %>%
        filter(!is.na(source)) %>%
        select(gene_name,KO,source) %>%
        distinct() %>%
        filter(!duplicated(gene_name)) %>%
        bind_rows(temp_combine) -> temp_merge_res           

    cat("The output contains",(temp_merge_res$gene_name %>% unique %>% length),"genes.\n",
"Each gene has only one annotation: ",(temp_merge_res %>% nrow()),"=", (temp_merge_res$gene_name %>% unique %>% length),"\n")

    vroom_write(
        temp_merge_res,
        file.path(out_path,paste(gsub("\\.KEGG_filter","",basename(diamond_file_path)),suffix,sep="_")),
        col_names = TRUE,
        delim = "\t",
        quote = "none")

    temp_multi %>%
        filter(!gene_name %in% temp_merge_res$gene_name) -> uncertain_gene_anno

    cat("There are",(uncertain_gene_anno$gene_name %>% unique %>% length),"genes with unresolved annotation results.\n")

    vroom_write(
        uncertain_gene_anno,
        file.path(out_path,paste(gsub("\\.KEGG_filter","",basename(diamond_file_path)),"uncertain.gene.anno",sep="_")),
        col_names = TRUE,
        delim = "\t",
        quote = "none")
}


file_path <- read.table(
    file_path,
    header = TRUE,
    colClasses = "character",
    sep = "\t"
)
for (i in 1:nrow(file_path)) {
    diamond_file_path=file_path[i,"KEGG_path"];
    hmm_file_path=file_path[i,"kofam_path"];
    hmm_diamond_merge(diamond_file_path, hmm_file_path,min_pident_thres,max_pident_thres,suffix)
}






