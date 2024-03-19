#!/share/home/weizhuang/miniconda3/envs/instrain/bin/Rscript

# use code: Rscript redun_kofam.r <work_path> <redun_kofam_file_name> <out_path>
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
              help="Working directory", metavar="character"),
  make_option(c("-f", "--file_path"), 
              type="character", default="KEGG_anno_name", dest="file_path",
              help="kofam detail-tsv file path", metavar="character"),
  make_option(c("-o", "--out_path"), 
              type="character", default=".",  dest="out_path",
              help="output directory", metavar="character"),
  make_option(c("-p", "--pident_thres"), 
              type="double", default="30", dest="pident_thres",
              help="The pident threshold used for filtering diamond KEGG annotation results.\n
              The default value is 30.", metavar="double"), 
  make_option(c("-e", "--evalue_thres"), 
              type="double", default="1e-5", dest="evalue_thres",
              help="The evalue threshold used for filtering diamond KEGG annotation results.\n
              The default value is 1e-5.", metavar="double"),                
  make_option(c("-s", "--suffix"), 
              type="character", default="filter",  dest="suffix",
              help="output file suffix", 
              metavar="character"),
  make_option(c("-g", "--GeneBank"), 
              type="logical", default="FALSE",  dest="GeneBank",
              help="Keep GeneBank information.\n
              The result will definitely contain duplicate 'gene_name'.", 
              metavar="character"),             
  make_option(c("-l", "--log"), 
              type="character", default="diamond.filter.log",  dest="log",
              help="filter log file name", 
              metavar="character")                
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
args
#args <- commandArgs(trailingOnly = TRUE)
work_dir <- args$work_dir
file_path <- args$file_path
out_path <- args$out_path
GeneBank  <- args$GeneBank
suffix  <- args$suffix
evalue_thres  <- args$evalue_thres
pident_thres  <- args$pident_thres
log  <- args$log
dir.create(out_path)


setwd(work_dir)
filenames <- read.table(file_path, header=FALSE)
colnames(filenames) <- "name"

only_KO  <- function(filename,...) {
    # import data
    temp_kegg <- vroom(
    filename,
    col_names = TRUE,
    delim = "\t");
    # filter data
    temp_kegg %>%
        filter(evalue <= evalue_thres,
               pident >= pident_thres) %>%
        group_by(qtitle) %>% 
        filter(bitscore  == max(bitscore)) %>%    
        filter(pident == max(pident)) %>%
        filter(evalue == min(evalue)) %>%
        filter(mismatch == min(mismatch)) %>%
        filter(gapopen == min(gapopen)) %>%
        ungroup() %>% 
        extract(col = qtitle, into = "gene_name", regex = "(\\S+)") %>%
        extract(col = stitle, into = c("anno_gene","definition"), regex = "(\\S+)\\s(.*)") %>%
        separate(col = definition, into = c("KEGG","GenBank"),sep = " \\| ") %>%
        mutate(GenBank  = str_remove(GenBank,"\\(GenBank\\) ")) %>%
        extract(col = KEGG,into = c("KO","KO_definition"),regex = "(\\S+)\\s(.*)") %>%
        filter( !KO %in% "no") %>% 
        #mutate(KO_definition = if_else(KO == "no", GenBank, KO_definition)) %>%
        #mutate(GenBank = if_else(KO == "no", GenBank, KO_definition)) %>%
        select("gene_name","KO","KO_definition",
                evalue, pident, bitscore, mismatch,
                gapopen,length) %>%
        distinct() %>%
        group_by(gene_name,KO) %>% 
        slice(1) %>%
        ungroup() -> filtered_data
    res <- list(filtered_data=filtered_data,
        redun=data.frame(
            filename = basename(filename),
            KO_num = nrow(filtered_data),
            gene_num = length(unique(filtered_data$gene_name))))
    return(res);       
    }


keep_genebank  <- function(filename,...) {
    # import data
    temp_kegg <- vroom(
    filename,
    col_names = TRUE,
    delim = "\t");
    # filter data
    temp_kegg %>%
        filter(evalue <= evalue_thres,
               pident >= pident_thres) %>%
        group_by(qtitle) %>% 
        filter(bitscore  == max(bitscore)) %>%    
        filter(pident == max(pident)) %>%
        filter(evalue == min(evalue)) %>%
        filter(mismatch == min(mismatch)) %>%
        filter(gapopen == min(gapopen)) %>%
        extract(col = qtitle, into = "gene_name", regex = "(\\S+)") %>%
        extract(col = stitle, into = c("anno_gene","definition"), regex = "(\\S+)\\s(.*)") %>%
        separate(col = definition, into = c("KEGG","GenBank"),sep = " \\| ") %>%
        mutate(GenBank  = str_remove(GenBank,"\\(GenBank\\) ")) %>%
        extract(col = KEGG,into = c("KO","KO_definition"),regex = "(\\S+)\\s(.*)") %>%
        mutate(KO_definition = if_else(KO == "no", GenBank, KO_definition)) %>%
        select("gene_name","KO","KO_definition",GenBank,
                evalue, pident, bitscore, mismatch,
                gapopen,length) -> filtered_data
    res <- list(filtered_data=filtered_data,
        redun=data.frame(
            filename = basename(filename),
            KO_num = nrow(filtered_data),
            gene_num = length(unique(filtered_data$gene_name))))
    return(res);   
    }


if (toupper(GeneBank) == "TRUE"){
    filter.log <- data.frame();
    for (filename in filenames$name) {
    results <- keep_genebank(filename);
    outpath = file.path(out_path,paste(basename(filename),suffix,sep="_"));
    vroom_write(
        x = results$filtered_data,
        file = outpath,
        col_names = TRUE,
        quote = "none",
        delim = "\t");
    filter.log <- rbind(filter.log,
                    results$redun);               
 }
 vroom_write(
        x = filter.log,
        file = log,
        col_names = TRUE,
        quote = "none",
        delim = "\t");
} else {
   filter.log <- data.frame();
   for (filename in filenames$name) {
    results <- only_KO(filename);
    outpath = file.path(out_path,paste(basename(filename),suffix,sep="_"));
    vroom_write(
        x = results$filtered_data,
        file = outpath,
        col_names = TRUE,
        quote = "none",
        delim = "\t");
    filter.log <- rbind(filter.log,
                    results$redun);               

 }
 vroom_write(
        x = filter.log,
        file = log,
        col_names = TRUE,
        quote = "none",
        delim = "\t");
}