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
              type="character", default="kofam_anno_name", dest="file_path",
              help="kofam detail-tsv file path", metavar="character"),
  make_option(c("-o", "--out_path"), 
              type="character", default=".",  dest="out_path",
              help="output directory", metavar="character"),
  make_option(c("-s", "--strict"), 
              type="logical", default="FALSE",  dest="strict",
              help="If TRUE,the score in the filtering results must be greater than the threshld", 
              metavar="logical"),
  make_option(c("-p", "--suffix"), 
              type="character", default="filter",  dest="suffix",
              help="output file suffix", 
              metavar="character"),
  make_option(c("-l", "--log"), 
              type="character", default="kofam.filter.log",  dest="log",
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
strict  <- args$strict
suffix  <- args$suffix
log  <- args$log
dir.create(out_path)


setwd(work_dir)
filenames <- read.table(file_path, header=FALSE)

selet_fun <- function(filename){
test <- vroom(filename, 
 col_names = FALSE,
 comment = "#",
 col_types =c("cccnnnc"),
 col_select = c(2:7),
 delim = "\t");
colnames(test) = c("gene_name","KO",
"thrshld","score","E_value","KO_definition");
if(strict %in% c("TRUE","true","True")){
filtered_data <- test %>%
  mutate(thrshld = if_else(is.na(thrshld),0,thrshld)) %>%
  mutate(score2 = score - thrshld) %>% 
  filter(E_value <= 1e-5) %>% 
  group_by(gene_name) %>%
  filter(score2 >= 0) %>%
  filter(E_value == min(E_value)) %>%
  filter(score2 == max(score2)) %>%
  ungroup() %>%
  select(gene_name, KO,"thrshld","score","E_value","KO_definition");
  } else {
  filtered_data <- test %>%
  mutate(thrshld = if_else(is.na(thrshld),0,thrshld)) %>%
  #mutate(score2 = score - thrshld) %>% 
  filter(E_value <= 1e-5) %>% 
  group_by(gene_name) %>%
  filter(score == max(score)) %>%
  filter(E_value == min(E_value)) %>%
  #filter(score2 == max(score2)) %>%
  ungroup() %>%
  select(gene_name, KO,"thrshld","score","E_value","KO_definition")
  }
res <- list(filtered_data=filtered_data,
        redun=data.frame(
            filename = filename,
            KO_num = nrow(filtered_data),
            gene_num = length(unique(filtered_data$gene_name))))
return(res);
}

filter.log <- data.frame();

for (filename in filenames[,1]) {
    results <- selet_fun(filename);
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
