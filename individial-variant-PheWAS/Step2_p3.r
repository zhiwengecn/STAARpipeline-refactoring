########### STARR PheWAS #############
###########################################################
# Individual analysis using STAARpipeline - Parallel Version
###########################################################
rm(list=ls())
gc()

## Load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelinePheWAS)
library(data.table)
library(stringr)


## Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if(length(idx) == 0) stop(paste("Missing argument", flag))
  if(idx == length(args)) stop(paste("No value provided for", flag))
  return(args[idx + 1])
}

jobs_path    <- get_arg("--jobs_path")
pheno_path   <- get_arg("--pheno_path")
step1_path   <- get_arg("--step1_path")
out_path     <- get_arg("--out_path")
step2.1_path <- get_arg("--step2.1_path")
step2.2_path <- get_arg("--step2.2_path")
## Input array ID from batch file (Harvard FAS RC cluster)
arrayid      <- as.integer(get_arg("--arrayid"))
begin        <- as.integer(get_arg("--begin"))
end          <- as.integer(get_arg("--end"))
stage        <- get_arg("--stage")
log_path     <- get_arg("--log_path")


current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
start_time <- Sys.time()
log_name <- paste0(log_path, 'p3', current_time, '.log')
cat(paste0("\n\nstart time: ", start_time, "\n"), file = log_name, append = TRUE)


## Number of jobs for each chromosome
jobs_num <- get(load(jobs_path))

## Load phenotype list in parallel-safe manner
list <- fread(paste0(pheno_path,'staar_fullsample_disease_sum.txt'))
selected_folders <- NULL
for (i in begin:end){
  name <- as.character(list[i,1])
  selected_folders <- c(selected_folders,name)
}  

print(paste0("####Line: ",begin,"_to_",end," ",stage," Begin####"))
print(selected_folders)

all_files <- list.files(path = paste0(step1_path,stage), pattern = ".Rdata$", full.names = TRUE, recursive = TRUE)
files <- grep(paste0("/(", paste(selected_folders, collapse = "|"), ")/.*_obj_nullmodel_match_v2[.]Rdata$"), 
              all_files, 
              value = TRUE,
              ignore.case = TRUE)
length_obj <- length(files)

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"


###########################################################
#           Main Function - Parallel Version
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1) {
  groupid <- arrayid
} else {
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr - 1]
}

start_loc <- (groupid - 1) * 10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc, jobs_num$end_loc[chr])

#@objnullmodel_add 
phenotype.id <- get(load(paste0(step2.1_path, "phenotypeid","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))

subset.num <- get(load(paste0(step2.1_path, "subsetnum","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))

results <- rep(list(c()), length_obj)

for (kk in 1:subset.num) 
{
  results_tmp <- get(load(paste0(step2.2_path, "results_list_kk", "_", kk, "_", begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))
  
  for(num in 1:length_obj) {
    # 先检查是否为NULL，再检查nrow
    if(!is.null(results_tmp[[num]])) {
      if(nrow(results_tmp[[num]]) > 0) {
        if(is.null(results[[num]])) {
          results[[num]] <- results_tmp[[num]]
        } else {
          results[[num]] <- rbind(results[[num]], results_tmp[[num]])
        }
      }
    }
  }
}

for (num in 1:length_obj) {
  if (!is.null(results[[num]])) {
    results[[num]] <- results[[num]][order(results[[num]][, 2]), ]
  }
}

save(results, file=paste0(out_path, "backup/", stage, "/from_", begin, "_to_", end, "_", arrayid, ".Rdata"))

for (n in 1:length(files)){
  file_name <- files[n]
  dis_name <- str_extract(file_name, "(?<=/)[^/]+(?=/[^/]+$)")
  ## file directory for the output file 
  out_path_dis=paste0(out_path,stage,"/",dis_name,"/")
  if (!dir.exists(out_path_dis)) {
    dir.create(out_path_dis, recursive = TRUE)
  } else { 
  }
  ## output file name
  output_name <- paste0(stage,"_",dis_name,"_Individual_Analysis")
  single_result <- results[n]
  save(single_result,file=paste0(out_path_dis,output_name,"_",arrayid,".Rdata"))
}
