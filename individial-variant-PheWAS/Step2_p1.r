########### STARR PheWAS #############
###########################################################
# Individual analysis using STAARpipeline - Parallel Version
# part1: process and save the genofile seperately
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
adgs         <- get_arg("--adgs")
log_path     <- get_arg("--log_path")
pheno_path   <- get_arg("--pheno_path")
step1_path   <- get_arg("--step1_path")
out_path     <- get_arg("--out_path")
step2.1_path <- get_arg("--step2.1_path")
begin        <- as.integer(get_arg("--begin"))
end          <- as.integer(get_arg("--end"))
stage        <- get_arg("--stage")
## Input array ID from batch file (Harvard FAS RC cluster)
arrayid      <- as.integer(get_arg("--arrayid"))

start_time <- Sys.time()
current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
log_name <- paste0(log_path, 'p1_', current_time, '.log')
cat(paste0("\n\nstart time: ", start_time, "\n"), file = log_name, append = TRUE)


## Number of jobs for each chromosome
jobs_num <- get(load(jobs_path))
## aGDS directory
agds_dir <- get(load(adgs))

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
obj_nullmodel_list <- lapply(files, function(f) {
  tmp_name <- load(f)
  get(tmp_name)
})
length_obj <- length(obj_nullmodel_list)

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

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path) #读一个几十G的文件？ 

start_loc <- (groupid - 1) * 10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc, jobs_num$end_loc[chr])

## move the phenotype.id here and remove obj_nullmodel_list to SAVE STORAGE
phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
  as.character(x$id_include)
}))
rm(obj_nullmodel_list)
gc()
save(phenotype.id, file=paste0(step2.1_path, "phenotypeid","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata"))

## parameters for genofile processing
mac_cutoff = 20
subset_variants_num=5e3
tol=.Machine$double.eps^0.25

## processing genofile

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant") {
  SNVlist <- filter == "PASS"
} else if(variant_type=="SNV") {
  SNVlist <- (filter == "PASS") & isSNV(genofile)
} else if(variant_type=="Indel"){
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

position <- as.numeric(seqGetData(genofile, "position"))

variant.id <- seqGetData(genofile, "variant.id")
is.in <- (SNVlist) & (position>=start_loc) & (position<=end_loc)
SNV.id <- variant.id[is.in]

if(length(SNV.id) == 0) {
  return(rep(list(NULL), length_obj))
}

## get AF, AC and perform initial filtering
seqSetFilter(genofile, variant.id = SNV.id, sample.id = phenotype.id)
AF_AC_Missing <- seqGetAF_AC_Missing(genofile, minor=FALSE, parallel=FALSE) ##?????
REF_AF <- AF_AC_Missing$af
REF_AC <- AF_AC_Missing$ac
Missing_rate <- AF_AC_Missing$miss
ALT_AC <- 2*round(length(phenotype.id)*(1-Missing_rate))-REF_AC
MAC <- ifelse(REF_AC>=ALT_AC, ALT_AC, REF_AC)

is.include <- !((MAC<mac_cutoff) | is.na(MAC))
SNV.id <- SNV.id[is.include]
REF_AF <- REF_AF[is.include]
Missing_rate <- Missing_rate[is.include]
rm(AF_AC_Missing, is.include)
gc()

seqResetFilter(genofile)

subset.num <- ceiling(length(SNV.id)/subset_variants_num)
save(subset.num, file=paste0(step2.1_path, "subsetnum","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata"))
write.table(subset.num, paste0(step2.1_path, "subsetnum","_",begin,"_to_",end,"_array", arrayid,"_", stage,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
for(kk in 1:subset.num) {
  print(paste0("======",Sys.time(),"Saving genofile_kk",kk,"_",begin,"_to_",end,"_array", arrayid,"_", stage,"======"))
  if(kk < subset.num) {
    is.in <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
  } else {
    is.in <- ((kk-1)*subset_variants_num+1):length(SNV.id)
  }
  
  ## extract the sparse matrix of genotypes
  REF_AF.in <- REF_AF[is.in]
  Missing_rate.in <- Missing_rate[is.in]
  Genotype_sp <- Genotype_sp_extraction(genofile, variant.id=SNV.id[is.in],
                                        sample.id=phenotype.id,
                                        REF_AF=REF_AF.in, Missing_rate=Missing_rate.in) 
  save(Genotype_sp, file=paste0(step2.1_path, "genofile_kk",kk,"_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata"))
  seqResetFilter(genofile)
}

end_time <- Sys.time()
run_time <- end_time - start_time
cat(paste0("end time: ", end_time, "\n"), file = log_name, append = TRUE)
cat(paste0("total cost: ", as.numeric(run_time, units = "mins"), " mins.\n"), file = log_name, append = TRUE)
