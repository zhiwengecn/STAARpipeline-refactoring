##################################################################################
# Individual (single-variant) analysis using STAARpipeline
# Modified: Save Rdata for each kk loop
##################################################################################
rm(list=ls())
gc()

## Load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(data.table)
library(stringr)
library(Matrix)

###########################################################
#           Parse command line arguments
###########################################################
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if(length(idx) == 0) stop(paste("Missing argument", flag))
  if(idx == length(args)) stop(paste("No value provided for", flag))
  return(args[idx + 1])
}

# Retrieve arguments
adgs_path   <- get_arg("--adgs")
null_path   <- get_arg("--null")
jobs_path   <- get_arg("--jobs")
p1_path     <- get_arg("--p1_path")
log_path    <- get_arg("--log_path")
array_start <- as.integer(get_arg("--array_start"))
array_end   <- as.integer(get_arg("--array_end"))
subset_variants_num <- as.numeric(get_arg("--subset_num")) 

cat("Command-line arguments:\n")
print(list(adgs_path=adgs_path, null_path=null_path, jobs_path=jobs_path,
  p1_path=p1_path, log_path=log_path, subset_variants_num=subset_variants_num,
  array_start=array_start, array_end=array_end))

# Check array range
if(array_start < 1 || array_start > 288 || array_end < 1 || array_end > 288) {
  stop("Error: --array_start and --array_end must be in the range 1..288")
}
if(array_end < array_start) stop("Error: --array_end must be >= --array_start")

current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
start_time <- Sys.time()
log_name <- paste0(log_path, 'p1_', array_start, '_', array_end, '_', current_time, '.log')
cat(paste0("\n\nstart time: ", start_time, "\n"), file = log_name, append = TRUE)

###########################################################
#           Load Input Files
###########################################################
agds_dir <- local({ e <- new.env(); load(adgs_path, envir=e); if(!is.null(e$agds_dir)) e$agds_dir else e$a })
obj_nullmodel <- local({ e <- new.env(); load(null_path, envir=e); e$obj_nullmodel })
jobs_num <- local({ e <- new.env(); load(jobs_path, envir=e); e$jobs_num })

## QC label and variant options
QC_label <- "annotation/info/QC_label"
variant_type <- "variant"
geno_missing_imputation <- "mean"


###########################################################
#           Helper Functions (Optimized)
###########################################################

# 1. Prepare Null Model Matrices (Run once per script execution)
prepare_null_matrices <- function(obj_nullmodel, SPA_p_filter=TRUE) {
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
  
  use_SPA <- if(!is.null(obj_nullmodel$use_SPA)) obj_nullmodel$use_SPA else FALSE
  
  P <- NULL; Sigma_i <- NULL; Sigma_iX <- NULL; cov <- NULL
  muhat <- NULL; XW <- NULL; XXWX_inv <- NULL
  
  if(SPA_p_filter || !use_SPA) {
    if(!obj_nullmodel$sparse_kins) {
      P <- obj_nullmodel$P
    } else {
      Sigma_i <- obj_nullmodel$Sigma_i
      Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
      cov <- obj_nullmodel$cov
    }
  }
  
  if(use_SPA) {
    muhat <- obj_nullmodel$fitted.values
    if(obj_nullmodel$relatedness) {
      if(!obj_nullmodel$sparse_kins) {
        XW <- obj_nullmodel$XW
        XXWX_inv <- obj_nullmodel$XXWX_inv
      } else {
        XW <- as.matrix(obj_nullmodel$XSigma_i)
        XXWX_inv <- as.matrix(obj_nullmodel$XXSigma_iX_inv)
      }
    } else {
      XW <- obj_nullmodel$XW
      XXWX_inv <- obj_nullmodel$XXWX_inv
    }
  }
  
  return(list(
    phenotype.id = phenotype.id, n_pheno = n_pheno, residuals = residuals.phenotype,
    P = P, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX, cov = cov,
    muhat = muhat, XW = XW, XXWX_inv = XXWX_inv, use_SPA = use_SPA,
    sparse_kins = obj_nullmodel$sparse_kins
  ))
}

# 2. Get Sample Mapping Indices (Run once per chromosome/file)
get_sample_mapping <- function(genofile, phenotype.id) {
  seqSetFilter(genofile, sample.id = phenotype.id, verbose = FALSE)
  id.genotype <- seqGetData(genofile, "sample.id")
  idx <- match(phenotype.id, id.genotype)

  if(any(is.na(idx))) stop("Error: Some phenotype IDs are missing in the GDS file!")
  seqResetFilter(genofile, verbose = FALSE)
  return(idx)
}

# 3. Select Target Variants (Run once per chunk)
select_target_variants <- function(genofile, start_loc, end_loc, 
                                   QC_label="annotation/filter", variant_type="variant") {
  
  # Read metadata ONLY once for the file or rely on current filter
  # Assuming genofile is open and we want to subset from the whole
  filter <- seqGetData(genofile, QC_label)
  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  
  pass_filter <- filter == "PASS"
  
  if(variant_type=="SNV") {
    is_snv <- isSNV(genofile)
    pass_filter <- pass_filter & is_snv
  }
  if(variant_type=="Indel") {
    is_snv <- isSNV(genofile)
    pass_filter <- pass_filter & !is_snv
  }
  
  is_in_region <- (position >= start_loc) & (position <= end_loc)
  
  target_ids <- variant.id[pass_filter & is_in_region]
  return(target_ids)
}

###########################################################
#           Main Function Loop
###########################################################

# [Pre-calc 1] Extract Null Model Matrices
cat("Pre-calculating Null Model matrices...\n")
precalculate_null_model <- prepare_null_matrices(obj_nullmodel, SPA_p_filter=TRUE)
save(precalculate_null_model, file=paste0(p1_path, "precalculate_null_model.Rdata"))
rm(obj_nullmodel); gc()
summary_info_file <- paste0(p1_path, "Analysis_Meta_Summary.txt")

for (arrayid in array_start:array_end) {
  
  cat("Processing arrayid:", arrayid, "\n")
  cat(paste0(Sys.time(), " Processing arrayid: ", arrayid, "\n"), file = log_name, append = TRUE)
  
  chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
  if (chr == 1) { groupid <- arrayid } else { groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr - 1] }
  
  agds.path <- agds_dir[chr]
  genofile <- seqOpen(agds.path)
  
  start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
  end_loc <- start_loc + 10e6 - 1
  end_loc <- min(end_loc, jobs_num$end_loc[chr])
  
  if (start_loc <= end_loc) {
    sample_indices <- get_sample_mapping(genofile, precalculate_null_model$phenotype.id)
    target_ids <- select_target_variants(genofile, start_loc, end_loc, QC_label, variant_type)
    n_targets <- length(target_ids)
    subset.num <- ceiling(n_targets/subset_variants_num)
    cat(paste0("Found ", n_targets, " variants.\n"))

    meta_line <- data.frame(
      Array_num = arrayid,
      Start_Loc = start_loc,
      End_Loc   = end_loc,
      target_ids_count = n_targets,
      subset_variants_num = subset_variants_num,
      kk_num = subset.num
    )
    Sys.sleep(runif(1, 0.1, 5.0))
    tryCatch({
      append_mode <- file.exists(summary_info_file)
      
      write.table(meta_line, file = summary_info_file, 
                  append = TRUE, quote = FALSE, sep = "\t", 
                  row.names = FALSE, 
                  col.names = !append_mode)
                  
      cat(paste0("Recorded meta info to: ", summary_info_file, "\n"))
      
    }, error = function(e) {
      cat(paste0("[WARNING] Failed to write meta info (File Locked?): ", e$message, "\n"))
      cat("[ACTION] You may need to manually check the summary file for this ArrayID.\n")
    })
    
    save(sample_indices, file=paste0(p1_path, arrayid, "_sample_indices.Rdata"))
    save(target_ids, file=paste0(p1_path, "target_ids_", arrayid, "_", start_loc, "_", end_loc, ".Rdata"))
  }
  seqClose(genofile)
  gc()
}

cat("======Finished!======\n")
end_time <- Sys.time()
run_time <- end_time - start_time
cat(paste0("end time: ", end_time, "\n"), file = log_name, append = TRUE)
cat(paste0("total cost: ", as.numeric(run_time, units = "mins"), " mins.\n"), file = log_name, append = TRUE)


