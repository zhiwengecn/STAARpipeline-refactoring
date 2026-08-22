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
adgs_path           <- get_arg("--adgs")
null_path           <- get_arg("--null")
jobs_path           <- get_arg("--jobs")
out_prefix          <- get_arg("--out")
p1_path             <- get_arg("--p1_path")
p2_path             <- get_arg("--p2_path")
log_path            <- get_arg("--log_path")
array_id            <- as.integer(get_arg("--array_id"))
subset_variants_num <- as.numeric(get_arg("--subset_num")) 
save_chunk_arg_idx  <- which(args == "--save_chunk")
save_chunk_size     <- if(length(save_chunk_arg_idx) > 0) as.integer(args[save_chunk_arg_idx + 1]) else 100
kk_start            <- as.integer(get_arg("--kk_start"))
kk_end              <- as.integer(get_arg("--kk_end"))

cat("Command-line arguments:\n")
print(list(adgs_path=adgs_path, null_path=null_path, jobs_path=jobs_path, out_prefix=out_prefix,
  p1_path=p1_path, p2_path=p2_path, log_path=log_path, array_id=array_id,
  subset_variants_num=subset_variants_num, save_chunk_size=save_chunk_size,
  kk_start=kk_start, kk_end=kk_end))

# Check array range
if(array_id < 1 || array_id > 288) {
  stop("Error: array_id must be in the range 1..288")
}

current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
start_time <- Sys.time()
log_name <- paste0(log_path, 'p2_', array_id, '_', kk_start, '_', kk_end, '_', current_time, '.log') 
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
#           Optimized Individual Analysis Function
###########################################################

Individual_Analysis <- function(
    chr,
    genofile, 
    target_variant_ids,
    sample_indices,
    precalc_model,
    start_loc, 
    end_loc,
    array_id,
    mac_cutoff=20,
    save_chunk_size=100,
    subset_variants_num=5e3,
    kk_start,
    kk_end,
    geno_missing_imputation=c("mean","minor"),
    tol=.Machine$double.eps^0.25,
    max_iter=1000,
    out_prefix_path,
    log_file_path,
    use_ancestry_informed=FALSE,
    p_filter_cutoff=0.05,
    SPA_p_filter=TRUE){
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    # Unpack precomputed model components
    phenotype.id <- precalc_model$phenotype.id
    residuals.phenotype <- precalc_model$residuals
    n_pheno <- precalc_model$n_pheno
    samplesize <- length(phenotype.id)
    use_SPA <- precalc_model$use_SPA
    
    P <- precalc_model$P
    Sigma_i <- precalc_model$Sigma_i
    Sigma_iX <- precalc_model$Sigma_iX
    cov <- precalc_model$cov
    muhat <- precalc_model$muhat
    XW <- precalc_model$XW
    XXWX_inv <- precalc_model$XXWX_inv
    sparse_kins <- precalc_model$sparse_kins
    
    if(use_ancestry_informed) return(NULL) # Skip AI for now
    results_kk <- NULL
    subset.num <- ceiling(length(target_variant_ids)/subset_variants_num)
    if(subset.num == 0) return(NULL)
    for(kk in kk_start:kk_end)
    {
        cat(paste0(Sys.time(), " kk: ", kk, "/", subset.num, "\n"), file = log_name, append = TRUE)

        start_idx <- (kk-1)*subset_variants_num + 1
        end_idx <- min(kk*subset_variants_num, length(target_variant_ids))
        current_ids <- target_variant_ids[start_idx:end_idx]
        
        # Skip AI for now
        seqSetFilter(genofile, variant.id=current_ids, sample.id=phenotype.id, verbose=FALSE)
        
        Geno_raw <- seqGetData(genofile, "$dosage")
        Geno <- Geno_raw[sample_indices,, drop=FALSE] # Align samples
        
        if(geno_missing_imputation=="mean") Geno <- matrix_flip_mean(Geno)
        if(geno_missing_imputation=="minor") Geno <- matrix_flip_minor(Geno)
        rm(Geno_raw)
        
        MAF <- Geno$MAF
        ALT_AF <- 1 - Geno$AF
        CHR <- as.numeric(seqGetData(genofile, "chromosome"))
        position <- as.numeric(seqGetData(genofile, "position"))
        REF <- as.character(seqGetData(genofile, "$ref"))
        ALT <- as.character(seqGetData(genofile, "$alt"))
        N <- rep(samplesize,length(CHR))

        if(!all(CHR==chr)) warning("chr mismatch!")

        # --- Analysis ---
        if((use_SPA)&!SPA_p_filter) {
          if(sum(MAF>(mac_cutoff-0.5)/samplesize/2)>=1)
          {
                Geno <- Geno$Geno
                Geno_common <- Geno[,(MAF>(mac_cutoff-0.5)/samplesize/2),drop=FALSE]
                # ... Extract Meta ...
                CHR_c <- CHR[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                POS_c <- position[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                REF_c <- REF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                ALT_c <- ALT[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                MAF_c <- MAF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                AF_c <- ALT_AF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                N_c <- N[(MAF>(mac_cutoff-0.5)/samplesize/2)]
                
                rm(Geno);
                gc()

                pvalue <- Individual_Score_Test_SPA(Geno_common,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
                results_temp <- data.frame(CHR=CHR_c,POS=POS_c,REF=REF_c,ALT=ALT_c,ALT_AF=AF_c,MAF=MAF_c,N=N_c, pvalue=pvalue)
                results_kk <- rbind(results_kk,results_temp)
          }
        } else
        {
            ## Common
            if(sum(MAF>=0.05)>=1) {
                Geno_common <- Geno$Geno[,MAF>=0.05, drop=FALSE]
                # ... Extract Meta ...
                CHR_c <- CHR[MAF>=0.05]; POS_c <- position[MAF>=0.05]; REF_c <- REF[MAF>=0.05]; 
                ALT_c <- ALT[MAF>=0.05]; MAF_c <- MAF[MAF>=0.05]; AF_c <- ALT_AF[MAF>=0.05]; N_c <- N[MAF>=0.05]

                if(sum(MAF>=0.05)==1) { Geno_common <- as.matrix(Geno_common,ncol=1) }

                # Test
                if(sparse_kins) {
                    if(n_pheno == 1) Score_test <- Individual_Score_Test(Geno_common, Sigma_i, Sigma_iX, cov, residuals.phenotype)
                    else { Geno_common <- Diagonal(n = n_pheno) %x% Geno_common; Score_test <- Individual_Score_Test_sp_multi(Geno_common, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno) }
                } else {
                    if(n_pheno == 1) Score_test <- Individual_Score_Test_denseGRM(Geno_common, P, residuals.phenotype)
                    else { Geno_common <- Diagonal(n = n_pheno) %x% Geno_common; Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_common, P, residuals.phenotype, n_pheno) }
                }

                # SPA
                if(use_SPA) {
                    pvalue <- exp(-Score_test$pvalue_log)
                    if(sum(pvalue < p_filter_cutoff)>=1) {
                        Geno_common_SPA <- Geno_common[,pvalue < p_filter_cutoff,drop=FALSE]
                        pvalue_SPA <- Individual_Score_Test_SPA(Geno_common_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
                        pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
                    }
                }
                
                # Store
                if(use_SPA) {
                    results_temp <- data.frame(CHR=CHR_c,POS=POS_c,REF=REF_c,ALT=ALT_c,ALT_AF=AF_c,MAF=MAF_c,N=N_c, pvalue=pvalue)
                } else {
                  if(n_pheno == 1) {
                    results_temp <- data.frame(CHR=CHR_c,POS=POS_c,REF=REF_c,ALT=ALT_c,ALT_AF=AF_c,MAF=MAF_c,N=N_c,
                                              pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
                                              Score=Score_test$Score,Score_se=Score_test$Score_se,
                                              Est=Score_test$Est,Est_se=Score_test$Est_se)
                  } else {
                    results_temp <- data.frame(CHR=CHR_c,POS=POS_c,REF=REF_c,ALT=ALT_c,ALT_AF=AF_c,MAF=MAF_c,N=N_c,
                                              pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
                    results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
                    colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
                  }
                }
                results_kk <- rbind(results_kk, results_temp)
            }

            ## Rare
            if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=1) 
            {
                Geno_rare <- Geno$Geno[,(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                # ... Extract Meta ...
                CHR_r <- CHR[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                POS_r <- position[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                REF_r <- REF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                ALT_r <- ALT[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                MAF_r <- MAF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                AF_r <- ALT_AF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
                N_r <- N[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]

                # Test logic (Simplified for brevity, assume similar to above)
                if(sparse_kins) {
                     if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2) {
                        if(n_pheno==1) { 
                          Geno_rare <- as(Geno_rare,"dgCMatrix"); Score_test <- Individual_Score_Test_sp(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype) 
                        } else { 
                          Geno_rare <- Diagonal(n=n_pheno) %x% Geno_rare; Score_test <- Individual_Score_Test_sp_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno) }
                     } else {
                        # Single rare variant case
                        if(n_pheno==1) {
                          Geno_rare <- as.matrix(Geno_rare,ncol=1); Score_test <- Individual_Score_Test(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
                        } else {
                          Geno_rare <- as.matrix(Diagonal(n = n_pheno) %x% Geno_rare); Score_test <- Individual_Score_Test_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno) }
                     }
                } else { # Dense
                     # Similar logic for dense...
                     if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2) {
                        if(n_pheno==1) {
                          Geno_rare <- as(Geno_rare,"dgCMatrix"); Score_test <- Individual_Score_Test_sp_denseGRM(Geno_rare, P, residuals.phenotype)
                        } else {
                          Geno_rare <- Diagonal(n=n_pheno) %x% Geno_rare; Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno)
                        }
                     } else {
                        if(n_pheno==1) {
                          Geno_rare <- as.matrix(Geno_rare,ncol=1); Score_test <- Individual_Score_Test_denseGRM(Geno_rare, P, residuals.phenotype)
                        } else { 
                          Geno_rare <- as.matrix(Diagonal(n = n_pheno) %x% Geno_rare); Score_test <- Individual_Score_Test_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno) }
                     }
                }

                # SPA Rare
                if(use_SPA) {
                    pvalue <- exp(-Score_test$pvalue_log)
                    if(sum(pvalue < p_filter_cutoff)>=1) {
                        Geno_rare_SPA <- as.matrix(Geno_rare)[,pvalue < p_filter_cutoff, drop=FALSE] # Handle drop dimensions
                        pvalue_SPA <- Individual_Score_Test_SPA(Geno_rare_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
                        pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
                    }
                }

                if(use_SPA) {
                    results_temp <- data.frame(CHR=CHR_r,POS=POS_r,REF=REF_r,ALT=ALT_r,ALT_AF=AF_r,MAF=MAF_r,N=N_r, pvalue=pvalue)
                } else {
                  if(n_pheno == 1)
                  {
                    results_temp <- data.frame(CHR=CHR_r,POS=POS_r,REF=REF_r,ALT=ALT_r,ALT_AF=AF_r,MAF=MAF_r,N=N_r,
                                              pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
                                              Score=Score_test$Score,Score_se=Score_test$Score_se,
                                              Est=Score_test$Est,Est_se=Score_test$Est_se)
                  }
                  else
                  {
                    results_temp <- data.frame(CHR=CHR_r,POS=POS_r,REF=REF_r,ALT=ALT_r,ALT_AF=AF_r,MAF=MAF_r,N=N_r,
                                              pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
                    results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
                    colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
                  }
                }
                results_kk <- rbind(results_kk, results_temp)
            }
        }
        if (kk %% save_chunk_size == 0 || kk == subset.num) {
          if(!is.null(results_kk)){
            chunk_filename <- paste0(p2_path, array_id, "_", start_loc, "_", end_loc, "_chunk_kk_", kk, ".Rdata")
            save(results_kk, file=chunk_filename)
            cat(paste0("Saved buffer to: ", chunk_filename, " (Current kk=", kk, ")\n"), file = log_name, append = TRUE)
          }
          results_kk <- NULL
        }
        rm(Geno, MAF, ALT_AF)
        seqResetFilter(genofile, verbose=FALSE)
        gc() 
    }
    
    return(NULL)
}

###########################################################
#           Main Function Loop
###########################################################

# [Pre-calc 1] Extract Null Model Matrices
cat("Loading Pre-calculate Null Model matrices...\n")
precalculate_null_model <- get(load(paste0(p1_path, "precalculate_null_model.Rdata")))
  
cat("Processing array_id: ", array_id, ", kk_start: ", kk_start, ", kk_end: ", kk_end, "...\n")
cat(paste0(Sys.time(), " Processing array_id: ", array_id, ". kk_start: ",
 kk_start, ", kk_end: ", kk_end, "\n"), file = log_name, append = TRUE)

chr <- which.max(array_id <= cumsum(jobs_num$individual_analysis_num))
if (chr == 1) { groupid <- array_id } else { groupid <- array_id - cumsum(jobs_num$individual_analysis_num)[chr - 1] }

agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc, jobs_num$end_loc[chr])

if (start_loc <= end_loc) {
  sample_indices <- get(load(paste0(p1_path, array_id, "_sample_indices.Rdata")))
  target_ids <- get(load(paste0(p1_path, "target_ids_", array_id, "_", start_loc, "_", end_loc, ".Rdata")))
  
  cat(paste0("Found ", length(target_ids), " variants.\n"))
  
  if(length(target_ids) > 0) {
      Individual_Analysis(
          chr = chr,
          genofile = genofile,
          target_variant_ids = target_ids,
          sample_indices = sample_indices,
          precalc_model = precalculate_null_model,
          start_loc = start_loc,
          end_loc = end_loc,
          array_id = array_id,
          out_prefix_path = out_prefix,
          log_file_path = log_name,
          mac_cutoff = 20,
          save_chunk_size = save_chunk_size,
          subset_variants_num = subset_variants_num, 
          kk_start = kk_start,
          kk_end = kk_end,
          geno_missing_imputation = geno_missing_imputation,
          SPA_p_filter = TRUE
      )
  }
}
seqClose(genofile)
gc()  

cat("======Finished!======\n")
end_time <- Sys.time()
run_time <- end_time - start_time
cat(paste0("end time: ", end_time, "\n"), file = log_name, append = TRUE)
cat(paste0("total cost: ", as.numeric(run_time, units = "mins"), " mins.\n"), file = log_name, append = TRUE)


