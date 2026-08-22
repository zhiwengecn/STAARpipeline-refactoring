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
log_path     <- get_arg("--log_path")
pheno_path   <- get_arg("--pheno_path")
step1_path   <- get_arg("--step1_path")
out_path     <- get_arg("--out_path")
step2.1_path <- get_arg("--step2.1_path")
step2.2_path <- get_arg("--step2.2_path")
arrayid      <- as.integer(get_arg("--arrayid"))
begin        <- as.integer(get_arg("--begin"))
end          <- as.integer(get_arg("--end"))
stage        <- get_arg("--stage")
kk_input     <- get_arg("--kk")

current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
start_time <- Sys.time()
log_name <- paste0(log_path, 'p2', '_', kk_input, '_', current_time, '.log')
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

phenotype.id <- get(load(paste0(step2.1_path, "phenotypeid","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))
subset.num <- get(load(paste0(step2.1_path, "subsetnum","_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))

staar_s2_PheWAS_serial <- function(chr, start_loc, end_loc, genofile,  
                                   mac_cutoff=20, subset_variants_num=5e3,
                                   QC_label="annotation/filter", 
                                   variant_type=c("variant", "SNV", "Indel"),
                                   geno_missing_imputation=c("mean", "minor"),
                                   tol=.Machine$double.eps^0.25, 
                                   max_iter=1000, SPA_p_filter=TRUE, 
                                   p_filter_cutoff=0.05, kk) {
  results_list <- rep(list(c()), length_obj)
  print(paste0("processing ",kk,"/",subset.num,"......"))
  if(subset.num == 0) {
    save(results_list, file = paste0(step2.2_path, 'results_list_kk', '_', kk, '_', begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata"))
    return(NULL)
  }
  
  # load the sparse matrix of genotypes
  Genotype_sp <- get(load(paste0(step2.1_path, "genofile_kk",kk,"_",begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata")))
  
  Geno <- Genotype_sp$Geno
  results_information <- Genotype_sp$results_information
  rm(Genotype_sp)
  gc()
  
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  for (num in 1:length_obj) {
    tmp <- str_extract(files[num], "(?<=/)[^/]+(?=/[^/]+$)")
    print(paste0("processing phenotype ",num,": ",tmp,"......"))
    
    # load obj_nullmodel
    obj_nullmodel <- get(load(files[num]))
    
    ## number of traits in analysis
    n_pheno <- obj_nullmodel$n.pheno
    
    ## SPA status
    if(!is.null(obj_nullmodel$use_SPA))
    {
      use_SPA <- obj_nullmodel$use_SPA
    }else
    {
      use_SPA <- FALSE
    }
    
    ## residuals and cov
    residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
    if(SPA_p_filter) {
      ### dense GRM
      if(!obj_nullmodel$sparse_kins) {
        P <- obj_nullmodel$P
      }
      
      ### sparse GRM
      if(obj_nullmodel$sparse_kins) {
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
      }
    }
    
    ## SPA
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
    } else {
      ### dense GRM
      if(!obj_nullmodel$sparse_kins) {
        P <- obj_nullmodel$P
      }
      
      ### sparse GRM
      if(obj_nullmodel$sparse_kins) {
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
      }
    }
    
    id.phenotype.num <- as.character(obj_nullmodel$id_include)
    id.phenotype.num.merge <- data.frame(id.phenotype.num)
    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge, phenotype.id.merge, 
                                                by = c("id.phenotype.num" = "phenotype.id"))
    phenotype.id.match <- id.phenotype.num.merge$index
    samplesize.num <- length(phenotype.id.match)
    
    ## Extract the genotype matrix and variant information of the current trait
    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
    MAC.in <- Matrix::colSums(Geno.num, na.rm = TRUE)
    Geno.num <- Geno.num[, MAC.in>=mac_cutoff, drop = FALSE]
    results_information.num <- results_information[MAC.in>=mac_cutoff, , drop = FALSE]
    MAC.in <- MAC.in[MAC.in>=mac_cutoff]
    
    Missing_num.in <- Missing_num.sp(Geno.num)
    Missing_rate.in <- Missing_num.in/samplesize.num
    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
    ALT_AF.in <- results_information.num$ALT_AF
    
    CHR <- results_information.num$CHR
    position <- results_information.num$position
    REF <- results_information.num$REF
    ALT <- results_information.num$ALT
    N <- rep(samplesize.num, length(CHR))
    
    if (geno_missing_imputation == "mean") {
      Geno.num <- na.replace.sp(Geno.num, m=2*MAF.in)
    }
    if (geno_missing_imputation == "minor") {
      Geno.num <- na.replace.sp(Geno.num, is_NA_to_Zero=TRUE)
      MAF.in <- MAC.in/(2*samplesize.num)
    }
    ALT_AF.in <- ifelse(ALT_AF.in>0.5, 1-MAF.in, MAF.in)
    
    if(!all(CHR==chr)) {
      warning("chr does not match the chromosome of genofile (the opened aGDS)!")
    }
    
    if((use_SPA) & !SPA_p_filter) {
      if(length(MAF.in)>=1) {
        Geno.num <- as.matrix(Geno.num)
        pvalue <- Individual_Score_Test_SPA(Geno.num, XW, XXWX_inv, residuals.phenotype, muhat, tol, max_iter)
        
        results_temp <- data.frame(CHR=CHR, POS=position, REF=REF, ALT=ALT, 
                                    ALT_AF=ALT_AF.in, MAF=MAF.in, N=N, pvalue=pvalue)
        
        results_list[[num]] <- rbind(results_list[[num]], results_temp)
      }
    } else {
      ## Common_variants or variants with relatively high missing rate
      is.common_highmissing <- (MAF.in>=0.01) | (Missing_rate.in>=0.01)
      if(sum(is.common_highmissing)>=1) {
        Geno_common <- Geno.num[, is.common_highmissing, drop=FALSE]
        
        CHR_common <- CHR[is.common_highmissing]
        position_common <- position[is.common_highmissing]
        REF_common <- REF[is.common_highmissing]
        ALT_common <- ALT[is.common_highmissing]
        MAF_common <- MAF.in[is.common_highmissing]
        ALT_AF_common <- ALT_AF.in[is.common_highmissing]
        N_common <- N[is.common_highmissing]
        
        ## Split into small chunks to run
        subset_variants_num_common <- 200
        subset.num_common <- ceiling(length(CHR_common)/subset_variants_num_common)
        
        for(kk_common in 1:subset.num_common) {
          if(kk_common < subset.num_common) {
            is.in_common_subset <- ((kk_common-1)*subset_variants_num_common+1):(kk_common*subset_variants_num_common)
          }
          if(kk_common == subset.num_common) {
            is.in_common_subset <- ((kk_common-1)*subset_variants_num_common+1):length(CHR_common)
          }
          
          Geno_common_subset <- Geno_common[, is.in_common_subset, drop=FALSE]
          
          CHR_common_subset <- CHR_common[is.in_common_subset]
          position_common_subset <- position_common[is.in_common_subset]
          REF_common_subset <- REF_common[is.in_common_subset]
          ALT_common_subset <- ALT_common[is.in_common_subset]
          MAF_common_subset <- MAF_common[is.in_common_subset]
          ALT_AF_common_subset <- ALT_AF_common[is.in_common_subset]
          N_common_subset <- N_common[is.in_common_subset]
          
          ## sparse GRM
          if(obj_nullmodel$sparse_kins) {
            if(n_pheno == 1) {
              Score_test <- Individual_Score_Test_sp(Geno_common_subset, Sigma_i, Sigma_iX, cov, residuals.phenotype)
            } else {
              Geno_common_subset <- Diagonal(n = n_pheno) %x% Geno_common_subset
              Score_test <- Individual_Score_Test_sp_multi(Geno_common_subset, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
            }
          }
          
          ## dense GRM
          if(!obj_nullmodel$sparse_kins) {
            if(n_pheno == 1) {
              Score_test <- Individual_Score_Test_sp_denseGRM(Geno_common_subset, P, residuals.phenotype)
            } else {
              Geno_common_subset <- Diagonal(n = n_pheno) %x% Geno_common_subset
              Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_common_subset, P, residuals.phenotype, n_pheno)
            }
          }
          
          ## SPA approximation for small p-values
          if(use_SPA) {
            pvalue <- exp(-Score_test$pvalue_log)
            
            if(any(is.infinite(Score_test$pvalue_log))) {
              print("pvalue_log has infinite value!!!")
            }
            if(any(is.na(Score_test$pvalue_log))) {
              print("pvalue_log has NA value!!!")
            }
            if(any(is.na(pvalue))) {
              print("pvalue NA value count")
              print(sum(is.na(pvalue)))
              pvalue[is.na(pvalue)] <- 999
            }
            if(any(is.infinite(pvalue))) {
              print("pvalue has infinite value!!!")
            }
            
            pvalue_1 <- na.omit(pvalue)
            
            if(sum(pvalue_1 < p_filter_cutoff)>=1) {
              is.common_subset_SPA <- as.vector(pvalue < p_filter_cutoff)
              Geno_common_subset_SPA <- Geno_common_subset[, is.common_subset_SPA, drop=FALSE]
              Geno_common_subset_SPA <- as.matrix(Geno_common_subset_SPA)
              
              pvalue_SPA <- Individual_Score_Test_SPA(Geno_common_subset_SPA, XW, XXWX_inv, residuals.phenotype, muhat, tol, max_iter)
              
              pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
            }
          }
          
          if(use_SPA) {
            results_temp <- data.frame(CHR=CHR_common_subset, POS=position_common_subset, 
                                        REF=REF_common_subset, ALT=ALT_common_subset,
                                        ALT_AF=ALT_AF_common_subset, MAF=MAF_common_subset,
                                        N=N_common_subset, pvalue=pvalue)
          } else {
            if(n_pheno == 1) {
              results_temp <- data.frame(CHR=CHR_common_subset, POS=position_common_subset, 
                                          REF=REF_common_subset, ALT=ALT_common_subset,
                                          ALT_AF=ALT_AF_common_subset, MAF=MAF_common_subset,
                                          N=N_common_subset, pvalue=exp(-Score_test$pvalue_log),
                                          pvalue_log10=Score_test$pvalue_log/log(10),
                                          Score=Score_test$Score, Score_se=Score_test$Score_se,
                                          Est=Score_test$Est, Est_se=Score_test$Est_se)
            } else {
              results_temp <- data.frame(CHR=CHR_common_subset, POS=position_common_subset, 
                                          REF=REF_common_subset, ALT=ALT_common_subset,
                                          ALT_AF=ALT_AF_common_subset, MAF=MAF_common_subset,
                                          N=N_common_subset, pvalue=exp(-Score_test$pvalue_log),
                                          pvalue_log10=Score_test$pvalue_log/log(10))
              results_temp <- cbind(results_temp, matrix(Score_test$Score, ncol=n_pheno))
              colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score", seq_len(n_pheno))
            }
          }
          results_list[[num]] <- rbind(results_list[[num]], results_temp)
        }
      }
      
      ## Rare_variants with relatively low missing rate
      is.rare_lowmissing <- (MAF.in<0.01) & (Missing_rate.in<0.01)
      if(sum(is.rare_lowmissing)>=1) {
        Geno_rare <- Geno.num[, is.rare_lowmissing, drop=FALSE]
        
        CHR_rare <- CHR[is.rare_lowmissing]
        position_rare <- position[is.rare_lowmissing]
        REF_rare <- REF[is.rare_lowmissing]
        ALT_rare <- ALT[is.rare_lowmissing]
        MAF_rare <- MAF.in[is.rare_lowmissing]
        ALT_AF_rare <- ALT_AF.in[is.rare_lowmissing]
        N_rare <- N[is.rare_lowmissing]
        
        ## sparse GRM
        if(obj_nullmodel$sparse_kins) {
          if(n_pheno == 1) {
            Score_test <- Individual_Score_Test_sp(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
          } else {
            Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
            Score_test <- Individual_Score_Test_sp_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
          }
        }
        
        ## dense GRM
        if(!obj_nullmodel$sparse_kins) {
          if(n_pheno == 1) {
            Score_test <- Individual_Score_Test_sp_denseGRM(Geno_rare, P, residuals.phenotype)
          } else {
            Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
            Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno)
          }
        }
        
        ## SPA approximation for small p-values
        if(use_SPA) {
          pvalue <- exp(-Score_test$pvalue_log)
          
          if(any(is.infinite(Score_test$pvalue_log))) {
            print("Rare pvalue_log has infinite value!!!")
          }
          if(any(is.na(Score_test$pvalue_log))) {
            print("Rare pvalue_log has NA value!!!")
          }
          if(any(is.na(pvalue))) {
            print("Rare pvalue NA value count:")
            print(sum(is.na(pvalue)))
            pvalue[is.na(pvalue)] <- 999
          }
          if(any(is.infinite(pvalue))) {
            print("Rare pvalue has infinite value!!!")
          }
          
          is.rare_SPA <- as.vector(pvalue < p_filter_cutoff)
          if(sum(is.rare_SPA)>=1) {
            pvalue_SPA <- c()
            Geno_rare_SPA <- Geno_rare[, is.rare_SPA, drop=FALSE]
            
            ## Split into small chunks to run
            subset_variants_num_SPA <- 50
            subset.num_SPA <- ceiling(sum(is.rare_SPA)/subset_variants_num_SPA)
            
            for(kk_SPA in 1:subset.num_SPA) {
              if(kk_SPA < subset.num_SPA) {
                is.in_SPA_subset <- ((kk_SPA-1)*subset_variants_num_SPA+1):(kk_SPA*subset_variants_num_SPA)
              }
              if(kk_SPA == subset.num_SPA) {
                is.in_SPA_subset <- ((kk_SPA-1)*subset_variants_num_SPA+1):sum(is.rare_SPA)
              }
              
              Geno_rare_SPA_subset <- Geno_rare_SPA[, is.in_SPA_subset, drop=FALSE]
              Geno_rare_SPA_subset <- as.matrix(Geno_rare_SPA_subset)
              
              if(length(is.in_SPA_subset)>=1) {
                pvalue_SPA_subset <- Individual_Score_Test_SPA(Geno_rare_SPA_subset, XW, XXWX_inv, 
                                                                residuals.phenotype, muhat, tol, max_iter)
                pvalue_SPA <- c(pvalue_SPA, pvalue_SPA_subset)
              }
            }
            pvalue[is.rare_SPA] <- pvalue_SPA
          }
        }
        
        if(use_SPA) {
          results_temp <- data.frame(CHR=CHR_rare, POS=position_rare, REF=REF_rare, ALT=ALT_rare,
                                      ALT_AF=ALT_AF_rare, MAF=MAF_rare, N=N_rare, pvalue=pvalue)
        } else {
          if(n_pheno == 1) {
            results_temp <- data.frame(CHR=CHR_rare, POS=position_rare, REF=REF_rare, ALT=ALT_rare,
                                        ALT_AF=ALT_AF_rare, MAF=MAF_rare, N=N_rare,
                                        pvalue=exp(-Score_test$pvalue_log), 
                                        pvalue_log10=Score_test$pvalue_log/log(10),
                                        Score=Score_test$Score, Score_se=Score_test$Score_se,
                                        Est=Score_test$Est, Est_se=Score_test$Est_se)
          } else {
            results_temp <- data.frame(CHR=CHR_rare, POS=position_rare, REF=REF_rare, ALT=ALT_rare,
                                        ALT_AF=ALT_AF_rare, MAF=MAF_rare, N=N_rare,
                                        pvalue=exp(-Score_test$pvalue_log), 
                                        pvalue_log10=Score_test$pvalue_log/log(10))
            results_temp <- cbind(results_temp, matrix(Score_test$Score, ncol=n_pheno))
            colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score", seq_len(n_pheno))
          }
        }
        
        results_list[[num]] <- rbind(results_list[[num]], results_temp)
      }
    }
  }

  save(results_list, file = paste0(step2.2_path, 'results_list_kk', '_', kk, '_', begin,"_to_",end,"_array", arrayid,"_", stage,".Rdata"))
}

## Execute analysis
a <- Sys.time()
if (start_loc <= end_loc) {
  staar_s2_PheWAS_serial(
    chr = chr, start_loc = start_loc, end_loc = end_loc,
    mac_cutoff = 20,
    QC_label = QC_label, variant_type = variant_type,
    geno_missing_imputation = geno_missing_imputation,
    max_iter = 1000,
    SPA_p_filter = TRUE,
    p_filter_cutoff = 0.05,
    kk=kk_input
  )
}
b <- Sys.time()
cat("Time taken:", b - a, "\n")
