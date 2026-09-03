##################################################################################
# Merge STAARpipeline Results
##################################################################################
rm(list=ls())
gc()

library(data.table)

# ================= 1. Argument Parsing =================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default=NULL) {
  idx <- which(args == flag)
  if(length(idx) > 0) return(args[idx + 1])
  if(!is.null(default)) return(default)
  stop(paste("Missing required argument:", flag))
}

out_path            <- get_arg("--out_path")
p2_path             <- get_arg("--p2_path")
array_id            <- as.numeric(get_arg("--array_id"))
log_path            <- get_arg("--log_path") 
summary_file        <- get_arg("--summary_file") 
subset_variants_num <- as.numeric(get_arg("--subset_num", default = 1000)) 
save_chunk_size     <- as.numeric(get_arg("--chunk_size", default = 100))

cat("Command-line arguments:\n")
print(list(array_id=array_id, summary_file=summary_file, subset_variants_num=subset_variants_num,
  save_chunk_size=save_chunk_size, p2_path=p2_path, out_path=out_path, log_path=log_path))

current_time <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
start_time <- Sys.time()
log_name <- paste0(log_path, 'p3_', array_id, '_', current_time, '.log') 
cat(paste0("\n\nstart time: ", start_time, "\n"), file = log_name, append = TRUE)

# ================= 2. Extract Meta Information =================

cat(paste0("Reading meta info from: ", summary_file, "\n"))

if (!file.exists(summary_file)) {
  stop("Summary file does not exist!")
}

meta_dt <- fread(summary_file, header = TRUE)
target_row <- meta_dt[Array_num == array_id]

if (nrow(target_row) == 0) {
  stop(paste0("ArrayID ", array_id, " not found in summary file!"))
} else if (nrow(target_row) > 1) {
  cat("Warning: Multiple entries found. Using the LAST entry.\n")
  target_row <- tail(target_row, 1)
}

# extract start_loc, end_loc, and count
n_variants <- as.numeric(target_row$target_ids_count)
start_loc  <- as.numeric(target_row$Start_Loc)
end_loc    <- as.numeric(target_row$End_Loc)

cat("================ Meta Info ================\n")
cat("Array ID:       ", array_id, "\n")
cat("Location:       ", start_loc, "-", end_loc, "\n")
cat("Total Variants: ", n_variants, "\n")
cat("===========================================\n")

# ================= 3. Merge =================

# 1. Calculate total number of subsets
subset.num <- ceiling(n_variants / subset_variants_num)

if (subset.num == 0) {
  warning("No variants indicated (subset.num = 0). Exiting.")
  q("no")
}

# 2. chunk indices
kk_indices <- c()
if (subset.num >= save_chunk_size) {
  kk_indices <- seq(from = save_chunk_size, to = subset.num, by = save_chunk_size)
}
if (!(subset.num %in% kk_indices)) {
  kk_indices <- c(kk_indices, subset.num)
}
kk_indices <- sort(unique(kk_indices))

cat(paste0("Expecting ", length(kk_indices), " chunk files to merge.\n"))

# 3. Construct file list
expected_files <- paste0(p2_path, array_id, "_", start_loc, "_", end_loc, 
                         "_chunk_kk_", kk_indices, ".Rdata")

# 4. Check for missing files
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  cat("\n[ERROR] Missing chunk files:\n")
  print(head(missing_files))
  stop("Merge aborted due to missing files.")
}

# 5. Load and merge
cat("\nStarting merge...\n")
data_list <- lapply(expected_files, function(f) {
  env <- new.env()
  load(f, envir = env)
  if (exists("results_kk", envir = env)) {
    return(get("results_kk", envir = env))
  } else {
    return(get(ls(envir = env)[1], envir = env))
  }
})

final_result <- rbindlist(data_list, fill = TRUE, use.names = TRUE)

# 6. Save
final_output_file <- paste0(out_path, "/ArrayID_", array_id, ".Rdata")
save(final_result, file = final_output_file)
cat(paste0("\nSaved merged file to: ", final_output_file, "\n"))

rm(data_list, final_result)
gc()

cat("======Finished!======\n")
end_time <- Sys.time()
run_time <- end_time - start_time
cat(paste0("end time: ", end_time, "\n"), file = log_name, append = TRUE)
cat(paste0("total cost: ", as.numeric(run_time, units = "mins"), " mins.\n"), file = log_name, append = TRUE)
