
library(data.table)
library(tidyverse)

chunk_size <- 3e6 # must be a multiple of 3
stopifnot(chunk_size %% 3 == 0)
skip <- 0
finished <- FALSE

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
total_lines <- R.utils::countLines(file_name)

while(!finished) {
  
  if(skip >= total_lines) {
    message("Done.")
    break
  }
  
  dt <- fread(file_name, skip=skip, nrows=chunk_size)
  print(paste(Sys.time(), "Processing rows", skip, "to", skip + nrow(dt)))
  
  # store prior to collapsing alternative mutations when averaging per position
  og_size <- nrow(dt)
  
  if(og_size == 0) {
    finished <- TRUE
    print("Done.")
    break
  }
  
  names(dt) <- c("chrom", "pos", "ref", "alt", "score")
  dt[, c("ref", "alt") := NULL] 
  
  dt[, benegas_score := mean(score), by=.(chrom, pos)]
  dt <- dt[!duplicated(pos)] %>% dplyr::select(., c(chrom, pos, benegas_score)) %>% setDT()
  chrom <- dt$chrom[1]
  dt[, benegas_score := round(benegas_score, digits=3)]
  
  fwrite(dt, paste0("split_scores/scores_chr", chrom, "_", min(dt$pos), "-", max(dt$pos), ".csv.gz"))
  
  skip <- skip + og_size
}

