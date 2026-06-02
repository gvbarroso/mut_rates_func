
library(data.table)
library(tidyverse)

chunk_size <- 3e6 # must be a multiple of 3
stopifnot(chunk_size %% 3 == 0)
skip <- 32 # number of header lines
finished <- FALSE

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]

file_name <- paste0(chr, "_rate_v5.2_TFBS_correction_all.vcf.bgz")
total_lines <- R.utils::countLines(file_name)
cat(paste0("total lines:", total_lines, "\n"))

while(!finished) {
  
  if(skip >= total_lines) {
    message("Done.")
    break
  }
  
  mut_map <- fread(file_name, skip=skip, nrows=chunk_size)
  print(paste(Sys.time(), "Processing rows", skip, "to", skip + nrow(mut_map)))
  # store prior to collapsing alternative mutations when averaging per position
  og_size <- nrow(mut_map)
  
  if(og_size == 0) {
    finished <- TRUE
    print("Done.")
    break
  }
  
  mut_map[, 3:7 := NULL] # remove by idx because of name irregularities between files (header)
  names(mut_map) <- c("chrom", "pos", "rates")
  
  mut_map[, roulette := as.numeric(sub(".*MR=([0-9.]+).*", "\\1", rates))]
  mut_map[, carlson := as.numeric(sub(".*MC=([0-9.]+).*", "\\1", rates))]
  mut_map[, gnomad := as.numeric(sub(".*MG=([0-9.]+).*", "\\1", rates))]
  mut_map[, triplet := sub(".*PN=([^;]+).*", "\\1", rates)] # 5-mer
  mut_map[, triplet := substr(triplet, 2, nchar(triplet) - 1)] # keeps inner 3-mer
  mut_map[, rates := NULL]
  
  # https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
  roulette_scale <- 1.015e-7 / 2 
  gnomad_scale <- roulette_scale
  carlson_scale <- 2.086e-9 / 2
  
  mut_map[, roulette := roulette * roulette_scale]
  mut_map[, carlson := carlson * carlson_scale]
  mut_map[, gnomad := gnomad * gnomad_scale]
  
  mut_map[, roulette := mean(roulette, na.rm=T), by=.(chrom, pos)]
  mut_map[, carlson := mean(carlson, na.rm=T), by=.(chrom, pos)]
  mut_map[, gnomad := mean(gnomad, na.rm=T), by=.(chrom, pos)]
  
  mut_map <- mut_map[!duplicated(pos)]
  
  fwrite(mut_map, paste0("split_muts/mut_map_chr", chr, "_", min(mut_map$pos), "-", max(mut_map$pos), ".csv.gz"))
  skip <- skip + og_size
}

