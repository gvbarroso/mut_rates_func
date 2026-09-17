library(arrow)
library(data.table)

# downloaded on 16/09/2026
# https://huggingface.co/datasets/songlab/gpn-star-scores

# Usage:
# Rscript parse_gpn-star.R llr/llr_chr1.parquet

args <- commandArgs(trailingOnly=TRUE)
f <- args[1]

ds <- open_dataset(f, format="parquet")

scanner <- Scanner$create(ds, columns=c("chrom", "pos", "ref", "alt", "llr_calibrated"), batch_size=100000)
reader <- scanner$ToRecordBatchReader()

chr <- sub("llr_(chr[^.]+)\\.parquet", "\\1", basename(f))
out <- sprintf("tsv/gpn-star-M.%s.tsv.gz", chr)

first <- TRUE

repeat {
  batch <- reader$read_next_batch()
  if(is.null(batch)) break
  
  x <- as.data.table(batch)[, .(chrom, pos, ref, alt, score=llr_calibrated)]
  x[, chrom := as.character(chrom)]
  
  fwrite(x, out, sep="\t", col.names=first, append=!first)
  
  first <- FALSE
}

cat("Wrote:", out, "\n")