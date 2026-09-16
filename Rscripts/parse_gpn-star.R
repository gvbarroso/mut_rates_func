library(arrow)
library(data.table)

# downloaded on 16/09/2026
# https://huggingface.co/datasets/songlab/gpn-star-scores

args <- commandArgs(trailingOnly=TRUE)
f <- args[1]
# e.g. f <- "M/3bd187238b5a3e66bf2f99a801d7f0a5a84e4522595a82be73e42e5e903e652a"
# Rscript M/3bd187238b5a3e66bf2f99a801d7f0a5a84e4522595a82be73e42e5e903e652a

ds <- open_dataset(f, format="parquet")

scanner <- Scanner$create(ds, columns=c("chrom", "pos", "ref", "alt", "llr_calibrated"), batch_size=100000)
reader <- scanner$ToRecordBatchReader()
files <- list()

repeat {
  batch <- reader$read_next_batch()
  if(is.null(batch)) break
  
  x <- as.data.table(batch)[, .(chrom, pos, ref, alt, llr_calibrated)]
  x[, chrom := as.character(chrom)]
  
  for(chr in unique(x$chrom)) {
    if(is.null(files[[chr]])) {
      files[[chr]] <- sprintf("gpn-star-M.chr%s.tsv.gz", chr)
      fwrite(x[chrom == chr], files[[chr]], sep="\t", col.names=FALSE)
    } else {
      fwrite(x[chrom == chr], files[[chr]], sep="\t", col.names=FALSE, append=TRUE)
    }
    
    cat(paste("wrote chr", chr, "\n"))
  }
}