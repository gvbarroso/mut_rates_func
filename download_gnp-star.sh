mkdir -p gpn-star-scores/{V,M,P}

aria2c -c -x 16 -s 16 -d gpn-star-scores/V \
  https://huggingface.co/datasets/songlab/gpn-star-scores/resolve/main/data/gpn-star-hg38-v100-200m/llr/llr_chr{1..22}.parquet

aria2c -c -x 16 -s 16 -d gpn-star-scores/M \
  https://huggingface.co/datasets/songlab/gpn-star-scores/resolve/main/data/gpn-star-hg38-m447-200m/llr/llr_chr{1..22}.parquet

aria2c -c -x 16 -s 16 -d gpn-star-scores/P \
  https://huggingface.co/datasets/songlab/gpn-star-scores/resolve/main/data/gpn-star-hg38-p243-200m/llr/llr_chr{1..22}.parquet 
