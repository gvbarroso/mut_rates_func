#!/bin/bash

cd /media/gvbarroso/extradrive1/gpn-star-scores

pip install -U huggingface_hub

hf download songlab/gpn-star-scores \
  --repo-type dataset \
  --include "data/gpn-star-hg38-m447-200m/llr/*" \
  --local-dir gpn-star-scores

hf download songlab/gpn-star-scores \
  --repo-type dataset \
  --include "data/gpn-star-hg38-v100-200m/llr/*" \
  --local-dir gpn-star-scores

hf download songlab/gpn-star-scores \
  --repo-type dataset \
  --include "data/gpn-star-hg38-p243-200m/llr/*" \
  --local-dir gpn-star-scores
