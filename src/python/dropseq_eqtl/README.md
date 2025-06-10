# python tools for tensorQTL workflow

## Installation

Requires python >= 3.12
```
pip install 'git+https://github.com/broadinstitute/Drop-seq.git#egg=dropseq_metadata&subdirectory=src/python/dropseq_eqtl'
```

## Command-line tools

- merge_parquet_files: Merge parquet files together to either another parquet file or a tab-delimited file.
- normalize_tensorqtl_expression: Normalize tensorQTL expression per gene with edgeR CPM then per donor with an inverse normal transformation.
- prepare_tensorqtl_data: Prepare and validate data by removing genes and covariates that would cause errors when running tensorQTL.
