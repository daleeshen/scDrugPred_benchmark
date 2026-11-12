```markdown
# scGPT

This directory contains the scripts used to run scGPT for single-cell drug response prediction.

## Run a single dataset

Use:

    bash ./run_single_example.sh MODEL_ROOT INPUT_FILE

where:

- `MODEL_ROOT` is the path to the `scGPT-main/scgpt` directory.
- `INPUT_FILE` is the input tab-separated TXT file.

Example:

    bash ./run_single_example.sh \
        /media/shenli/Expansion/DRMREF/scDrugBenchmark/foundation_models/scGPT-main/scgpt \
        /media/shenli/Expansion/DRMREF/scDrugBenchmark/foundation_models/scGPT-main/scdata/GSE108397.txt

The pretrained scGPT model files should be located under:

    MODEL_ROOT/save/

The single-run workflow consists of two steps:

    stp1_get_ebd_auto.py
        Input TXT -> scGPT cell embeddings

    stp2_benchmarking_main_EBD_auto.py
        Cell embeddings -> cross-validation results and cell-level predictions

Each run creates timestamped output directories:

    MODEL_ROOT/samples/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/logs/single_run_YYYYMMDD_HHMMSS/

The main outputs are:

    MODEL_ROOT/samples/single_run_YYYYMMDD_HHMMSS/<dataset>_samples.pkl
    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/<dataset>.csv
    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/<dataset>_predictions.csv
    MODEL_ROOT/logs/single_run_YYYYMMDD_HHMMSS/<dataset>.log
```
