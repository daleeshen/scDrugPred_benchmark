```markdown
# scFoundation

This directory contains the scripts used to run scFoundation for single-cell drug response prediction.

## Run a single dataset

Use:

    bash run_single_example.sh MODEL_ROOT INPUT_FILE

where:

- `MODEL_ROOT` is the path to the `scFoundation-main/model` directory.
- `INPUT_FILE` is the input tab-separated TXT file.

Example:

    bash run_single_example.sh \
        /path/to/scFoundation-main/model \
        /path/to/GSE108397.txt

An example data is here: scFoundation-main/model/scdata/GSE108397.txt

To specify a Python executable:

    PYTHON_BIN=/path/to/python bash run_single_example.sh \
        /path/to/scFoundation-main/model \
        /path/to/GSE108397.txt

The required scFoundation files should be located at:

    MODEL_ROOT/OS_scRNA_gene_index.19264.tsv
    MODEL_ROOT/models/models.ckpt

Each run creates timestamped output directories such as:

    MODEL_ROOT/samples/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/labels/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/embeds/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/
    MODEL_ROOT/logs/single_run_YYYYMMDD_HHMMSS/

The final evaluation results are saved in:

    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/<dataset>.csv
    MODEL_ROOT/results/single_run_YYYYMMDD_HHMMSS/<dataset>_predictions.csv
```
