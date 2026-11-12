# SCAD

This directory contains the SCAD implementation used in the benchmark. The pipeline supports both non-SMOTE and SMOTE training and evaluates three feature settings: `norm`, `ppi`, and `tp4k`.

## Required SCAD Data

The SCAD `data/` directory is provided separately through the `scDrugBenchmark_env` archive because of its large size.

Download `scDrugBenchmark_env` from Zenodo:

**Zenodo:** https://doi.org/10.5281/zenodo.23260161

After extraction, locate:

```text
scDrugBenchmark_env/
└── SCAD/
    └── SCAD/
        └── data/
```

Copy the entire `data/` directory into the corresponding SCAD directory in the GitHub repository:

```bash
cp -r scDrugBenchmark_env/methods/SCAD/SCAD/data \
      scDrugBenchmark_github/methods/SCAD/SCAD/
```

After copying, the directory should contain:

```text
scDrugBenchmark_github/
└── methods/
    └── SCAD/
        └── SCAD/
            └── data/
```

This step must be completed before setting up the Docker environment and running SCAD.

## Input

SCAD requires:

- a raw count matrix in CSV format, with genes as rows and cells as columns;
- a metadata CSV containing the response label column.

Example input:

```text
data/example/balanced_1000/scdeal_input/
├── balanced_sampled_GSE108397_1round_seurat_afterAnno_raw_count_matrix.csv
└── balanced_sampled_GSE108397_1round_seurat_afterAnno_meta_data.csv
```

## Docker Environment

SCAD uses the shared Python Docker environment in:

```text
methods/shared_main_env/
```

Build the image from `methods/shared_main_env/`:

```bash
sudo docker build -t scdrugbenchmark-main .
```

Start the container from the root directory of the `scDrugBenchmark_github` repository:

```bash
sudo docker run --rm -it \
-v "$(pwd)":/workspace \
scdrugbenchmark-main
```

Inside the container, navigate to SCAD:

```bash
cd /workspace/methods/SCAD
```

Then follow the commands in the Run section below.

## Run

From `methods/SCAD/`:

```bash
chmod +x run_scad_example.sh

./run_scad_example.sh \
--id GSE108397 \
--drug PLX4720 \
--count ../../data/example/balanced_1000/scdeal_input/balanced_sampled_GSE108397_1round_seurat_afterAnno_raw_count_matrix.csv \
--meta ../../data/example/balanced_1000/scdeal_input/balanced_sampled_GSE108397_1round_seurat_afterAnno_meta_data.csv \
--output-dir ./output/GSE108397_PLX4720
```

Default labels are:

```text
condition
resistant
sensitive
```

Alternative labels can be specified as:

```bash
./run_scad_example.sh \
--id DATASET_ID \
--drug DRUG \
--count counts.csv \
--meta metadata.csv \
--label-col condition \
--res-label pre-resistant \
--sens-label pre-sensitive
```

Main optional arguments:

```text
--solid-tumor 0|1
--label-col NAME
--res-label LABEL
--sens-label LABEL
--max-iter N
--gpu 0|1
--overwrite
```

Use the following command to see all available options:

```bash
./run_scad_example.sh --help
```

## Output

Each run is stored in the directory specified by `--output-dir`.

Final per-cell predictions are saved as:

```text
*_predicted_y.tsv
```

Both non-SMOTE and SMOTE predictions are generated separately for the `norm`, `ppi`, and `tp4k` feature settings.