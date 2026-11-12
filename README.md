# scDrugBenchmark

This repository contains the code used in the study:

**Predicting Pre-treatment Resistance or Post-treatment Effect? A Systematic Benchmarking of Single-Cell Drug Response Models**

The repository provides scripts for benchmark dataset construction, model execution, performance evaluation, lineage-tracing analyses, training-label substitution, scDEAL ablation analysis, PDAC analysis, complementary evaluation of single-cell foundation models, and figure generation.

Large benchmark datasets, pretrained model files, reference data, and runtime environments are distributed separately through Zenodo to keep the GitHub repository lightweight and easier to navigate.

## Repository Structure

```text
scDrugPred_benchmark/
├── analyses/
│   ├── label_substitution/
│   ├── lineage_tracing/
│   │   ├── ReSisTrace/
│   │   └── Watermelon/
│   ├── PDAC/
│   └── scDEAL_ablation/
├── benchmark_construction/
├── data/
├── evaluation/
├── Figure/
├── foundation_models/
│   ├── scfoundation-main/
│   └── scGPT-main/
└── methods/
    ├── Beyondcell/
    ├── CaDRReS-sc/
    ├── DREEP/
    ├── DrugFormer/
    ├── Precily/
    ├── SCAD/
    ├── scDEAL/
    ├── scDr/
    └── scIDUC/
```

### `benchmark_construction/`

Scripts used to generate the balanced and imbalanced benchmark datasets.

### `methods/`

Benchmark implementations and execution scripts for the nine task-specific single-cell drug response prediction methods:

- Beyondcell
- CaDRReS-sc
- DREEP
- DrugFormer
- Precily
- SCAD
- scDEAL
- scDr
- scIDUC

Method-specific instructions, required files, environment setup, input format, and example commands are provided in the corresponding `README.md` files.

### `foundation_models/`

Scripts used for the complementary evaluation of:

- scGPT
- scFoundation

The benchmark includes within-dataset and cross-dataset evaluations. Large pretrained model files and other required resources are distributed separately through Zenodo.

### `analyses/`

Scripts for additional analyses performed in the study:

- `label_substitution/` – training-label comparison and substitution between SCAD and scDEAL
- `lineage_tracing/` – ReSisTrace and Watermelon analyses
- `scDEAL_ablation/` – scDEAL component ablation analysis
- `PDAC/` – downstream analysis of the in-house PDAC organoid dataset

### `evaluation/`

Scripts used to summarize and evaluate prediction results from the balanced and imbalanced benchmark settings.

### `data/`

Scripts used to prepare method-specific example inputs.

The example datasets required to reproduce the demonstration workflows are distributed through Zenodo. See `data/README.md` for details.

### `Figure/`

Scripts used to generate the main and supplementary figures presented in the manuscript.

## Data and Runtime Environments

Large files are not stored directly in this GitHub repository.

Supporting data and runtime environments are archived on Zenodo:

**Zenodo:** https://doi.org/10.5281/zenodo.23260161

The Zenodo record contains:

- example benchmark datasets and processed data;
- large reference files;
- pretrained model files;
- method-specific resources excluded from GitHub because of file size;
- shared Python and R runtime environments;
- method-specific runtime environments required for reproduction.

The main files include:

```text
scDrugBenchmark_data.tar.gz
scDrugBenchmark_env_CaDRReS-sc.tar.gz
scDrugBenchmark_env_DrugFormer.tar.gz
scDrugBenchmark_env_Precily.tar.gz
scDrugBenchmark_env_SCAD.tar.gz
scDrugBenchmark_env_scDEAL.tar.gz
scDrugBenchmark_env_shared_main_env.tar.gz
scDrugBenchmark_env_shared_r440_env.tar.gz
```

Some methods require files from the corresponding Zenodo environment archive to be copied into the GitHub method directory before execution. Detailed instructions are provided in each method-specific `README.md`.

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/daleeshen/scDrugPred_benchmark.git
cd scDrugPred_benchmark
```

### 2. Download the required data and environment files

Download the required archives from:

**https://doi.org/10.5281/zenodo.23260161**

Only the environment archives required for the method of interest need to be downloaded.

After extraction, a typical local organization may look like:

```text
project/
├── scDrugPred_benchmark/
├── data/
└── scDrugBenchmark_env/
```

### 3. Prepare method-specific resources

Some methods require large reference data, pretrained models, or method-specific files distributed through Zenodo.

Follow the instructions in the corresponding method directory before building the environment or running the example:

```text
methods/<METHOD>/README.md
```

### 4. Prepare example input data

Method-specific example inputs can be generated using the scripts in:

```text
data/
```

For example:

```bash
cd data
Rscript prepare_tpm_example.R
Rscript prepare_drugformer_example.R
Rscript prepare_precily_example.R
Rscript prepare_scdeal_example.R
```

Only the preparation script required for the method being tested needs to be executed.

### 5. Run a benchmark method

Each supported method directory contains its own execution instructions and example command.

For example:

```text
methods/CaDRReS-sc/README.md
methods/DrugFormer/README.md
methods/Precily/README.md
methods/SCAD/README.md
methods/scDEAL/README.md
```

Please follow the method-specific README because the required runtime environment and input format differ among methods.

## Benchmark Settings

The benchmark includes:

- balanced datasets with different numbers of sampled cells;
- imbalanced datasets with different resistant/sensitive cell ratios;
- lineage-tracing evaluation using ReSisTrace and Watermelon;
- complementary within-dataset and cross-dataset evaluation using scGPT and scFoundation;
- training-label substitution and scDEAL component ablation analyses.

## Reproducibility

To improve reproducibility and traceability:

- benchmark construction scripts are provided in `benchmark_construction/`;
- method-specific execution scripts are provided in `methods/`;
- example input preparation scripts are provided in `data/`;
- performance evaluation scripts are provided in `evaluation/`;
- downstream and additional analyses are provided in `analyses/`;
- figure-generation scripts are provided in `Figure/`;
- large supporting data and runtime resources are archived on Zenodo.

The method-specific README files describe the files, environments, and execution steps required for each workflow.

## Data Availability

The public datasets used in the benchmark are described in the manuscript.

Processed benchmark data, example datasets, pretrained model files, and additional resources required for reproduction are available from Zenodo:

**https://doi.org/10.5281/zenodo.23260161**

## Code Availability

The benchmarking and analysis code is available in this GitHub repository:

**https://github.com/daleeshen/scDrugPred_benchmark**

Large runtime environments and supporting data are archived separately on Zenodo:

**https://doi.org/10.5281/zenodo.23260161**

## Citation

If you use this benchmark or the accompanying resources, please cite:

> Shen L, et al. *Predicting Pre-treatment Resistance or Post-treatment Effect? A Systematic Benchmarking of Single-Cell Drug Response Models.*

Citation information will be updated upon publication.
