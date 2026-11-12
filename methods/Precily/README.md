# Precily

This directory contains the benchmark-specific scripts used to run Precily.

The GitHub repository contains:

```text
methods/Precily/
├── run_precily.R
└── run_precily_example.sh
```

The complete Precily method files are provided separately through Zenodo because of their large size.

## Download the Precily Files

Download the `scDrugBenchmark_env` archive from Zenodo:

**Zenodo:** https://doi.org/10.5281/zenodo.23260161

After extraction, locate:

```text
scDrugBenchmark_env/
└── Precily/
    └── precily/
```

Copy the entire `precily/` directory into:

```text
scDrugBenchmark_github/
└── methods/
    └── Precily/
```

For example:

```bash
cp -r scDrugBenchmark_env/Precily/precily \
      scDrugBenchmark_github/methods/Precily/
```

After copying, the directory should contain:

```text
methods/Precily/
├── precily/
├── run_precily.R
└── run_precily_example.sh
```

This step must be completed before building the Precily Docker environment and running the benchmark.

## Docker Environment

Precily requires two Conda environments inside the Docker container:

- `r440` – R 4.4.0 and the required R/Bioconductor packages
- `precily` – Python and TensorFlow required by Precily

The `r440` environment is used as the main environment, while the Precily R script accesses the Python environment through:

```text
/opt/conda/envs/precily/bin/python
```

Precily is built on top of the shared R 4.4.0 Docker image located in:

```text
methods/shared_r440_env/
```

Therefore, the shared R image must be built before building the Precily image.

### 1. Build the shared R environment

From:

```text
methods/shared_r440_env/
```

run:

```bash
sudo docker build -t scdrugbenchmark-r440 .
```

This creates the base image:

```text
scdrugbenchmark-r440
```

### 2. Build the Precily environment

The Precily Docker image inherits the shared R environment and additionally includes the Precily Python environment.

From:

```text
methods/Precily/
```

run:

```bash
sudo docker build -t scdrugbenchmark-precily .
```

The resulting image contains:

```text
/opt/conda/envs/r440
/opt/conda/envs/precily
```

### 3. Start the Precily container

From the root directory of the `scDrugBenchmark_github` repository, run:

```bash
sudo docker run --rm -it \
-v "$(pwd)":/workspace \
scdrugbenchmark-precily
```

The container automatically activates the `r440` environment.

Inside the container, the prompt should therefore start with:

```text
(r440)
```

The R environment can be checked with:

```bash
which R
```

which should return:

```text
/opt/conda/envs/r440/bin/R
```

The Precily Python environment can be checked with:

```bash
/opt/conda/envs/precily/bin/python --version
```

Inside the container, move to:

```bash
cd /workspace/methods/Precily
```

## Prepare Example Inputs

Precily uses log2(TPM + 1) expression matrices as input.

From the `data/` directory, run:

```bash
Rscript prepare_precily_example.R
```

The prepared inputs are saved under:

```text
data/example/
├── balanced_1000/precily_input/
├── imbalanced_1_3/precily_input/
└── imbalanced_1_30/precily_input/
```

## Run the Example

From the `methods/Precily/` directory, run:

```bash
./run_precily_example.sh
```

The script automatically runs all three rounds for the balanced, 1:3 imbalanced, and 1:30 imbalanced example datasets.

Predictions are saved to:

```text
example_output/
├── balanced_1000/
├── imbalanced_1_3/
└── imbalanced_1_30/
```

To run Precily on an individual dataset:

```bash
Rscript run_precily.R \
<logTPM_file> \
<metadata_file> \
<drug> \
<cancer_type> \
<output_dir>
```