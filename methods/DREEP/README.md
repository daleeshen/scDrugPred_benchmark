# DREEP

This directory contains the scripts used to run DREEP in the benchmark.

## Docker Environment

DREEP uses the shared R 4.4.0 Docker environment located in:

    methods/shared_r440_env/

Before running DREEP, build the shared Docker image from the `methods/shared_r440_env/` directory:

    sudo docker build -t scdrugbenchmark-r440 .

Then start the container from the root directory of the `scDrugBenchmark` repository:

    sudo docker run --rm -it \
      -v "$(pwd)":/workspace \
      scdrugbenchmark-r440

The container automatically activates the `r440` Conda environment.

Inside the container, move to the DREEP directory:

    cd /workspace/methods/DREEP

You can verify that DREEP is available with:

    R -e 'library(DREEP); packageVersion("DREEP")'

## Run the Example

The example workflow includes three balanced rounds, three 1:3 imbalanced rounds, and three 1:30 imbalanced rounds for PLX-4720.

From the `methods/DREEP/` directory, run:

```bash
./run_dreep_example.sh
```

Predictions are saved to:

```text
example_output/
├── balanced_1000/
├── imbalanced_1_3/
└── imbalanced_1_30/
```

To run DREEP on an individual dataset:

```bash
Rscript run_dreep.R <input_file> <drug> <output_dir> <reference_rds>
```

