# scDr

This directory contains the scripts used to run scDr in the benchmark.

## Docker Environment

scDr uses the shared R 4.4.0 Docker environment located in:

    methods/shared_r440_env/

Build the shared Docker image from the `methods/shared_r440_env/` directory:

    sudo docker build -t scdrugbenchmark-r440 .

Then start the container from the root directory of the `scDrugBenchmark` repository:

    sudo docker run --rm -it \
      -v "$(pwd)":/workspace \
      scdrugbenchmark-r440

The container automatically activates the `r440` Conda environment.

Inside the container, move to the scDr directory:

    cd /workspace/methods/scDr

## Run the Example

The example workflow includes three balanced rounds, three 1:3 imbalanced rounds, and three 1:30 imbalanced rounds.

From the `methods/scDr/` directory, run:

    ./run_scdr_example.sh

Predictions are saved to:

    example_output/
    ├── balanced_1000/
    ├── imbalanced_1_3/
    └── imbalanced_1_30/

The drug signatures required by scDr are provided in:

    scDr/drug_up_signature.RDS
    scDr/drug_down_signature.RDS

To run scDr on an individual dataset:

    Rscript run_scDr.R <input_file> <output_dir> <up_signature> <down_signature>
