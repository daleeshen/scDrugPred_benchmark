# Beyondcell

This directory contains the scripts used to run Beyondcell in the benchmark.

## Docker

A Docker environment is provided for Beyondcell to reproduce the software environment used in the benchmark.

From the `methods/Beyondcell/` directory, build the Docker image with:

    sudo docker build -t scdrugbenchmark-beyondcell .

Then start the container from the repository root:

    sudo docker run --rm -it \
      -v "$(pwd)":/workspace \
      scdrugbenchmark-beyondcell

The container starts directly in the `beyondcell` Conda environment.

You can verify the environment with:

    which R
    R --version
    R -e "library(beyondcell); packageVersion('beyondcell')"

After entering the container, go to:

    cd /workspace/methods/Beyondcell

and run the example workflow as described below.

The Docker image is built from:

    Dockerfile
    beyondcell_environment.yml

The Docker image itself is not included in the repository and should be built locally from these files.

## Run the Example

The example workflow includes three balanced rounds, three 1:3 imbalanced rounds, and three 1:30 imbalanced rounds for PLX-4720.

From the `methods/Beyondcell/` directory, run:

```bash
./run_beyondcell_example.sh
```

The script automatically reads the example Seurat objects from:

```text
../../data/example/
```

Predictions are saved to:

```text
example_output/
├── balanced_1000/
├── imbalanced_1_3/
└── imbalanced_1_30/
```

To run Beyondcell on an individual dataset:

```bash
Rscript run_beyondcell.R <input_file> <drug> <output_dir> <PSc|SSc>
```

For example:

```bash
Rscript run_beyondcell.R \
    ../../data/example/balanced_1000/balanced_sampled_GSE108397_1round_seurat_afterAnno.RData \
    PLX-4720 \
    ./example_output/balanced_1000 \
    PSc
```
