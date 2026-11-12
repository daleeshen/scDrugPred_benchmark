# scDEAL Ablation Analysis

This folder contains the downstream analysis code for the scDEAL ablation experiments used in Figure 5E and Figure 5F.

The ablation models are run from the scDEAL method directory, while the resulting prediction files are summarized and visualized using the R script in this folder.

## Files

analyses/ablation/
- ablation_analysis.R
- README.md

methods/scDEAL/scDEAL/
- run_ablation.sh

## Ablation settings

The full scDEAL model uses:

- integrated bulk data
- DAE pretraining
- cluster regularization
- transfer learning

The following ablations are evaluated:

- no_cluster: remove cluster regularization
- no_dae: replace DAE with AE
- gdsc_only: use only GDSC bulk data
- ccle_only: use only CCLE bulk data
- no_transfer: use the predictions generated before transfer learning

The no-transfer condition does not require a separate model run. It is obtained from the pre-transfer prediction columns:

sens_preds_pret
sens_label_pret

## Step 1. Run the ablation models

Go to the scDEAL method directory:

    cd methods/scDEAL/scDEAL

Run:

    bash run_ablation.sh

The script first identifies the best bottleneck and sampling parameters selected for the full scDEAL model for each dataset, drug, round, and evaluation setting.

These parameters are then fixed when running the ablation models.

The following settings are evaluated:

- balanced
- 1:3 imbalance
- 1:30 imbalance

Each setting is evaluated across three sampling rounds.

The script runs:

- no_cluster
- no_dae
- gdsc_only
- ccle_only

and saves the corresponding prediction results.

## Step 2. Analyze the ablation results

After all ablation runs are completed, go to:

    cd analyses/ablation

Run:

    Rscript ablation_analysis.R

The R script:

- reads all ablation prediction results
- calculates performance metrics
- checks missing and duplicated results
- combines the ablation results with the full scDEAL benchmark results
- calculates paired performance changes as:

    Delta = Full - Ablated

- summarizes component contributions
- generates the Figure 5E ablation plots
- generates the component-importance heatmap
- generates the mean-rank/lollipop plot
- generates the comprehensive importance ranking matrix

The comprehensive component analysis also incorporates the label-definition effect obtained from the label-substitution analysis.

## Important: update paths before running

The scripts currently contain absolute paths from the original analysis environment.

Before running the analysis on another system, update the paths in both:

    methods/scDEAL/scDEAL/run_ablation.sh

and:

    analyses/ablation/ablation_analysis.R

In particular, check and update paths for:

- the scDEAL source-code directory
- the Python executable
- bulk expression and label files
- balanced and imbalanced single-cell datasets
- dataset-drug mapping file
- full-model benchmark result files
- ablation output directories
- label-substitution result files used in the comprehensive ranking analysis

The expected directory structure may differ across systems, so all absolute paths should be checked before execution.

## Recommended execution order

    1. Run the main scDEAL benchmark and obtain the full-model results.
    2. Run the label-substitution analysis if the comprehensive importance ranking is required.
    3. Run methods/scDEAL/scDEAL/run_ablation.sh.
    4. Run analyses/ablation/ablation_analysis.R.
    5. Collect the generated summary tables and figures.
