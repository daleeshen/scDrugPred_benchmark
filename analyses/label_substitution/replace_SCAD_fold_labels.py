import os
import shutil
import pandas as pd
from pathlib import Path


# ============================================================
# Paths
# ============================================================

PAIR_FILE = (
    "/media/shenli/Expansion/DRMREF/New_fig5/"
    "Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
)

SCAD_ROOT = Path(
    "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source"
)

SCDEAL_ROOT = Path(
    "/media/shenli/Expansion/DRMREF/New_fig5/"
    "SCAD_common_source_scDEAL_label"
)


# ============================================================
# Drug name conversion
# ============================================================

def convert_drug(drug):
    mapping = {
        "PLX-4720": "PLX4720",
        "Palbociclib": "PD-0332991",
    }
    return mapping.get(drug, drug)


# ============================================================
# Scenarios
# ============================================================

scenarios = [
    "balanced",
    "imbalance_1_3",
    "imbalance_1_30"
]

rounds = [
    "round_1",
    "round_2",
    "round_3"
]

genesets = {
    "all": {
        "source_suffix": "",
        "strat_dir": "stratified_norm"
    },
    "ppi": {
        "source_suffix": "_ppi",
        "strat_dir": "stratified_ppi"
    },
    "tp4k": {
        "source_suffix": "_tp4k",
        "strat_dir": "stratified_tp4k"
    }
}


# ============================================================
# Dataset-drug pairs
# ============================================================

pairs = pd.read_csv(PAIR_FILE)

dataset_col = pairs.columns[0]
drug_col = pairs.columns[1]


# ============================================================
# Main
# ============================================================

summary = []

for scenario in scenarios:

    for round_name in rounds:

        for _, row in pairs.iterrows():

            dataset = str(row[dataset_col]).strip()
            drug_original = str(row[drug_col]).strip()
            drug = convert_drug(drug_original)

            folder_name = f"{dataset}_{drug}"

            scad_folder = (
                SCAD_ROOT /
                scenario /
                round_name /
                folder_name
            )

            scdeal_folder = (
                SCDEAL_ROOT /
                scenario /
                round_name /
                folder_name
            )

            print("\n==============================================")
            print(f"{scenario} | {round_name}")
            print(f"{dataset} | {drug}")
            print("==============================================")


            for geneset_name, info in genesets.items():

                suffix = info["source_suffix"]
                strat_name = info["strat_dir"]

                # ------------------------------------------------
                # Source file containing scDEAL labels
                # ------------------------------------------------

                source_scdeal_file = (
                    scdeal_folder /
                    f"Source_exprs_resp_z.{drug}{suffix}.tsv"
                )

                if not source_scdeal_file.exists():
                    raise FileNotFoundError(
                        f"Missing scDEAL Source file:\n"
                        f"{source_scdeal_file}"
                    )


                # ------------------------------------------------
                # Original SCAD stratified folds
                # ------------------------------------------------

                src_strat = scad_folder / strat_name
                dst_strat = scdeal_folder / strat_name

                if not src_strat.exists():
                    raise FileNotFoundError(
                        f"Missing SCAD stratified folder:\n"
                        f"{src_strat}"
                    )


                # ------------------------------------------------
                # Copy whole stratified folder
                # ------------------------------------------------

                shutil.copytree(
                    src_strat,
                    dst_strat,
                    dirs_exist_ok=True
                )


                # ------------------------------------------------
                # Read scDEAL-label source
                # ------------------------------------------------

                source_scdeal = pd.read_csv(
                    source_scdeal_file,
                    sep="\t",
                    index_col=0
                )

                source_scdeal.index = (
                    source_scdeal.index.astype(str)
                )

                if "response" not in source_scdeal.columns:
                    raise ValueError(
                        f"'response' column missing:\n"
                        f"{source_scdeal_file}"
                    )


                # ------------------------------------------------
                # Find source fold directory
                #
                # Handles names such as:
                # source_5_folds_norm
                # source_5_folds__ppi
                # source_5_folds__tp4k
                # ------------------------------------------------

                source_fold_dirs = [
                    p for p in dst_strat.iterdir()
                    if p.is_dir()
                    and p.name.startswith("source_5_folds")
                ]

                if len(source_fold_dirs) != 1:
                    raise RuntimeError(
                        f"Expected exactly one source fold dir in:\n"
                        f"{dst_strat}\n"
                        f"Found: {source_fold_dirs}"
                    )

                source_fold_dir = source_fold_dirs[0]


                # ------------------------------------------------
                # Replace Y only
                # ------------------------------------------------

                for split_num in range(1, 6):

                    split_dir = (
                        source_fold_dir /
                        f"split{split_num}"
                    )

                    for phase in ["train", "val"]:

                        x_file = (
                            split_dir /
                            f"X_{phase}_source.tsv"
                        )

                        y_file = (
                            split_dir /
                            f"Y_{phase}_source.tsv"
                        )

                        X = pd.read_csv(
                            x_file,
                            sep="\t",
                            index_col=0
                        )

                        X.index = X.index.astype(str)

                        missing_ids = (
                            set(X.index) -
                            set(source_scdeal.index)
                        )

                        if missing_ids:
                            raise ValueError(
                                f"Missing {len(missing_ids)} IDs "
                                f"in scDEAL source labels:\n"
                                f"{dataset} | {drug} | "
                                f"{scenario} | {round_name} | "
                                f"{geneset_name} | "
                                f"split{split_num} | {phase}"
                            )


                        # exact same IDs and order as X
                        y_new = source_scdeal.loc[
                            X.index,
                            "response"
                        ].copy()

                        y_new.name = "response"

                        y_new.to_csv(
                            y_file,
                            sep="\t",
                            header=True
                        )


                        summary.append({
                            "scenario": scenario,
                            "round": round_name,
                            "dataset": dataset,
                            "drug": drug_original,
                            "geneset": geneset_name,
                            "split": split_num,
                            "phase": phase,
                            "n_samples": len(y_new),
                            "n_R": int((y_new == 0).sum()),
                            "n_S": int((y_new == 1).sum())
                        })


                print(
                    f"OK: {geneset_name}"
                )


# ============================================================
# Save summary
# ============================================================

summary_df = pd.DataFrame(summary)

summary_file = (
    "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source_scDEAL_label/"
    "SCAD_scDEAL_label_fold_summary.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)

print("\n==============================================")
print("ALL DONE")
print(f"Summary: {summary_file}")
print("==============================================")