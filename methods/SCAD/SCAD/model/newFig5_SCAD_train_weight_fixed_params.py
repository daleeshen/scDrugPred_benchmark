import os
import re
import random
import argparse
from pathlib import Path
from collections import Counter
from itertools import cycle

import numpy as np
import pandas as pd
import torch
import sklearn.preprocessing as sk

from sklearn.metrics import roc_auc_score, average_precision_score
from torch.utils.data.sampler import WeightedRandomSampler

from SCADmodules import *


# ============================================================
# Basic paths
# ============================================================

PAIR_FILE = Path(
    "/media/shenli/Expansion/DRMREF/New_fig5/"
    "Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
)

DATA_ROOTS = {
    "SCAD_label": Path(
        "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source"
    ),
    "scDEAL_label": Path(
        "/media/shenli/Expansion/DRMREF/New_fig5/"
        "SCAD_common_source_scDEAL_label"
    )
}

OUTPUT_ROOT = Path(
    "/media/shenli/Expansion/DRMREF/New_fig5/"
    "SCAD_fixed_original_params_results"
)


# ============================================================
# Experiment settings
# ============================================================

SCENARIOS = {
    "balanced": None,
    "imbalance_1_3": "1_3",
    "imbalance_1_30": "1_30"
}

ROUNDS = [1, 2, 3]

FEATURES = {

    "norm": {
        "stratified": "stratified_norm",
        "source_dir": "source_5_folds_norm",
        "target_dir": "target_5_folds_norm",
        "source_suffix": ""
    },

    "ppi": {
        "stratified": "stratified_ppi",
        "source_dir": "source_5_folds__ppi",
        "target_dir": "target_5_folds__ppi",
        "source_suffix": "_ppi"
    },

    "tp4k": {
        "stratified": "stratified_tp4k",
        "source_dir": "source_5_folds__tp4k",
        "target_dir": "target_5_folds__tp4k",
        "source_suffix": "_tp4k"
    }
}


# ============================================================
# Parameter filename pattern
# ============================================================

PARAM_PATTERN = re.compile(
    r"^hdim(?P<hdim>\d+)"
    r"_zdim(?P<zdim>\d+)"
    r"_lr(?P<lr>[\d.eE+-]+)"
    r"_epoch(?P<epoch>\d+)"
    r"_lambda1(?P<lambda1>[\d.eE+-]+)"
    r"_mbS(?P<mbS>\d+)"
    r"_mbT(?P<mbT>\d+)"
    r"_seed(?P<seed>\d+)\.tsv$"
)


# ============================================================
# Args
# ============================================================

def get_args():

    parser = argparse.ArgumentParser(
        description="Run SCAD using original best parameters."
    )

    parser.add_argument(
        "--gpu",
        type=int,
        default=1,
        choices=[0, 1]
    )

    parser.add_argument(
        "--label_mode",
        type=str,
        default="both",
        choices=["both", "SCAD_label", "scDEAL_label"]
    )

    return parser.parse_args()


# ============================================================
# Utilities
# ============================================================

def convert_drug(drug):

    mapping = {
        "PLX-4720": "PLX4720",
        "Palbociclib": "PD-0332991"
    }

    return mapping.get(drug, drug)


def set_seed(seed):

    random.seed(seed)
    np.random.seed(seed)

    os.environ["PYTHONHASHSEED"] = str(seed)

    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# ============================================================
# Locate original parameter result directory
# ============================================================

def get_param_root(scenario, round_num, feature):

    # --------------------------------------------------------
    # Balanced
    # --------------------------------------------------------

    if scenario == "balanced":

        # Round 1
        if round_num == 1:

            return Path(
                "/media/shenli/Expansion/SCAD/"
                "1st_round_result_history/"
                "best_param_res_weight/"
                f"best_param_results_{feature}"
            )

        # Round 2 / 3
        else:

            return Path(
                f"/media/shenli/Expansion/DRMREF/"
                f"{round_num}_round_1000cells_res/"
                f"1000cells/SCAD/"
                f"best_param_results_{feature}"
            )

    # --------------------------------------------------------
    # Imbalanced
    # --------------------------------------------------------

    ratio = SCENARIOS[scenario]

    return Path(
        "/media/shenli/Expansion/DRMREF/"
        "unbalanced_test_results/"
        f"{round_num}_round_unbalanced_results/"
        f"{ratio}/SCAD/"
        "best_param_res_weight/"
        f"best_param_res_weight_{feature}"
    )


# ============================================================
# Find parameter TSV
# ============================================================

def find_parameter_files(
    scenario,
    round_num,
    feature,
    dataset,
    drug
):

    root = get_param_root(
        scenario,
        round_num,
        feature
    )

    dataset_folder = root / f"{dataset}_{drug}"

    # fallback: dataset_*drug
    if not dataset_folder.exists():

        folders = list(
            root.glob(f"{dataset}_*{drug}")
        )

        if len(folders) != 1:
            raise FileNotFoundError(
                f"\nCannot uniquely find dataset folder:\n"
                f"root = {root}\n"
                f"dataset = {dataset}\n"
                f"drug = {drug}\n"
                f"found = {folders}"
            )

        dataset_folder = folders[0]


    test_root = (
        dataset_folder /
        "bimodal" /
        "best_models" /
        "FX" /
        "test"
    )

    if not test_root.exists():
        raise FileNotFoundError(
            f"Missing test directory:\n{test_root}"
        )


    candidates = []

    for file in test_root.glob("mbS*_mbT*/*.tsv"):

        match = PARAM_PATTERN.match(file.name)

        if not match:
            continue

        p = match.groupdict()

        params = {
            "hdim": int(p["hdim"]),
            "zdim": int(p["zdim"]),
            "lr": float(p["lr"]),
            "epoch": int(p["epoch"]),
            "lambda1": float(p["lambda1"]),
            "mbS": int(p["mbS"]),
            "mbT": int(p["mbT"]),
            "seed": int(p["seed"])
        }

        candidates.append(
            (file, params)
        )


    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No parameter TSV found:\n{test_root}"
        )


    # normally use original seed42 results
    seed42 = [
        x for x in candidates
        if x[1]["seed"] == 42
    ]

    if len(seed42) > 0:
        candidates = seed42


    print(
        f"Found {len(candidates)} parameter set(s) "
        f"for {dataset}_{drug} | {feature}"
    )

    for file, params in candidates:
        print(f"  {file.name}")
        print(f"  {params}")


    return candidates


# ============================================================
# Prediction / evaluation
# ============================================================

def predict_label(
    XTestCells,
    gen_model,
    map_model,
    device
):

    gen_model.eval()
    map_model.eval()

    gen_model.to(device)
    map_model.to(device)

    with torch.no_grad():

        F_xt_test = gen_model(
            XTestCells.to(device)
        )

        yhatt_test = map_model(
            F_xt_test
        )

    return yhatt_test


def evaluate_model(
    XTestCells,
    YTestCells,
    gen_model,
    map_model,
    device
):

    XTestCells = XTestCells.to(device)
    YTestCells = YTestCells.to(device)

    y_predicted = predict_label(
        XTestCells,
        gen_model,
        map_model,
        device
    )

    yt = (
        YTestCells
        .view(-1)
        .cpu()
        .numpy()
    )

    yp = (
        y_predicted
        .view(-1)
        .cpu()
        .numpy()
    )

    auc = roc_auc_score(
        yt,
        yp
    )

    apr = average_precision_score(
        yt,
        yp
    )

    return auc, apr, yp


# ============================================================
# Run one dataset
# ============================================================

def run_one(
    device,
    label_mode,
    scenario,
    round_num,
    dataset,
    drug_original,
    feature,
    params,
    parameter_file
):

    drug = convert_drug(
        drug_original
    )

    folder_name = f"{dataset}_{drug}"

    feature_info = FEATURES[feature]

    data_root = (
        DATA_ROOTS[label_mode] /
        scenario /
        f"round_{round_num}" /
        folder_name
    )

    # --------------------------------------------------------
    # Parameters
    # --------------------------------------------------------

    h_dim = params["hdim"]
    z_dim = params["zdim"]

    lr = params["lr"]
    epoch = params["epoch"]
    lam1 = params["lambda1"]

    mbS = params["mbS"]
    mbT = params["mbT"]

    seed = params["seed"]


    # --------------------------------------------------------
    # Skip if already finished
    # --------------------------------------------------------

    param_str = (
        f"hdim{h_dim}"
        f"_zdim{z_dim}"
        f"_lr{lr}"
        f"_epoch{epoch}"
        f"_lambda1{lam1}"
        f"_mbS{mbS}"
        f"_mbT{mbT}"
        f"_seed{seed}"
    )

    output_dir = (
            OUTPUT_ROOT /
            label_mode /
            scenario /
            f"round_{round_num}" /
            folder_name /
            feature
    )

    metrics_file = output_dir / f"{param_str}_metrics.tsv"
    pred_file = output_dir / f"{param_str}_predicted_y.tsv"

    if metrics_file.exists() and pred_file.exists():
        print(
            f"SKIP: already finished | "
            f"{label_mode} | {scenario} | round_{round_num} | "
            f"{folder_name} | {feature} | {param_str}"
        )
        return None


    # --------------------------------------------------------
    # Fold paths
    # --------------------------------------------------------

    stratified_path = (
        data_root /
        feature_info["stratified"]
    )

    source_fold_root = (
        stratified_path /
        feature_info["source_dir"]
    )

    target_fold_root = (
        stratified_path /
        feature_info["target_dir"]
    )


    # --------------------------------------------------------
    # Complete source data for weight calculation
    # --------------------------------------------------------

    source_file = (
        data_root /
        (
            f"Source_exprs_resp_z."
            f"{drug}"
            f"{feature_info['source_suffix']}.tsv"
        )
    )

    if not source_file.exists():

        raise FileNotFoundError(
            f"Missing source file:\n{source_file}"
        )


    YGDSC = pd.read_csv(
        source_file,
        sep="\t",
        index_col=0
    )

    counts = Counter(
        YGDSC["response"].astype(int)
    )

    if 0 not in counts or 1 not in counts:

        raise ValueError(
            f"Source contains only one class:\n"
            f"{source_file}\n"
            f"{counts}"
        )


    class_sample_count = np.array([
        counts[0] / len(YGDSC),
        counts[1] / len(YGDSC)
    ])

    print("\n--------------------------------------------------")
    print(
        f"{label_mode} | "
        f"{scenario} | "
        f"round_{round_num} | "
        f"{folder_name} | "
        f"{feature}"
    )
    print(f"Parameter file: {parameter_file}")
    print(f"Parameters: {params}")
    print("--------------------------------------------------")


    # IMPORTANT:
    # Reset seed separately for each label condition.
    set_seed(seed)


    split_metrics = []
    predictions = []


    # ========================================================
    # 5 folds
    # ========================================================

    for split_num in range(1, 6):

        split = f"split{split_num}"

        print(f"Processing {split}...")


        # ----------------------------------------------------
        # Source
        # ----------------------------------------------------

        source_split = (
            source_fold_root /
            split
        )

        xtr_s = pd.read_csv(
            source_split /
            "X_train_source.tsv",
            sep="\t",
            index_col=0
        )

        ytr_s = pd.read_csv(
            source_split /
            "Y_train_source.tsv",
            sep="\t",
            index_col=0
        )


        # ----------------------------------------------------
        # Target
        # ----------------------------------------------------

        target_split = (
            target_fold_root /
            split
        )

        xtr_t = pd.read_csv(
            target_split /
            "X_train_target.tsv",
            sep="\t",
            index_col=0
        )

        ytr_t = pd.read_csv(
            target_split /
            "Y_train_target.tsv",
            sep="\t",
            index_col=0
        )

        xte_t = pd.read_csv(
            target_split /
            "X_test_target.tsv",
            sep="\t",
            index_col=0
        )

        yte_t = pd.read_csv(
            target_split /
            "Y_test_target.tsv",
            sep="\t",
            index_col=0
        )


        # ----------------------------------------------------
        # Standardization
        # ----------------------------------------------------

        scaler = sk.StandardScaler().fit(
            pd.concat([
                xtr_s,
                xtr_t
            ]).values
        )

        xtr_s_n = scaler.transform(
            xtr_s.values
        )

        xtr_t_n = scaler.transform(
            xtr_t.values
        )

        xte_t_n = scaler.transform(
            xte_t.values
        )


        # ----------------------------------------------------
        # Weighted sampler
        # ----------------------------------------------------

        source_labels = (
            ytr_s.iloc[:, 0]
            .astype(int)
            .values
        )

        sample_weights = np.array([
            1.0 / class_sample_count[label]
            for label in source_labels
        ])

        sampler = WeightedRandomSampler(
            torch.DoubleTensor(sample_weights),
            len(sample_weights),
            replacement=True
        )


        loader_s = torch.utils.data.DataLoader(

            torch.utils.data.TensorDataset(
                torch.FloatTensor(xtr_s_n),
                torch.FloatTensor(ytr_s.values)
            ),

            batch_size=mbS,
            sampler=sampler
        )


        loader_t = torch.utils.data.DataLoader(

            torch.utils.data.TensorDataset(
                torch.FloatTensor(xtr_t_n),
                torch.FloatTensor(ytr_t.values)
            ),

            batch_size=mbT,
            shuffle=True
        )


        # ----------------------------------------------------
        # Model
        # ----------------------------------------------------

        gen = FX(
            0.5,
            xtr_s_n.shape[1],
            h_dim,
            z_dim
        ).to(device)


        map_net = MTLP(
            0.5,
            h_dim,
            z_dim
        ).to(device)


        dg = Discriminator(
            0.5,
            h_dim,
            z_dim
        ).to(device)


        opt = torch.optim.Adagrad(

            list(gen.parameters()) +
            list(map_net.parameters()) +
            list(dg.parameters()),

            lr=lr
        )


        criterion = torch.nn.BCELoss()


        # ----------------------------------------------------
        # Training
        # ----------------------------------------------------

        for it in range(epoch):

            gen.train()
            map_net.train()
            dg.train()


            for ds, dt in zip(
                cycle(loader_s),
                loader_t
            ):

                xs = ds[0].to(device)

                ys = (
                    ds[1]
                    .view(-1, 1)
                    .to(device)
                )

                xt = dt[0].to(device)


                if (
                    xs.size(0) < 5
                    or xt.size(0) < 5
                ):
                    continue


                f_xs = gen(xs)
                f_xt = gen(xt)


                # classification loss
                loss_c = criterion(
                    map_net(f_xs),
                    ys
                )


                # domain loss
                domain_label = torch.cat([

                    torch.ones(
                        f_xs.size(0),
                        1
                    ),

                    torch.zeros(
                        f_xt.size(0),
                        1
                    )

                ]).to(device)


                loss_dg = criterion(

                    dg(
                        torch.cat(
                            [f_xs, f_xt],
                            dim=0
                        )
                    ),

                    domain_label
                )


                loss = (
                    loss_c +
                    lam1 * loss_dg
                )


                opt.zero_grad()

                loss.backward()

                opt.step()


        # ----------------------------------------------------
        # Test
        # ----------------------------------------------------

        auc, apr, pred_y = evaluate_model(

            torch.FloatTensor(
                xte_t_n
            ),

            torch.FloatTensor(
                yte_t.values
            ),

            gen,
            map_net,
            device
        )


        split_metrics.append({

            "split": split,
            "AUC": auc,
            "APR": apr
        })


        for cell_id, true_y, pred in zip(

            yte_t.index,

            yte_t.iloc[:, 0].values,

            pred_y
        ):

            predictions.append({

                "split": split,
                "cell_id": cell_id,
                "true_label": true_y,
                "prediction": pred
            })


        print(
            f"{split}: "
            f"AUC={auc:.4f}, "
            f"APR={apr:.4f}"
        )


    # ========================================================
    # Save
    # ========================================================

    param_str = (

        f"hdim{h_dim}"
        f"_zdim{z_dim}"
        f"_lr{lr}"
        f"_epoch{epoch}"
        f"_lambda1{lam1}"
        f"_mbS{mbS}"
        f"_mbT{mbT}"
        f"_seed{seed}"
    )


    output_dir = (

        OUTPUT_ROOT /
        label_mode /
        scenario /
        f"round_{round_num}" /
        folder_name /
        feature
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    metrics_df = pd.DataFrame(
        split_metrics
    )

    predictions_df = pd.DataFrame(
        predictions
    )


    metrics_df.to_csv(

        output_dir /
        f"{param_str}_metrics.tsv",

        sep="\t",
        index=False
    )


    predictions_df.to_csv(

        output_dir /
        f"{param_str}_predicted_y.tsv",

        sep="\t",
        index=False
    )


    mean_auc = metrics_df["AUC"].mean()
    std_auc = metrics_df["AUC"].std(ddof=0)

    mean_apr = metrics_df["APR"].mean()
    std_apr = metrics_df["APR"].std(ddof=0)


    return {

        "label_mode": label_mode,
        "scenario": scenario,
        "round": round_num,
        "dataset": dataset,
        "drug": drug_original,
        "drug_SCAD": drug,
        "feature": feature,

        "hdim": h_dim,
        "zdim": z_dim,
        "lr": lr,
        "epoch": epoch,
        "lambda1": lam1,
        "mbS": mbS,
        "mbT": mbT,
        "seed": seed,

        "mean_AUC": mean_auc,
        "sd_AUC": std_auc,

        "mean_APR": mean_apr,
        "sd_APR": std_apr,

        "original_parameter_file":
            str(parameter_file)
    }


# ============================================================
# Main
# ============================================================

def main():

    args = get_args()

    device = (
        "cuda"
        if args.gpu == 1
        and torch.cuda.is_available()
        else "cpu"
    )

    print(
        f"Using device: {device}"
    )


    # --------------------------------------------------------
    # label conditions
    # --------------------------------------------------------

    if args.label_mode == "both":

        label_modes = [
            "SCAD_label",
            "scDEAL_label"
        ]

    else:

        label_modes = [
            args.label_mode
        ]


    # --------------------------------------------------------
    # Dataset-drug pairs
    # --------------------------------------------------------

    pairs = pd.read_csv(
        PAIR_FILE
    )

    dataset_col = pairs.columns[0]
    drug_col = pairs.columns[1]


    all_results = []


    # ========================================================
    # Batch
    # ========================================================

    for scenario in SCENARIOS:

        for round_num in ROUNDS:

            for _, row in pairs.iterrows():

                dataset = str(
                    row[dataset_col]
                ).strip()

                drug_original = str(
                    row[drug_col]
                ).strip()

                drug = convert_drug(
                    drug_original
                )


                for feature in FEATURES:

                    # -----------------------------------------
                    # Get original best parameters ONCE
                    # -----------------------------------------

                    parameter_candidates = find_parameter_files(
                        scenario,
                        round_num,
                        feature,
                        dataset,
                        drug
                    )

                    for parameter_file, params in parameter_candidates:

                        print(
                            "\n=========================================="
                        )

                        print("Original best parameter:")
                        print(parameter_file)
                        print(params)

                        print(
                            "=========================================="
                        )

                        # 同一套参数分别跑两种 label
                        for label_mode in label_modes:
                            result = run_one(

                                device=device,

                                label_mode=label_mode,

                                scenario=scenario,

                                round_num=round_num,

                                dataset=dataset,

                                drug_original=drug_original,

                                feature=feature,

                                params=params,

                                parameter_file=parameter_file
                            )

                            all_results.append(result)

                        # progressive save
                        summary_df = pd.DataFrame(
                            all_results
                        )

                        OUTPUT_ROOT.mkdir(
                            parents=True,
                            exist_ok=True
                        )

                        summary_df.to_csv(
                            OUTPUT_ROOT /
                            "SCAD_fixed_original_params_summary.csv",
                            index=False
                        )


    print(
        "\n=========================================="
    )

    print("ALL FINISHED")

    print(
        OUTPUT_ROOT /
        "SCAD_fixed_original_params_summary.csv"
    )

    print(
        "=========================================="
    )


if __name__ == "__main__":
    main()