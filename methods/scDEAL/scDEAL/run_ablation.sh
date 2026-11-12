#!/bin/bash

# ============================================================
# scDEAL ablation experiments for Figure 5E
#
# Full model:
#   integrate + DAE + mod=new
#
# Re-run four ablations:
#   1. no_cluster = integrate + DAE + mod=ori
#   2. no_dae     = integrate + AE  + mod=new
#   3. gdsc_only  = old       + DAE + mod=new
#   4. ccle_only  = new       + DAE + mod=new
#
# No-transfer:
#   DO NOT re-run.
#   Use sens_preds_pret / sens_label_pret from Full model output.
# ============================================================


# ================= 实验场景 =================
round_list=(1 2 3)
imbalance_list=("1_3" "1_30")


# ================= 固定路径 =================
SCDEAL_CODE_DIR="/media/shenli/Expansion/scDrug_benchmark/scDEAL"
PYTHON_EXEC="/home/shenli/py_env/main_env/bin/python"

# Figure 5中使用的9个dataset-drug pairs
mapping_csv="/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"

# ------------------------------------------------------------
# Full scDEAL已经筛选好的best-performance结果
# ------------------------------------------------------------

# Balanced, 1000 cells
# columns:
# dataset drug bottleneck sampling_method AUROC ACC F1 round
balanced_best_csv="/media/shenli/Expansion/DRMREF/1_round_1000cells_res/1000cells_scDEAL.csv"

# Imbalanced
# columns:
# round ratio bottleneck sampling_method dataset AUROC AUPRC
# NormAUPRC ACC Recall Precision F1 drug
imbalanced_best_csv="/media/shenli/Expansion/DRMREF/unb_scDEAL_best_performance.csv"


# ------------------------------------------------------------
# Figure 5E输出目录
# ------------------------------------------------------------
FIGURE5E_ROOT="/media/shenli/Expansion/DRMREF/New_fig5/Figure5E"
ABLAT_OUTPUT_ROOT="${FIGURE5E_ROOT}/scDEAL_ablation"

mkdir -p "$FIGURE5E_ROOT"
mkdir -p "$ABLAT_OUTPUT_ROOT"


# ------------------------------------------------------------
# scDEAL bulk data
# ------------------------------------------------------------
bulk_expression_file="${SCDEAL_CODE_DIR}/data/ALL_expression.csv"
bulk_label_file="${SCDEAL_CODE_DIR}/data/ALL_label_binary_wf.csv"


# ------------------------------------------------------------
# Balanced 1000-cell数据
# ------------------------------------------------------------
balanced_data_dir_list=(
    "/media/shenli/Expansion/DRMREF/1_round_raw_1000cells/data_for_SCAD"
    "/media/shenli/Expansion/DRMREF/2_round_raw_1000cells/data_for_SCAD"
    "/media/shenli/Expansion/DRMREF/3_round_raw_1000cells/data_for_SCAD"
)


cd "$SCDEAL_CODE_DIR" || exit 1


# ============================================================
# Step 0
# 从两个best-performance CSV生成统一的parameter table
#
# 最终格式：
# setting    round    dataset    drug    bottleneck    sampling_method
#
# setting:
# balanced
# 1_3
# 1_30
# ============================================================

PARAM_TABLE="${FIGURE5E_ROOT}/scDEAL_ablation_best_params.tsv"


echo "============================================================"
echo "Preparing Full-model best parameter table..."
echo "============================================================"


"$PYTHON_EXEC" - \
    "$mapping_csv" \
    "$balanced_best_csv" \
    "$imbalanced_best_csv" \
    "$PARAM_TABLE" <<'PY'

import sys
import re
import pandas as pd


mapping_file = sys.argv[1]
balanced_file = sys.argv[2]
imbalanced_file = sys.argv[3]
output_file = sys.argv[4]


def read_table(path):
    """
    Automatically detect comma/tab separator.
    """
    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        encoding="utf-8-sig"
    )


def clean_string(x):
    return str(x).strip().strip('"').strip("'")


def make_key(x):
    return clean_string(x).lower()


def normalize_ratio(x):
    x = clean_string(x)
    x = x.replace(":", "_")
    x = x.replace("/", "_")
    x = re.sub(r"\s+", "", x)
    return x


# ============================================================
# 1. Read dataset-drug mapping
# ============================================================

mapping_raw = read_table(mapping_file)

if mapping_raw.shape[1] < 2:
    raise ValueError(
        "mapping_csv must contain at least two columns: dataset and drug."
    )

mapping = mapping_raw.iloc[:, :2].copy()
mapping.columns = ["dataset", "drug"]

mapping["dataset"] = mapping["dataset"].map(clean_string)
mapping["drug"] = mapping["drug"].map(clean_string)

mapping["_dataset_key"] = mapping["dataset"].map(make_key)
mapping["_drug_key"] = mapping["drug"].map(make_key)

mapping = mapping.drop_duplicates(
    subset=["_dataset_key", "_drug_key"]
).reset_index(drop=True)

print(f"Number of dataset-drug pairs in mapping: {len(mapping)}")


# ============================================================
# 2. Balanced best parameters
# ============================================================

balanced = read_table(balanced_file)
balanced.columns = [str(x).strip() for x in balanced.columns]

required_balanced = {
    "dataset",
    "drug",
    "bottleneck",
    "sampling_method",
    "round"
}

missing = required_balanced - set(balanced.columns)

if missing:
    raise ValueError(
        "Balanced best-performance CSV missing columns: "
        + ", ".join(sorted(missing))
    )


balanced["dataset"] = balanced["dataset"].map(clean_string)
balanced["drug"] = balanced["drug"].map(clean_string)

balanced["_dataset_key"] = balanced["dataset"].map(make_key)
balanced["_drug_key"] = balanced["drug"].map(make_key)

balanced["round"] = pd.to_numeric(
    balanced["round"],
    errors="raise"
).astype(int)

balanced["bottleneck"] = pd.to_numeric(
    balanced["bottleneck"],
    errors="raise"
).astype(int)

balanced["sampling_method"] = (
    balanced["sampling_method"]
    .astype(str)
    .str.strip()
)


balanced_merge = mapping.merge(
    balanced,
    on=["_dataset_key", "_drug_key"],
    how="inner",
    suffixes=("_map", "_best")
)


balanced_out = pd.DataFrame({
    "setting": "balanced",
    "round": balanced_merge["round"],
    "dataset": balanced_merge["dataset_map"],
    "drug": balanced_merge["drug_map"],
    "bottleneck": balanced_merge["bottleneck"],
    "sampling_method": balanced_merge["sampling_method"]
})


balanced_out = balanced_out[
    balanced_out["round"].isin([1, 2, 3])
].copy()


# ============================================================
# 3. Imbalanced best parameters
# ============================================================

imbalanced = read_table(imbalanced_file)
imbalanced.columns = [str(x).strip() for x in imbalanced.columns]

required_imbalanced = {
    "round",
    "ratio",
    "bottleneck",
    "sampling_method",
    "dataset",
    "drug"
}

missing = required_imbalanced - set(imbalanced.columns)

if missing:
    raise ValueError(
        "Imbalanced best-performance CSV missing columns: "
        + ", ".join(sorted(missing))
    )


imbalanced["dataset"] = imbalanced["dataset"].map(clean_string)
imbalanced["drug"] = imbalanced["drug"].map(clean_string)

imbalanced["_dataset_key"] = imbalanced["dataset"].map(make_key)
imbalanced["_drug_key"] = imbalanced["drug"].map(make_key)

imbalanced["round"] = pd.to_numeric(
    imbalanced["round"],
    errors="raise"
).astype(int)

imbalanced["bottleneck"] = pd.to_numeric(
    imbalanced["bottleneck"],
    errors="raise"
).astype(int)

imbalanced["sampling_method"] = (
    imbalanced["sampling_method"]
    .astype(str)
    .str.strip()
)

imbalanced["setting"] = imbalanced["ratio"].map(normalize_ratio)


# Figure 5只使用1:3和1:30
imbalanced = imbalanced[
    imbalanced["setting"].isin(["1_3", "1_30"])
].copy()


imbalanced_merge = mapping.merge(
    imbalanced,
    on=["_dataset_key", "_drug_key"],
    how="inner",
    suffixes=("_map", "_best")
)


imbalanced_out = pd.DataFrame({
    "setting": imbalanced_merge["setting"],
    "round": imbalanced_merge["round"],
    "dataset": imbalanced_merge["dataset_map"],
    "drug": imbalanced_merge["drug_map"],
    "bottleneck": imbalanced_merge["bottleneck"],
    "sampling_method": imbalanced_merge["sampling_method"]
})


imbalanced_out = imbalanced_out[
    imbalanced_out["round"].isin([1, 2, 3])
].copy()


# ============================================================
# 4. Combine
# ============================================================

result = pd.concat(
    [balanced_out, imbalanced_out],
    ignore_index=True
)


# Check duplicate best parameter records
duplicate_cols = [
    "setting",
    "round",
    "dataset",
    "drug"
]

duplicates = result[
    result.duplicated(
        duplicate_cols,
        keep=False
    )
]

if len(duplicates) > 0:
    print("\nERROR: duplicated best parameter records found:")
    print(duplicates.to_string(index=False))
    raise SystemExit(1)


# ============================================================
# 5. Check whether every expected combination exists
# ============================================================

expected = []

for setting in ["balanced", "1_3", "1_30"]:
    for round_id in [1, 2, 3]:
        for _, row in mapping.iterrows():
            expected.append(
                (
                    setting,
                    round_id,
                    row["dataset"],
                    row["drug"]
                )
            )


actual = set(
    zip(
        result["setting"],
        result["round"],
        result["dataset"],
        result["drug"]
    )
)

missing_records = [
    x for x in expected
    if x not in actual
]


if missing_records:

    print("\nWARNING: some best-parameter records are missing and will be skipped:")
    for x in missing_records:
        print(x)


# ============================================================
# 6. Sort and save
# ============================================================

setting_order = {
    "balanced": 0,
    "1_3": 1,
    "1_30": 2
}

result["_order"] = result["setting"].map(setting_order)

result = result.sort_values(
    [
        "_order",
        "round",
        "dataset",
        "drug"
    ]
)

result = result.drop(columns="_order")


result.to_csv(
    output_file,
    sep="\t",
    index=False
)


print("\nParameter table successfully created:")
print(output_file)

print("\nNumber of runs:")
print(
    result.groupby(["setting", "round"])
    .size()
    .to_string()
)

print("\nPreview:")
print(result.head(10).to_string(index=False))

PY


if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Failed to prepare best parameter table."
    exit 1
fi


echo ""
echo "Best parameter table successfully prepared:"
echo "$PARAM_TABLE"
echo ""


# ============================================================
# Function:
# run one ablation × one scenario × one round
# ============================================================

run_scdeal_group() {

    setting="$1"
    scenario_name="$2"
    round_id="$3"
    data_dir="$4"
    run_output="$5"
    file_prefix="$6"

    ablation_name="$7"
    bulk_mode="$8"
    dimreduce="$9"
    mod="${10}"


    log_dir="${run_output}/logs"
    prediction_dir="${run_output}/prediction_results"


    mkdir -p \
        "$log_dir" \
        "$prediction_dir" \
        "${run_output}/save/bulk_pre" \
        "${run_output}/save/bulk_encoder" \
        "${run_output}/save/sc_pre" \
        "${run_output}/save/sc_encoder" \
        "${run_output}/save/logs"


    found_any=0


    # ========================================================
    # Read Full-model best parameters
    # ========================================================

    while IFS=$'\t' read -r \
        setting_csv \
        round_csv \
        DATASET \
        DRUG_NAME \
        bottleneck_value \
        sampling_method

    do

        [ -z "$DATASET" ] && continue

        found_any=1


        # ----------------------------------------------------
        # Input single-cell files
        # ----------------------------------------------------

        sc_expression_file="${data_dir}/${file_prefix}${DATASET}_raw_count_matrix.csv"
        sc_meta_file="${data_dir}/${file_prefix}${DATASET}_meta_data.csv"


        if [ ! -f "$sc_expression_file" ]; then

            echo "WARNING: expression file missing:"
            echo "$sc_expression_file"
            continue

        fi


        if [ ! -f "$sc_meta_file" ]; then

            echo "WARNING: metadata file missing:"
            echo "$sc_meta_file"
            continue

        fi


        DATA_TAG="$DATASET"


        safe_dataset="${DATASET//\//_}"
        safe_dataset="${safe_dataset// /_}"

        safe_drug="${DRUG_NAME//\//_}"
        safe_drug="${safe_drug// /_}"


        # ----------------------------------------------------
        # Output names
        # ----------------------------------------------------

        LOG_FILE="${log_dir}/${safe_dataset}_${safe_drug}_${ablation_name}_bottle${bottleneck_value}_${sampling_method}.log"


        final_result="${prediction_dir}/${safe_dataset}_${safe_drug}_${ablation_name}_${dimreduce}_bottle${bottleneck_value}_mod${mod}_bulk${bulk_mode}_${sampling_method}_prediction_results.csv"


        # ----------------------------------------------------
        # Allow the script to resume without re-running
        # completed jobs
        # ----------------------------------------------------

        if [ -f "$final_result" ]; then

            echo "SKIP: result already exists:"
            echo "$final_result"
            continue

        fi


        echo "======================================================================================" | tee -a "$LOG_FILE"

        echo "Ablation: ${ablation_name}" | tee -a "$LOG_FILE"
        echo "Scenario: ${scenario_name}" | tee -a "$LOG_FILE"
        echo "Round: ${round_id}" | tee -a "$LOG_FILE"

        echo "Dataset: ${DATASET}" | tee -a "$LOG_FILE"
        echo "Drug: ${DRUG_NAME}" | tee -a "$LOG_FILE"

        echo "Bottleneck: ${bottleneck_value}" | tee -a "$LOG_FILE"
        echo "Sampling: ${sampling_method}" | tee -a "$LOG_FILE"

        echo "Bulk mode: ${bulk_mode}" | tee -a "$LOG_FILE"
        echo "Dimreduce: ${dimreduce}" | tee -a "$LOG_FILE"
        echo "Mod: ${mod}" | tee -a "$LOG_FILE"

        echo "Timestamp: $(date)" | tee -a "$LOG_FILE"

        echo "======================================================================================" | tee -a "$LOG_FILE"



        # ====================================================
        # Step 1
        # Train bulk model
        # ====================================================

        echo ">> [Step 1/2] Running bulkmodel.py..." \
            | tee -a "$LOG_FILE"


        "$PYTHON_EXEC" bulkmodel.py \
            --root_path "${run_output}/" \
            --data "$bulk_expression_file" \
            --label "$bulk_label_file" \
            --drug "$DRUG_NAME" \
            --data_name "$DATA_TAG" \
            --bottleneck "$bottleneck_value" \
            --sampling "$sampling_method" \
            --dimreduce "$dimreduce" \
            --encoder_h_dims "512,256" \
            --predictor_h_dims "256,128" \
            --dropout 0.1 \
            --lr 0.5 \
            --mod "$mod" \
            --bulk "$bulk_mode" \
            --printgene "F" \
            --epochs 500 \
            --batch_size 200 \
            --bulk_model "${run_output}/save/bulk_pre/" \
            --bulk_encoder "${run_output}/save/bulk_encoder/" \
            --log "${run_output}/save/logs/" \
            >> "$LOG_FILE" 2>&1


        if [ $? -ne 0 ]; then

            echo "ERROR: bulkmodel.py failed. Skipping scmodel." \
                | tee -a "$LOG_FILE"

            continue

        fi


        echo "OK: bulkmodel.py finished successfully." \
            | tee -a "$LOG_FILE"



        # ====================================================
        # Step 2
        # Transfer + prediction
        # ====================================================

        echo ">> [Step 2/2] Running scmodel_auto.py..." \
            | tee -a "$LOG_FILE"


        "$PYTHON_EXEC" scmodel_auto.py \
            --root_path "${run_output}/" \
            --raw_data_path "$data_dir" \
            --sc_expression_file "$sc_expression_file" \
            --meta_file "$sc_meta_file" \
            --sc_data "$DATA_TAG" \
            --drug "$DRUG_NAME" \
            --bottleneck "$bottleneck_value" \
            --sampling "$sampling_method" \
            --dimreduce "$dimreduce" \
            --bulk_h_dims "512,256" \
            --predictor_h_dims "256,128" \
            --dropout 0.1 \
            --lr 0.5 \
            --mod "$mod" \
            --bulk "$bulk_mode" \
            --bulk_data "$bulk_expression_file" \
            --label "$bulk_label_file" \
            --printgene "F" \
            --checkpoint "False" \
            --bulk_model_path "${run_output}/save/bulk_pre/" \
            --sc_model_path "${run_output}/save/sc_pre/" \
            --sc_encoder_path "${run_output}/save/sc_encoder/" \
            --logging_file "${run_output}/save/logs/" \
            >> "$LOG_FILE" 2>&1



        # ====================================================
        # Save prediction result
        # ====================================================

        if [ $? -eq 0 ]; then

            echo "OK: scmodel_auto.py finished successfully." \
                | tee -a "$LOG_FILE"


            default_result="${run_output}/shuyu_data_Ola_${dimreduce}_bottle${bottleneck_value}_mod${mod}_${sampling_method}_prediction_results.csv"


            if [ -f "$default_result" ]; then

                mv "$default_result" "$final_result"

                echo "Result saved:" \
                    | tee -a "$LOG_FILE"

                echo "$final_result" \
                    | tee -a "$LOG_FILE"

            else

                echo "WARNING: scmodel completed, but prediction result was not found:" \
                    | tee -a "$LOG_FILE"

                echo "$default_result" \
                    | tee -a "$LOG_FILE"

            fi


        else

            echo "ERROR: scmodel_auto.py failed." \
                | tee -a "$LOG_FILE"

        fi


        echo "" | tee -a "$LOG_FILE"


    done < <(
        awk -F '\t' \
            -v target_setting="$setting" \
            -v target_round="$round_id" \
            'NR > 1 && $1 == target_setting && $2 == target_round {print $0}' \
            "$PARAM_TABLE"
    )


    if [ "$found_any" -eq 0 ]; then

        echo "WARNING: No parameter records found:"
        echo "setting=${setting}, round=${round_id}"

    fi
}



# ============================================================
# Four ablation configurations
#
# format:
# name | bulk | dimreduce | mod
# ============================================================

ablation_configs=(

    "no_cluster|integrate|DAE|ori"

    "no_dae|integrate|AE|new"

    "gdsc_only|old|DAE|new"

    "ccle_only|new|DAE|new"

)



# ============================================================
# Run all ablations
# ============================================================

for config in "${ablation_configs[@]}"
do


    IFS='|' read -r \
        ablation_name \
        bulk_mode \
        dimreduce \
        mod \
        <<< "$config"


    echo ""
    echo "################################################################################"
    echo "Starting ablation:"
    echo ""
    echo "  ${ablation_name}"
    echo ""
    echo "  bulk      = ${bulk_mode}"
    echo "  dimreduce = ${dimreduce}"
    echo "  mod       = ${mod}"
    echo ""
    echo "################################################################################"
    echo ""



    # ========================================================
    # Part 1
    # Imbalanced:
    #
    # 1:3
    # 1:30
    #
    # rounds 1-3
    # ========================================================

    for imbalance_ratio in "${imbalance_list[@]}"
    do


        for round_id in "${round_list[@]}"
        do


            imbalanced_data_dir="/media/shenli/Expansion/DRMREF/unbalanced_original/${round_id}_round_unbalanced_sampling/${imbalance_ratio}/metadata_count"


            imbalanced_output="${ABLAT_OUTPUT_ROOT}/${ablation_name}/imbalanced/${imbalance_ratio}/round${round_id}"


            run_scdeal_group \
                "$imbalance_ratio" \
                "imbalanced_${imbalance_ratio}_round${round_id}" \
                "$round_id" \
                "$imbalanced_data_dir" \
                "$imbalanced_output" \
                "" \
                "$ablation_name" \
                "$bulk_mode" \
                "$dimreduce" \
                "$mod"


        done

    done



    # ========================================================
    # Part 2
    # Balanced:
    #
    # 1000 cells
    # rounds 1-3
    # ========================================================

    for i in "${!round_list[@]}"
    do


        round_id="${round_list[$i]}"

        balanced_data_dir="${balanced_data_dir_list[$i]}"


        balanced_output="${ABLAT_OUTPUT_ROOT}/${ablation_name}/balanced/round${round_id}"


        run_scdeal_group \
            "balanced" \
            "balanced_1000cells_round${round_id}" \
            "$round_id" \
            "$balanced_data_dir" \
            "$balanced_output" \
            "balanced_sampled_" \
            "$ablation_name" \
            "$bulk_mode" \
            "$dimreduce" \
            "$mod"


    done


done


echo ""
echo "################################################################################"
echo "All scDEAL Figure 5E ablation experiments finished."
echo ""
echo "Results:"
echo "${ABLAT_OUTPUT_ROOT}"
echo "################################################################################"