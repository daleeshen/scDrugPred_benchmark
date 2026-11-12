#!/bin/bash

# ================= Parameter lists =================
bottleneck_list=(8 16 32 64 128 256 512)
sampling_list=("upsampling" "downsampling" "SMOTE")
round_list=(1 2 3)
imbalance_list=("1_3" "1_30")

# ================= Fixed paths =================
SCDEAL_CODE_DIR="/media/shenli/Expansion/scDrug_benchmark/scDEAL"
PYTHON_EXEC="/home/shenli/py_env/main_env/bin/python"

mapping_csv="/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
balanced_output_root="/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/scDEAL_balanced"
imbalanced_output_root="/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/scDEAL_imbalanced"

bulk_expression_file="${SCDEAL_CODE_DIR}/data/ALL_expression.csv"
bulk_label_file="${SCDEAL_CODE_DIR}/data/ALL_label_binary_wf.csv"

balanced_data_dir_list=(
    "/media/shenli/Expansion/DRMREF/1_round_raw_1000cells/1st_round_data_for_SCAD"
    "/media/shenli/Expansion/DRMREF/2_round_raw_1000cells/data_for_SCAD"
    "/media/shenli/Expansion/DRMREF/3_round_raw_1000cells/data_for_SCAD"
)

cd "$SCDEAL_CODE_DIR" || exit 1

# Arguments: scenario name, data directory, output directory, file prefix
run_scdeal_group() {
    scenario_name="$1"
    data_dir="$2"
    run_output="$3"
    file_prefix="$4"

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

    tail -n +2 "$mapping_csv" | while IFS=',' read -r DATASET DRUG_NAME; do
        DATASET=$(echo "$DATASET" | sed 's/\r//g; s/^"//; s/"$//; s/^[[:space:]]*//; s/[[:space:]]*$//')
        DRUG_NAME=$(echo "$DRUG_NAME" | sed 's/\r//g; s/^"//; s/"$//; s/^[[:space:]]*//; s/[[:space:]]*$//')

        sc_expression_file="${data_dir}/${file_prefix}${DATASET}_raw_count_matrix.csv"
        sc_meta_file="${data_dir}/${file_prefix}${DATASET}_meta_data.csv"

        [ -f "$sc_expression_file" ] && [ -f "$sc_meta_file" ] || continue

        DATA_TAG="$DATASET"
        safe_dataset="${DATASET//\//_}"
        safe_dataset="${safe_dataset// /_}"
        safe_drug="${DRUG_NAME//\//_}"
        safe_drug="${safe_drug// /_}"

        for bottleneck_value in "${bottleneck_list[@]}"; do
            for sampling_method in "${sampling_list[@]}"; do

                LOG_FILE="${log_dir}/${safe_dataset}_${safe_drug}_bottleneck${bottleneck_value}_sampling${sampling_method}.log"

                echo "======================================================================================" | tee -a "$LOG_FILE"
                echo "Scenario: ${scenario_name} | Dataset: ${DATASET} | Drug: ${DRUG_NAME} | Bottleneck: ${bottleneck_value} | Sampling: ${sampling_method}" | tee -a "$LOG_FILE"
                echo "Timestamp: $(date)" | tee -a "$LOG_FILE"
                echo "======================================================================================" | tee -a "$LOG_FILE"

                # ---------------------------------------------------------
                # Step 1: Bulk model
                # ---------------------------------------------------------
                echo ">> [Step 1/2] Running bulkmodel.py..." | tee -a "$LOG_FILE"

                "$PYTHON_EXEC" bulkmodel.py \
                    --root_path "${run_output}/" \
                    --data "$bulk_expression_file" \
                    --label "$bulk_label_file" \
                    --drug "$DRUG_NAME" \
                    --data_name "$DATA_TAG" \
                    --bottleneck "$bottleneck_value" \
                    --sampling "$sampling_method" \
                    --dimreduce "DAE" \
                    --encoder_h_dims "512,256" \
                    --predictor_h_dims "256,128" \
                    --dropout 0.1 \
                    --lr 0.5 \
                    --mod "new" \
                    --printgene "F" \
                    --epochs 500 \
                    --batch_size 200 \
                    --bulk_model "${run_output}/save/bulk_pre/" \
                    --bulk_encoder "${run_output}/save/bulk_encoder/" \
                    --log "${run_output}/save/logs/" >> "$LOG_FILE" 2>&1

                if [ $? -ne 0 ]; then
                    echo "❌ bulkmodel.py failed. Skipping scmodel." | tee -a "$LOG_FILE"
                    continue
                fi

                echo "✅ bulkmodel.py finished successfully." | tee -a "$LOG_FILE"

                # ---------------------------------------------------------
                # Step 2: Single-cell model
                # ---------------------------------------------------------
                echo ">> [Step 2/2] Running scmodel_auto.py..." | tee -a "$LOG_FILE"

                "$PYTHON_EXEC" scmodel_auto.py \
                    --root_path "${run_output}/" \
                    --raw_data_path "$data_dir" \
                    --sc_expression_file "$sc_expression_file" \
                    --meta_file "$sc_meta_file" \
                    --sc_data "$DATA_TAG" \
                    --drug "$DRUG_NAME" \
                    --bottleneck "$bottleneck_value" \
                    --sampling "$sampling_method" \
                    --dimreduce "DAE" \
                    --bulk_h_dims "512,256" \
                    --predictor_h_dims "256,128" \
                    --dropout 0.1 \
                    --lr 0.5 \
                    --mod "new" \
                    --printgene "F" \
                    --checkpoint "False" \
                    --bulk_model_path "${run_output}/save/bulk_pre/" \
                    --sc_model_path "${run_output}/save/sc_pre/" \
                    --sc_encoder_path "${run_output}/save/sc_encoder/" \
                    --logging_file "${run_output}/save/logs/" >> "$LOG_FILE" 2>&1

                if [ $? -eq 0 ]; then
                    echo "✅ scmodel_auto.py finished successfully." | tee -a "$LOG_FILE"

                    default_result="${run_output}/shuyu_data_Ola_DAE_bottle${bottleneck_value}_modnew_${sampling_method}_prediction_results.csv"
                    if [ -f "$default_result" ]; then
                        mv "$default_result" \
                           "${prediction_dir}/${safe_dataset}_${safe_drug}_DAE_bottle${bottleneck_value}_modnew_${sampling_method}_prediction_results.csv"
                    fi
                else
                    echo "❌ scmodel_auto.py failed." | tee -a "$LOG_FILE"
                fi

                echo "" | tee -a "$LOG_FILE"
            done
        done
    done
}

# ======================================================================
# Part 1: Run imbalanced datasets first
#         (1_3 and 1_30; rounds 1, 2, and 3 for each ratio)
# ======================================================================
for imbalance_ratio in "${imbalance_list[@]}"; do
    for round_id in "${round_list[@]}"; do
        imbalanced_data_dir="/media/shenli/Expansion/DRMREF/unbalanced_original/${round_id}_round_unbalanced_sampling/${imbalance_ratio}/metadata_count"
        imbalanced_output="${imbalanced_output_root}/${imbalance_ratio}/${round_id}"

        run_scdeal_group \
            "imbalanced_${imbalance_ratio}_round${round_id}" \
            "$imbalanced_data_dir" \
            "$imbalanced_output" \
            ""
    done
done

# ======================================================================
# Part 2: Run balanced datasets
#         (rounds 1, 2, and 3)
# ======================================================================
for i in "${!round_list[@]}"; do
    round_id="${round_list[$i]}"
    balanced_data_dir="${balanced_data_dir_list[$i]}"
    balanced_output="${balanced_output_root}/${round_id}"

    run_scdeal_group \
        "balanced_round${round_id}" \
        "$balanced_data_dir" \
        "$balanced_output" \
        "balanced_sampled_"
done
