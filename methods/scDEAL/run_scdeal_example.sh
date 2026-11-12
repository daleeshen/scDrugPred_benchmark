#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CODE_DIR="${SCRIPT_DIR}/scDEAL"
PYTHON_EXEC="${PYTHON_EXEC:-python}"

DATA_DIR="${ROOT_DIR}/data/example"
BULK_EXPR="${CODE_DIR}/data/ALL_expression.csv"
BULK_LABEL="${CODE_DIR}/data/ALL_label_binary_wf.csv"

OUT_DIR="${SCRIPT_DIR}/example_output"
MODEL_DIR="${SCRIPT_DIR}/example_models"
LOG_DIR="${SCRIPT_DIR}/example_logs"

DRUG="PLX-4720"
DATA_TAG="GSE108397"

BOTTLENECKS=(8 16 32 64 128 256 512)
SAMPLINGS=("upsampling" "downsampling" "SMOTE")
SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

mkdir -p "${OUT_DIR}" "${MODEL_DIR}" "${LOG_DIR}"

cd "${CODE_DIR}"

for bottleneck in "${BOTTLENECKS[@]}"; do
    for sampling in "${SAMPLINGS[@]}"; do

        combo="bottle${bottleneck}_${sampling}"
        bulk_root="${MODEL_DIR}/bulk/${combo}"
        bulk_log_dir="${LOG_DIR}/bulk"
        bulk_log="${bulk_log_dir}/${combo}.log"

        mkdir -p \
            "${bulk_root}/save/bulk_pre" \
            "${bulk_root}/save/bulk_encoder" \
            "${bulk_root}/save/logs" \
            "${bulk_log_dir}"

        echo "================================================================================"
        echo "Bulk model | Bottleneck: ${bottleneck} | Sampling: ${sampling}"
        echo "================================================================================"

        "${PYTHON_EXEC}" bulkmodel.py \
            --root_path "${bulk_root}/" \
            --data "${BULK_EXPR}" \
            --label "${BULK_LABEL}" \
            --drug "${DRUG}" \
            --data_name "${DATA_TAG}" \
            --bottleneck "${bottleneck}" \
            --sampling "${sampling}" \
            --dimreduce "DAE" \
            --encoder_h_dims "512,256" \
            --predictor_h_dims "256,128" \
            --dropout 0.1 \
            --lr 0.5 \
            --mod "new" \
            --printgene "F" \
            --epochs 500 \
            --batch_size 200 \
            --checkpoint "False" \
            --bulk_model "${bulk_root}/save/bulk_pre/" \
            --bulk_encoder "${bulk_root}/save/bulk_encoder/" \
            --log "${bulk_root}/save/logs/" \
            > "${bulk_log}" 2>&1

        echo "Bulk model completed."

        for scenario in "${SCENARIOS[@]}"; do
            input_dir="${DATA_DIR}/${scenario}/scdeal_input"
            scenario_out="${OUT_DIR}/${scenario}/${combo}"
            scenario_log="${LOG_DIR}/${scenario}/${combo}"

            mkdir -p "${scenario_out}" "${scenario_log}"

            shopt -s nullglob
            expression_files=("${input_dir}"/*_raw_count_matrix.csv)

            if [ ${#expression_files[@]} -eq 0 ]; then
                echo "No scDEAL input files found in ${input_dir}"
                continue
            fi

            echo "--------------------------------------------------------------------------------"
            echo "Scenario: ${scenario} | Bottleneck: ${bottleneck} | Sampling: ${sampling}"
            echo "--------------------------------------------------------------------------------"

            for sc_expression_file in "${expression_files[@]}"; do
                filename="$(basename "${sc_expression_file}")"
                sample_id="${filename%_raw_count_matrix.csv}"
                sc_meta_file="${input_dir}/${sample_id}_meta_data.csv"

                if [ ! -f "${sc_meta_file}" ]; then
                    echo "Missing metadata file: ${sc_meta_file}"
                    exit 1
                fi

                sample_root="${MODEL_DIR}/sc/${scenario}/${combo}/${sample_id}"
                log_file="${scenario_log}/${sample_id}.log"

                mkdir -p \
                    "${sample_root}/save/sc_pre" \
                    "${sample_root}/save/sc_encoder" \
                    "${sample_root}/save/logs"

                echo "Running: ${sample_id}"

                "${PYTHON_EXEC}" scmodel_auto.py \
                    --root_path "${sample_root}/" \
                    --raw_data_path "${input_dir}" \
                    --sc_expression_file "${sc_expression_file}" \
                    --meta_file "${sc_meta_file}" \
                    --sc_data "${DATA_TAG}" \
                    --drug "${DRUG}" \
                    --bottleneck "${bottleneck}" \
                    --sampling "${sampling}" \
                    --dimreduce "DAE" \
                    --bulk_h_dims "512,256" \
                    --predictor_h_dims "256,128" \
                    --dropout 0.1 \
                    --lr 0.5 \
                    --mod "new" \
                    --printgene "F" \
                    --checkpoint "False" \
                    --bulk_data "${BULK_EXPR}" \
                    --label "${BULK_LABEL}" \
                    --bulk_model_path "${bulk_root}/save/bulk_pre/" \
                    --sc_model_path "${sample_root}/save/sc_pre/" \
                    --sc_encoder_path "${sample_root}/save/sc_encoder/" \
                    --logging_file "${sample_root}/save/logs/" \
                    > "${log_file}" 2>&1

                result_file="${sample_root}/${sample_id}_DAE_bottle${bottleneck}_modnew_${sampling}_prediction_results.csv"

                if [ -f "${result_file}" ]; then
                    mv "${result_file}" "${scenario_out}/"
                    echo "Prediction saved: ${scenario_out}/$(basename "${result_file}")"
                else
                    echo "Prediction result not found: ${result_file}"
                    exit 1
                fi
            done
        done
    done
done

echo "================================================================================"
echo "scDEAL example completed."
echo "Predictions: ${OUT_DIR}"
echo "Logs: ${LOG_DIR}"
echo "Models: ${MODEL_DIR}"
echo "================================================================================"
