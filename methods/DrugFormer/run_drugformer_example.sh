#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON_EXEC="${PYTHON_EXEC:-python}"

CODE_DIR="${SCRIPT_DIR}/DrugFormer"
DATASET_SCRIPT="${CODE_DIR}/Dataset_auto.py"
MAIN_SCRIPT="${CODE_DIR}/main_auto.py"
DICT_FILE="${CODE_DIR}/token_gene_dictionary.pickle"

DATA_DIR="${ROOT_DIR}/data/example"
CELL_DT_DIR="${SCRIPT_DIR}/cell_dt"
OUT_DIR="${SCRIPT_DIR}/example_output"
LOG_DIR="${SCRIPT_DIR}/example_logs"

SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

mkdir -p "${CELL_DT_DIR}" "${OUT_DIR}" "${LOG_DIR}"
shopt -s nullglob

for scenario in "${SCENARIOS[@]}"; do
    input_dir="${DATA_DIR}/${scenario}/drugformer_input"
    dataset_root="${CELL_DT_DIR}/${scenario}"
    output_root="${OUT_DIR}/${scenario}"
    log_root="${LOG_DIR}/${scenario}"

    if [ "${scenario}" = "balanced_1000" ]; then
        EP_NUM=3
        TRAIN_BATCH=12
        TEST_BATCH=24
    else
        EP_NUM=5
        TRAIN_BATCH=8
        TEST_BATCH=8
    fi

    mkdir -p "${dataset_root}" "${output_root}" "${log_root}"

    echo "============================================================"
    echo "Running DrugFormer: ${scenario}"
    echo "Epochs: ${EP_NUM} | Train batch: ${TRAIN_BATCH} | Test batch: ${TEST_BATCH}"
    echo "============================================================"

    input_files=("${input_dir}"/*_DrugFormer.txt)

    if [ ${#input_files[@]} -eq 0 ]; then
        echo "No DrugFormer input files found in ${input_dir}"
        continue
    fi

    for input_file in "${input_files[@]}"; do
        filename="$(basename "${input_file}")"
        project="${filename%_DrugFormer.txt}"

        dataset_dir="${dataset_root}/${project}"
        save_dir="${output_root}/${project}"
        project_log_dir="${log_root}/${project}"

        rm -rf "${dataset_dir}" "${save_dir}" "${project_log_dir}"
        mkdir -p "${save_dir}" "${project_log_dir}"

        echo "------------------------------------------------------------"
        echo "Project: ${project}"
        echo "Input: ${filename}"
        echo "------------------------------------------------------------"

        echo ">> Step 1/2: Preparing dataset"
        "${PYTHON_EXEC}" "${DATASET_SCRIPT}" \
            --input_file "${input_file}" \
            --dict_file "${DICT_FILE}" \
            --output_dir "${dataset_dir}" \
            --pos_keyword "sensitive" \
            --neg_keyword "resistant"

        echo ">> Step 2/2: Running DrugFormer"
        "${PYTHON_EXEC}" "${MAIN_SCRIPT}" \
            --project_name "${project}" \
            --data_path "${dataset_root}" \
            --save_path "${save_dir}" \
            --log_folder "${project_log_dir}" \
            --ep_num "${EP_NUM}" \
            --train_batch_size "${TRAIN_BATCH}" \
            --test_batch_size "${TEST_BATCH}"

        echo "${project} completed."
    done
done

echo "============================================================"
echo "DrugFormer example completed."
echo "Logs: ${LOG_DIR}"
echo "Models: ${OUT_DIR}"
echo "============================================================"
