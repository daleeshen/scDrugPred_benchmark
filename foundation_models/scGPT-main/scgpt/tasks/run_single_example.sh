#!/usr/bin/env bash

# Run scGPT on one dataset:
#   1. stp1_get_ebd_auto.py
#   2. stp2_benchmarking_main_EBD_auto.py

set -uo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "  bash run_single_example.sh MODEL_ROOT INPUT_FILE"
    exit 1
fi

MODEL_ROOT="$1"
INPUT_FILE="$2"

EMBED_SCRIPT="${MODEL_ROOT}/tasks/stp1_get_ebd_auto.py"
BENCHMARK_SCRIPT="${MODEL_ROOT}/tasks/stp2_benchmarking_main_EBD_auto.py"

MODEL_DIR="${MODEL_ROOT}/save"

SAMPLES_ROOT="${MODEL_ROOT}/samples"
RESULTS_ROOT="${MODEL_ROOT}/results"
LOGS_ROOT="${MODEL_ROOT}/logs"

PYTHON_BIN="${PYTHON_BIN:-python}"
N_SPLITS="${N_SPLITS:-3}"
BATCH_SIZE="${BATCH_SIZE:-64}"
MAX_LENGTH="${MAX_LENGTH:-1200}"

for required_path in \
    "${EMBED_SCRIPT}" \
    "${BENCHMARK_SCRIPT}" \
    "${MODEL_DIR}" \
    "${INPUT_FILE}"
do
    if [[ ! -e "${required_path}" ]]; then
        echo "Required path not found: ${required_path}"
        exit 1
    fi
done

if ! cd "${MODEL_ROOT}"; then
    echo "Cannot enter model directory: ${MODEL_ROOT}"
    exit 1
fi

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_NAME="single_run_${RUN_ID}"
stem="$(basename "${INPUT_FILE}" .txt)"

sample_dir="${SAMPLES_ROOT}/${RUN_NAME}"
result_dir="${RESULTS_ROOT}/${RUN_NAME}"
log_dir="${LOGS_ROOT}/${RUN_NAME}"

embed_output="${sample_dir}/${stem}_samples.pkl"
metrics_output="${result_dir}/${stem}.csv"
predictions_output="${result_dir}/${stem}_predictions.csv"
log_file="${log_dir}/${stem}.log"

input_root="$(dirname "${INPUT_FILE}")"

if ! mkdir -p \
    "${sample_dir}" \
    "${result_dir}" \
    "${log_dir}"
then
    echo "ERROR: failed to create output directories"
    exit 1
fi

{
    echo "============================================================"
    echo "scGPT single-dataset run"
    echo "Dataset: ${INPUT_FILE}"
    echo "Run ID : ${RUN_NAME}"
    echo "Started: $(date '+%F %T')"
    echo "============================================================"

    echo "[1/2] Running stp1_get_ebd_auto.py"

    "${PYTHON_BIN}" "${EMBED_SCRIPT}" \
        --input_file "${INPUT_FILE}" \
        --input_root "${input_root}" \
        --output_root "${sample_dir}" \
        --model_dir "${MODEL_DIR}" \
        --batch_size "${BATCH_SIZE}" \
        --max_length "${MAX_LENGTH}"

    if [[ ! -s "${embed_output}" ]]; then
        echo "ERROR: embedding output was not created: ${embed_output}"
        exit 1
    fi

    echo "Embedding output: ${embed_output}"

    echo "[2/2] Running stp2_benchmarking_main_EBD_auto.py"

    "${PYTHON_BIN}" "${BENCHMARK_SCRIPT}" \
        --data_file "${embed_output}" \
        --n_splits "${N_SPLITS}" \
        --output_file "${metrics_output}"

    if [[ ! -s "${metrics_output}" ]]; then
        echo "ERROR: metrics file was not created: ${metrics_output}"
        exit 1
    fi

    if [[ ! -s "${predictions_output}" ]]; then
        echo "ERROR: predictions file was not created: ${predictions_output}"
        exit 1
    fi

    expected_cells="$(awk 'END {print NR - 1}' "${INPUT_FILE}")"
    prediction_cells="$(awk 'END {print NR - 1}' "${predictions_output}")"

    if [[ "${prediction_cells}" -ne "${expected_cells}" ]]; then
        echo "ERROR: prediction row count does not match input cell count"
        echo "Input cells     : ${expected_cells}"
        echo "Prediction rows : ${prediction_cells}"
        exit 1
    fi

    echo "============================================================"
    echo "Completed: $(date '+%F %T')"
    echo "Embedding output  : ${embed_output}"
    echo "Metrics output    : ${metrics_output}"
    echo "Predictions output: ${predictions_output}"
    echo "Prediction rows   : ${prediction_cells}"
    echo "============================================================"
} 2>&1 | tee "${log_file}"

