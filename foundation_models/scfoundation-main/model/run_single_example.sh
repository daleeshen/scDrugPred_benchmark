#!/usr/bin/env bash

# Run scFoundation on one dataset:
#   1. dataset_making.py
#   2. get_ebd.py
#   3. benchmarking_main_EBD_auto.py

set -uo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "  bash run_scfoundation_example.sh MODEL_ROOT INPUT_FILE"
    exit 1
fi

MODEL_ROOT="$1"
INPUT_FILE="$2"

DATASET_SCRIPT="${MODEL_ROOT}/dataset_making.py"
EMBED_SCRIPT="${MODEL_ROOT}/get_ebd.py"
BENCHMARK_SCRIPT="${MODEL_ROOT}/benchmarking_main_EBD_auto.py"

GENE_INDEX="${MODEL_ROOT}/OS_scRNA_gene_index.19264.tsv"
MODEL_CKPT="${MODEL_ROOT}/models/models.ckpt"

SAMPLES_ROOT="${MODEL_ROOT}/samples"
LABELS_ROOT="${MODEL_ROOT}/labels"
EMBEDS_ROOT="${MODEL_ROOT}/embeds"
RESULTS_ROOT="${MODEL_ROOT}/results"
LOGS_ROOT="${MODEL_ROOT}/logs"

PYTHON_BIN="${PYTHON_BIN:-python}"
N_SPLITS="${N_SPLITS:-3}"


for required_path in \
    "${DATASET_SCRIPT}" \
    "${EMBED_SCRIPT}" \
    "${BENCHMARK_SCRIPT}" \
    "${GENE_INDEX}" \
    "${MODEL_CKPT}" \
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
label_dir="${LABELS_ROOT}/${RUN_NAME}"
embed_dir="${EMBEDS_ROOT}/${RUN_NAME}"
result_dir="${RESULTS_ROOT}/${RUN_NAME}"
log_dir="${LOGS_ROOT}/${RUN_NAME}"

sample_output="${sample_dir}/${stem}_samples.pkl"
label_output="${label_dir}/${stem}_labels.npy"
embed_output="${embed_dir}/${stem}_embeds.npy"
metrics_output="${result_dir}/${stem}.csv"
predictions_output="${result_dir}/${stem}_predictions.csv"
log_file="${log_dir}/${stem}.log"

if ! mkdir -p \
    "${sample_dir}" \
    "${label_dir}" \
    "${embed_dir}" \
    "${result_dir}" \
    "${log_dir}"
then
    echo "ERROR: failed to create output directories"
    exit 1
fi

{
    echo "============================================================"
    echo "scFoundation single-dataset run"
    echo "Dataset: ${INPUT_FILE}"
    echo "Run ID : ${RUN_NAME}"
    echo "Started: $(date '+%F %T')"
    echo "============================================================"

    echo "[1/3] Running dataset_making.py"

    "${PYTHON_BIN}" "${DATASET_SCRIPT}" \
        --data_path "${INPUT_FILE}" \
        --gene_index "${GENE_INDEX}" \
        --model_path "${MODEL_CKPT}" \
        --sample_output "${sample_output}" \
        --label_output "${label_output}"

    if [[ ! -s "${sample_output}" ]]; then
        echo "ERROR: sample output was not created: ${sample_output}"
        exit 1
    fi

    if [[ ! -s "${label_output}" ]]; then
        echo "ERROR: label output was not created: ${label_output}"
        exit 1
    fi

    echo "Sample output: ${sample_output}"
    echo "Label output : ${label_output}"

    echo "[2/3] Running get_ebd.py"

    "${PYTHON_BIN}" "${EMBED_SCRIPT}" \
        --data_path "${INPUT_FILE}" \
        --gene_index "${GENE_INDEX}" \
        --model_path "${MODEL_CKPT}" \
        --output_file "${embed_output}"

    if [[ ! -s "${embed_output}" ]]; then
        echo "ERROR: embedding output was not created: ${embed_output}"
        exit 1
    fi

    echo "Embedding output: ${embed_output}"

    echo "[3/3] Running benchmarking_main_EBD_auto.py"

    "${PYTHON_BIN}" "${BENCHMARK_SCRIPT}" \
        --data_file "${embed_output}" \
        --label_file "${label_output}" \
        --input_file "${INPUT_FILE}" \
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
    echo "Sample output     : ${sample_output}"
    echo "Label output      : ${label_output}"
    echo "Embedding output  : ${embed_output}"
    echo "Metrics output    : ${metrics_output}"
    echo "Predictions output: ${predictions_output}"
    echo "Prediction rows   : ${prediction_cells}"
    echo "============================================================"
} 2>&1 | tee "${log_file}"