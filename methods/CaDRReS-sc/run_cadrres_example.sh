#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CADRRES_DIR="${SCRIPT_DIR}/CaDRReS-sc"
PYTHON_EXEC="${PYTHON_EXEC:-python}"
RUN_SCRIPT="${CADRRES_DIR}/auto_script.py"
MODEL_DIR="${CADRRES_DIR}/pretrained_models"
GDSC_EXP="${CADRRES_DIR}/data/GDSC/GDSC_exp.tsv"
ESSENTIAL_GENES="${CADRRES_DIR}/data/essential_genes.txt"

DATA_DIR="${ROOT_DIR}/data/example"
OUT_DIR="${SCRIPT_DIR}/example_output"
SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

mkdir -p "${OUT_DIR}"
shopt -s nullglob

for scenario in "${SCENARIOS[@]}"; do
    input_dir="${DATA_DIR}/${scenario}/TPM_metadata"
    output_dir="${OUT_DIR}/${scenario}"
    mkdir -p "${output_dir}"

    echo "============================================================"
    echo "Running CaDRReS-sc: ${scenario}"
    echo "============================================================"

    tpm_files=("${input_dir}"/*_TPM.csv)
    if [ ${#tpm_files[@]} -eq 0 ]; then
        echo "No TPM files found in ${input_dir}"
        continue
    fi

    for tpm_file in "${tpm_files[@]}"; do
        echo "Input: $(basename "${tpm_file}")"

        "${PYTHON_EXEC}" "${RUN_SCRIPT}" \
            --tpm "${tpm_file}" \
            --outdir "${output_dir}" \
            --model-dir "${MODEL_DIR}" \
            --gdsc-exp "${GDSC_EXP}" \
            --essential-genes "${ESSENTIAL_GENES}"
    done
done

echo "============================================================"
echo "CaDRReS-sc example completed."
echo "Results: ${OUT_DIR}"
echo "============================================================"
