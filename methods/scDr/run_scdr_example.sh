#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RSCRIPT_EXEC="${RSCRIPT_EXEC:-Rscript}"

RUN_SCRIPT="${SCRIPT_DIR}/run_scdr.R"
SCDR_DIR="${SCRIPT_DIR}/scDr"
UP_SIGNATURE="${SCDR_DIR}/drug_up_signature.RDS"
DOWN_SIGNATURE="${SCDR_DIR}/drug_down_signature.RDS"

DATA_DIR="${ROOT_DIR}/data/example"
OUT_DIR="${SCRIPT_DIR}/example_output"

SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

mkdir -p "${OUT_DIR}"
shopt -s nullglob

for scenario in "${SCENARIOS[@]}"; do
    input_dir="${DATA_DIR}/${scenario}"
    output_dir="${OUT_DIR}/${scenario}"

    mkdir -p "${output_dir}"

    echo "============================================================"
    echo "Running scDr: ${scenario}"
    echo "============================================================"

    input_files=("${input_dir}"/*.RData "${input_dir}"/*.rds)

    if [ ${#input_files[@]} -eq 0 ]; then
        echo "No Seurat files found in ${input_dir}"
        continue
    fi

    for input_file in "${input_files[@]}"; do
        echo "------------------------------------------------------------"
        echo "Input: $(basename "${input_file}")"
        echo "------------------------------------------------------------"

        "${RSCRIPT_EXEC}" "${RUN_SCRIPT}" \
            "${input_file}" \
            "${output_dir}" \
            "${UP_SIGNATURE}" \
            "${DOWN_SIGNATURE}"

        echo "$(basename "${input_file}") completed."
    done
done

echo "============================================================"
echo "scDr example completed."
echo "Results: ${OUT_DIR}"
echo "============================================================"
