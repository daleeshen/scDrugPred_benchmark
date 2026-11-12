#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RSCRIPT_EXEC="${RSCRIPT_EXEC:-Rscript}"

RUN_SCRIPT="${SCRIPT_DIR}/run_beyondcell.R"
DATA_DIR="${ROOT_DIR}/data/example"
OUT_DIR="${SCRIPT_DIR}/example_output"

TARGET_DRUG="PLX-4720"
COLLECTION="PSc"
SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

mkdir -p "${OUT_DIR}"

for scenario in "${SCENARIOS[@]}"; do
    input_dir="${DATA_DIR}/${scenario}"
    output_dir="${OUT_DIR}/${scenario}"
    mkdir -p "${output_dir}"

    echo "============================================================"
    echo "Running Beyondcell: ${scenario}"
    echo "============================================================"

    found=false
    while IFS= read -r -d '' input_file; do
        found=true
        echo "Input: $(basename "${input_file}")"

        "${RSCRIPT_EXEC}" "${RUN_SCRIPT}" \
            "${input_file}" \
            "${TARGET_DRUG}" \
            "${output_dir}" \
            "${COLLECTION}"
    done < <(
        find "${input_dir}" -maxdepth 1 -type f \
            \( -name "*.RData" -o -name "*.rds" \) \
            -print0 | sort -z
    )

    if [ "${found}" = false ]; then
        echo "No Seurat files found in ${input_dir}"
    fi
done

echo "============================================================"
echo "Beyondcell example completed."
echo "Results: ${OUT_DIR}"
echo "============================================================"
