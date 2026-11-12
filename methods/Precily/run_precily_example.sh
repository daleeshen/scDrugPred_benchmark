#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RSCRIPT_EXEC="${RSCRIPT_EXEC:-Rscript}"

RUN_SCRIPT="${SCRIPT_DIR}/run_precily.R"
DATA_DIR="${ROOT_DIR}/data/example"
OUT_DIR="${SCRIPT_DIR}/example_output"

TARGET_DRUG="PLX-4720"
CANCER_TYPE="SKCM"

SCENARIOS=("balanced_1000" "imbalanced_1_3" "imbalanced_1_30")

export LD_LIBRARY_PATH="/opt/conda/envs/precily/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "${OUT_DIR}"
shopt -s nullglob

for scenario in "${SCENARIOS[@]}"; do
    input_dir="${DATA_DIR}/${scenario}/precily_input"
    output_dir="${OUT_DIR}/${scenario}"

    mkdir -p "${output_dir}"

    echo "============================================================"
    echo "Running Precily: ${scenario}"
    echo "Drug: ${TARGET_DRUG} | Cancer type: ${CANCER_TYPE}"
    echo "============================================================"

    logtpm_files=("${input_dir}"/*_logTPM.csv)

    if [ ${#logtpm_files[@]} -eq 0 ]; then
        echo "No Precily input files found in ${input_dir}"
        continue
    fi

    for logtpm_file in "${logtpm_files[@]}"; do
        filename="$(basename "${logtpm_file}")"
        sample_id="${filename%_logTPM.csv}"
        metadata_file="${input_dir}/${sample_id}_meta_data.csv"

        if [ ! -f "${metadata_file}" ]; then
            echo "Metadata file not found: ${metadata_file}"
            exit 1
        fi

        echo "------------------------------------------------------------"
        echo "Sample: ${sample_id}"
        echo "------------------------------------------------------------"

        "${RSCRIPT_EXEC}" "${RUN_SCRIPT}" \
            "${logtpm_file}" \
            "${metadata_file}" \
            "${TARGET_DRUG}" \
            "${CANCER_TYPE}" \
            "${output_dir}"

        echo "${sample_id} completed."
    done
done

echo "============================================================"
echo "Precily example completed."
echo "Results: ${OUT_DIR}"
echo "============================================================"
