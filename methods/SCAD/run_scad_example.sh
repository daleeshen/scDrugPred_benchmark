#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="${SCRIPT_DIR}/SCAD"
PYTHON_EXEC="${PYTHON_EXEC:-/home/shenli/py_env/main_env/bin/python}"

ID=""
DRUG=""
COUNT_FILE=""
META_FILE=""
OUTPUT_DIR=""
SOLID_TUMOR=1
LABEL_COL="condition"
RES_LABEL="resistant"
SENS_LABEL="sensitive"
MAX_ITER=20
GPU=1
OVERWRITE=0

usage() {
    echo "Usage:"
    echo "  $0 --id ID --drug DRUG --count COUNT.csv --meta META.csv [options]"
    echo ""
    echo "Required:"
    echo "  --id ID"
    echo "  --drug DRUG"
    echo "  --count FILE"
    echo "  --meta FILE"
    echo ""
    echo "Options:"
    echo "  --output-dir DIR       Output directory"
    echo "  --solid-tumor 0|1      Default: 1"
    echo "  --label-col NAME       Default: condition"
    echo "  --res-label LABEL      Default: resistant"
    echo "  --sens-label LABEL     Default: sensitive"
    echo "  --max-iter N           Hyperparameter search iterations, default: 20"
    echo "  --gpu 0|1              Default: 1"
    echo "  --overwrite            Remove existing output directory before running"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --id) ID="$2"; shift 2 ;;
        --drug) DRUG="$2"; shift 2 ;;
        --count) COUNT_FILE="$2"; shift 2 ;;
        --meta) META_FILE="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --solid-tumor) SOLID_TUMOR="$2"; shift 2 ;;
        --label-col) LABEL_COL="$2"; shift 2 ;;
        --res-label) RES_LABEL="$2"; shift 2 ;;
        --sens-label) SENS_LABEL="$2"; shift 2 ;;
        --max-iter) MAX_ITER="$2"; shift 2 ;;
        --gpu) GPU="$2"; shift 2 ;;
        --overwrite) OVERWRITE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1"; usage; exit 1 ;;
    esac
done

[[ -n "$ID" ]] || { echo "ERROR: --id is required"; exit 1; }
[[ -n "$DRUG" ]] || { echo "ERROR: --drug is required"; exit 1; }
[[ -n "$COUNT_FILE" ]] || { echo "ERROR: --count is required"; exit 1; }
[[ -n "$META_FILE" ]] || { echo "ERROR: --meta is required"; exit 1; }
[[ -f "$COUNT_FILE" ]] || { echo "ERROR: Count matrix not found: $COUNT_FILE"; exit 1; }
[[ -f "$META_FILE" ]] || { echo "ERROR: Metadata not found: $META_FILE"; exit 1; }

case "$DRUG" in
    PLX4720|PLX-4720)
        PREP_DRUG="PLX-4720"
        MODEL_DRUG="PLX4720"
        ;;
    Palbociclib|PD-0332991)
        PREP_DRUG="Palbociclib"
        MODEL_DRUG="PD-0332991"
        ;;
    TAE684|NVP-TAE684)
        PREP_DRUG="TAE684"
        MODEL_DRUG="NVP-TAE684"
        ;;
    *)
        PREP_DRUG="$DRUG"
        MODEL_DRUG="$DRUG"
        ;;
esac

if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${SCRIPT_DIR}/output/${ID}_${MODEL_DRUG}"
fi

if [[ "$OVERWRITE" -eq 1 && -d "$OUTPUT_DIR" ]]; then
    rm -rf "$OUTPUT_DIR"
fi

if [[ -d "$OUTPUT_DIR" && -n "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory is not empty: $OUTPUT_DIR"
    echo "Use --overwrite to replace it."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

if [[ "$MODEL_DRUG" == "Osimertinib" && -f "${CODE_DIR}/data/response/GDSC_response.all_drugs_with_Osimertinib.tsv" ]]; then
    GDSC_RESPONSE="${GDSC_RESPONSE:-${CODE_DIR}/data/response/GDSC_response.all_drugs_with_Osimertinib.tsv}"
    GDSC_IC50="${GDSC_IC50:-${CODE_DIR}/data/response/GDSC_response.logIC50.all_drugs_with_Osimertinib.tsv}"
else
    GDSC_RESPONSE="${GDSC_RESPONSE:-${CODE_DIR}/data/response/GDSC_response.all_drugs.tsv}"
    GDSC_IC50="${GDSC_IC50:-${CODE_DIR}/data/response/GDSC_response.logIC50.all_drugs.tsv}"
fi

HGNC="${HGNC:-${CODE_DIR}/data/original/HGNC_symbol_all_genes.tsv}"
GDSC_RMA="${GDSC_RMA:-${CODE_DIR}/data/original/gdsc-rma_gene-expression.csv}"
COSMIC="${COSMIC:-${CODE_DIR}/data/response/TableS1E.xlsx}"

for f in "$GDSC_RESPONSE" "$GDSC_IC50" "$HGNC" "$GDSC_RMA" "$COSMIC"; do
    [[ -f "$f" ]] || { echo "ERROR: Required SCAD reference file not found: $f"; exit 1; }
done

echo "================================================================================"
echo "Running SCAD"
echo "Dataset:        ${ID}"
echo "Drug:           ${MODEL_DRUG}"
echo "Count matrix:   ${COUNT_FILE}"
echo "Metadata:       ${META_FILE}"
echo "Output:         ${OUTPUT_DIR}"
echo "================================================================================"

echo ">> Step 1/5: single-cell preprocessing"
"$PYTHON_EXEC" "${CODE_DIR}/Scanpy_data_preprocess_auto.py" \
    --input "$COUNT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --id "$ID" \
    --drug "$PREP_DRUG"

echo ">> Step 2/5: source/target preprocessing"
"$PYTHON_EXEC" "${CODE_DIR}/data/processing/preprocessing_GSE_drugs_norm_auto.py" \
    --id "$ID" \
    --drug "$PREP_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --meta_path "$META_FILE" \
    --data_path "$OUTPUT_DIR" \
    --gdsc_response_path "$GDSC_RESPONSE" \
    --gdsc_ic50_path "$GDSC_IC50" \
    --hgnc_path "$HGNC" \
    --gdsc_rma_path "$GDSC_RMA" \
    --cosmic_path "$COSMIC" \
    --label_col "$LABEL_COL" \
    --res_label "$RES_LABEL" \
    --sens_label "$SENS_LABEL"

echo ">> Step 3/5: 5-fold splitting"
"$PYTHON_EXEC" "${CODE_DIR}/data/split_norm/split_data_SCAD_5fold_norm_auto.py" \
    --id "$ID" \
    --drug "$MODEL_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --response_file "$GDSC_RESPONSE" \
    --data_path "${OUTPUT_DIR}/data/split_norm/" \
    --n_splits 5 \
    --seed 42

echo ">> Step 4a/5: hyperparameter search - non-SMOTE"
"$PYTHON_EXEC" "${CODE_DIR}/model/SCAD_5foldsCV_on_source_arg_auto.py" \
    --id "$ID" \
    --drug "$MODEL_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --base_path "$OUTPUT_DIR" \
    --gpu "$GPU" \
    --max_iter "$MAX_ITER"

echo ">> Step 4b/5: hyperparameter search - SMOTE"
"$PYTHON_EXEC" "${CODE_DIR}/model/SCAD_5foldsCV_on_source_arg_smote_auto.py" \
    --id "$ID" \
    --drug "$MODEL_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --base_path "$OUTPUT_DIR" \
    --gpu "$GPU" \
    --max_iter "$MAX_ITER"

echo ">> Step 5a/5: final model - non-SMOTE"
"$PYTHON_EXEC" "${CODE_DIR}/model/SCAD_train_binarized_5folds_auto.py" \
    --id "$ID" \
    --drug "$MODEL_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --base_path "$OUTPUT_DIR" \
    --gpu "$GPU"

echo ">> Step 5b/5: final model - SMOTE"
"$PYTHON_EXEC" "${CODE_DIR}/model/SCAD_train_binarized_5folds_smote_auto.py" \
    --id "$ID" \
    --drug "$MODEL_DRUG" \
    --solid_tumor "$SOLID_TUMOR" \
    --base_path "$OUTPUT_DIR" \
    --gpu "$GPU"

echo "================================================================================"
echo "SCAD completed: ${ID} | ${MODEL_DRUG}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Prediction files:"
find "$OUTPUT_DIR" -type f -name "*predicted_y.tsv" -print
echo "================================================================================"
