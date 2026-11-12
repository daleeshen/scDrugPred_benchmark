#!/usr/bin/env bash
set -euo pipefail

SCAD_CODE_DIR="/media/shenli/Expansion/SCAD"
SPLIT_SCRIPT="${SCAD_CODE_DIR}/data/split_norm/split_data_SCAD_5fold_norm_auto.py"
LABEL_REPLACE_SCRIPT="./replace_SCAD_fold_labels.py"

PAIR_FILE="/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
BASE="/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source"
PYTHON_EXEC="${PYTHON_EXEC:-python}"

SCENARIOS=(
  "balanced round_1"
  "balanced round_2"
  "balanced round_3"
  "imbalance_1_3 round_1"
  "imbalance_1_3 round_2"
  "imbalance_1_3 round_3"
  "imbalance_1_30 round_1"
  "imbalance_1_30 round_2"
  "imbalance_1_30 round_3"
)

cd "$SCAD_CODE_DIR"

tail -n +2 "$PAIR_FILE" | while IFS=',' read -r DATASET DRUG REST; do
  DATASET=$(echo "$DATASET" | tr -d '"\r' | xargs)
  DRUG=$(echo "$DRUG" | tr -d '"\r' | xargs)

  for ITEM in "${SCENARIOS[@]}"; do
    SCENARIO=$(echo "$ITEM" | awk '{print $1}')
    ROUND=$(echo "$ITEM" | awk '{print $2}')
    DATA_PATH="${BASE}/${SCENARIO}/${ROUND}/"

    echo "=================================================="
    echo "Dataset  : $DATASET"
    echo "Drug     : $DRUG"
    echo "Scenario : $SCENARIO"
    echo "Round    : $ROUND"
    echo "=================================================="

    "$PYTHON_EXEC" "$SPLIT_SCRIPT" \
      --id "$DATASET" \
      --drug "$DRUG" \
      --solid_tumor 1 \
      --data_path "$DATA_PATH" \
      --n_splits 5 \
      --seed 42
  done
done

echo "SCAD original-label folds generated."

"$PYTHON_EXEC" "$LABEL_REPLACE_SCRIPT"

echo "scDEAL-label folds generated using identical SCAD splits."
