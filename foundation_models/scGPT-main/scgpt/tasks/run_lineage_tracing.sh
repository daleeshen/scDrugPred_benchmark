#!/bin/bash

DATA_ROOT="/media/shenli/Expansion/DRMREF/LLM/scGPT-main/scdata"
MODEL_DIR="/media/shenli/Expansion/DRMREF/LLM/scGPT-main/scgpt/save"
OUTPUT_ROOT="/media/shenli/Expansion/DRMREF/LLM/scGPT-main/scgpt/samples"
RESULT_ROOT="/media/shenli/Expansion/DRMREF/LLM/scGPT-main/scgpt/results"

FILES=(
    "$DATA_ROOT/lineage_tracing/ReSisTrace_ola_input.txt"
    "$DATA_ROOT/lineage_tracing/Watermelon_round1_input.txt"
    "$DATA_ROOT/lineage_tracing/Watermelon_round2_input.txt"
    "$DATA_ROOT/lineage_tracing/Watermelon_round3_input.txt"
)

for file in "${FILES[@]}"
do

    echo ""
    echo "============================================================"
    echo "Processing:"
    echo "$file"
    echo "============================================================"

    relative_path="${file#$DATA_ROOT/}"
    relative_dir=$(dirname "$relative_path")
    filename=$(basename "$file")
    dataset_name="${filename%.txt}"

    EBD_FILE="$OUTPUT_ROOT/$relative_dir/${dataset_name}_samples.pkl"
    RESULT_FILE="$RESULT_ROOT/$relative_dir/${dataset_name}.csv"

    # ========================================================
    # Step 1: Generate scGPT embedding
    # ========================================================
    echo ""
    echo "-------------------- Step 1: Embedding --------------------"
    echo "Input : $file"
    echo "Output: $EBD_FILE"
    echo "------------------------------------------------------------"

    python stp1_get_ebd_auto.py \
        --input_file "$file" \
        --input_root "$DATA_ROOT" \
        --output_root "$OUTPUT_ROOT" \
        --model_dir "$MODEL_DIR" \
        --batch_size 64 \
        --max_length 1200

    if [ ! -f "$EBD_FILE" ]; then
        echo "ERROR: Embedding file was not generated:"
        echo "$EBD_FILE"
        echo "Skip benchmark for this dataset."
        continue
    fi

    # ========================================================
    # Step 2: Benchmark embedding
    # ========================================================
    echo ""
    echo "-------------------- Step 2: Benchmark --------------------"
    echo "Embedding: $EBD_FILE"
    echo "Result   : $RESULT_FILE"
    echo "------------------------------------------------------------"

    python stp2_benchmarking_main_EBD_auto.py \
        --data_file "$EBD_FILE" \
        --output_file "$RESULT_FILE"

    echo ""
    echo "============================================================"
    echo "Finished:"
    echo "$file"
    echo "============================================================"

done