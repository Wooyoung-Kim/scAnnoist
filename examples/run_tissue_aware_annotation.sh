#!/bin/bash
#
# Example script: Tissue-Aware CellTypist Annotation
# This script demonstrates how to use the tissue-aware annotation feature via CLI
#

set -e  # Exit on error

echo "=================================================="
echo "Tissue-Aware CellTypist Annotation Example"
echo "=================================================="

# Configuration
DATA_FILE="your_data.h5ad"
OUTPUT_DIR="./celltypist_output"
TISSUE_MAP="tissue_mapping.csv"
SAMPLE_TYPE_MAP="sample_type_mapping.csv"

# ============================================
# Example 1: Single tissue annotation
# ============================================
example_1() {
    echo ""
    echo "Example 1: Single Tissue Annotation"
    echo "------------------------------------"

    scrna-agent annotate \
        --file ${DATA_FILE} \
        --tissue "Spleen" \
        --species mouse \
        --sample-col sample \
        --output ${OUTPUT_DIR}/single_tissue.h5ad

    echo "✓ Single tissue annotation completed"
    echo "  Output: ${OUTPUT_DIR}/single_tissue.h5ad"
}

# ============================================
# Example 2: Multi-tissue with mapping file
# ============================================
example_2() {
    echo ""
    echo "Example 2: Multi-Tissue Annotation"
    echo "------------------------------------"

    # Create tissue mapping file
    cat > ${TISSUE_MAP} << EOF
sample,tissue
sample1,Spleen
sample2,Blood
sample3,Liver
sample4,Lung
EOF

    echo "Created tissue mapping:"
    cat ${TISSUE_MAP}

    # Run annotation
    scrna-agent annotate \
        --file ${DATA_FILE} \
        --tissue-map ${TISSUE_MAP} \
        --species mouse \
        --sample-col sample \
        --output ${OUTPUT_DIR}/multi_tissue.h5ad

    echo "✓ Multi-tissue annotation completed"
    echo "  Output: ${OUTPUT_DIR}/multi_tissue.h5ad"
}

# ============================================
# Example 3: With sample type information
# ============================================
example_3() {
    echo ""
    echo "Example 3: Annotation with Sample Type"
    echo "----------------------------------------"

    # Create tissue mapping
    cat > ${TISSUE_MAP} << EOF
sample,tissue
sample1,Blood
sample2,Blood
sample3,Lung
sample4,Spleen
EOF

    # Create sample type mapping
    cat > ${SAMPLE_TYPE_MAP} << EOF
sample,sample_type
sample1,PBMC
sample2,tumor
sample3,normal
sample4,normal
EOF

    echo "Created mappings:"
    echo "Tissue map:"
    cat ${TISSUE_MAP}
    echo ""
    echo "Sample type map:"
    cat ${SAMPLE_TYPE_MAP}

    # Run annotation
    scrna-agent annotate \
        --file ${DATA_FILE} \
        --tissue-map ${TISSUE_MAP} \
        --sample-type-map ${SAMPLE_TYPE_MAP} \
        --species mouse \
        --sample-col sample \
        --output ${OUTPUT_DIR}/with_sample_type.h5ad

    echo "✓ Annotation with sample type completed"
    echo "  Output: ${OUTPUT_DIR}/with_sample_type.h5ad"
}

# ============================================
# Example 4: Download models first
# ============================================
example_4() {
    echo ""
    echo "Example 4: Download Models"
    echo "--------------------------"

    # List available models
    echo "Available models in registry:"
    scrna-agent models list

    # Download models for specific tissue
    echo ""
    echo "Downloading models for Spleen..."
    scrna-agent models pull --tissue Spleen

    # Or download all models
    # echo "Downloading all models..."
    # scrna-agent models pull --all

    echo "✓ Model download completed"
}

# ============================================
# Example 5: Allow missing tissue data
# ============================================
example_5() {
    echo ""
    echo "Example 5: Allow Missing Tissue"
    echo "--------------------------------"

    # Create incomplete tissue mapping (missing some samples)
    cat > ${TISSUE_MAP} << EOF
sample,tissue
sample1,Spleen
sample2,Blood
EOF
    # Note: sample3 and sample4 are not included

    # Run annotation with --allow-missing-tissue flag
    scrna-agent annotate \
        --file ${DATA_FILE} \
        --tissue-map ${TISSUE_MAP} \
        --allow-missing-tissue \
        --species mouse \
        --sample-col sample \
        --output ${OUTPUT_DIR}/allow_missing.h5ad

    echo "✓ Annotation with missing tissue completed"
    echo "  Missing samples were labeled as 'unknown'"
}

# ============================================
# Main menu
# ============================================
show_menu() {
    echo ""
    echo "Select an example to run:"
    echo "  1) Single tissue annotation"
    echo "  2) Multi-tissue annotation"
    echo "  3) Annotation with sample type"
    echo "  4) Download models"
    echo "  5) Allow missing tissue"
    echo "  6) Run all examples"
    echo "  0) Exit"
    echo ""
}

run_all() {
    example_4  # Download models first
    example_1
    example_2
    example_3
    example_5
}

# ============================================
# Script execution
# ============================================

# Check if data file exists
if [ ! -f "${DATA_FILE}" ]; then
    echo "Warning: Data file '${DATA_FILE}' not found!"
    echo "Please update DATA_FILE variable in this script to point to your h5ad file."
    echo ""
    echo "For demonstration purposes, you can still run:"
    echo "  - Example 4 (Download models)"
    echo ""
fi

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Interactive mode if no arguments
if [ $# -eq 0 ]; then
    while true; do
        show_menu
        read -p "Enter choice [0-6]: " choice

        case $choice in
            1) example_1 ;;
            2) example_2 ;;
            3) example_3 ;;
            4) example_4 ;;
            5) example_5 ;;
            6) run_all ;;
            0) echo "Exiting..."; exit 0 ;;
            *) echo "Invalid option. Please try again." ;;
        esac
    done
else
    # Command-line mode
    case $1 in
        1) example_1 ;;
        2) example_2 ;;
        3) example_3 ;;
        4) example_4 ;;
        5) example_5 ;;
        all) run_all ;;
        *) echo "Usage: $0 [1|2|3|4|5|all]"; exit 1 ;;
    esac
fi

echo ""
echo "=================================================="
echo "All done! Check the output directory for results:"
echo "  ${OUTPUT_DIR}/"
echo "=================================================="
