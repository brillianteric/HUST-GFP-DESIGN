#!/bin/bash
# NOTE:
# This is an example SLURM template used for the 4EUL workflow.
# Edit partition names, module names, bind paths, and external tool paths
# before running on a different HPC cluster.
set -euo pipefail

STAGE="${1:?Usage: $0 stage1|stage2 JSON_LIST OUTPUT_DIR [START_ID] [END_ID]}"
JSON_LIST="${2:?JSON list file is required}"
OUTPUT_DIR="${3:?Output directory is required}"
START_ID="${4:-1}"
END_ID="${5:-999999}"

case "${STAGE}" in
  stage1) SLURM_SCRIPT="hpc/alphafold3_stage1_data_pipeline.slurm" ;;
  stage2) SLURM_SCRIPT="hpc/alphafold3_stage2_inference.slurm" ;;
  *) echo "Unknown stage: ${STAGE}" >&2; exit 2 ;;
esac

mkdir -p logs "${OUTPUT_DIR}"
mapfile -t JSON_FILES < "${JSON_LIST}"

for ((idx=START_ID; idx<=END_ID && idx<=${#JSON_FILES[@]}; idx++)); do
  json_path="${JSON_FILES[$((idx-1))]}"
  name="$(basename "${json_path}" .json)"
  out_dir="${OUTPUT_DIR}/${name}"
  log_file="logs/${STAGE}_${name}.out"

  if [[ -s "${out_dir}/${name}_data.json" || -s "${out_dir}/summary_confidences.json" ]]; then
    echo "[SKIP] existing output: ${out_dir}"
    continue
  fi

  echo "[SUBMIT] ${STAGE} ${idx}: ${json_path}"
  sbatch -o "${log_file}" "${SLURM_SCRIPT}" "${json_path}" "${OUTPUT_DIR}"
  sleep 2
done
