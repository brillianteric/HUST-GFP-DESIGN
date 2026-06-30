#!/bin/bash
# NOTE:
# This is an example SLURM template used for the 4EUL workflow.
# Edit partition names, module names, bind paths, and external tool paths
# before running on a different HPC cluster.
#SBATCH -J af3-stage1
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

set -euo pipefail

JSON_PATH="$1"
OUTPUT_DIR="$2"

: "${AF3_BASEDIR:?Set AF3_BASEDIR to the AlphaFold 3 Apptainer installation}"
: "${AF3_DB_DIR:?Set AF3_DB_DIR to the AlphaFold 3 database directory}"
: "${AF3_MODEL_DIR:?Set AF3_MODEL_DIR to the AlphaFold 3 model weights directory}"
: "${AF3_RUNTIME:=apptainer}"
: "${AF3_EXTRA_BINDS:=}"

source "${AF3_BASEDIR}/scripts/env.sh"

BIND_ARGS=()
for bind_path in ${AF3_EXTRA_BINDS}; do
  BIND_ARGS+=(--bind "${bind_path}")
done

"${AF3_RUNTIME}" exec \
  "${BIND_ARGS[@]}" \
  --bind "${AF3_BASEDIR}/case" \
  --bind "${AF3_DB_DIR}" \
  "${AF3_BASEDIR}/app/AlphaFold-v3.0.0-Apptainer.sif" \
  python /app/alphafold/run_alphafold.py \
    --norun_inference \
    --json_path="${JSON_PATH}" \
    --model_dir="${AF3_MODEL_DIR}" \
    --db_dir="${AF3_DB_DIR}" \
    --output_dir="${OUTPUT_DIR}"
