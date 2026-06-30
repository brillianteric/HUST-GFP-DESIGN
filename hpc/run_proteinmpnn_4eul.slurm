#!/bin/bash
#SBATCH -J mpnn-4eul
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --output=logs/proteinmpnn_4eul.out

set -euo pipefail

module purge
: "${CUDA_MODULE:=nvidia/cuda/11.3}"
module load "${CUDA_MODULE}"

: "${CONDA_SH:?Set CONDA_SH, for example ~/miniconda3/etc/profile.d/conda.sh}"
: "${PROTEINMPNN_ENV:=ProteinMPNN}"
: "${PROTEINMPNN_ROOT:?Set PROTEINMPNN_ROOT to the ProteinMPNN repository}"

source "${CONDA_SH}"
conda activate "${PROTEINMPNN_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
CHROMO_PARSE_SCRIPT="${CHROMO_PARSE_SCRIPT:-${REPO_ROOT}/scripts/parse_multiple_chains_chromophore.py}"
PDB_DIR="${REPO_ROOT}/data/raw"
OUTPUT_DIR="${REPO_ROOT}/outputs/proteinmpnn_4eul"
PARSED_JSONL="${OUTPUT_DIR}/parsed_pdbs_4EUL.jsonl"
FIXED_JSONL="${OUTPUT_DIR}/fixed_4EUL_positions.jsonl"
FIXED_SUMMARY="${OUTPUT_DIR}/fixed_4EUL_positions_summary.csv"

mkdir -p "${OUTPUT_DIR}" "${REPO_ROOT}/logs"

echo "========== ENV =========="
hostname
which python
python --version
nvidia-smi || true

echo "========== STEP 1: chromophore-aware PDB parsing =========="
python "${CHROMO_PARSE_SCRIPT}" \
  --input_path "${PDB_DIR}" \
  --output_path "${PARSED_JSONL}"

echo "========== STEP 2: 4EUL fixed positions =========="
python "${REPO_ROOT}/scripts/make_4eul_fixed_positions.py" \
  --pdb "${PDB_DIR}/4EUL.pdb" \
  --parsed-jsonl "${PARSED_JSONL}" \
  --output-jsonl "${FIXED_JSONL}" \
  --summary-csv "${FIXED_SUMMARY}" \
  --cutoff 5.0

echo "========== STEP 3: ProteinMPNN sequence design =========="
python "${PROTEINMPNN_ROOT}/protein_mpnn_run.py" \
  --jsonl_path "${PARSED_JSONL}" \
  --fixed_positions_jsonl "${FIXED_JSONL}" \
  --out_folder "${OUTPUT_DIR}" \
  --num_seq_per_target "${NUM_SEQ_PER_TARGET:-20}" \
  --sampling_temp "${SAMPLING_TEMP:-0.1}" \
  --seed "${MPNN_SEED:-37}" \
  --batch_size "${BATCH_SIZE:-20}"

echo "========== DONE =========="
echo "ProteinMPNN output: ${OUTPUT_DIR}"
