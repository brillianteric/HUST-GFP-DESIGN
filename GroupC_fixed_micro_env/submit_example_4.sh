#!/bin/bash
#SBATCH -p qdagexclu11
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --output=example_4_GroupC_fixed.out

set -euo pipefail

# =========================
# 环境
# =========================
module purge
module load nvidia/cuda/11.3

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ProteinMPNN
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

echo "========== ENV =========="
hostname
which python
python --version
nvidia-smi || true


# =========================
# 路径设置
# =========================
WORK_DIR="/work/home/mayongze/ProteinMPNN_project/Graduation_Design/GroupC_fixed_micro_env"
cd "$WORK_DIR"

MPNN_ROOT="/work/share/acssag9hf5/AI_protein_generation/Models/ProteinMPNN"

folder_with_pdbs="only_4EUL"

output_dir="only_4EUL_seq"

PARSE_SCRIPT="/work/home/mayongze/ProteinMPNN_project/Graduation_Design/parse_multiple_chains_chromophore.py"
MAKE_GROUPC_SCRIPT="make_groupC_fixed_positions_1.py"

mkdir -p "$output_dir"

path_for_parsed_chains="${output_dir}/parsed_pdbs_GroupC.jsonl"
path_for_fixed_positions="${output_dir}/fixed_GroupC_positions.jsonl"
path_for_fixed_summary="${output_dir}/fixed_GroupC_positions_summary.csv"

echo "========== PATHS =========="
echo "WORK_DIR=$WORK_DIR"
echo "MPNN_ROOT=$MPNN_ROOT"
echo "folder_with_pdbs=$folder_with_pdbs"
echo "output_dir=$output_dir"
echo "PARSE_SCRIPT=$PARSE_SCRIPT"
echo "MAKE_GROUPC_SCRIPT=$MAKE_GROUPC_SCRIPT"
echo "path_for_parsed_chains=$path_for_parsed_chains"
echo "path_for_fixed_positions=$path_for_fixed_positions"
echo "path_for_fixed_summary=$path_for_fixed_summary"


# =========================
# Step 1. chromophore-aware parse
# =========================
echo
echo "========== STEP 1: parse PDBs with chromophore-aware parser =========="

python "$PARSE_SCRIPT" \
  --input_path "$folder_with_pdbs" \
  --output_path "$path_for_parsed_chains"

echo "[OK] parsed jsonl: $path_for_parsed_chains"


# =========================
# Step 2. 生成 GroupC fixed positions
# GroupC = GroupB + 发色团 5 Å 壳层残基
# =========================
echo
echo "========== STEP 2: make GroupC fixed positions =========="

python "$MAKE_GROUPC_SCRIPT" \
  --pdb_dir "$folder_with_pdbs" \
  --parsed_jsonl "$path_for_parsed_chains" \
  --output_jsonl "$path_for_fixed_positions" \
  --summary_csv "$path_for_fixed_summary" \
  --cutoff 5.0 \
  --chain A

echo "[OK] fixed positions jsonl: $path_for_fixed_positions"
cat "$path_for_fixed_positions"

echo
echo "[OK] fixed positions summary:"
cat "$path_for_fixed_summary"


# =========================
# Step 3. 简单检查 fixed positions 是否覆盖所有 target
# =========================
echo
echo "========== STEP 3: validate fixed positions =========="

python - <<PY
import json
from pathlib import Path

parsed_path = Path("$path_for_parsed_chains")
fixed_path = Path("$path_for_fixed_positions")

parsed_names = []
for line in parsed_path.open():
    if line.strip():
        parsed_names.append(json.loads(line)["name"])

fixed = json.loads(fixed_path.read_text().strip())
fixed_names = list(fixed.keys())

print("[INFO] parsed names:", parsed_names)
print("[INFO] fixed names :", fixed_names)

missing = sorted(set(parsed_names) - set(fixed_names))
extra = sorted(set(fixed_names) - set(parsed_names))

if missing:
    raise SystemExit(f"[ERROR] missing fixed positions for: {missing}")

if extra:
    print(f"[WARN] extra fixed names: {extra}")

for name in parsed_names:
    chain_dict = fixed[name]
    for chain, positions in chain_dict.items():
        print(f"[OK] {name} chain {chain}: {len(positions)} fixed positions")

print("[OK] fixed positions validation passed")
PY


# =========================
# Step 4. Run ProteinMPNN
# =========================
echo
echo "========== STEP 4: run ProteinMPNN GroupC =========="

python "$MPNN_ROOT/protein_mpnn_run.py" \
  --jsonl_path "$path_for_parsed_chains" \
  --fixed_positions_jsonl "$path_for_fixed_positions" \
  --out_folder "$output_dir" \
  --num_seq_per_target 20 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 20

echo
echo "========== FINISHED =========="
date