# HUST GFP Design

This repository contains a 4EUL-centered GFP design workflow for chromophore-aware fixed-position sequence design, structure prediction, and final candidate filtering.

The workflow includes:

1. ProteinMPNN sequence design with fixed functionally important residues.
2. AlphaFold 3 stage 1 data-pipeline search for MSA/template features.
3. AlphaFold 3 stage 2 structure prediction from augmented JSON files.
4. Structure-level and chromophore-microenvironment filtering for final candidate selection.

---

## Repository Layout

```text
data/
  raw/                         Prepared single-chain 4EUL PDB
  proteinmpnn/                 ProteinMPNN FASTA output for 4EUL designs
  fixed_positions/             Fixed-position JSONL and audit summary
  af3_inputs/                  Top-10 AlphaFold 3 input JSON files

hpc/
  run_proteinmpnn_4eul.slurm
  alphafold3_stage1_data_pipeline.slurm
  alphafold3_stage2_inference.slurm
  submit_af3_stage_jobs.sh

scripts/
  parse_multiple_chains_chromophore.py
  make_4eul_fixed_positions.py
  generate_af3_jsons_from_mpnn.py
  convert_all_af3_cif_to_pdb.py
  calc_structure_metrics_exclude_chromophore_1.py
  compute_microenv_metrics.py
  select_final_4eul_candidates.py

results/
  metrics/                     4EUL-only structure and microenvironment metrics
  selected_candidates/         Final selected candidates and audit table
  example_models/              Representative AF3-predicted PDB examples

outputs/                       Runtime outputs generated during reproduction;
                               ignored by git unless explicitly curated
```

---

## External Dependencies Not Included

This repository only contains project-specific inputs, scripts, processed metrics, and selected outputs. Large third-party runtime assets are intentionally not included.

The following external dependencies are **not** redistributed in this repository:

* ProteinMPNN official source code, runtime environment, and model weights;
* AlphaFold 3 official source code, model parameters, genetic/template databases, and Docker/Singularity/Apptainer images;
* large raw AlphaFold 3 runtime outputs;
* cluster-specific software modules, conda environments, and absolute HPC paths.

Users should install ProteinMPNN and AlphaFold 3 from their official repositories and provide local paths to the corresponding environments, model weights, databases, and containers in the HPC wrapper scripts.

Official repositories:

```text
ProteinMPNN:
https://github.com/dauparas/ProteinMPNN

AlphaFold 3:
https://github.com/google-deepmind/alphafold3
```

---

## Design Strategy

The design target is the 4EUL GFP scaffold. The main goal is to generate sequence variants while preserving the GFP fold, the chromophore-forming region, and the local chromophore microenvironment.

ProteinMPNN sampling was constrained by a chromophore-aware fixed-position strategy. The final fixed-position set combines:

1. literature-supported functionally important residues;
2. residues located within 5 Angstrom of the chromophore.

This strategy is intended to allow sequence diversification while keeping the structural and functional core of the fluorescent protein stable.

For 4EUL, the fixed-position files are recorded in:

```text
data/fixed_positions/fixed_4EUL_positions.jsonl
data/fixed_positions/fixed_4EUL_positions_summary.csv
```

The project-specific ProteinMPNN preprocessing scripts are:

```text
scripts/parse_multiple_chains_chromophore.py
scripts/make_4eul_fixed_positions.py
```

The chromophore-aware parser converts fluorescent-protein PDB files into ProteinMPNN-compatible parsed JSONL format while handling common GFP chromophore records such as `CRO`, `CR2`, and `NRQ`. The fixed-position script then uses the parsed JSONL and the original 4EUL PDB to define the final chromophore-aware fixed residues.

---

## Installation

Install the lightweight Python dependencies used by the project-specific preprocessing, JSON generation, and filtering scripts:

```bash
pip install -r requirements.txt
```

This does **not** install ProteinMPNN or AlphaFold 3. Those tools should be installed separately following their official repositories.

---

## Prepare 4EUL Input

Place the prepared single-chain 4EUL PDB file at:

```text
data/raw/4EUL.pdb
```

If redistribution of the PDB file is not allowed in a public release, provide a preparation note instead, including:

* PDB ID: `4EUL`;
* selected chain;
* any residue renumbering or cleanup procedure;
* whether non-protein records were retained for chromophore-aware parsing.

---

## Generate ProteinMPNN Fixed Positions

First, parse the prepared 4EUL PDB using the chromophore-aware parser:

```bash
python scripts/parse_multiple_chains_chromophore.py \
  --input_path data/raw \
  --output_path data/fixed_positions/parsed_pdbs_4EUL.jsonl
```

Then regenerate the chromophore-aware fixed-position file:

```bash
python scripts/make_4eul_fixed_positions.py \
  --pdb data/raw/4EUL.pdb \
  --parsed-jsonl data/fixed_positions/parsed_pdbs_4EUL.jsonl \
  --output-jsonl data/fixed_positions/fixed_4EUL_positions.jsonl \
  --summary-csv data/fixed_positions/fixed_4EUL_positions_summary.csv \
  --cutoff 5.0
```

Expected outputs:

```text
data/fixed_positions/parsed_pdbs_4EUL.jsonl
data/fixed_positions/fixed_4EUL_positions.jsonl
data/fixed_positions/fixed_4EUL_positions_summary.csv
```

---

## Run ProteinMPNN Sequence Design

ProteinMPNN is treated as an external dependency. The official ProteinMPNN source code, runtime environment, and model weights are not duplicated in this repository.

Set the local ProteinMPNN paths before submitting the job:

```bash
export CONDA_SH=~/miniconda3/etc/profile.d/conda.sh
export PROTEINMPNN_ENV=ProteinMPNN
export PROTEINMPNN_ROOT=/path/to/ProteinMPNN
```

Submit the ProteinMPNN design job:

```bash
sbatch hpc/run_proteinmpnn_4eul.slurm
```

The wrapper script performs three steps:

1. parse `data/raw/4EUL.pdb` with the chromophore-aware parser;
2. regenerate fixed positions using the 5 Angstrom chromophore shell;
3. run the official `protein_mpnn_run.py` entry point with the project-specific fixed-position constraints.

The design settings used for the included results were:

```text
num_seq_per_target = 20
sampling_temp      = 0.1
seed               = 37
batch_size         = 20
```

The expected ProteinMPNN FASTA output is:

```text
data/proteinmpnn/4EUL_designs.fa
```

---

## Generate AlphaFold 3 Input JSON Files

Generate AlphaFold 3 input JSON files from the ProteinMPNN FASTA output:

```bash
python scripts/generate_af3_jsons_from_mpnn.py \
  --input-fasta data/proteinmpnn/4EUL_designs.fa \
  --output-dir data/af3_inputs \
  --top-n 10 \
  --model-seed 1 \
  --allow-x
```

This produces the top-10 AF3 input JSON files and a file list for batch submission.

Expected outputs:

```text
data/af3_inputs/
data/af3_inputs/af3_json_files.list
```

---

## Run AlphaFold 3

AlphaFold 3 is treated as an external dependency. This repository does not redistribute the official AlphaFold 3 source code, model parameters, genetic/template databases, or container images.

Users should follow the official Google DeepMind AlphaFold 3 repository to:

1. obtain access to model parameters;
2. prepare the required genetic/template databases;
3. build or use a compatible Docker/Singularity/Apptainer runtime;
4. expose local paths through the environment variables below.

Set the local AF3 paths:

```bash
export AF3_BASEDIR=/path/to/AlphaFold-v3.0.0-Apptainer
export AF3_DB_DIR=/path/to/alphafold3/database
export AF3_MODEL_DIR=/path/to/alphafold3/models
```

The helper scripts in `hpc/` only wrap the standard AF3 two-stage execution pattern. They do not provide AF3 itself.

---

### AlphaFold 3 Stage 1: Data Pipeline

Stage 1 performs sequence/template database search and generates augmented `*_data.json` files containing MSA/template features.

Submit stage 1 jobs:

```bash
bash hpc/submit_af3_stage_jobs.sh \
  stage1 \
  data/af3_inputs/af3_json_files.list \
  outputs/af3_stage1 \
  1 10
```

Expected stage 1 outputs:

```text
outputs/af3_stage1/
```

The generated `*_data.json` files are used as input for stage 2.

---

### AlphaFold 3 Stage 2: Structure Inference

After stage 1 finishes, create a list of augmented JSON files:

```bash
find outputs/af3_stage1 -name "*_data.json" | sort > outputs/stage2_input_files.list
```

Submit stage 2 jobs:

```bash
bash hpc/submit_af3_stage_jobs.sh \
  stage2 \
  outputs/stage2_input_files.list \
  outputs/af3_stage2 \
  1 10
```

Stage 2 runs AF3 inference from the augmented JSON files, typically with `--norun_data_pipeline`, and therefore does not repeat the database search.

Expected stage 2 outputs:

```text
outputs/af3_stage2/
```

Large raw AF3 runtime outputs are ignored by git unless manually curated into `results/`.

---

## Recompute Metrics

After AF3 stage 2 finishes, convert predicted mmCIF models to PDB files:

```bash
python scripts/convert_all_af3_cif_to_pdb.py
```

Compute 4EUL structure-level metrics and merge AF3 confidence metrics:

```bash
python scripts/calc_structure_metrics_exclude_chromophore_1.py
```

The `--confidence-dir` should point to raw AF3 stage 2 outputs containing
`*_confidences.json`. If those files are absent, the script still computes RMSD
and approximate TM-score, but pLDDT/PAE columns are left blank and final
filtering cannot be fully reproduced from scratch.

Compute chromophore microenvironment metrics:

```bash
python scripts/compute_microenv_metrics.py
```

Expected metric outputs:

```text
results/metrics/4EUL_structure_metrics_with_plddt_pae.csv
results/metrics/4EUL_microenv_metrics.csv
```

---

## Candidate Selection

The final selection step combines structure-level metrics and chromophore-microenvironment metrics.

Run:

```bash
python scripts/select_final_4eul_candidates.py \
  --structure results/metrics/4EUL_structure_metrics_with_plddt_pae.csv \
  --microenv results/metrics/4EUL_microenv_metrics.csv \
  --out results/selected_candidates/selected_4EUL_candidate_ids.csv \
  --audit-out results/selected_candidates/selected_4EUL_audit.csv
```

Expected outputs:

```text
results/selected_candidates/selected_4EUL_candidate_ids.csv
results/selected_candidates/selected_4EUL_audit.csv
```

---

## Filtering Criteria

The final selection script applies conservative structure and microenvironment cutoffs:

| Metric                                                         |      Cutoff |
| -------------------------------------------------------------- | ----------: |
| CA RMSD excluding chromophore-forming residues                 |  `<= 0.325` |
| Approximate CA TM-score excluding chromophore-forming residues | `>= 0.9965` |
| Mean atom pLDDT                                                |   `>= 94.8` |
| Mean PAE                                                       |   `<= 2.17` |
| 5 Angstrom shell CA RMSD                                       |  `<= 0.325` |
| 8 Angstrom shell CA RMSD                                       |   `<= 0.29` |
| Chromophore contact F1                                         |   `>= 0.84` |
| Chromophore clash count within 2 Angstrom                      |      `<= 5` |
| 5 Angstrom shell mutation count                                |       `= 0` |

These criteria are intended to retain candidates with high predicted structural confidence, minimal scaffold distortion, preserved chromophore-shell geometry, and no mutations in the immediate 5 Angstrom chromophore shell.

---

## Selected Candidates

The included audit currently selects six final 4EUL candidates:

```text
4EUL_rank01_score0.5902
4EUL_rank02_score0.5915
4EUL_rank05_score0.6056
4EUL_rank06_score0.6123
4EUL_rank07_score0.6176
4EUL_rank09_score0.6225
```

The selected candidate list and audit table are provided in:

```text
results/selected_candidates/selected_4EUL_candidate_ids.csv
results/selected_candidates/selected_4EUL_audit.csv
```

Representative predicted structures can be placed under:

```text
results/example_models/
```

---

## Reproducibility Notes

For full reproduction, users should provide or verify the following external resources and version information:

1. **4EUL structure source and preparation**

   * PDB ID and chain used;
   * single-chain extraction procedure;
   * whether chromophore records were retained;
   * any residue cleanup or renumbering.

2. **ProteinMPNN setup**

   * official ProteinMPNN commit or release;
   * model checkpoint/version;
   * conda or Python environment;
   * command-line parameters used for design.

3. **AlphaFold 3 setup**

   * official AlphaFold 3 commit or release;
   * Docker/Singularity/Apptainer container version;
   * local model parameter directory;
   * database release or HPC database path;
   * stage 1 and stage 2 command-line parameters.

4. **HPC-specific settings**

   * SLURM partition names;
   * GPU type;
   * module names;
   * local path variables;
   * output storage location.

5. **Metric generation**

   * post-AF3 scripts or documented commands used to regenerate:

     * `results/metrics/4EUL_structure_metrics_with_plddt_pae.csv`;
     * `results/metrics/4EUL_microenv_metrics.csv`.

The included `hpc/` scripts are templates and may need small edits for a specific cluster environment.

---

## Files Intentionally Excluded

The following files are intentionally excluded from version control:

```text
outputs/
*.sif
*.pt
*.pth
*.params
*.cif.gz
large raw AF3 output directories
third-party model checkpoints
third-party genetic databases
```

This keeps the repository lightweight and avoids redistributing third-party assets that should be obtained from their official sources.

---

## Citation and License Notes

Please cite the relevant original resources when using this workflow or derived results:

* ProteinMPNN official repository and publication;
* AlphaFold 3 official repository and publication;
* Protein Data Bank;
* the original 4EUL structure source.

A project license file should be added before public release. Third-party tools, model parameters, databases, and structures remain subject to their own licenses and terms of use.

---

## Minimal Reproduction Summary

A minimal reproduction of the project-specific workflow is:

```bash
# 1. Install project script dependencies
pip install -r requirements.txt

# 2. Parse 4EUL with the chromophore-aware parser
python scripts/parse_multiple_chains_chromophore.py \
  --input_path data/raw \
  --output_path data/fixed_positions/parsed_pdbs_4EUL.jsonl

# 3. Generate chromophore-aware fixed positions
python scripts/make_4eul_fixed_positions.py \
  --pdb data/raw/4EUL.pdb \
  --parsed-jsonl data/fixed_positions/parsed_pdbs_4EUL.jsonl \
  --output-jsonl data/fixed_positions/fixed_4EUL_positions.jsonl \
  --summary-csv data/fixed_positions/fixed_4EUL_positions_summary.csv \
  --cutoff 5.0

# 4. Run ProteinMPNN through the HPC wrapper
sbatch hpc/run_proteinmpnn_4eul.slurm

# 5. Generate top-10 AF3 input JSON files
python scripts/generate_af3_jsons_from_mpnn.py \
  --input-fasta data/proteinmpnn/4EUL_designs.fa \
  --output-dir data/af3_inputs \
  --top-n 10 \
  --model-seed 1 \
  --allow-x

# 6. Run AF3 stage 1
bash hpc/submit_af3_stage_jobs.sh \
  stage1 \
  data/af3_inputs/af3_json_files.list \
  outputs/af3_stage1 \
  1 10

# 7. Run AF3 stage 2
find outputs/af3_stage1 -name "*_data.json" | sort > outputs/stage2_input_files.list

bash hpc/submit_af3_stage_jobs.sh \
  stage2 \
  outputs/stage2_input_files.list \
  outputs/af3_stage2 \
  1 10

# 8. Select final candidates from included metrics
python scripts/select_final_4eul_candidates.py \
  --structure results/metrics/4EUL_structure_metrics_with_plddt_pae.csv \
  --microenv results/metrics/4EUL_microenv_metrics.csv \
  --out results/selected_candidates/selected_4EUL_candidate_ids.csv \
  --audit-out results/selected_candidates/selected_4EUL_audit.csv
```
