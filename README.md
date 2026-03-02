## MPNN Affinity Maturation

### Overview

This repository contains a small toolkit for protein interface design and affinity maturation built around:

- **ProteinMPNN** (vendored in `ProteinMPNN/`) for sequence design given a fixed backbone.
- **FlowPacker**-style side-chain sampling for generating all-atom models from designed sequences (`flowpacker_sampler_pipe.py`, `sampler_pdb_colab.py`).
- **Interface energy analysis** tools using Rosetta log files (`get_interface_energy.py`) or direct PyRosetta scoring (`per_residue_energy_pyrosetta.py`).

The typical use case is to:

1. Design sequences on a protein–protein complex backbone with ProteinMPNN.
2. Build side-chain conformations with FlowPacker.
3. Quantify per-residue and pairwise interface energies to guide affinity maturation.

For details on core ProteinMPNN functionality, see the upstream README in `ProteinMPNN/README.md`.

### Interactive notebook

You can explore a small example interactively in Google Colab:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nedru004/MPNN_affinity_maturation/blob/main/test_notebook.ipynb)

### Repository layout

- `ProteinMPNN/` – Upstream ProteinMPNN code, license, and README.
  - `protein_mpnn_run.py` – Main ProteinMPNN design script (see its README for flags and examples).
  - `protein_mpnn_utils.py` and `helper_scripts/` – Utilities for parsing PDBs, defining designed/fixed chains, tying positions, biasing AAs, etc.
- `flowpacker_sampler_pipe.py` – FlowPacker-based sampling pipeline that:
  - Pre-processes PDBs based on a CSV of sequences.
  - Runs a trained FlowPacker model to generate side-chain conformations.
  - Writes sampled PDBs and basic side-chain/chi metrics.
- `sampler_pdb_colab.py` – Earlier/alternative script for running the same FlowPacker sampler, with support for:
  - Supplying either a CSV with sequences **or** a directory of PDBs.
  - Preparing batch PDBs for the FlowPacker dataloader.
- `get_interface_energy.py` – Post-processing script for **Rosetta** output logs:
  - Parses per–residue–pair energies from Rosetta `.out` files.
  - Identifies interface residue pairs based on Cβ/Cα distances.
  - Aggregates and plots per-residue interface energies on the binder chain.
- `per_residue_energy_pyrosetta.py` – Pure **PyRosetta** version of the interface energy analysis:
  - Runs a FastRelax-like protocol around the interface.
  - Computes pairwise residue energies directly via PyRosetta scoring.
  - Produces CSVs and plots similar to `get_interface_energy.py`.
- `test/` – Example PDBs (e.g. `4dn4_ccl2_8_4.pdb` and processed variants) for testing the pipeline.
- `test_notebook.ipynb` – Scratch notebook for interactive experimentation.
- `requirements.txt` – Minimal Python requirements used during development (currently only `torch` is pinned; see below for additional dependencies).

### Dependencies and installation

You will typically want a recent Python (≥3.9) and CUDA-enabled PyTorch for non-trivial models.

- **Base Python packages** (install via `pip` or `conda`), in addition to those in `requirements.txt`:
  - `pandas`, `numpy`, `tqdm`, `pyyaml`, `easydict`
  - `biopython`, `matplotlib`, `seaborn`
  - `torch` (if not already installed to match your CUDA stack)
- **External toolkits**:
  - **ProteinMPNN** – Uses PyTorch; follow instructions in `ProteinMPNN/README.md` for model weights.
  - **FlowPacker** – This repo expects a separate FlowPacker installation.
  - **PyRosetta** – Required only for `per_residue_energy_pyrosetta.py` (licensed distribution; install from the official PyRosetta source).
  - **Rosetta** (optional) – Only needed if you want to consume existing Rosetta `.out` logs with `get_interface_energy.py`.

A minimal setup might look like:

```bash
python -m venv .venv
source .venv/bin/activate

# Core dependencies
pip install -r requirements.txt
pip install pandas numpy tqdm pyyaml easydict biopython matplotlib seaborn

# Install PyTorch to match your CUDA / CPU setup, e.g.:
pip install torch
```

Then install or point to:

- **ProteinMPNN weights** (see `ProteinMPNN/README.md`).
- **FlowPacker** (cloned elsewhere, see below).
- **PyRosetta** if you want to run the PyRosetta-based scripts.

### Using ProteinMPNN

All standard ProteinMPNN usage is unchanged; examples and command-line flags are documented in `ProteinMPNN/README.md`. At a high level:

```bash
cd ProteinMPNN
python protein_mpnn_run.py \
  --path_to_model_weights /path/to/proteinmpnn_weights \
  --model_name v_48_020 \
  --pdb_path /path/to/complex.pdb \
  --pdb_path_chains A,B \
  --out_folder /path/to/mpnn_outputs \
  --num_seq_per_target 32
```

This produces designed sequences (e.g. in FASTA/JSONL form) that can be fed into the FlowPacker sampling step.

### FlowPacker side-chain sampling (`flowpacker_sampler_pipe.py`)

The main FlowPacker pipeline script is `flowpacker_sampler_pipe.py`. It expects:

- A **FlowPacker installation** and its Python dependencies.
- A **CSV file** describing which sequences to graft onto which PDBs.
- A **config** file for the FlowPacker model (YAML or experiment name, depending on your setup).

#### Environment

Set the `FLOWPACKER_REPO` environment variable so the script can import FlowPacker:

```bash
export FLOWPACKER_REPO=/path/to/flowpacker
```

#### CSV format

`flowpacker_sampler_pipe.py` uses a CSV with at least the following columns:

- `link_name` – File name of the template PDB in the input directory (e.g. `4dn4_ccl2_8_4.pdb`).
- `seq` – Designed amino-acid sequence(s). For multichain designs, you can provide:
  - Slash-separated segments (e.g. `HEAVY/LIGHT`);
  - Concatenated sequences that are split by chain length; or
  - Chain-labeled forms (e.g. `A:HEAVY/C:LIGHT`).
- `seq_idx` – Index used to disambiguate multiple sequences per backbone; appended to output PDB names.

The script will:

- Parse the CSV.
- For each row, read the corresponding PDB from `before_pdb_dir` (derived from the FlowPacker config).
- Replace residue names on the specified binder chain(s) with the one-letter sequences.
- Write new PDBs with names like `link_name_seq_idx.pdb` into an `after_pdbs_batch` directory.

#### Running the sampler

Basic usage:

```bash
python flowpacker_sampler_pipe.py \
  path/to/config.yaml \
  --save_dir path/to/output_dir \
  --csv_file path/to/designs.csv \
  --binder_chain A
```

Key arguments:

- `config` (positional) – FlowPacker config file or experiment identifier.
- `--save_dir` – Where to write sampled PDBs and metrics.
- `--csv_file` – CSV of designs (see above).
- `--binder_chain` – Comma-separated list of binder chain IDs (default `"A"`; also configurable via `FLOWPACKER_BINDER_CHAIN`).
- `--save_traj` – If `True`, also save full sampling trajectories.
- `--use_gt_masks` – Use ground-truth chi and atom masks from the dataset.

Outputs include:

- Sampled PDBs under `<save_dir>/run_1/`.
- A `output_dict.pth` file with per-sample metrics (chi MAE/accuracy, RMSD, etc.).
- A copy of the processed PDBs under an `after_pdbs` directory for downstream steps.

### Alternative FlowPacker script (`sampler_pdb_colab.py`)

`sampler_pdb_colab.py` is an earlier or more notebook-oriented version of the FlowPacker sampler. It supports two modes:

- **CSV mode** (similar to `flowpacker_sampler_pipe.py`): grafts sequences from a CSV onto PDBs.
- **Directory mode**: if no CSV is provided, it copies all `.pdb` files from `--pdb_dir` into an `after_pdbs_batch` directory with normalized names.

Example:

```bash
python sampler_pdb_colab.py \
  path/to/config.yaml \
  --pdb_dir path/to/input_pdbs \
  --save_dir path/to/output_dir
```

The script writes samples into `save_dir`, plus a flattened `after_pdbs` directory with processed PDBs.

### Interface energy analysis with Rosetta logs (`get_interface_energy.py`)

`get_interface_energy.py` consumes Rosetta `.out` files and produces per-residue interface energy summaries and heatmaps.

Inputs:

- Either:
  - `--input_csv` – CSV with a `pdbpath` column pointing to structures; or
  - `--input_pdbdir` – Directory containing `.pdb` files.
- `--rosetta_dir` – Directory that contains one subdirectory per PDB, each with `out/<pdbname>.out` Rosetta log files.
- `--binder_id` / `--target_id` – Chain IDs for binder and target.
- `--output_dir` – Directory where results (CSV + plots) will be written.
- `--interface_dist` – Cβ/Cα distance cutoff (Å) to define interface pairs (default ~12 Å).

Example:

```bash
python get_interface_energy.py \
  --input_pdbdir path/to/pdbs \
  --rosetta_dir path/to/rosetta_results \
  --binder_id L \
  --target_id R \
  --output_dir path/to/interface_results \
  --interface_dist 12.0
```

The script will:

- Parse `ResResE` lines in Rosetta logs into a per-pair dataframe.
- Identify interface residue pairs based on the distance cutoff.
- Sum interface energies per binder residue.
- Write `residue_energy.csv` and a KDE plot of binder residue energy distribution.

### Interface energy analysis with PyRosetta (`per_residue_energy_pyrosetta.py`)

`per_residue_energy_pyrosetta.py` replicates the interface energy computation using **only PyRosetta**, without relying on external Rosetta `.out` logs.

High-level behavior:

- Initializes PyRosetta with `init_pyrosetta`.
- Applies a FastRelax-like protocol focused on residues within 20 Å of the binder/target chains.
- Uses the scorefunction to compute per-pair energies for residues in the interface.
- Filters to favorable (negative) pair energies and aggregates per-residue contributions.

Usage:

```bash
python per_residue_energy_pyrosetta.py \
  --pdb path/to/complex.pdb \
  --binder_id L \
  --target_id R \
  --interface_dist 12.0 \
  --output_dir path/to/pyrosetta_interface \
  --extra_flags "-database /path/to/rosetta_database"
```

Outputs:

- `<pdbname>_pair_energies_pyrosetta.csv` – Per-pair interface energies (including per-term columns).
- `residue_energy_pyrosetta.csv` – Per-binder-residue summed interface energies.
- `<pdbname>_interface_binder_residues_score_heatmap.png` – Heatmap of interface energies.

### Example end-to-end workflow

An example affinity-maturation workflow using this repository might be:

1. **Design sequences with ProteinMPNN** on a complex backbone (`ProteinMPNN/protein_mpnn_run.py`), writing sequences to a CSV compatible with `flowpacker_sampler_pipe.py`.
2. **Sample side chains with FlowPacker** using `flowpacker_sampler_pipe.py`, generating all-atom models for each design.
3. **Compute interface energies** for each sampled structure:
   - Either from Rosetta `.out` logs with `get_interface_energy.py`, or
   - Directly from PDBs with `per_residue_energy_pyrosetta.py`.
4. **Rank designs** by per-residue or total interface energy and inspect top candidates.

### Citing ProteinMPNN

If you use this repository in academic work, please cite the original ProteinMPNN paper (see also `ProteinMPNN/README.md`):

> Dauparas, J. *et al.* Robust deep learning–based protein sequence design using ProteinMPNN. *Science* **378**, 49–56 (2022).

