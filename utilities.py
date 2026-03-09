import pandas as pd
import os
import ast
from tqdm import tqdm
import csv
from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import json
from pathlib import Path
import shutil
from typing import Optional, Dict, List


def _build_pdb_index(pdb_root: Path) -> Dict[str, List[Path]]:
    """
    Recursively index all PDB files under `pdb_root` by stem.
    """
    index: Dict[str, List[Path]] = {}
    for pdb_file in pdb_root.rglob("*.pdb"):
        index.setdefault(pdb_file.stem, []).append(pdb_file)
    return index


def _build_file_index(root: Path, patterns: List[str]) -> Dict[str, List[Path]]:
    """
    Recursively index files under `root` by stem for multiple glob patterns.
    """
    index: Dict[str, List[Path]] = {}
    for pattern in patterns:
        for file_path in root.rglob(pattern):
            index.setdefault(file_path.stem, []).append(file_path)
    return index


def _candidate_names_from_path(json_file: Path, extra_candidates: Optional[List[str]] = None) -> List[str]:
    """
    Derive candidate names for matching a JSON file to an original PDB
    based on its filename and containing folders.
    """
    candidates: List[str] = []

    # Start with the JSON stem
    candidates.append(json_file.stem)

    # Add up to a few parent folder names (closest first)
    parent = json_file.parent
    for _ in range(3):
        if parent is None or parent == parent.parent:
            break
        candidates.append(parent.name)
        parent = parent.parent

    # Add any extra supplied candidates (e.g., stripped prefixes)
    if extra_candidates:
        candidates.extend(extra_candidates)

    # Deduplicate while preserving order
    seen = set()
    ordered: List[str] = []
    for name in candidates:
        if name not in seen:
            seen.add(name)
            ordered.append(name)

    return ordered


def _find_matching_pdb(
    json_file: Path,
    pdb_index: Dict[str, List[Path]],
    extra_candidates: Optional[List[str]] = None,
) -> Optional[Path]:
    """
    Try to find a PDB in `pdb_index` that matches this JSON file.

    Matching strategy:
    - First try direct stem matches against candidate names.
    - Then fall back to prefix/substring heuristics.
    """
    if not pdb_index:
        return None

    candidates = _candidate_names_from_path(json_file, extra_candidates)
    stems = list(pdb_index.keys())

    # Direct equality on stem
    for name in candidates:
        if name in pdb_index:
            return pdb_index[name][0]

    # Heuristics: prefix / substring matches
    for name in candidates:
        for stem in stems:
            if stem.startswith(name) or name.startswith(stem) or name in stem:
                return pdb_index[stem][0]

    return None


def _find_matching_file(
    source_path: Path,
    file_index: Dict[str, List[Path]],
    extra_candidates: Optional[List[str]] = None,
) -> Optional[Path]:
    """
    Find a matching file (PDB/CIF/etc.) from an indexed file set.
    """
    if not file_index:
        return None

    candidates = _candidate_names_from_path(source_path, extra_candidates)
    stems = list(file_index.keys())

    for name in candidates:
        if name in file_index:
            return file_index[name][0]

    for name in candidates:
        for stem in stems:
            if stem.startswith(name) or name.startswith(stem) or name in stem:
                return file_index[stem][0]

    return None


def _extract_ca_atoms(structure) -> Dict[tuple, "Atom"]:
    """
    Collect CA atoms keyed by (chain_id, residue_number) from the first model.
    """
    atoms: Dict[tuple, "Atom"] = {}
    # Use the first model only for consistency
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    # residue.id is a tuple like (' ', resseq, icode)
                    resseq = residue.id[1]
                    key = (chain.id, resseq)
                    atoms[key] = residue["CA"]
        break
    return atoms


def _compute_aligned_rmsd(original_pdb: Path, predicted_pdb: Path) -> Optional[float]:
    """
    Compute CA-aligned RMSD between original and predicted PDB structures.

    Returns None if RMSD cannot be computed (e.g., too few matching residues).
    """
    try:
        if original_pdb.suffix.lower() == ".cif":
            orig_parser = MMCIFParser(QUIET=True)
        else:
            orig_parser = PDBParser(QUIET=True)

        if predicted_pdb.suffix.lower() == ".cif":
            pred_parser = MMCIFParser(QUIET=True)
        else:
            pred_parser = PDBParser(QUIET=True)

        orig_struct = orig_parser.get_structure("orig", str(original_pdb))
        pred_struct = pred_parser.get_structure("pred", str(predicted_pdb))
    except Exception as e:
        print(f"Error parsing PDBs for RMSD ({original_pdb}, {predicted_pdb}): {e}")
        return None

    orig_atoms = _extract_ca_atoms(orig_struct)
    pred_atoms = _extract_ca_atoms(pred_struct)

    common_keys = sorted(set(orig_atoms.keys()) & set(pred_atoms.keys()))
    if len(common_keys) < 3:
        # Need at least 3 points to define a stable superposition
        return None

    fixed_atoms = [orig_atoms[k] for k in common_keys]
    moving_atoms = [pred_atoms[k] for k in common_keys]

    sup = Superimposer()
    try:
        sup.set_atoms(fixed_atoms, moving_atoms)
    except Exception as e:
        print(f"Error during RMSD superposition ({original_pdb}, {predicted_pdb}): {e}")
        return None

    return float(sup.rms)

# merge key residues into one pdb file
def process_pdb_mutation_and_renumber(csv, pdb_output_dir, 
                                      renumber_chain='B',
                                      start_index=1):
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    
    df = pd.read_csv(csv)
    os.makedirs(pdb_output_dir, exist_ok=True)

    print(f"Start processing: Renaming Chain A & Renumbering Chain {renumber_chain} (start from {start_index})...")

    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        
        if not os.path.exists(row["pdb_path"]):
            print(f"Warning: {row['pdb_path']} not found, skipping...")
            continue

        resi_to_resname = {}
        if "key_res" in row and pd.notna(row["key_res"]):
            try:
                for resi, aa in ast.literal_eval(row["key_res"]).items():
                    resi_to_resname[int(resi)] = one_to_three.get(aa[1], "UNK")
            except Exception as e:
                print(f"Error parsing key_res for {row['pdb_name']}: {e}")

        with open(row["pdb_path"], "r") as f:
            pdb_lines = f.readlines()

        new_lines = []
        
        current_renumber_idx = start_index - 1
        last_seen_resi_id = None

        for line in pdb_lines:
            if not line.startswith("ATOM"):
                new_lines.append(line)
                continue

            chain_id = line[21]
            

            if chain_id == "A":
                resi = int(line[22:26])
                if resi in resi_to_resname:
                    new_resname = resi_to_resname[resi]

                    line = line[:17] + new_resname.ljust(3) + line[20:]

            elif chain_id == renumber_chain:
                current_resi_id_str = line[22:27]                
                if current_resi_id_str != last_seen_resi_id:
                    current_renumber_idx += 1
                    last_seen_resi_id = current_resi_id_str

                new_resi_str = f"{current_renumber_idx:>4}"
                line = line[:22] + new_resi_str + " " + line[27:]
            new_lines.append(line)
        output_pdb = os.path.join(pdb_output_dir, os.path.basename(row["pdb_name"])+".pdb")
        with open(output_pdb, "w") as f:
            f.writelines(new_lines)
    print("All done!")


def direct_fasta_to_csv(input_dirs: list, output_csv: str, suffix: str = ".pdb"):

    seen_seqs = set()

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)

    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["link_name", "seq", "seq_idx"])
        for folder in input_dirs:
            if not os.path.exists(folder): continue
            files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(('.fasta', '.fa'))]

            for file_path in files:
                base_name = os.path.splitext(os.path.basename(file_path))[0]

                for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
                    if i == 0:
                        continue

                    seq_str = str(record.seq)
                    if seq_str in seen_seqs:
                        continue
                    seen_seqs.add(seq_str)
                    link_name = f"{base_name}{suffix}"
                    seq_idx = str(i)

                    writer.writerow([link_name, seq_str, seq_idx])

    print(f"✅ Processing complete! {len(seen_seqs)} unique sequences have been written to: {output_csv}")

def filter_boltz_scores(
    input_dir: str,
    output_csv: str,
    output_dir: str,
    original_pdb_dir: Optional[str] = None,
    threshold_pTM: float = 0.8,
    threshold_ipTM: float = 0.8,
    threshold_rmsd: float = 2.0,
):
    """
    Extracts scores from Boltz JSON files in the given directory (including subdirectories)
    and outputs them to a CSV file.
    """
    score_keys = [
        "confidence_score",
        "ptm",
        "iptm",
        "ligand_iptm",
        "protein_iptm",
        "complex_plddt",
        "complex_iplddt",
        "complex_pde",
        "complex_ipde"
    ]
    
    data_list = []
    input_path = Path(input_dir)
    json_files = list(input_path.rglob("*.json"))
    
    # Ensure output directory for selected PDBs exists
    os.makedirs(output_dir, exist_ok=True)

    # If an original PDB root is provided, index its contents once
    pdb_root: Optional[Path] = None
    pdb_index: Dict[str, List[Path]] = {}
    if original_pdb_dir is not None:
        pdb_root = Path(original_pdb_dir)
        if not pdb_root.exists():
            print(f"Warning: original PDB directory {original_pdb_dir} does not exist; falling back to local JSON directory search.")
            pdb_root = None
        else:
            pdb_index = _build_pdb_index(pdb_root)

    # Index predicted structures under input_dir (Boltz outputs may be nested and/or CIF)
    predicted_index = _build_file_index(input_path, ["*.pdb", "*.cif"])

    for json_file in tqdm(json_files, desc="Parsing JSON files"):
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                
            if "confidence_score" in data:
                row = {"filename": json_file.name, "filepath": str(json_file)}
                for key in score_keys:
                    row[key] = data.get(key, None)
                
                # Also extract chain specific scores if present (optional but could be useful)
                if "chains_ptm" in data:
                    for chain_idx, val in data["chains_ptm"].items():
                        row[f"chain_{chain_idx}_ptm"] = val
                # Determine original and predicted PDB paths for RMSD computation and copying
                original_pdb_path: Optional[Path] = None
                predicted_pdb_path: Optional[Path] = None

                # Locate original PDB in the provided root (if any)
                if pdb_root is not None:
                    extra_candidates: List[str] = []
                    json_stem = json_file.stem
                    if json_stem.startswith("confidence_"):
                        extra_candidates.append(json_stem[len("confidence_"):])

                    original_pdb_path = _find_matching_pdb(json_file, pdb_index, extra_candidates)

                # Locate predicted structure from indexed outputs
                json_stem = json_file.stem
                extra_candidates: List[str] = []
                if json_stem.startswith("confidence_"):
                    extra_candidates.append(json_stem[len("confidence_"):])
                predicted_pdb_path = _find_matching_file(json_file, predicted_index, extra_candidates)

                # Compute RMSD if we have both original and predicted structures
                rmsd_val: Optional[float] = None
                if original_pdb_path is not None and predicted_pdb_path is not None:
                    rmsd_val = _compute_aligned_rmsd(original_pdb_path, predicted_pdb_path)
                row["rmsd"] = rmsd_val
                row["original_pdb_path"] = str(original_pdb_path) if original_pdb_path is not None else None
                row["predicted_structure_path"] = str(predicted_pdb_path) if predicted_pdb_path is not None else None

                # Record scores and determine if this prediction passes thresholds
                ptm_val = data.get("ptm")
                iptm_val = data.get("iptm")
                passes_ptm_ipTM = (
                    ptm_val is not None
                    and iptm_val is not None
                    and ptm_val >= threshold_pTM
                    and iptm_val >= threshold_ipTM
                )

                passes_rmsd = True
                # Apply RMSD filter only when we have an original PDB root and a positive threshold
                if pdb_root is not None and threshold_rmsd is not None and threshold_rmsd > 0:
                    passes_rmsd = rmsd_val is not None and rmsd_val <= threshold_rmsd

                row["passes_threshold"] = passes_ptm_ipTM and passes_rmsd
                data_list.append(row)

                # If this prediction passes all thresholds, copy the appropriate PDB
                if passes_ptm_ipTM and passes_rmsd:
                    pdb_path_to_copy: Optional[Path] = None
                    dest_path: Optional[Path] = None

                    if pdb_root is not None:
                        # Copy the original PDB, preserving folder structure under output_dir
                        pdb_path_to_copy = original_pdb_path
                        if pdb_path_to_copy is not None and pdb_path_to_copy.is_file():
                            try:
                                rel_path = pdb_path_to_copy.relative_to(pdb_root)
                            except ValueError:
                                # If for some reason it's not under pdb_root, just flatten
                                rel_path = pdb_path_to_copy.name
                            dest_path = Path(output_dir) / rel_path
                    else:
                        # No original root provided: copy the predicted PDB itself
                        pdb_path_to_copy = predicted_pdb_path
                        if pdb_path_to_copy is not None and pdb_path_to_copy.is_file():
                            dest_path = Path(output_dir) / pdb_path_to_copy.name

                    if pdb_path_to_copy is not None and dest_path is not None:
                        dest_path.parent.mkdir(parents=True, exist_ok=True)
                        try:
                            shutil.copy2(pdb_path_to_copy, dest_path)
                        except Exception as copy_err:
                            print(f"Error copying PDB for {json_file}: {copy_err}")
                    else:
                        print(f"Warning: No PDB file found for {json_file}")
                        
        except Exception as e:
            print(f"Error parsing {json_file}: {e}")
            
    if not data_list:
        print(f"No valid Boltz score JSON files found in {input_dir}")
        return
        
    df = pd.DataFrame(data_list)
    
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv, index=False)
    n_passed = df["passes_threshold"].sum()
    n_rmsd = df["rmsd"].notna().sum() if "rmsd" in df.columns else 0
    print(f"✅ Extracted scores from {len(df)} files to {output_csv}")
    print(f"   {int(n_passed)} PDB(s) passed thresholds (pTM ≥ {threshold_pTM}, ipTM ≥ {threshold_ipTM}, RMSD ≤ {threshold_rmsd})")
    print(f"   RMSD computed for {int(n_rmsd)}/{len(df)} entries")


def filter_protenix_scores(
    input_dir: str,
    output_csv: str,
    output_dir: str,
    original_pdb_dir: Optional[str] = None,
    threshold_pTM: float = 0.8,
    threshold_ipTM: float = 0.8,
):
    """
    Recursively finds Protenix `summary_confidence.json` files under `input_dir`,
    extracts key scores, and writes them to `output_csv`.

    Adds a `passes_threshold` column indicating whether
    ptm >= threshold_pTM and iptm >= threshold_ipTM.
    """
    simple_keys = [
        "plddt",
        "gpde",
        "ptm",
        "iptm",
        "ranking_score",
        "disorder",
        "num_recycles",
        "ipsae_pae_cutoff",
        "ipsae_max",
        "ipsae_interface_max",
        "has_clash",
    ]
    chain_keys = [
        "chain_gpde",
        "chain_ptm",
        "chain_iptm",
        "chain_plddt",
    ]

    data_list = []
    input_path = Path(input_dir)
    json_files = list(input_path.rglob("summary_confidence.json"))

    # Ensure output directory for selected PDBs exists
    os.makedirs(output_dir, exist_ok=True)

    # If an original PDB root is provided, index its contents once
    pdb_root: Optional[Path] = None
    pdb_index: Dict[str, List[Path]] = {}
    if original_pdb_dir is not None:
        pdb_root = Path(original_pdb_dir)
        if not pdb_root.exists():
            print(f"Warning: original PDB directory {original_pdb_dir} does not exist; PDBs will not be copied.")
            pdb_root = None
        else:
            pdb_index = _build_pdb_index(pdb_root)

    for json_file in tqdm(json_files, desc="Parsing Protenix summary_confidence.json files"):
        try:
            with open(json_file, "r", encoding="utf-8") as f:
                data = json.load(f)

            row = {"filename": json_file.name, "filepath": str(json_file)}

            # Top-level scalar metrics
            for key in simple_keys:
                row[key] = data.get(key, None)

            # Per-chain metrics (flatten lists)
            for key in chain_keys:
                values = data.get(key)
                if isinstance(values, list):
                    for idx, val in enumerate(values):
                        row[f"{key}_{idx}"] = val

            # Threshold flag based on ptm / iptm
            ptm_val = data.get("ptm")
            iptm_val = data.get("iptm")
            passes = (
                ptm_val is not None
                and iptm_val is not None
                and ptm_val >= threshold_pTM
                and iptm_val >= threshold_ipTM
            )
            row["passes_threshold"] = passes

            # If this prediction passes the thresholds, try to copy the matching PDB
            if passes and pdb_root is not None:
                pdb_path = _find_matching_pdb(json_file, pdb_index)

                if pdb_path is not None and pdb_path.is_file():
                    try:
                        rel_path = pdb_path.relative_to(pdb_root)
                    except ValueError:
                        rel_path = pdb_path.name

                    dest_path = Path(output_dir) / rel_path
                    dest_path.parent.mkdir(parents=True, exist_ok=True)

                    try:
                        shutil.copy2(pdb_path, dest_path)
                    except Exception as copy_err:
                        print(f"Error copying PDB for {json_file}: {copy_err}")
                else:
                    print(f"Warning: No original PDB file found for {json_file}")

            data_list.append(row)

        except Exception as e:
            print(f"Error parsing {json_file}: {e}")

    if not data_list:
        print(f"No Protenix summary_confidence.json files found in {input_dir}")
        return

    df = pd.DataFrame(data_list)
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv, index=False)
    n_passed = df["passes_threshold"].sum()
    print(f"✅ Extracted Protenix scores from {len(df)} files to {output_csv}")
    print(f"   {int(n_passed)} PDB(s) passed thresholds (pTM ≥ {threshold_pTM}, ipTM ≥ {threshold_ipTM})")

def generate_boltz_yamls_from_pdbs(input_dir: str, output_dir: str = None, use_template: bool = False):
    """
    Reads all PDB files from a directory, extracts the sequence for each chain,
    and outputs a .yaml file formatted for Boltz prediction for each PDB.
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    input_path = Path(input_dir)
    if output_dir is None:
        output_dir = input_dir
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_files = list(input_path.glob("*.pdb"))
    
    for pdb_file in tqdm(pdb_files, desc="Generating Boltz YAMLs"):
        chains = {}
        last_resi_id = None
        kept_residues = set()        # (chain_id, resi_id) pairs that map to standard AAs
        cleaned_atom_lines = []      # lines to write into the cleaned template PDB
        
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    chain_id = line[21]
                    resi_id = line[22:27].strip()  # include insertion code if present
                    altloc = line[16]
                    
                    # keep only first altloc (blank, 'A', or '1')
                    if altloc not in [' ', 'A', '1']:
                        continue
                    
                    if chain_id not in chains:
                        chains[chain_id] = []
                    
                    current_res = (chain_id, resi_id)
                    aa = three_to_one.get(res_name, 'X')

                    # avoid sequence duplication if atom for same residue appears consecutively
                    if last_resi_id != current_res:
                        # Only include residues that map to a standard amino acid
                        if aa != 'X':
                            chains[chain_id].append(aa)
                            kept_residues.add(current_res)
                        last_resi_id = current_res

                    # If this residue is part of the cleaned sequence, keep its atoms
                    if current_res in kept_residues:
                        cleaned_atom_lines.append(line)
        
        if not chains:
            print(f"Warning: No ATOM lines found in {pdb_file}")
            continue
            
        yaml_lines = ["version: 1\nsequences:\n"]
        # Ensure deterministic ordering of chains
        for chain_id in sorted(chains.keys()):
            seq_list = chains[chain_id]
            if not seq_list:
                continue
            seq_str = "".join(seq_list)
            yaml_lines.append("  - protein:\n")
            yaml_lines.append(f"      id: {chain_id}\n")
            yaml_lines.append(f"      sequence: {seq_str}\n")
            if use_template:
                # when using templates, request an empty MSA so Boltz
                # uses only the provided template information
                yaml_lines.append("      msa: empty\n")

        cleaned_pdb_path = None
        if use_template:
            # Write a cleaned template PDB containing only the residues
            # that appear in the sequences above and only protein ATOM records.
            cleaned_pdb_path = Path(output_dir) / f"{pdb_file.stem}_template.pdb"
            with open(cleaned_pdb_path, "w") as tpl_f:
                for line in cleaned_atom_lines:
                    tpl_f.write(line)
                tpl_f.write("END\n")

            # Add a single templates block that references the cleaned PDB.
            # Boltz will assign template chain IDs incrementally (A1, A2, B1, etc)
            # based on the chain names present in the PDB file.
            chain_ids = sorted(chains.keys())
            chain_list_str = ", ".join(chain_ids)
            yaml_lines.append("templates:\n")
            # Use full path so Boltz can find the template regardless of CWD
            yaml_lines.append(f"  - pdb: {cleaned_pdb_path.resolve()}\n")
            yaml_lines.append(f"    chain_id: [{chain_list_str}]\n")
            yaml_lines.append(f"    template_id: [{chain_list_str}]\n")

        output_yaml = Path(output_dir) / f"{pdb_file.stem}.yaml"
        with open(output_yaml, "w") as f:
            f.writelines(yaml_lines)
            
    print(f"✅ Generated {len(pdb_files)} Boltz YAML files in {output_dir}")