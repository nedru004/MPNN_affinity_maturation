import pandas as pd
import os
import ast
from tqdm import tqdm
import csv
from Bio import SeqIO
import json
from pathlib import Path

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

def extract_boltz_scores_to_csv(input_dir: str, output_csv: str):
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
                        
                data_list.append(row)
        except Exception as e:
            print(f"Error parsing {json_file}: {e}")
            
    if not data_list:
        print(f"No valid Boltz score JSON files found in {input_dir}")
        return
        
    df = pd.DataFrame(data_list)
    
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"✅ Extracted scores from {len(df)} files to {output_csv}")

def generate_boltz_yamls_from_pdbs(input_dir: str, output_dir: str = None):
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
        
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    chain_id = line[21]
                    resi_id = line[22:27].strip() # include insertion code if present
                    altloc = line[16]
                    
                    # keep only first altloc (blank, 'A', or '1')
                    if altloc not in [' ', 'A', '1']:
                        continue
                        
                    if chain_id not in chains:
                        chains[chain_id] = []
                        
                    current_res = (chain_id, resi_id)
                    # avoid sequence duplication if atom for same residue appears consecutively
                    if last_resi_id != current_res:
                        chains[chain_id].append(three_to_one.get(res_name, 'X'))
                        last_resi_id = current_res
                        
        if not chains:
            print(f"Warning: No ATOM lines found in {pdb_file}")
            continue
            
        yaml_lines = ["version: 1\nsequences:\n"]
        for chain_id, seq_list in chains.items():
            # remove unrecognised amino acids ('X') that might act as weird gaps
            seq_str = "".join(seq_list).replace('X', '')
            if seq_str:
                yaml_lines.append("  - protein:\n")
                yaml_lines.append(f"      id: {chain_id}\n")
                yaml_lines.append(f"      sequence: {seq_str}\n")
                
        output_yaml = Path(output_dir) / f"{pdb_file.stem}.yaml"
        with open(output_yaml, "w") as f:
            f.writelines(yaml_lines)
            
    print(f"✅ Generated {len(pdb_files)} Boltz YAML files in {output_dir}")