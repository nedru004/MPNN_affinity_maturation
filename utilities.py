import pandas as pd
import os
import ast
from tqdm import tqdm


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