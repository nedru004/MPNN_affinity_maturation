import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List

import pandas as pd

from utilities import direct_fasta_to_csv, process_pdb_mutation_and_renumber


def run_mpnn_sequence_design(
    working_dir: Path,
    input_pdb_dir: Path,
    chains_to_design: str = "A",
    num_seq_per_target: int = 8,
    sampling_temp: float = 0.5,
    seed: int = 37,
    fixed_positions: str = "11 12 13 14 15",
    mpnn_repo: Path | None = None,
    bias_AA_file: Path | None = None,
) -> None:
    """
    Step 1: Run ProteinMPNN sequence design (with and without AA bias).
    """
    if mpnn_repo is None:
        mpnn_repo = Path("ProteinMPNN")
    mpnn_script_path = mpnn_repo / "protein_mpnn_run.py"

    test_dir = working_dir / "test"
    test_dir.mkdir(parents=True, exist_ok=True)

    folder_with_pdbs = input_pdb_dir
    path_for_parsed_chains = test_dir / "parsed_chains.jsonl"
    path_for_assigned_chains = test_dir / "assigned_chains.jsonl"
    path_for_fixed_positions = test_dir / "fixed_positions.jsonl"
    path_for_fixed_positions_csv = test_dir / "fixed_positions.csv"
    out_dir = test_dir / "mpnn_out"
    out_dir_bias = test_dir / "mpnn_out_bias"

    if bias_AA_file is None:
        bias_AA_file = working_dir / "bias_AA.jsonl"

    # write fixed positions CSV used by helper script
    path_for_fixed_positions_csv.write_text(fixed_positions + "\n")

    # Helper scripts shipped with ProteinMPNN (assumed present)
    subprocess.run(
        [
            sys.executable,
            str(mpnn_repo / "helper_scripts" / "parse_multiple_chains.py"),
            "--input_path",
            str(folder_with_pdbs),
            "--output_path",
            str(path_for_parsed_chains),
        ],
        check=True,
    )

    subprocess.run(
        [
            sys.executable,
            str(mpnn_repo / "helper_scripts" / "assign_fixed_chains.py"),
            "--input_path",
            str(path_for_parsed_chains),
            "--output_path",
            str(path_for_assigned_chains),
            "--chain_list",
            chains_to_design,
        ],
        check=True,
    )

    subprocess.run(
        [
            sys.executable,
            str(mpnn_repo / "helper_scripts" / "make_fixed_positions_dict.py"),
            "--input_path",
            str(path_for_parsed_chains),
            "--output_path",
            str(path_for_fixed_positions),
            "--chain_list",
            chains_to_design,
            "--position_list",
            fixed_positions,
        ],
        check=True,
    )

    # Main ProteinMPNN run (unbiased)
    args_first = [
        sys.executable,
        str(mpnn_script_path),
        "--jsonl_path",
        str(path_for_parsed_chains),
        "--chain_id_jsonl",
        str(path_for_assigned_chains),
        "--fixed_positions_jsonl",
        str(path_for_fixed_positions),
        "--out_folder",
        str(out_dir),
        "--num_seq_per_target",
        str(num_seq_per_target),
        "--sampling_temp",
        str(sampling_temp),
        "--seed",
        str(seed),
        "--batch_size",
        str(num_seq_per_target),
        "--omit_AAs",
        "C",
    ]
    subprocess.run(args_first, check=True)

    # ProteinMPNN run with AA bias
    args_second = [
        sys.executable,
        str(mpnn_script_path),
        "--jsonl_path",
        str(path_for_parsed_chains),
        "--chain_id_jsonl",
        str(path_for_assigned_chains),
        "--fixed_positions_jsonl",
        str(path_for_fixed_positions),
        "--out_folder",
        str(out_dir_bias),
        "--num_seq_per_target",
        str(num_seq_per_target),
        "--sampling_temp",
        str(sampling_temp),
        "--seed",
        str(seed),
        "--batch_size",
        str(num_seq_per_target),
        "--omit_AAs",
        "C",
        "--bias_AA_jsonl",
        str(bias_AA_file),
    ]
    subprocess.run(args_second, check=True)


def run_flowpacker(
    working_dir: Path,
    input_pdb_dir: Path,
    use_gt_masks: bool,
) -> None:
    """
    Step 2: Run flowpacker to generate backbone designs from MPNN sequences.
    """
    test_dir = working_dir / "test"
    mpnn_seq_dirs: List[Path] = [
        test_dir / "mpnn_out" / "seqs",
        test_dir / "mpnn_out_bias" / "seqs",
    ]
    mpnn_final_csv = test_dir / "mpnn_final_result.csv"

    # Convert MPNN FASTA outputs into CSV required by sampler
    direct_fasta_to_csv(
        input_dirs=[str(p) for p in mpnn_seq_dirs],
        output_csv=str(mpnn_final_csv),
        suffix=".pdb",
    )

    # Patch flowpacker loader config to point at this repo's base.yaml
    loader_path = Path("/content/flowpacker/utils/loader.py")
    if loader_path.exists():
        line_number_to_edit = 45
        new_content = f"    config_dir = '{working_dir}/base.yaml'"
        with loader_path.open("r") as f:
            lines = f.readlines()
        if 0 < line_number_to_edit <= len(lines):
            lines[line_number_to_edit - 1] = new_content + "\n"
        with loader_path.open("w") as f:
            f.writelines(lines)

    flow_out_dir = test_dir / "flowpacker_out"

    subprocess.run(
        [
            sys.executable,
            str(working_dir / "sampler_pdb_colab.py"),
            "base",
            "--use_gt_masks",
            str(use_gt_masks),
            "--pdb_dir",
            str(input_pdb_dir),
            "--save_dir",
            str(flow_out_dir),
            "--csv_file",
            str(mpnn_final_csv),
        ],
        check=True,
    )


def run_structure_scoring_and_filtering(
    working_dir: Path,
    run_protenix: bool,
    run_boltz: bool,
    threshold_pTM: float = 0.8,
    threshold_ipTM: float = 0.8,
    rmsd_threshold: float = 3.0,
) -> None:
    """
    Step 3: Run structure scoring (Protenix and/or Boltz) and filter PDBs.
    """
    from utilities import filter_protenix_scores, filter_boltz_scores, generate_boltz_yamls_from_pdbs

    test_dir = working_dir / "test"
    flow_after_pdbs = test_dir / "flowpacker_out" / "after_pdbs"
    filtered_pdb_dir = test_dir / "filtered_pdb"
    structure_score_csv = test_dir / "structure_score_summary.csv"

    if run_protenix:
        score_py = Path("/content/ProtenixScore/score.py")
        cli_py = Path("/content/ProtenixScore/cli.py")

        if score_py.exists():
            with score_py.open("r") as f:
                lines = f.readlines()
            lines[19 - 1] = 'import sys\nsys.path.insert(0, "/content/ProtenixScore/Protenix_fork")\n'
            lines[1071 - 1] = (
                "    from ProtenixScore.msa_colabfold import ColabFoldMSAConfig, ensure_msa_dir\n"
            )
            with score_py.open("w") as f:
                f.writelines(lines)

        if cli_py.exists():
            with cli_py.open("r") as f:
                lines = f.readlines()
            lines[5 - 1] = "from ProtenixScore.score import run_score\n"
            with cli_py.open("w") as f:
                f.writelines(lines)

        subprocess.run(
            [
                sys.executable,
                "ProtenixScore/cli.py",
                "score",
                "--input",
                str(flow_after_pdbs),
                "--output",
                str(test_dir / "protenix_scores"),
                "--recursive",
            ],
            check=True,
        )

        filter_protenix_scores(
            input_dir=str(test_dir / "protenix_scores"),
            output_csv=str(structure_score_csv),
            output_dir=str(filtered_pdb_dir),
            original_pdb_dir=str(flow_after_pdbs),
            threshold_pTM=threshold_pTM,
            threshold_ipTM=threshold_ipTM,
        )

    if run_boltz:
        boltz_input_dir = test_dir / "boltz_input"
        boltz_out_dir = test_dir / "boltz_out"

        generate_boltz_yamls_from_pdbs(str(flow_after_pdbs), str(boltz_input_dir), use_template=False)

        # Run Boltz prediction over all YAMLs
        for yaml_path in boltz_input_dir.glob("*.yaml"):
            subprocess.run(
                [
                    "bash",
                    "-lc",
                    f"source /content/boltz_env/bin/activate; boltz predict {yaml_path} --out_dir {boltz_out_dir} --use_msa_server",
                ],
                check=True,
            )

        filter_boltz_scores(
            input_dir=str(boltz_out_dir),
            output_csv=str(structure_score_csv),
            output_dir=str(filtered_pdb_dir),
            original_pdb_dir=str(flow_after_pdbs),
            threshold_pTM=threshold_pTM,
            threshold_ipTM=threshold_ipTM,
            threshold_rmsd=rmsd_threshold,
        )


def run_pyrosetta_interface_energy(
    working_dir: Path,
    binder_chain: str = "A",
    target_chain: str = "M",
    interface_dist: float = 10.0,
    energy_threshold: float = -5.0,
) -> Path:
    """
    Step 4: Run PyRosetta interface energy on filtered PDBs and aggregate
    key residues into a single CSV used for motif merging.
    """
    test_dir = working_dir / "test"
    filtered_pdb_dir = test_dir / "filtered_pdb"
    xml_protocol = working_dir / "update.xml"

    filtered_pdbs = sorted(filtered_pdb_dir.glob("*.pdb"))
    if not filtered_pdbs:
        raise FileNotFoundError(f"No PDBs found in {filtered_pdb_dir} for PyRosetta scoring.")

    rows = []
    for pdb_path in filtered_pdbs:
        pdb_name = pdb_path.stem
        out_dir = test_dir / "af3score" / pdb_name
        out_dir.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                sys.executable,
                str(working_dir / "per_residue_energy_pyrosetta.py"),
                "--pdb",
                str(pdb_path),
                "--binder_id",
                binder_chain,
                "--target_id",
                target_chain,
                "--interface_dist",
                str(interface_dist),
                "--output_dir",
                str(out_dir),
                "--xml_protocol",
                str(xml_protocol),
                "--energy_threshold",
                str(energy_threshold),
            ],
            check=True,
        )

        csv_path = out_dir / "fixed_residue_pyrosetta.csv"
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            rows.append(df)

    if not rows:
        raise RuntimeError("No fixed_residue_pyrosetta.csv files were produced by PyRosetta step.")

    merged_df = pd.concat(rows, ignore_index=True)
    merged_csv = test_dir / "fixed_residue_pyrosetta.csv"
    merged_df.to_csv(merged_csv, index=False)
    return merged_csv


def run_merge_motif_pdb(
    fixed_residue_csv: Path,
    working_dir: Path,
    target_chain: str = "M",
    start_index: int = 9,
) -> Path:
    """
    Step 5: Merge key residues into motif PDBs and renumber target chain.
    """
    test_dir = working_dir / "test"
    merge_dir = test_dir / "merge_motif_pdb"
    merge_dir.mkdir(parents=True, exist_ok=True)

    process_pdb_mutation_and_renumber(
        csv=str(fixed_residue_csv),
        pdb_output_dir=str(merge_dir),
        renumber_chain=target_chain,
        start_index=start_index,
    )
    return merge_dir


def run_pipeline(
    working_dir: str | Path,
    input_pdb_dir: str | Path,
    chains_to_design: str = "A",
    num_seq_per_target: int = 8,
    sampling_temp: float = 0.5,
    seed: int = 37,
    fixed_positions: str = "11 12 13 14 15",
    run_protenix: bool = False,
    run_boltz: bool = True,
    binder_chain: str = "A",
    target_chain: str = "M",
    interface_dist: float = 10.0,
    energy_threshold: float = -5.0,
) -> Path:
    """
    Run the full design pipeline end-to-end.

    This mirrors the steps in `test_notebook.ipynb`, and is designed so that
    the final `merge_motif_pdb` directory can be reused as the input_pdb_dir
    for subsequent affinity-maturation iterations.
    """
    working_dir = Path(working_dir).resolve()
    input_pdb_dir = Path(input_pdb_dir).resolve()

    run_mpnn_sequence_design(
        working_dir=working_dir,
        input_pdb_dir=input_pdb_dir,
        chains_to_design=chains_to_design,
        num_seq_per_target=num_seq_per_target,
        sampling_temp=sampling_temp,
        seed=seed,
        fixed_positions=fixed_positions,
    )

    run_flowpacker(
        working_dir=working_dir,
        input_pdb_dir=input_pdb_dir,
        use_gt_masks=True,
    )

    run_structure_scoring_and_filtering(
        working_dir=working_dir,
        run_protenix=run_protenix,
        run_boltz=run_boltz,
    )

    fixed_csv = run_pyrosetta_interface_energy(
        working_dir=working_dir,
        binder_chain=binder_chain,
        target_chain=target_chain,
        interface_dist=interface_dist,
        energy_threshold=energy_threshold,
    )

    merge_dir = run_merge_motif_pdb(
        fixed_residue_csv=fixed_csv,
        working_dir=working_dir,
        target_chain=target_chain,
        start_index=9,
    )

    return merge_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the full MPNN affinity-maturation pipeline.")
    parser.add_argument(
        "--working_dir",
        type=str,
        default=os.getcwd(),
        help="Path to the MPNN_affinity_maturation repo (default: current directory).",
    )
    parser.add_argument(
        "--input_pdb_dir",
        type=str,
        required=True,
        help="Directory containing input PDBs for MPNN / flowpacker.",
    )
    parser.add_argument(
        "--chains_to_design",
        type=str,
        default="A",
        help='Chain IDs to design (ProteinMPNN --chain_list, default: "A").',
    )
    parser.add_argument(
        "--num_seq_per_target",
        type=int,
        default=8,
        help="Number of sequences per design target for ProteinMPNN.",
    )
    parser.add_argument(
        "--sampling_temp",
        type=float,
        default=0.5,
        help="Sampling temperature for ProteinMPNN.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=37,
        help="Random seed for ProteinMPNN.",
    )
    parser.add_argument(
        "--fixed_positions",
        type=str,
        default="11 12 13 14 15",
        help="Space-separated list of residue positions to keep fixed in ProteinMPNN.",
    )
    parser.add_argument(
        "--run_protenix",
        action="store_true",
        help="If set, run Protenix-based scoring/filtering.",
    )
    parser.add_argument(
        "--run_boltz",
        action="store_true",
        help="If set, run Boltz-based scoring/filtering.",
    )
    parser.add_argument(
        "--binder_chain",
        type=str,
        default="A",
        help="Binder chain ID for PyRosetta interface energy.",
    )
    parser.add_argument(
        "--target_chain",
        type=str,
        default="M",
        help="Target chain ID for PyRosetta interface energy.",
    )
    parser.add_argument(
        "--interface_dist",
        type=float,
        default=10.0,
        help="Interface distance cutoff (Å) for defining interface residues.",
    )
    parser.add_argument(
        "--energy_threshold",
        type=float,
        default=-5.0,
        help="Energy threshold (kcal/mol) for defining key residues in PyRosetta.",
    )

    args = parser.parse_args()

    merge_dir = run_pipeline(
        working_dir=args.working_dir,
        input_pdb_dir=args.input_pdb_dir,
        chains_to_design=args.chains_to_design,
        num_seq_per_target=args.num_seq_per_target,
        sampling_temp=args.sampling_temp,
        seed=args.seed,
        fixed_positions=args.fixed_positions,
        run_protenix=args.run_protenix,
        run_boltz=args.run_boltz,
        binder_chain=args.binder_chain,
        target_chain=args.target_chain,
        interface_dist=args.interface_dist,
        energy_threshold=args.energy_threshold,
    )

    print(f"✅ Pipeline complete. Merged motif PDBs are in: {merge_dir}")


if __name__ == "__main__":
    main()

