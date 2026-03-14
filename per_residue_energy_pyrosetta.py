import argparse
import os
import glob
import pandas as pd
import pyrosetta
from pyrosetta import rosetta

# residue selectors and task operations to mimic the XML protocol
from pyrosetta.rosetta.core.select import residue_selector as rs
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking,
    OperateOnResidueSubset,
    PreventRepackingRLT,
)

# FastRelax mover
from pyrosetta.rosetta.protocols.relax import FastRelax

# reuse your existing logic for interface detection + aggregation (same directory)
from get_interface_energy import (
    get_residue_pairs_within_distance,
    get_interface_energy,
)


def init_pyrosetta(extra_flags: str = ""):
    # add whatever flags you normally need (e.g. database path, silent flags, etc.)
    pyrosetta.init(f"-ignore_zero_occupancy false {extra_flags}")


def relax_pose(
    pdb_file: str,
    binder_chain: str,
    target_chain: str,
    xml_protocol: str | None = None,
):
    """
    Load a pose from pdb_file, run relaxation, and return the relaxed pose
    and scorefunction.

    If xml_protocol is provided, use RosettaScripts to run exactly that XML
    (e.g. your update.xml with FastRelax + ScoreCutoffFilter). Otherwise,
    fall back to the hand-translated FastRelax setup.
    """
    pose = rosetta.core.import_pose.pose_from_file(pdb_file)
    scorefxn = pyrosetta.get_fa_scorefxn()

    if xml_protocol is not None:
        # Use XmlObjects helper to load and run the RosettaScripts XML protocol.
        # This mirrors how RosettaScripts is used in compiled Rosetta.
        with open(xml_protocol, "r") as f:
            xml_string = f.read()
        xml_objs = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
            xml_string
        )
        # Top-level protocol mover is conventionally named "ParsedProtocol"
        rs_mover = xml_objs.get_mover("ParsedProtocol")
        rs_mover.apply(pose)
        # Ensure energies are populated under the standard full-atom scorefunction
        scorefxn(pose)
        return pose, scorefxn

    # Fallback: approximate XML behavior using explicit FastRelax + selectors
    apply_fast_relax_with_task(pose, scorefxn, binder_chain, target_chain)
    scorefxn(pose)  # populate pose.energies() after relax
    return pose, scorefxn


def apply_fast_relax_with_task(
    pose: rosetta.core.pose.Pose,
    scorefxn: rosetta.core.scoring.ScoreFunction,
    binder_chain: str,
    target_chain: str,
):
    """
    Approximate the XML:
      - Neighborhood around binder/target chains within 20 Å
      - Prevent repacking outside that neighborhood
      - RestrictToRepacking globally
      - FastRelax (2 repeats) with sidechain moves only (no backbone)
    """
    # Residue selectors
    binder_sel = rs.ChainSelector(binder_chain)
    target_sel = rs.ChainSelector(target_chain)
    focus_sel = rs.OrResidueSelector(binder_sel, target_sel)

    nbr_sel = rs.NeighborhoodResidueSelector()
    nbr_sel.set_focus_selector(focus_sel)
    nbr_sel.set_distance(20.0)
    nbr_sel.set_include_focus_in_subset(True)

    others_sel = rs.NotResidueSelector(nbr_sel)

    # TaskFactory: Restrict to repacking, and prevent repacking on "others"
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), others_sel))

    # MoveMap: no backbone moves, allow sidechain (chi) moves
    mm = rosetta.core.kinematics.MoveMap()
    #mm.set_bb(False)
    mm.set_chi(True)

    relax = FastRelax(scorefxn, 2)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.apply(pose)


def compute_interface_pair_energies_with_pyrosetta(
    pose: rosetta.core.pose.Pose,
    pdb_file: str,
    binder_chain: str,
    target_chain: str,
    scorefxn: rosetta.core.scoring.ScoreFunction,
    distance_threshold: float = 10.0,
):
    """
    Compute pairwise energies for residue pairs at the interface between
    binder_chain and target_chain using a (typically relaxed) pose.
    Returns a DataFrame with columns:
        binder_id, target_id, binder_res, target_res, total, [per-term columns...]
    """
    weights = scorefxn.weights()
    pdb_info = pose.pdb_info()

    # Allow multiple target chains, e.g. "A,B,C"
    if isinstance(target_chain, str) and "," in target_chain:
        target_chains = [t.strip() for t in target_chain.split(",") if t.strip()]
    else:
        target_chains = [target_chain]

    rows = []
    for tgt_chain in target_chains:
        # use your existing Biopython-based distance selector to define interface pairs
        interface_pairs = get_residue_pairs_within_distance(
            pdb_file, binder_chain, tgt_chain, distance_threshold=distance_threshold
        )

        for binder_resnum, target_resnum in interface_pairs:
            # map PDB numbering (chain, resnum) to pose index
            i = pdb_info.pdb2pose(binder_chain, binder_resnum)
            j = pdb_info.pdb2pose(tgt_chain, target_resnum)
            if i == 0 or j == 0:
                continue

            # compute pairwise energy for this residue pair using the scorefunction
            emap = rosetta.core.scoring.EMapVector()
            # context-independent and context-dependent 2-body energies
            scorefxn.eval_ci_2b(pose.residue(i), pose.residue(j), pose, emap)
            scorefxn.eval_cd_2b(pose.residue(i), pose.residue(j), pose, emap)

            # weighted total and per-term contributions over all score types
            total = 0.0
            per_term = {}
            end_enum = int(rosetta.core.scoring.end_of_score_type_enumeration)
            for st_int in range(1, end_enum):
                st = rosetta.core.scoring.ScoreType(st_int)
                w = weights[st]
                if w == 0.0:
                    continue
                val = emap[st]
                if val == 0.0:
                    continue
                total += w * val
                term_name = rosetta.core.scoring.name_from_score_type(st)
                per_term[term_name] = per_term.get(term_name, 0.0) + val

            row = {
                "binder_id": binder_chain,
                "target_id": tgt_chain,
                "binder_res": str(binder_resnum),
                "target_res": str(target_resnum),
                "total": total,
            }
            # add raw (unweighted) per-term energies as separate columns
            row.update(per_term)
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["binder_id", "target_id", "binder_res", "target_res", "total"])

    df = pd.DataFrame(rows)
    # match original pipeline: keep only favorable (negative) pair energies
    df = df[df["total"] < 0].reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--binder_id", required=True, help="Binder chain ID (e.g. L)")
    parser.add_argument("--target_id", required=True, help="Target chain ID (e.g. R)")
    parser.add_argument(
        "--interface_dist",
        type=float,
        default=12.0,
        help="Interface distance threshold between target and binder",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where residue_energy.csv and plot will be written",
    )
    parser.add_argument(
        "--extra_flags",
        default="",
        help="Extra PyRosetta init flags, if needed (e.g. database path)",
    )
    parser.add_argument(
        "--xml_protocol",
        default=None,
        help="Optional RosettaScripts XML file to run for relax (e.g. update.xml)",
    )

    parser.add_argument(
        "--energy_threshold",
        type=float,
        default=-5.0,
        help="Energy threshold for identifying strongly favorable residues",
    )

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    init_pyrosetta(args.extra_flags)

    # step 1: relax the input structure
    pose, scorefxn = relax_pose(
        pdb_file=args.pdb,
        binder_chain=args.binder_id,
        target_chain=args.target_id,
        xml_protocol=args.xml_protocol,
    )

    # optionally save relaxed PDB for inspection
    pdb_basename = os.path.splitext(os.path.basename(args.pdb))[0]
    relaxed_pdb_path = os.path.join(args.output_dir, f"{pdb_basename}_relaxed.pdb")
    pose.dump_pdb(relaxed_pdb_path)

    # step 2: compute per-residue-pair energies at the interface with PyRosetta
    interchain_score_df = compute_interface_pair_energies_with_pyrosetta(
        pose=pose,
        pdb_file=args.pdb,
        binder_chain=args.binder_id,
        target_chain=args.target_id,
        scorefxn=scorefxn,
        distance_threshold=args.interface_dist,
    )

    # step 3: save per-pair energies (including per-term columns) for inspection
    per_pair_csv = os.path.join(
        args.output_dir, f"{pdb_basename}_pair_energies_pyrosetta.csv"
    )
    interchain_score_df.to_csv(per_pair_csv, index=False)

    # step 4: reuse your existing aggregation + plotting code
    plot_path = os.path.join(
        args.output_dir, f"{pdb_basename}_interface_binder_residues_score_heatmap.png"
    )

    summed_dict = get_interface_energy(
        interchain_score_df,
        # get_interface_energy expects the original interface_pair set
        interface_pair={
            (int(row["binder_res"]), int(row["target_res"]))
            for _, row in interchain_score_df.iterrows()
        },
        binder_id=args.binder_id,
        plot_path=plot_path,
    )

    # save per-binder-residue summed energies as a small one-row CSV
    out_df = pd.DataFrame(
        [{"pdbpath": args.pdb, "binder_id": args.binder_id, "binder_energy": summed_dict}]
    )
    out_csv = os.path.join(args.output_dir, f"{pdb_basename}_residue_energy_pyrosetta.csv")
    out_df.to_csv(out_csv, index=False)

    # step 5: identify strongly favorable residues and record their amino acids
    # filter for any residue with energy < -5.0 and capture the one-letter AA
    pdb_info = pose.pdb_info()
    key_res = {}
    for resnum_str, energy in summed_dict.items():
        energy_val = float(energy)
        if energy_val >= args.energy_threshold:
            continue
        pdb_resnum = int(resnum_str)
        pose_idx = pdb_info.pdb2pose(args.binder_id, pdb_resnum)
        if pose_idx == 0:
            continue
        aa_one_letter = pose.residue(pose_idx).name1()
        key_res[resnum_str] = [energy_val, aa_one_letter]

    fixed_residue_df = pd.DataFrame(
        [
            {
                "pdb_path": args.pdb,
                "pdb_name": pdb_basename,
                "target_id": args.target_id,
                "binder_id": args.binder_id,
                "key_res": key_res,
                "num_fixed_residues": len(key_res),
            }
        ]
    )
    fixed_csv = os.path.join(args.output_dir, f"{pdb_basename}_fixed_residue_pyrosetta.csv")
    fixed_residue_df.to_csv(fixed_csv, index=False)

    print(f"Saved per-residue interface energies to {out_csv}")
    print(f"Saved filtered key residues (E < {args.energy_threshold}) with amino acids to {fixed_csv}")
    print(f"Saved heatmap to {plot_path}")


if __name__ == "__main__":
    main()