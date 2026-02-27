import argparse
import os

from Bio import PDB
import numpy as np
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




def get_residue_pairs_within_distance(pdb_file, binder_id, target_id, distance_threshold=10.0):

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    binder_chain = model[binder_id]
    target_chain = model[target_id]

    selected_pairs = set()  
    # selected_target = set()
    # selected_binder = set()
    for res1 in binder_chain:
        coord1 = get_cb_or_ca(res1)
        for res2 in target_chain:
            coord2 = get_cb_or_ca(res2)
            if coord1 is not None and coord2 is not None:
                distance = np.linalg.norm(coord1 - coord2)
                if distance <= distance_threshold:
                    selected_pairs.add((res1.id[1], res2.id[1])) 
                    # selected_target.add(res1.id[1])
                    # selected_binder.add(res2.id[1])

    return selected_pairs


def get_interface_energy(interchain_score_df, interface_pair, binder_id, plot_path):
    # interchain_score_df = pd.read_csv(interchain_score_path)
    #('L', 'R', l_res, r_res): (pair_idx1, pair_idx2, l_res, r_res)

    interchain_score_df['in_interface'] = interchain_score_df.apply(lambda row: (int(row['binder_res']), int(row['target_res'])) in interface_pair, axis=1)
    interface_score_df = interchain_score_df.loc[interchain_score_df['in_interface']==True]
    summed_df = interface_score_df.groupby('binder_res')['total'].sum().reset_index()
    summed_dict = summed_df.set_index('binder_res')['total'].to_dict()
    print(interface_score_df.head(2))
    print(interface_score_df.columns)
    plot_score(interchain_score_df, plot_path)
    return summed_dict


def init_pyrosetta(extra_flags: str = ""):
    # add whatever flags you normally need (e.g. database path, silent flags, etc.)
    pyrosetta.init(f"-ignore_zero_occupancy false {extra_flags}")


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
    mm.set_bb(False)
    mm.set_chi(True)

    relax = FastRelax(scorefxn, 2)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.apply(pose)


def compute_interface_pair_energies_with_pyrosetta(
    pdb_file: str,
    binder_chain: str,
    target_chain: str,
    distance_threshold: float = 10.0,
):
    """
    Compute pairwise energies for residue pairs at the interface between
    binder_chain and target_chain using PyRosetta.
    Returns a DataFrame with columns:
        binder_id, target_id, binder_res, target_res, total, [per-term columns...]
    """
    pose = rosetta.core.import_pose.pose_from_file(pdb_file)

    # standard full-atom scorefunction (similar to what FastRelax uses)
    scorefxn = pyrosetta.get_fa_scorefxn()
    # mimic the XML FastRelax protocol before scoring
    apply_fast_relax_with_task(pose, scorefxn, binder_chain, target_chain)
    scorefxn(pose)  # populate pose.energies() after relax

    weights = scorefxn.weights()
    pdb_info = pose.pdb_info()

    # use your existing Biopython-based distance selector to define interface pairs
    interface_pairs = get_residue_pairs_within_distance(
        pdb_file, binder_chain, target_chain, distance_threshold=distance_threshold
    )

    rows = []
    for binder_resnum, target_resnum in interface_pairs:
        # map PDB numbering (chain, resnum) to pose index
        i = pdb_info.pdb2pose(binder_chain, binder_resnum)
        j = pdb_info.pdb2pose(target_chain, target_resnum)
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
            "target_id": target_chain,
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

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    init_pyrosetta(args.extra_flags)

    # step 1: compute per-residue-pair energies at the interface with PyRosetta
    interchain_score_df = compute_interface_pair_energies_with_pyrosetta(
        pdb_file=args.pdb,
        binder_chain=args.binder_id,
        target_chain=args.target_id,
        distance_threshold=args.interface_dist,
    )

    # step 2: save per-pair energies (including per-term columns) for inspection
    pdb_basename = os.path.splitext(os.path.basename(args.pdb))[0]
    per_pair_csv = os.path.join(
        args.output_dir, f"{pdb_basename}_pair_energies_pyrosetta.csv"
    )
    interchain_score_df.to_csv(per_pair_csv, index=False)

    # step 3: reuse your existing aggregation + plotting code
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

    # save per-binder-residue summed energies as a small one-row CSV, or extend as needed
    out_df = pd.DataFrame(
        [{"pdbpath": args.pdb, "binder_id": args.binder_id, "binder_energy": summed_dict}]
    )
    out_csv = os.path.join(args.output_dir, "residue_energy_pyrosetta.csv")
    out_df.to_csv(out_csv, index=False)

    print(f"Saved per-residue interface energies to {out_csv}")
    print(f"Saved heatmap to {plot_path}")


if __name__ == "__main__":
    main()