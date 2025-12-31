# --*-- conding:utf-8 --*--
# @time:4/24/25 08:55
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:create_D_R_info.py

"""
Assemble RMSD and docking experiment data into per-protein JSON files,
naming output files as '{pdb_id}_RMSD_docking_result.json', placing
the overall average affinity as the first field under "docking".

For each pdb_id folder under BASE_DIR:
  1. Read RMSD from rmsd_info.txt
  2. Locate the docking directory under DOCK_ROOT matching pattern '{pdb_id}_*'
  3. Parse summary.txt for per-run averages and final average affinity
  4. Parse each run's log for seed and detailed mode results (9 modes)
  5. Write a JSON file into QDockBank/{pdb_id}

Assumes:
  - rmsd_info.txt is in the current working directory
  - quantum_docking_result/ and QDockBank/ are sibling directories
"""

import re
import json
from pathlib import Path

# Input files and directories
RMSD_FILE = Path("rmsd_info.txt")
DOCK_ROOT = Path("quantum_docking_result")
BASE_DIR = Path("QDockBank")

def parse_rmsd(path):
    rmsd = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pdb, val = line.split()
            rmsd[pdb] = float(val)
    return rmsd

def parse_summary(summary_path):
    per_run = {}
    avg_affinity = None
    with summary_path.open() as f:
        for line in f:
            line = line.strip()
            # match per-run averages
            m_run = re.match(
                r"run(\d+)\s*=>\s*affinity=([\-\d.]+),\s*rmsd_l\.b\.=([\-\d.]+),\s*rmsd_u\.b\.=([\-\d.]+)",
                line
            )
            if m_run:
                i, aff, lb, ub = m_run.groups()
                per_run[int(i)] = {
                    "affinity": float(aff),
                    "rmsd_l_b": float(lb),
                    "rmsd_u_b": float(ub)
                }
            # match final average line
            elif line.startswith("affinity="):
                parts = dict(re.findall(r"(\w+)=([\-\d.]+)", line))
                avg_affinity = float(parts.get("affinity", 0))
    return per_run, avg_affinity


def parse_run_log(run_dir):
    log_file = next(run_dir.glob("*ligand_docking_log.txt"), None)
    if not log_file:
        raise FileNotFoundError(f"No log file in {run_dir}")
    lines = [l.rstrip() for l in log_file.open()]
    seed = None
    # find seed
    for line in lines:
        m = re.search(r"Using random seed:\s*(\d+)", line)
        if m:
            seed = int(m.group(1))
            break
    # find mode table
    start = next((i for i,l in enumerate(lines) if l.strip().startswith("mode")), None)
    sep = next((i for i,l in enumerate(lines[start:]) if l.strip().startswith("-----")), None)
    if start is None or sep is None:
        raise ValueError(f"Mode table missing in {log_file}")
    modes = []
    idx = start + sep + 1
    for line in lines[idx:idx+9]:
        parts = line.split()
        if len(parts) >= 4:
            modes.append({
                "mode": int(parts[0]),
                "affinity": float(parts[1]),
                "rmsd_l_b": float(parts[2]),
                "rmsd_u_b": float(parts[3])
            })
    return seed, modes


def main():
    rmsd_map = parse_rmsd(RMSD_FILE)

    for pdb_dir in BASE_DIR.iterdir():
        if not pdb_dir.is_dir():
            continue
        pdb_id = pdb_dir.name
        # locate docking folder by pattern
        docking_dir = next(DOCK_ROOT.rglob(f"{pdb_id}_*"), None)
        if not docking_dir or not docking_dir.is_dir():
            print(f"Warning: docking folder for {pdb_id} not found; skipping")
            continue

        # parse summary.txt
        per_run, avg_affinity = parse_summary(docking_dir / "summary.txt")

        # parse each run log
        runs = []
        i = 1
        while True:
            run_dir = docking_dir / f"run{i}"
            if not run_dir.exists():
                break
            seed, modes = parse_run_log(run_dir)
            runs.append({
                "run": i,
                "seed": seed,
                "average": per_run.get(i, {}),
                "modes": modes
            })
            i += 1

        # assemble JSON
        record = {
            "rmsd": rmsd_map.get(pdb_id),
            "rmsd_method": "C-alpha coordinates",
            "rmsd_tool": "BioPython",
            "docking": {
                "average_affinity": avg_affinity,
                "runs": runs
            }
        }

        out_file = pdb_dir / f"{pdb_id}_RMSD_docking_result.json"
        with out_file.open("w") as f:
            json.dump(record, f, indent=2)
        print(f"Created {out_file}")

if __name__ == "__main__":
    main()




