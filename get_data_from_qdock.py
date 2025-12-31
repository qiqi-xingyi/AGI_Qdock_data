# --*-- conding:utf-8 --*--
# @time:12/30/25 19:41
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:get_data_from_qdock.py


import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def find_metadata_json(folder: Path) -> Optional[Path]:
    """
    Find metadata json in a pdb folder.
    Priority: *metadata.json (any), if multiple, pick the shortest name (deterministic).
    """
    candidates = list(folder.glob("*metadata.json"))
    if not candidates:
        return None
    # Deterministic selection
    candidates.sort(key=lambda p: (len(p.name), p.name))
    return candidates[0]


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def extract_protein_info(meta: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    pi = meta.get("protein_information")
    if not isinstance(pi, dict):
        return None

    pdb_id = pi.get("pdb_id")
    seq = pi.get("sequence")
    seq_len = pi.get("sequence_length")

    if not isinstance(pdb_id, str) or not isinstance(seq, str):
        return None
    if not isinstance(seq_len, int):
        # fallback: compute from sequence
        seq_len = len(seq)

    chain = pi.get("chain", "")
    residues = pi.get("residues", "")

    return {
        "pdb_id": pdb_id,
        "sequence": seq.strip(),
        "sequence_length": int(seq_len),
        "chain": chain if isinstance(chain, str) else "",
        "residues": residues if isinstance(residues, str) else "",
    }


def choose_unique_sequence(
    sequence: str,
    selected_set: set,
) -> Tuple[Optional[str], str]:
    """
    Return (selected_sequence, reason_if_failed).
    - If len<=5: select the whole seq if unique.
    - If len>5: select a unique 5-mer by left-to-right scanning.
    """
    seq = sequence.strip()
    n = len(seq)
    if n == 0:
        return None, "empty_sequence"

    if n <= 5:
        if seq in selected_set:
            return None, "duplicate_short_sequence"
        return seq, ""

    # n > 5: find a unique 5-mer
    k = 5
    for i in range(0, n - k + 1):
        sub = seq[i : i + k]
        if sub not in selected_set:
            return sub, ""
    return None, "no_unique_5mer"


def main():
    parser = argparse.ArgumentParser(
        description="Traverse QDockBank and extract globally-unique sequences (<=5 or unique 5-mer) into CSV."
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Path to QDockBank root folder (contains pdbid subfolders).",
    )
    parser.add_argument(
        "--out_csv",
        required=True,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--out_log",
        default="",
        help="Optional log file path for skipped entries.",
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    out_csv = Path(args.out_csv).expanduser().resolve()
    out_log = Path(args.out_log).expanduser().resolve() if args.out_log else None

    if not root.exists() or not root.is_dir():
        raise SystemExit(f"[ERROR] root not found or not a directory: {root}")

    selected_set = set()
    rows: List[Dict[str, Any]] = []
    skipped: List[Dict[str, Any]] = []

    # Iterate subfolders deterministically
    subfolders = [p for p in root.iterdir() if p.is_dir()]
    subfolders.sort(key=lambda p: p.name)

    for folder in subfolders:
        meta_path = find_metadata_json(folder)
        if meta_path is None:
            skipped.append(
                {"folder": folder.name, "reason": "metadata_not_found", "meta_path": ""}
            )
            continue

        try:
            meta = load_json(meta_path)
        except Exception as e:
            skipped.append(
                {
                    "folder": folder.name,
                    "reason": f"metadata_read_error:{type(e).__name__}",
                    "meta_path": str(meta_path),
                }
            )
            continue

        info = extract_protein_info(meta)
        if info is None:
            skipped.append(
                {
                    "folder": folder.name,
                    "reason": "protein_information_missing_or_invalid",
                    "meta_path": str(meta_path),
                }
            )
            continue

        original_seq = info["sequence"]
        original_len = info["sequence_length"]

        selected_seq, fail_reason = choose_unique_sequence(original_seq, selected_set)
        if selected_seq is None:
            skipped.append(
                {
                    "folder": folder.name,
                    "pdb_id": info["pdb_id"],
                    "reason": fail_reason,
                    "meta_path": str(meta_path),
                    "original_sequence": original_seq,
                    "original_length": original_len,
                }
            )
            continue

        # Global uniqueness guaranteed by selected_set
        selected_set.add(selected_seq)

        row = {
            "pdb_id": info["pdb_id"],
            "chain": info["chain"],
            "residues": info["residues"],
            "original_sequence": original_seq,
            "original_length": original_len,
            "selected_sequence": selected_seq,
            "selected_length": len(selected_seq),
            "metadata_json": str(meta_path),
            "folder": folder.name,
        }
        rows.append(row)

    # Write CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id",
        "chain",
        "residues",
        "original_sequence",
        "original_length",
        "selected_sequence",
        "selected_length",
        "metadata_json",
        "folder",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Optional log
    if out_log is not None:
        out_log.parent.mkdir(parents=True, exist_ok=True)
        with out_log.open("w", newline="", encoding="utf-8") as f:
            # write as CSV for easy inspection
            log_fields = sorted({k for item in skipped for k in item.keys()})
            writer = csv.DictWriter(f, fieldnames=log_fields)
            writer.writeheader()
            writer.writerows(skipped)

    print(f"[OK] Root: {root}")
    print(f"[OK] Selected unique sequences: {len(rows)}")
    print(f"[OK] Skipped: {len(skipped)}")
    print(f"[OK] CSV written to: {out_csv}")
    if out_log is not None:
        print(f"[OK] Log written to: {out_log}")


if __name__ == "__main__":
    main()
