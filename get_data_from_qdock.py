# --*-- conding:utf-8 --*--
# @time:12/30/25 19:41
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:get_data_from_qdock.py

from pathlib import Path
import csv
import json
from typing import Any, Dict, List, Optional, Tuple

# =========================
# Hardcoded config (edit here)
# =========================
SCRIPT_DIR = Path(__file__).resolve().parent   # directory of this script

ROOT_REL = Path("QDockBank")                  # relative to SCRIPT_DIR
OUT_CSV_REL = Path("outputs/qdockbank_selected_sequences.csv")
OUT_LOG_REL = Path("outputs/qdockbank_selected_sequences_skipped.csv")

K = 5  # target length for long sequences

# =========================
# Helpers
# =========================
def find_metadata_json(folder: Path) -> Optional[Path]:
    """Find *metadata.json in a pdb folder. If multiple, pick deterministically."""
    candidates = list(folder.glob("*metadata.json"))
    if not candidates:
        return None
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

    seq = seq.strip()
    if not isinstance(seq_len, int):
        seq_len = len(seq)

    chain = pi.get("chain", "")
    residues = pi.get("residues", "")

    return {
        "pdb_id": pdb_id,
        "sequence": seq,
        "sequence_length": int(seq_len),
        "chain": chain if isinstance(chain, str) else "",
        "residues": residues if isinstance(residues, str) else "",
    }

def choose_unique_sequence(sequence: str, selected_set: set) -> Tuple[Optional[str], str]:
    """
    全局唯一规则：
    - len<=K：只能选原序列；若重复 -> 跳过
    - len>K：从左到右滑窗找第一个不重复的 K-mer；若都重复 -> 跳过
    """
    seq = sequence.strip()
    n = len(seq)
    if n == 0:
        return None, "empty_sequence"

    if n <= K:
        if seq in selected_set:
            return None, "duplicate_short_sequence"
        return seq, ""

    for i in range(0, n - K + 1):
        sub = seq[i:i + K]
        if sub not in selected_set:
            return sub, ""
    return None, "no_unique_kmer"

# =========================
# Main
# =========================
def main():
    root = ROOT_DIR.expanduser().resolve()
    out_csv = OUT_CSV.expanduser().resolve()
    out_log = OUT_LOG.expanduser().resolve()

    if not root.exists() or not root.is_dir():
        raise RuntimeError(f"[ERROR] ROOT_DIR not found or not a directory: {root}")

    selected_set = set()
    rows: List[Dict[str, Any]] = []
    skipped: List[Dict[str, Any]] = []

    subfolders = [p for p in root.iterdir() if p.is_dir()]
    subfolders.sort(key=lambda p: p.name)

    for folder in subfolders:
        meta_path = find_metadata_json(folder)
        if meta_path is None:
            skipped.append({"folder": folder.name, "reason": "metadata_not_found", "meta_path": ""})
            continue

        try:
            meta = load_json(meta_path)
        except Exception as e:
            skipped.append({
                "folder": folder.name,
                "reason": f"metadata_read_error:{type(e).__name__}",
                "meta_path": str(meta_path),
            })
            continue

        info = extract_protein_info(meta)
        if info is None:
            skipped.append({
                "folder": folder.name,
                "reason": "protein_information_missing_or_invalid",
                "meta_path": str(meta_path),
            })
            continue

        original_seq = info["sequence"]
        original_len = info["sequence_length"]

        selected_seq, fail_reason = choose_unique_sequence(original_seq, selected_set)
        if selected_seq is None:
            skipped.append({
                "folder": folder.name,
                "pdb_id": info["pdb_id"],
                "reason": fail_reason,
                "meta_path": str(meta_path),
                "original_sequence": original_seq,
                "original_length": original_len,
            })
            continue

        selected_set.add(selected_seq)

        rows.append({
            "pdb_id": info["pdb_id"],
            "chain": info["chain"],
            "residues": info["residues"],
            "original_sequence": original_seq,
            "original_length": original_len,
            "selected_sequence": selected_seq,
            "selected_length": len(selected_seq),
            "metadata_json": str(meta_path),
            "folder": folder.name,
        })

    # write outputs
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
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    out_log.parent.mkdir(parents=True, exist_ok=True)
    if skipped:
        log_fields = sorted({k for item in skipped for k in item.keys()})
        with out_log.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=log_fields)
            w.writeheader()
            w.writerows(skipped)
    else:
        # still create an empty file with minimal header
        with out_log.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["folder", "reason", "meta_path"])
            w.writeheader()

    print(f"[OK] Root: {root}")
    print(f"[OK] Selected (global-unique): {len(rows)}")
    print(f"[OK] Skipped: {len(skipped)}")
    print(f"[OK] CSV: {out_csv}")
    print(f"[OK] Log: {out_log}")

if __name__ == "__main__":
    main()
