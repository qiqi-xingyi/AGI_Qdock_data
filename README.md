# QDockBank Dataset

**QDockBank** is a benchmark dataset of protein fragment structures predicted using utility‐level quantum computers, specifically designed for protein–ligand docking tasks. This README provides an overview of the dataset’s contents, directory structure, and usage instructions.

---

## 📦 Dataset Contents

1. **Predicted Structures**\
   Final protein fragment models in PDB format with PDB‐compliant headers containing dataset metadata.

2. **Quantum Metadata**\
   Per‐fragment JSON files (`{pdb_id}_metadata.json`) including:

   - **protein\_information**:
     - `pdb_id`
     - `sequence`
     - `sequence_length`
     - `chain` (always "A")
     - `residues` (start–end indices)
   - **quantum\_metadata**:
     - `number_of_qubits`
     - `circuit_depth`
     - `lowest_energy`
     - `highest_energy`
     - `energy_range`
     - `execution_time_s`

3. **Docking & RMSD Results**\
   Per‐fragment JSON files (`{pdb_id}_RMSD_docking_result.json`) containing:

   - `rmsd`: root‐mean‐square deviation vs. experimental structure
   - `rmsd_method` / `rmsd_tool`
   - **docking**:
     - `average_affinity` (top‐pose average across 20 runs)
     - `runs`: list of 20 objects, each with:
       - `run`: run index
       - `seed`: random seed used
       - `average`: affinity & RMSD bounds
       - `modes`: 9 binding modes (`affinity`, `rmsd_l_b`, `rmsd_u_b`)

4. **Index File**\
   `index.csv` or `index.json` summarizing all fragments: group (L/M/S), sequence, length, residue range.

---

## 🗂 Directory Structure

```
QDockBank/
├─ index.csv
├─ quantum_metadata/
│    ├─ 1xyz/
│    │    ├─ 1xyz_metadata.json
│    │    ├─ 1xyz.pdb
│    │    └─ 1xyz_RMSD_docking_result.json
│    └─ …
└─ quantum_docking_result/  (raw docking logs & outputs)
```

> Each fragment folder (`{pdb_id}`) contains:
>
> - `{pdb_id}.pdb`
> - `{pdb_id}_metadata.json`
> - `{pdb_id}_RMSD_docking_result.json`

---

## 🚀 Usage

1. **Visualize structures**\
   Open PDB files in PyMOL, ChimeraX, or other structural viewers.

2. **Programmatic access**\
   Load JSON metadata in Python:

   ```python
   import json
   from pathlib import Path

   folder = Path('QDockBank/1xyz')
   meta = json.load((folder/'1xyz_metadata.json').open())
   docking = json.load((folder/'1xyz_RMSD_docking_result.json').open())
   ```

---

For questions or issues, please open an issue in the project repository.

