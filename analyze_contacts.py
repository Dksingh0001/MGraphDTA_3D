import csv
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import PDBParser
    HAS_BIOPDB = True
except ImportError:
    HAS_BIOPDB = False

BASE = Path(r"D:\MGraphDTA_3D")
PDB_FILE = BASE / "3EFL.pdb"
OUTPUT_FILE = BASE / "ligand_residue_contacts.csv"

if not PDB_FILE.exists():
    raise FileNotFoundError(f"Structure not found: {PDB_FILE}")

def is_heavy_atom(element):
    return (element or "").strip().upper() not in {"H", "D", ""}

if HAS_BIOPDB:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3EFL", str(PDB_FILE))
    model = structure[0]
    chain = model["A"]
    ligand_residues = [
        residue for residue in chain
        if residue.resname == "706"
    ]
    if not ligand_residues:
        raise ValueError("Ligand 706 was not found in chain A.")

    ligand_atoms = []
    for residue in ligand_residues:
        for atom in residue:
            elem = (atom.element or "").strip().upper()
            if is_heavy_atom(elem):
                ligand_atoms.append(atom)

    protein_atoms = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        for atom in residue:
            elem = (atom.element or "").strip().upper()
            if is_heavy_atom(elem):
                protein_atoms.append((residue, atom))

    ligand_coords = np.array(
        [atom.coord for atom in ligand_atoms],
        dtype=float
    )
    ligand_names = [atom.name for atom in ligand_atoms]
    protein_atom_records = [
        (residue.id[1], residue.resname, atom.name, np.asarray(atom.coord, dtype=float))
        for residue, atom in protein_atoms
    ]
else:
    # Standalone PDB parser (no biopython dependency required)
    ligand_atoms_list = []
    protein_atom_records = []
    with open(PDB_FILE, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            is_het = line.startswith("HETATM")
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain_id = line[21].strip()
            resseq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            elem = line[76:78].strip().upper() if len(line) >= 78 else atom_name[0].upper()
            if not is_heavy_atom(elem):
                continue
            if chain_id == "A" and resname == "706":
                ligand_atoms_list.append((atom_name, np.array([x, y, z], dtype=float)))
            elif chain_id == "A" and not is_het and resname != "HOH":
                protein_atom_records.append((resseq, resname, atom_name, np.array([x, y, z], dtype=float)))

    if not ligand_atoms_list or not protein_atom_records:
        raise ValueError("Could not find ligand or protein heavy atoms.")

    ligand_coords = np.array([coord for _, coord in ligand_atoms_list], dtype=float)
    ligand_names = [name for name, _ in ligand_atoms_list]


# Calculate all ligand-atom to protein-atom distances.
results = {}

for resseq, resname, atom_name, atom_coord in protein_atom_records:
    distances = np.linalg.norm(ligand_coords - atom_coord, axis=1)

    ligand_index = int(np.argmin(distances))
    distance = float(distances[ligand_index])

    key = (resseq, resname)
    record = {
        "chain": "A",
        "residue_number": resseq,
        "residue_name": resname,
        "minimum_distance_A": distance,
        "protein_atom": atom_name,
        "ligand_atom": ligand_names[ligand_index],
    }

    if (
        key not in results
        or distance < results[key]["minimum_distance_A"]
    ):
        results[key] = record

contacts = sorted(
    results.values(),
    key=lambda item: item["minimum_distance_A"]
)

with OUTPUT_FILE.open("w", newline="", encoding="utf-8") as f:
    fieldnames = [
        "chain",
        "residue_number",
        "residue_name",
        "minimum_distance_A",
        "protein_atom",
        "ligand_atom",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    for contact in contacts:
        writer.writerow({
            **contact,
            "minimum_distance_A":
                f'{contact["minimum_distance_A"]:.3f}',
        })

print("=== EXPERIMENTAL STRUCTURE CONTACT ANALYSIS ===")
print("PDB:", PDB_FILE.name)
print("Ligand heavy atoms:", len(ligand_coords))
print("Protein heavy atoms:", len(protein_atom_records))
print()
print("Nearest protein residues:")
print(
    f'{"Residue":<12} {"Min distance (A)":>16} '
    f'{"Protein atom":>14} {"Ligand atom":>12}'
)

for contact in contacts[:20]:
    residue_label = (
        f'{contact["residue_name"]}'
        f'{contact["residue_number"]}'
    )
    print(
        f'{residue_label:<12} '
        f'{contact["minimum_distance_A"]:>16.3f} '
        f'{contact["protein_atom"]:>14} '
        f'{contact["ligand_atom"]:>12}'
    )

within_4 = sum(
    contact["minimum_distance_A"] <= 4.0
    for contact in contacts
)

print()
print("Residues with a heavy-atom distance <= 4.0 A:", within_4)
print("Full residue table saved to:", OUTPUT_FILE)
print("These are geometric contacts, not proof of binding affinity.")
