# MGraphDTA 3D: Interactive Multi-Complex Molecular Viewer & Contact Analyzer

An interactive WebGL-powered 3D molecular viewer and geometric binding pocket analyzer for drug–target complexes studied in **MGraphDTA** (Graph Neural Network for Explainable Drug–Target Binding Affinity Prediction).

Built with [3Dmol.js](https://3dmol.org/), this project visualizes experimental X-ray crystallographic structures of diverse human kinase targets in complex with FDA-approved and investigational kinase inhibitors from the Davis kinome benchmark.

---

## 🌟 Key Features

### 1. ⊞ Multi-Complex Gallery (Simultaneous 3D Grid)
* Live, simultaneous rendering of all 5 kinase–inhibitor complexes in an interactive card grid.
* Independent camera controls, 3D rotation, auto-spin toggles, and live contact metric badges for each complex.
* Direct one-click jump to detailed inspection for any structure.

### 2. ⇄ Dual Complex Comparison
* Side-by-side split comparison of any two complexes.
* Compare binding pocket conformations, ligand poses, and heavy-atom contact distributions between targets.
* View synchronization to align orientations.

### 3. 🔍 Detailed Pocket Inspector & Contact Analysis
* High-resolution viewer with toggles for:
  * **Protein Chain A** (spectrum rainbow cartoon)
  * **Binding Pocket** (orange stick representation within $\le 4.0\text{ \AA}$)
  * **Ligand** (stick rendering)
  * **Auto-Spin** and **Full Chain A** camera reset
* **Real-Time Geometric Contact Table**: Calculates minimum heavy-atom distances between ligand and protein residues. Clicking any residue row highlights it in **cyan sticks** and centers the camera on the residue.

### 4. ⚡ Zero-CORS Instant Loading
* Includes pre-indexed structure coordinate cache (`structures/structures_data.js`), allowing the application to work out-of-the-box both over HTTP and when opened directly via `file://` in any browser.

---

## 🧬 Structures Included

All structures are experimental X-ray crystallography complexes from the Protein Data Bank (PDB):

| Target | Gene | PDB ID | Ligand / Drug | Description & Medical Context |
| :--- | :---: | :---: | :---: | :--- |
| **ABL1** | `ABL1` | **2HYY** | Imatinib (`STI`) | BCR-ABL tyrosine kinase inhibitor (Gleevec) for Chronic Myeloid Leukemia |
| **BRAF** | `BRAF` | **1UWH** | Sorafenib (`BAX`) | Oncogenic RAF serine/threonine kinase inhibitor (Nexavar) for melanoma |
| **VEGFR2** | `KDR` | **3EFL** | Motesanib (`706`) | Anti-angiogenesis receptor tyrosine kinase inhibitor (AMG 706) |
| **EGFR** | `EGFR` | **1M17** | Erlotinib (`AQ4`) | Epidermal growth factor receptor inhibitor (Tarceva) for NSCLC |
| **SRC** | `SRC` | **1YOL** | CGP77675 (`S03`) | Proto-oncogene tyrosine kinase Src selective inhibitor |

---

## 🚀 Getting Started

### Option 1: Run with Local HTTP Server (Recommended)
```bash
# Clone the repository
git clone https://github.com/Dksingh0001/MGraphDTA_3D.git
cd MGraphDTA_3D

# Start local HTTP server using Python
python -m http.server 8080
```
Open **[http://localhost:8080](http://localhost:8080)** in your browser.

### Option 2: Open Directly in Browser
Because structure data is bundled in `structures/structures_data.js`, you can directly double-click `index.html` or `index_multi_structure.html` to open it in Chrome, Edge, Firefox, or Safari without CORS fetch errors.

---

## 📊 Offline Python Contact Analysis

To compute and export heavy-atom contacts offline:
```bash
python analyze_contacts.py
```
This parses `3EFL.pdb`, calculates all ligand-to-protein heavy-atom Euclidean distances, prints the nearest contact residues within $4.0\text{ \AA}$, and saves the full table to `ligand_residue_contacts.csv`.

*Note: The script runs with standard Python libraries and NumPy (Biopython is optional).*

---

## 📁 Repository Structure

```
MGraphDTA_3D/
├── index.html                   # Main interactive 3D multi-structure viewer application
├── index_multi_structure.html   # Multi-structure viewer & comparison interface
├── analyze_contacts.py          # Standalone Python script for geometric contact analysis
├── ligand_residue_contacts.csv  # Pre-computed contact table for VEGFR2–Motesanib (3EFL)
├── 3EFL.pdb                     # PDB structure file for VEGFR2 complex
├── 3EFL.fasta                   # Target protein sequence
├── structures/                  # PDB structures directory
│   ├── 1M17.pdb                 # EGFR + Erlotinib
│   ├── 1UWH.pdb                 # BRAF + Sorafenib
│   ├── 1YOL.pdb                 # SRC + CGP77675
│   ├── 2HYY.pdb                 # ABL1 + Imatinib
│   ├── 3EFL.pdb                 # VEGFR2 + Motesanib
│   └── structures_data.js       # Pre-indexed structure data bundle for zero-CORS loading
├── .gitignore                   # Git ignore patterns
└── README.md                    # Project documentation
```

---

## 🔗 Related Projects
* [MGraphDTA](https://github.com/Dksingh0001/MGraphDTA): Graph Neural Network implementation for explainable drug–target binding affinity prediction on Davis and KI datasets.

---

## 📄 License
This project is open-source under the MIT License.
