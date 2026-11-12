# Exercise week 04: Protein Contact Maps & Structural Analysis

This repository contains the code and results for **Exercise 4** of the bioinformatics project. The goal of this exercise is to parse protein PDB structures, extract $C_\alpha$ atom coordinates, compute pairwise Euclidean distance matrices, and generate binary contact maps using **Biopython** and **Matplotlib**.

---

## Objectives
- **PDB Parsing:** Automatically extract $C_\alpha$ atoms from structural files using Biopython's parser.
- **Distance Matrix & Contact Mapping:** Compute inter-residue distances and apply an $8.0\,\text{\AA}$ distance threshold to build binary contact maps ($1$ for contact, $0$ otherwise).
- **Comparative Visualization:** Process a dataset of 5 distinct protein structures and plot their contact maps side-by-side in a high-resolution comparative figure.

---

## Dataset (Proteins Analyzed)
The script processes the following 5 PDB files:
1. `1BOL.pdb`
2. `1L60.pdb`
3. `2H6R.pdb`
4. `5B61.pdb`
5. `101M.pdb`

---

## Repository Structure

```text
exercise-4-contact-maps/
├── src/
│   └── contact_maps.py        # Main Python script with the analysis pipeline
├── outputs/
│   └── contact_maps_comparison.png  # Generated 5-panel comparative plot (300 DPI)
├── data/                      # PDB input files
├── requirements.txt           # Project dependencies
└── README.md                  # Project documentation

