#!/usr/bin/env python
# coding: utf-8

# In[18]:


"""
Protein Contact Map Generator & Structural Comparison
Author: Alessandro Mieli
Description: Extracts C-alpha atoms using Biopython, computes 8.0 Å contact maps, 
and visualizes structural classes side-by-side.
"""


# In[19]:


from Bio import PDB
import matplotlib.pyplot as plt
import os

def create_biopython_map(file_path, threshold=8.0):
    """
    Parses a PDB file, extracts C-alpha atoms, calculates the distance matrix 
    using Biopython, and returns a binary contact map (0 or 1).
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("my_protein", file_path)

    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca_atoms.append(residue["CA"])

    residues_num = len(ca_atoms)
    print(f"Biopython extracted {residues_num} alpha carbons (CA).")

    # Calculating the distance matrix using Biopython's interatomic subtraction
    distance_matrix = []
    for i in range(residues_num):
        row = []
        for j in range(residues_num):
            distance = ca_atoms[i] - ca_atoms[j]
            row.append(distance)
        distance_matrix.append(row)

    # Conversion to a binary matrix (1 if <= threshold, 0 otherwise)
    contact_map = []
    for row in distance_matrix:
        new_row = [1 if d <= threshold else 0 for d in row]
        contact_map.append(new_row)

    return contact_map, ca_atoms, distance_matrix


# In[12]:


# --- MAIN EXECUTION AND MULTIPLE DISPLAYS ---
if __name__ == "__main__":
    # Make sure the ‘outputs’ folder exists to save the results
    os.makedirs("outputs", exist_ok=True)

    # Lista delle 5 proteine di classi strutturali diverse
    proteins_list = ["1BOL.pdb", "1L60.pdb", "2H6R.pdb", "5B61.pdb", "101M.pdb"]

    fig, axes = plt.subplots(1, len(proteins_list), figsize=(15, 3))

    for idx, file_path in enumerate(proteins_list):
        print(f"\nProcessing: {file_path}...")
        try:
            my_map, ca_list, distance_matrix = create_biopython_map(file_path, threshold=8.0)

            # Plot della mappa di contatto
            axes[idx].imshow(my_map, cmap="binary", origin="lower")
            axes[idx].set_title(file_path.split(".")[0].upper(), fontsize=10)
            axes[idx].set_xlabel("Residue", fontsize=8)
            if idx == 0:
                axes[idx].set_ylabel("Residue", fontsize=8)

        except FileNotFoundError:
            print(f"Warning: The file '{file_path}' was not found in the folder.")
            axes[idx].text(0.5, 0.5, "File not found", ha='center', va='center', fontsize=8)
            axes[idx].axis('off')

    plt.tight_layout()

    # Save in high resolution to the “outputs” folder
    output_path = "outputs/contact_maps_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nComparison chart successfully saved in: {output_path}")

    plt.show()

