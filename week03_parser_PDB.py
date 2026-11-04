#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Writing a PDB file parser from scratch
#1. Open the file in read mode and read the first row
file_path = "spermine_synthase.pdb"
with open(file_path, 'r') as open_file:
    first_row = open_file.readline()

    print(first_row)



# In[4]:


#2. Filtering the atom rows (ATOM and HETATM)
with open(file_path, "r") as open_file:
    for row in open_file:
        if row.startswith("ATOM") or row.startswith("HETATM"):
            print(row.strip())
            break


# In[5]:


#3. Slicing the row to extract the data
with open(file_path, "r") as open_file:
    for row in open_file:
        if row.startswith("ATOM") or row.startswith("HETATM"):
            atom_name = row[12:16].strip()
            residue_name = row[17:20].strip()
            chain = row[21].strip()

            x_coordinate = float(row[30:38])
            y_coordinate = float(row[38:46])
            z_coordinate = float(row[46:54])

            b_factor = float(row[60:66])

            print(f"Atom: {atom_name}, Residue: {residue_name}, Chain: {chain}, X: {x_coordinate} B-factor: {b_factor}")
            break


# In[9]:


#4. Put it all together and test it with BioPython
from Bio import PDB
def my_parser_pdb(file_path):
    collected_atoms = []

    with open(file_path, "r") as open_file:
        for row in open_file:
            if row.startswith("ATOM") or row.startswith("HETATM"):
                atom_data = { 
                    "name": row[12:16].strip(), 
                    "residue": row[17:20].strip(), 
                    "chain": row[21].strip(), 
                    "x": float(row[30:38]), 
                    "y": float(row[38:46]), 
                    "z": float(row[46:54]), 
                    "b_factor": float(row[60:66]) 
                }
                collected_atoms.append(atom_data)
    return collected_atoms



# In[10]:


target_file = "spermine_synthase.pdb" 
my_atoms = my_parser_pdb(target_file)
print(f"My parser found {len(my_atoms)} atoms.")

official_parser = PDB.PDBParser(QUIET=True)
official_structure = official_parser.get_structure("protein_1", target_file)

official_atoms = list(official_structure.get_atoms())
print(f"Biopython found {len(official_atoms)} atoms.")

assert len(my_atoms) == len(official_atoms), "Note: The number of atoms does not match"
print("Test passed successfully!")


# In[ ]:




