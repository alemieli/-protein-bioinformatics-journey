{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "7f3243fd-a407-41c5-81ec-4f893b11add2",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "HEADER    TRANSFERASE                             04-FEB-08   3C6K              \n",
      "\n"
     ]
    }
   ],
   "source": [
    "# Writing a PDB file parser from scratch\n",
    "#1. Open the file in read mode and read the first row\n",
    "file_path = \"spermine_synthase.pdb\"\n",
    "with open(file_path, 'r') as open_file:\n",
    "    first_row = open_file.readline()\n",
    "\n",
    "    print(first_row)\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "ea395db7-9f1f-4a28-86bb-53d6157ebc0e",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "ATOM      1  N   GLY A   5      -9.688  76.291  22.007  1.00 28.08           N\n"
     ]
    }
   ],
   "source": [
    "#2. Filtering the atom rows (ATOM and HETATM)\n",
    "with open(file_path, \"r\") as open_file:\n",
    "    for row in open_file:\n",
    "        if row.startswith(\"ATOM\") or row.startswith(\"HETATM\"):\n",
    "            print(row.strip())\n",
    "            break"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "3e8d63dc-ba5d-42ee-b5b6-9327f5f4ea6f",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Atom: N, Residue: GLY, Chain: A, X: -9.688 B-factor: 28.08\n"
     ]
    }
   ],
   "source": [
    "#3. Slicing the row to extract the data\n",
    "with open(file_path, \"r\") as open_file:\n",
    "    for row in open_file:\n",
    "        if row.startswith(\"ATOM\") or row.startswith(\"HETATM\"):\n",
    "            atom_name = row[12:16].strip()\n",
    "            residue_name = row[17:20].strip()\n",
    "            chain = row[21].strip()\n",
    "\n",
    "            x_coordinate = float(row[30:38])\n",
    "            y_coordinate = float(row[38:46])\n",
    "            z_coordinate = float(row[46:54])\n",
    "\n",
    "            b_factor = float(row[60:66])\n",
    "\n",
    "            print(f\"Atom: {atom_name}, Residue: {residue_name}, Chain: {chain}, X: {x_coordinate} B-factor: {b_factor}\")\n",
    "            break\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "78ac987b-add8-41ee-9ca6-18bf7b98e1d8",
   "metadata": {},
   "outputs": [],
   "source": [
    "#4. Put it all together and test it with BioPython\n",
    "from Bio import PDB\n",
    "def my_parser_pdb(file_path):\n",
    "    collected_atoms = []\n",
    "\n",
    "    with open(file_path, \"r\") as open_file:\n",
    "        for row in open_file:\n",
    "            if row.startswith(\"ATOM\") or row.startswith(\"HETATM\"):\n",
    "                atom_data = { \n",
    "                    \"name\": row[12:16].strip(), \n",
    "                    \"residue\": row[17:20].strip(), \n",
    "                    \"chain\": row[21].strip(), \n",
    "                    \"x\": float(row[30:38]), \n",
    "                    \"y\": float(row[38:46]), \n",
    "                    \"z\": float(row[46:54]), \n",
    "                    \"b_factor\": float(row[60:66]) \n",
    "                }\n",
    "                collected_atoms.append(atom_data)\n",
    "    return collected_atoms\n",
    "            \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "96ab7007-6da1-4af2-8202-037f8b102d07",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "My parser found 12240 atoms.\n",
      "Biopython found 12240 atoms.\n",
      "Test passed successfully!\n"
     ]
    }
   ],
   "source": [
    "target_file = \"spermine_synthase.pdb\" \n",
    "my_atoms = my_parser_pdb(target_file)\n",
    "print(f\"My parser found {len(my_atoms)} atoms.\")\n",
    "\n",
    "official_parser = PDB.PDBParser(QUIET=True)\n",
    "official_structure = official_parser.get_structure(\"protein_1\", target_file)\n",
    "\n",
    "official_atoms = list(official_structure.get_atoms())\n",
    "print(f\"Biopython found {len(official_atoms)} atoms.\")\n",
    "\n",
    "assert len(my_atoms) == len(official_atoms), \"Note: The number of atoms does not match\"\n",
    "print(\"Test passed successfully!\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2b0f47c2-2dca-45cb-96ee-53042bba7e53",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.17"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
