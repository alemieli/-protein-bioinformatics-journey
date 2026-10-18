{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "3c820fc6-ba2f-4b9a-9af3-f766df682d51",
   "metadata": {},
   "outputs": [],
   "source": [
    "from Bio import SeqIO\n",
    "p53 = SeqIO.read (\"p53.fasta\",\"fasta\")\n",
    "pdl1 = SeqIO.read (\"pdl1.fasta\",\"fasta\")\n",
    "prion = SeqIO.read (\"prion.fasta\",\"fasta\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "af9ce93e-6ea5-46d5-9890-e4aaa275d703",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "p53: 393 aa\n",
      "PD-L1: 290 aa\n",
      "PrP: 253 aa\n"
     ]
    }
   ],
   "source": [
    "proteins = {\n",
    "    'p53' : p53 , \n",
    "    'PD-L1' : pdl1 , \n",
    "    'PrP' : prion}\n",
    "\n",
    "for name, protein in proteins.items():\n",
    "    print(f\"{name}: {len(protein.seq)} aa\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "37cfaa91-a5cb-4785-b8c1-50a58a39200b",
   "metadata": {},
   "outputs": [],
   "source": [
    "from collections import Counter"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "cc88307e-b93c-4e23-acf2-b03f68fc782c",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "=== p53 ===\n",
      "A: 24 (6.11%)\n",
      "C: 10 (2.54%)\n",
      "D: 20 (5.09%)\n",
      "E: 30 (7.63%)\n",
      "F: 11 (2.80%)\n",
      "G: 23 (5.85%)\n",
      "H: 13 (3.31%)\n",
      "I: 8 (2.04%)\n",
      "K: 20 (5.09%)\n",
      "L: 32 (8.14%)\n",
      "M: 12 (3.05%)\n",
      "N: 14 (3.56%)\n",
      "P: 43 (10.94%)\n",
      "Q: 15 (3.82%)\n",
      "R: 26 (6.62%)\n",
      "S: 39 (9.92%)\n",
      "T: 22 (5.60%)\n",
      "V: 18 (4.58%)\n",
      "W: 4 (1.02%)\n",
      "Y: 9 (2.29%)\n",
      "\n",
      "=== PD-L1 ===\n",
      "A: 16 (5.52%)\n",
      "C: 6 (2.07%)\n",
      "D: 14 (4.83%)\n",
      "E: 21 (7.24%)\n",
      "F: 11 (3.79%)\n",
      "G: 12 (4.14%)\n",
      "H: 9 (3.10%)\n",
      "I: 19 (6.55%)\n",
      "K: 19 (6.55%)\n",
      "L: 30 (10.34%)\n",
      "M: 7 (2.41%)\n",
      "N: 14 (4.83%)\n",
      "P: 10 (3.45%)\n",
      "Q: 12 (4.14%)\n",
      "R: 15 (5.17%)\n",
      "S: 13 (4.48%)\n",
      "T: 25 (8.62%)\n",
      "V: 23 (7.93%)\n",
      "W: 3 (1.03%)\n",
      "Y: 11 (3.79%)\n",
      "\n",
      "=== PrP ===\n",
      "A: 10 (3.95%)\n",
      "C: 4 (1.58%)\n",
      "D: 6 (2.37%)\n",
      "E: 9 (3.56%)\n",
      "F: 7 (2.77%)\n",
      "G: 45 (17.79%)\n",
      "H: 10 (3.95%)\n",
      "I: 9 (3.56%)\n",
      "K: 11 (4.35%)\n",
      "L: 12 (4.74%)\n",
      "M: 11 (4.35%)\n",
      "N: 12 (4.74%)\n",
      "P: 17 (6.72%)\n",
      "Q: 14 (5.53%)\n",
      "R: 11 (4.35%)\n",
      "S: 15 (5.93%)\n",
      "T: 13 (5.14%)\n",
      "V: 15 (5.93%)\n",
      "W: 9 (3.56%)\n",
      "Y: 13 (5.14%)\n"
     ]
    }
   ],
   "source": [
    "for name, record in proteins.items():\n",
    "    seq = str(record.seq)\n",
    "    composition = Counter(seq)\n",
    "    print(f'\\n=== {name} ===')\n",
    "    for aa, n in sorted(composition.items()):\n",
    "        percentage = (n / len(seq)) * 100\n",
    "        print(f'{aa}: {n} ({percentage:.2f}%)')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "7c22d3a6-4cbf-4ce6-b7b2-3ace0626b82f",
   "metadata": {},
   "outputs": [],
   "source": [
    "kd_scale = { 'A' : +1.8 , 'R' : -4.5 , 'N' : -3.5 , 'D' : -3.5 , \n",
    "             'C' : +2.5 , 'Q' : -3.5 , 'E' : -3.5 , 'G' : -0.4 , \n",
    "             'H' : -3.2 , 'I' : +4.5 , 'L' : +3.8 , 'K' : -3.9 ,\n",
    "             'M' : +1.9 , 'F' : +2.8 , 'P' : -1.6 , 'S' : -0.8 , \n",
    "             'T' : -0.7 , 'W' : -0.9 , 'Y' : -1.3 , 'V' : +4.2 }"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "9c3dff9b-b895-496c-bd8e-377a26ae1b27",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "p53: -0.758\n",
      "PD-L1: -0.226\n",
      "PrP: -0.560\n"
     ]
    }
   ],
   "source": [
    "for name , record in proteins.items():\n",
    "    seq = str(record.seq)\n",
    "    average = sum(kd_scale[aa] for aa in seq) / len(seq)\n",
    "    print(f'{name}: {average:.3f}')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "49345349-4929-4079-96cc-1d66314564d9",
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
