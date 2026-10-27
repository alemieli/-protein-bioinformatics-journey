{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "a5bdfdd3-c909-4994-ac61-a9e7fcd21aca",
   "metadata": {},
   "outputs": [],
   "source": [
    "#1. Defining the two sequences to align\n",
    "seq1 = \"ADCNYRQCLCRPM\"\n",
    "seq2 = \"AYCYNRCKCRDP\"\n",
    "\n",
    "#2. Definition of simple scores\n",
    "match = 1\n",
    "mismatch = 0\n",
    "gap_penalty = -1\n",
    "\n",
    "#3. Let's set the dimensions of the matrix (+1 for the row/column of gaps)\n",
    "rows = len(seq1) + 1\n",
    "columns = len(seq2) + 1\n",
    "\n",
    "#4. Creating a matrix filled with zeros (0), this line creates a table of 8 lists, each containing 8 zeros\n",
    "matrix = [[0 for j in range(columns)] for i in range(rows)]\n",
    "\n",
    "#5. Initialization of the progressive penalty for gaps in the first column and the first row\n",
    "for i in range(rows):\n",
    "    matrix[i][0] = i * gap_penalty\n",
    "\n",
    "for j in range(columns):\n",
    "    matrix[0][j] = j * gap_penalty"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "a23454225dea345a",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "final score matrix:\n",
      "[0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12]\n",
      "[-1, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10]\n",
      "[-2, 0, 1, 0, -1, -2, -3, -4, -5, -6, -7, -7, -8]\n",
      "[-3, -1, 0, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7]\n",
      "[-4, -2, -1, 1, 2, 2, 1, 0, -1, -2, -3, -4, -5]\n",
      "[-5, -3, -1, 0, 2, 2, 2, 1, 0, -1, -2, -3, -4]\n",
      "[-6, -4, -2, -1, 1, 2, 3, 2, 1, 0, 0, -1, -2]\n",
      "[-7, -5, -3, -2, 0, 1, 2, 3, 2, 1, 0, 0, -1]\n",
      "[-8, -6, -4, -2, -1, 0, 1, 3, 3, 3, 2, 1, 0]\n",
      "[-9, -7, -5, -3, -2, -1, 0, 2, 3, 3, 3, 2, 1]\n",
      "[-10, -8, -6, -4, -3, -2, -1, 1, 2, 4, 3, 3, 2]\n",
      "[-11, -9, -7, -5, -4, -3, -1, 0, 1, 3, 5, 4, 3]\n",
      "[-12, -10, -8, -6, -5, -4, -2, -1, 0, 2, 4, 5, 5]\n",
      "[-13, -11, -9, -7, -6, -5, -3, -2, -1, 1, 3, 4, 5]\n"
     ]
    }
   ],
   "source": [
    "#6. Filling the Matrix (Dynamic Programming)\n",
    "for i in range(1, rows):\n",
    "    for j in range(1, columns):\n",
    "        \n",
    "        #If the amino acids match\n",
    "        if seq1[i-1] == seq2[j-1]:\n",
    "            diagonal_score = matrix[i-1][j-1] + match\n",
    "        else:\n",
    "            diagonal_score = matrix[i-1][j-1] + mismatch\n",
    "        \n",
    "        #Score for the indels\n",
    "        up_score = matrix[i-1][j] + gap_penalty\n",
    "        left_score = matrix[i][j-1] + gap_penalty\n",
    "\n",
    "        #Let's choose the maximum value among the three paths\n",
    "        matrix[i][j] = max(diagonal_score, up_score, left_score)\n",
    "        \n",
    "print(\"final score matrix:\")\n",
    "for row in matrix:\n",
    "    print(row)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "5e134c96b438c528",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Alignment Result\n",
      "ADCNYRQCLCRPM\n",
      "AYC-YNRCKCRDP\n"
     ]
    }
   ],
   "source": [
    "#7. Traceback: Backward Alignment Reconstruction\n",
    "i = len(seq1)\n",
    "j = len(seq2)\n",
    "\n",
    "alignment_seq1 = \"\"\n",
    "alignment_seq2 = \"\"\n",
    "\n",
    "while i > 0 or j > 0:\n",
    "    \n",
    "    #if it comes from diagonal\n",
    "    if i > 0 and j > 0:\n",
    "        if seq1[i-1] == seq2[j-1]:\n",
    "            diagonal_point = match\n",
    "        else:\n",
    "            diagonal_point = mismatch\n",
    "        if matrix[i][j] == matrix[i-1][j-1] + diagonal_point:\n",
    "            alignment_seq1 += seq1[i-1]\n",
    "            alignment_seq2 += seq2[j-1]\n",
    "            i -= 1\n",
    "            j -= 1\n",
    "            continue\n",
    "            \n",
    "    #if it comes from up\n",
    "    if i > 0 and matrix[i][j] == matrix[i-1][j] + gap_penalty:\n",
    "        alignment_seq1 += seq1[i-1]\n",
    "        alignment_seq2 += \"-\"\n",
    "        i -=  1\n",
    "        continue\n",
    "\n",
    "    #if it comes from left \n",
    "    if i > 0 and matrix[i][j] == matrix[i][j-1] + gap_penalty:\n",
    "        alignment_seq1 += \"-\"\n",
    "        alignment_seq2 += seq2[j-1]\n",
    "        j -=  1\n",
    "        continue\n",
    "\n",
    "alignment_seq1 = alignment_seq1[::-1]\n",
    "alignment_seq2 = alignment_seq2[::-1]\n",
    "\n",
    "print(\"\\nAlignment Result\")\n",
    "print(alignment_seq1)\n",
    "print(alignment_seq2)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "ed6bc0be-5998-4720-9544-5788f70c2022",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Real alignment score: 5\n",
      "Random scores mean: 1.58\n",
      "Random standard deviation: 1.02\n",
      "Final Z-score: 3.36\n"
     ]
    }
   ],
   "source": [
    "#8. Calculation Z-score\n",
    "import random\n",
    "import statistics\n",
    "\n",
    "#A function that implements the Needleman-Wunsch algorithm to return the final score\n",
    "def calculate_nw_score(s1, s2):\n",
    "    r = len(s1) + 1\n",
    "    c = len(s2) + 1\n",
    "    mat = [[0 for _ in range(c)] for _ in range(r)]\n",
    "\n",
    "    for i in range(r):\n",
    "        mat[i][0] = i * gap_penalty\n",
    "    for j in range(c):\n",
    "        mat[0][j] = j * gap_penalty\n",
    "\n",
    "    for i in range(1, r):\n",
    "        for j in range(1, c):\n",
    "            diag = mat[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch)\n",
    "            up = mat[i-1][j] + gap_penalty\n",
    "            left = mat[i][j-1] + gap_penalty\n",
    "            mat[i][j] = max(diag, up, left)\n",
    "    \n",
    "    return mat[-1][-1]\n",
    "\n",
    "real_score = calculate_nw_score(seq1, seq2)\n",
    "num_iterations = 100\n",
    "random_scores = []\n",
    "\n",
    "seq2_list = list(seq2)\n",
    "\n",
    "for _ in range(num_iterations):\n",
    "    random.shuffle(seq2_list)\n",
    "    seq2_random = \"\".join(seq2_list)\n",
    "\n",
    "    random_score = calculate_nw_score(seq1, seq2_random)\n",
    "    random_scores.append(random_score)\n",
    "\n",
    "random_mean = statistics.mean(random_scores)\n",
    "st_dev = statistics.stdev(random_scores) if len(set(random_scores)) > 1 else 1.0\n",
    "\n",
    "if st_dev == 0:\n",
    "    z_score = 0.0\n",
    "else:\n",
    "    z_score = (real_score - random_mean) / st_dev\n",
    "\n",
    "print(f\"\\nReal alignment score: {real_score}\")\n",
    "print(f\"Random scores mean: {random_mean:.2f}\")\n",
    "print(f\"Random standard deviation: {st_dev:.2f}\")\n",
    "print(f\"Final Z-score: {z_score:.2f}\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "ce271676-528c-4ed8-a69b-0d66a9593442",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Biopython score: 5.0\n",
      "\n",
      "Alignment visualization:\n",
      "target            0 ADCNYRQCLCRPM 13\n",
      "                  0 |.|..|-|.||.. 13\n",
      "query             0 AYCYNR-CKCRDP 12\n",
      "\n"
     ]
    }
   ],
   "source": [
    "#9. Code for comparison with Biopython\n",
    "from Bio.Align import PairwiseAligner\n",
    "\n",
    "aligner = PairwiseAligner()\n",
    "aligner.mode = 'global'\n",
    "\n",
    "aligner.match_score = 1\n",
    "aligner.mismatch_score = 0\n",
    "aligner.open_gap_score = -1\n",
    "aligner.extend_gap_score = -1\n",
    "\n",
    "alignments = aligner.align(seq1, seq2)\n",
    "\n",
    "best_alignment = alignments[0]\n",
    "\n",
    "print(f\"Biopython score: {best_alignment.score}\")\n",
    "print(\"\\nAlignment visualization:\")\n",
    "print(best_alignment)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1f6c46b9-c8a4-49e3-866c-4403854fa537",
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
