#!/usr/bin/env python
# coding: utf-8

# In[14]:


#1. Defining the two sequences to align
seq1 = "ADCNYRQCLCRPM"
seq2 = "AYCYNRCKCRDP"

#2. Definition of simple scores
match = 1
mismatch = 0
gap_penalty = -1

#3. Let's set the dimensions of the matrix (+1 for the row/column of gaps)
rows = len(seq1) + 1
columns = len(seq2) + 1

#4. Creating a matrix filled with zeros (0), this line creates a table of 8 lists, each containing 8 zeros
matrix = [[0 for j in range(columns)] for i in range(rows)]

#5. Initialization of the progressive penalty for gaps in the first column and the first row
for i in range(rows):
    matrix[i][0] = i * gap_penalty

for j in range(columns):
    matrix[0][j] = j * gap_penalty


# In[15]:


#6. Filling the Matrix (Dynamic Programming)
for i in range(1, rows):
    for j in range(1, columns):

        #If the amino acids match
        if seq1[i-1] == seq2[j-1]:
            diagonal_score = matrix[i-1][j-1] + match
        else:
            diagonal_score = matrix[i-1][j-1] + mismatch

        #Score for the indels
        up_score = matrix[i-1][j] + gap_penalty
        left_score = matrix[i][j-1] + gap_penalty

        #Let's choose the maximum value among the three paths
        matrix[i][j] = max(diagonal_score, up_score, left_score)

print("final score matrix:")
for row in matrix:
    print(row)


# In[17]:


#7. Traceback: Backward Alignment Reconstruction
i = len(seq1)
j = len(seq2)

alignment_seq1 = ""
alignment_seq2 = ""

while i > 0 or j > 0:

    #if it comes from diagonal
    if i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            diagonal_point = match
        else:
            diagonal_point = mismatch
        if matrix[i][j] == matrix[i-1][j-1] + diagonal_point:
            alignment_seq1 += seq1[i-1]
            alignment_seq2 += seq2[j-1]
            i -= 1
            j -= 1
            continue

    #if it comes from up
    if i > 0 and matrix[i][j] == matrix[i-1][j] + gap_penalty:
        alignment_seq1 += seq1[i-1]
        alignment_seq2 += "-"
        i -=  1
        continue

    #if it comes from left 
    if i > 0 and matrix[i][j] == matrix[i][j-1] + gap_penalty:
        alignment_seq1 += "-"
        alignment_seq2 += seq2[j-1]
        j -=  1
        continue

alignment_seq1 = alignment_seq1[::-1]
alignment_seq2 = alignment_seq2[::-1]

print("\nAlignment Result")
print(alignment_seq1)
print(alignment_seq2)


# In[18]:


#8. Calculation Z-score
import random
import statistics

#A function that implements the Needleman-Wunsch algorithm to return the final score
def calculate_nw_score(s1, s2):
    r = len(s1) + 1
    c = len(s2) + 1
    mat = [[0 for _ in range(c)] for _ in range(r)]

    for i in range(r):
        mat[i][0] = i * gap_penalty
    for j in range(c):
        mat[0][j] = j * gap_penalty

    for i in range(1, r):
        for j in range(1, c):
            diag = mat[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch)
            up = mat[i-1][j] + gap_penalty
            left = mat[i][j-1] + gap_penalty
            mat[i][j] = max(diag, up, left)

    return mat[-1][-1]

real_score = calculate_nw_score(seq1, seq2)
num_iterations = 100
random_scores = []

seq2_list = list(seq2)

for _ in range(num_iterations):
    random.shuffle(seq2_list)
    seq2_random = "".join(seq2_list)

    random_score = calculate_nw_score(seq1, seq2_random)
    random_scores.append(random_score)

random_mean = statistics.mean(random_scores)
st_dev = statistics.stdev(random_scores) if len(set(random_scores)) > 1 else 1.0

if st_dev == 0:
    z_score = 0.0
else:
    z_score = (real_score - random_mean) / st_dev

print(f"\nReal alignment score: {real_score}")
print(f"Random scores mean: {random_mean:.2f}")
print(f"Random standard deviation: {st_dev:.2f}")
print(f"Final Z-score: {z_score:.2f}")


# In[19]:


#9. Code for comparison with Biopython
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = 'global'

aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -1
aligner.extend_gap_score = -1

alignments = aligner.align(seq1, seq2)

best_alignment = alignments[0]

print(f"Biopython score: {best_alignment.score}")
print("\nAlignment visualization:")
print(best_alignment)


# In[ ]:

 },
 "nbformat": 4,
 "nbformat_minor": 5
}
