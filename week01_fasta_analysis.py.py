#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
p53 = SeqIO.read ("p53.fasta","fasta")
pdl1 = SeqIO.read ("pdl1.fasta","fasta")
prion = SeqIO.read ("prion.fasta","fasta")


# In[5]:


proteins = {
    'p53' : p53 , 
    'PD-L1' : pdl1 , 
    'PrP' : prion}

for name, protein in proteins.items():
    print(f"{name}: {len(protein.seq)} aa")


# In[6]:


from collections import Counter


# In[7]:


for name, record in proteins.items():
    seq = str(record.seq)
    composition = Counter(seq)
    print(f'\n=== {name} ===')
    for aa, n in sorted(composition.items()):
        percentage = (n / len(seq)) * 100
        print(f'{aa}: {n} ({percentage:.2f}%)')


# In[9]:


kd_scale = { 'A' : +1.8 , 'R' : -4.5 , 'N' : -3.5 , 'D' : -3.5 , 
             'C' : +2.5 , 'Q' : -3.5 , 'E' : -3.5 , 'G' : -0.4 , 
             'H' : -3.2 , 'I' : +4.5 , 'L' : +3.8 , 'K' : -3.9 ,
             'M' : +1.9 , 'F' : +2.8 , 'P' : -1.6 , 'S' : -0.8 , 
             'T' : -0.7 , 'W' : -0.9 , 'Y' : -1.3 , 'V' : +4.2 }


# In[13]:


for name , record in proteins.items():
    seq = str(record.seq)
    average = sum(kd_scale[aa] for aa in seq) / len(seq)
    print(f'{name}: {average:.3f}')


# In[ ]:
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
