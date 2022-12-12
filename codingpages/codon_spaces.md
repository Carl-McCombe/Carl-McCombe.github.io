---
layout: default
title: CarlMcCombe
description: PhD Candidate at The Australian National University
---

## Adding or Removing Spaces between Codons

[Download the notebook](/assets/ipython_notebooks/codon_spaces.ipynb) or read below. 

This notebook will take a nucleic acid sequence and ask whether the user wants to add sequences between codons in the sequence or remove them.

### Imports

This notebook requires the re library, which comes pre-installed if installing Python via [Anaconda](https://www.anaconda.com/products/individual). 

```python
import re
regex = re.compile("[^a-zA-Z]")
```

### Functions

The first function adds codon spaces whereas the second function removes codon spaces, the program will ask the user which function they want to use. 

```python
def addcodonspaces(sequence):
    sequence = regex.sub('', sequence)
    new_sequence = " ".join(sequence[i:i+3] for i in range(0, len(sequence), 3))
    j = len(sequence)
    print('Length of sequence is '+str(j)+' nucleotides')
    print('Here is the sequence with spaces between the codons: '+new_sequence)
    return

def removecodonspaces(sequence):
    new_sequence = re.sub(' ', '', sequence)
    j = len(sequence)
    print('Length of sequence is '+str(j)+' nucleotides')
    print('Here is the sequence without spaces between the codons: '+new_sequence)
    return
```

### User Input

The user must choose to either add or remove spaces and supply the sequence. The output from this cell is a print out of the new sequence that can be copied by the user. 

```python
choice = input("Enter 'a' if you want to add codon spaces or 'r' if you want to remove codon spaces (lowercase only):")
if choice == 'a':
    addcodonspaces(input("Enter your sequence to add codon spaces:"))
elif choice == 'r':
    removecodonspaces(input("Enter your sequence to remove codon spaces:"))
else:
    print("Incorrect input, type either 'a' to add codon spaces or 'r' to remove codon spaces from a sequence.")
    choice = input("Press 'a' if you want to add codon spaces or 'r' if you want to remove codon spaces (lowercase only):")
```
