---
layout: default
title: CarlMcCombe
description: PhD Candidate at The Australian National University
---

## Plotting SEC-MALS Data 

[Download the notebook](/assets/ipython_notebooks/MALS.ipynb) or read below. 

Aim: To produce easy to interpret plots detailing the most important findings from a MALS experiment. The following code uses pandas (to read the data) and matplotlib (to create the plots). Specifically this code has been designed for a MALS experiment looking to identify monomeric and homodimeric forms of the same protein.

### Imports

```python
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.lines as mlines
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
```

### Reading & Plotting Data
The cell below imports the data in the form of a csv file. In the csv file should include the elution volume in one column, the refractive index trace for the protein solution in another column, the elution volume corresponding to the region used for the mass estimate in another column, and the mass estimate in another column.
```python
data = pd.read_csv('data.csv')
Figure, axes = plt.subplots(1, 2, figsize =(5.5,3.5),sharey='row')
A = Figure.axes[0]
B = Figure.axes[1]
A.plot(data.Volume, data.A, color=(0,0,0), label='Refractive Index', linewidth = 2)
A.plot(data.Volume_B, data.B, color='tomato', label='Mass Estimate', linewidth = 4)
B.plot(data.Volume, data.C, color=(0,0,0), label='Refractive Index', linewidth = 2)
B.plot(data.Volume_C, data.D, color='tomato', label='Mass Estimate', linewidth = 4)
A.set_ylabel('Molecular Mass (Da)', fontsize=12)
Figure.text(0.5, 0, 'Volume (mL)', ha='center', fontsize=12)
A.set_title('Predicted Monomer', fontsize=12)
B.set_title('Predicted Homodimer', fontsize=12)
A.spines["top"].set_visible(False)
B.spines["top"].set_visible(False)
A.spines["right"].set_visible(False)
B.spines["right"].set_visible(False)
A.grid(which='major', axis='y')
B.grid(which='major', axis='y')
A.tick_params(axis='both', direction='in', labelsize=12, which='both')
B.tick_params(axis='both', direction='in', labelsize=12, which='both')
A.axhline(y=16619.75, xmin=0, xmax=1, linestyle = 'dashed', linewidth=4, alpha=0.75, color=(0.41, 0.67, 0.4))
B.axhline(y=16619.75, xmin=0, xmax=1, linestyle = 'dashed', linewidth=4, alpha=0.75, color=(0.41, 0.67, 0.4))
A.axhline(y=33239.5, xmin=0, xmax=1, linestyle = 'dashed', linewidth=4, alpha=0.75, color=(0.78, 0.60, 1))
B.axhline(y=33239.5, xmin=0, xmax=1, linestyle = 'dashed', linewidth=4, alpha=0.75, color=(0.78, 0.60, 1))
A.minorticks_on()
B.minorticks_on()
DimerLegend = mlines.Line2D([], [],color=(0.78, 0.60, 1), linewidth=4, linestyle = 'dashed', label='Dimer Mass')
MonomerLegend = mlines.Line2D([], [],color=(0.41, 0.67, 0.40), linewidth=4, linestyle = 'dashed', label='Monomer Mass')
MassLegend = mlines.Line2D([], [],color='tomato', linewidth=4, label='Mass Estimate')
RefractiveIndex = mlines.Line2D([], [],color=(0,0,0), linewidth=2, label='Refractive Index')
legend = B.legend(handles=[RefractiveIndex, MassLegend, MonomerLegend, DimerLegend],loc='center right', fontsize=12, borderaxespad=-11)
plt.savefig('savelocation/filename', dpi=600, bbox_inches='tight')
```
### Example plot 
I used a slightly modified version of this code to produce the SEC-MALS plot below. 


![SEC-MALS Plot](/assets/images/MALS.png)
