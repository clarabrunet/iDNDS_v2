# iDNDS
Authors: Clara Brunet Villarejo, Miho Nakamura, Pere Puigbo

Year: 2024

Title: Indel-aware analysis of positive selection

**Description**:

This program calculates dNdS (modified from  Nei-Gojobori 1986) taking into account indels in the pairwise alignment. Indels are treated as non-synonymous substitutions in the pair-wise alignment.

**Reference**: 

Saeka Shimochi, Clara Brunet, Margalida Fontcuberta-Rigo, Katja Hrovat, Pere Puigbò, Miho Nakamura (2024) **Bone mechano-response is driven by transitions in vertebrate evolution**. In preparation. 

**Structure**:

- libs/: Contains the Python libraries.
- script/: Contains the main launcher script.


**Requirements**:

Python

Biopython

**Usage**:

The main launcher script takes a FASTA file as input containing two aligned codon sequences.

```bash
python iDNDS.py <fasta_file>
```
**Output**: 

- Nd: Total observed non-synonymous substitutions.
- Sd: Total observed synonymous substitutions.
- N: Total expected non-synonymous sites.
- S: Total expected synonymous sites.
- pn and ps: Proportion of observed non-synonymous and synonymous changes.
- dn and ds: Adjusted rates of non-synonymous and synonymous changes.
-    dN/dS: Positive selection ratio.
---
**No Resposibility Disclaimer**: 

The authors and contributors of this repository take no responsibility for any harm or damages that may arise from the use of this software. This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement. In no event shall the authors or contributors be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.
