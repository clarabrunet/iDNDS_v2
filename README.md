# iDNDS
Authors: Clara Brunet Villarejo, Miho Nakamura, Pere Puigbo

Year: 2024

Title: Indel-aware analysis of positive selection

**Description**:

This program calculates non-synonymous (dN) and synonymous (dS) substitution ratio (dN/dS), using an algorithem based on Nei-Gojobori (1986), that takes into account indels in the pairwise alignment. Indels are treated as non-synonymous substitutions in the pair-wise alignment.

**Reference**: 

Saeka Shimochi, Clara Brunet, Margalida Fontcuberta-Rigo, Katja Hrovat, Pere Puigbò, Miho Nakamura (2024) **Bone mechano-response is driven by transitions in vertebrate evolution**. In preparation. 

**Structure**:

- **src/** – Python source code  
  - `main.py` – Main calculation script  

- **scripts/** – Execution scripts  
  - `run_demo.sh` – Example demo launcher  

- **data/** – Input and output data  
  - `sample/example.fasta` – Demo FASTA input file  
  - `expected_output/` – Output folder  

- **LICENSE** – License file (CC BY 4.0)  
- **CITATION.cff** – Citation metadata    
- **requirements.txt** – Python dependencies  
- **README.md** – Documentation


## 4. Requirements & Installation

This software requires **Python 3.9+** and the following libraries:

- Biopython  
- NumPy  
- Pandas  
- Matplotlib  
- Scikit-learn  

To install all dependencies:

```bash
pip install -r requirements.txt

```
- Tested on Python 3.11.4 under Windows 10 and Ubuntu 22.04.
- Typical install time: less than 2 minutes on a standard desktop computer

**Usage**:

The main launcher script takes a FASTA file as input containing two aligned codon sequences.

```bash
python iDNDS.py <fasta_file>
```
The input FASTA file used in the demo is located at:
```bash
data/sample/example.fasta
```
Expected runtime: less than 5 seconds on a normal desktop computer.

To reproduce the demo results, run:
```bash
bash scripts/run_demo.sh
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
