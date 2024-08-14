#from mutations_analysis_amb_gaps import analyze_mutations_amb_gaps
from mutations_analysis_amb_gaps_NEW import analyze_mutations_amb_gaps_NEW
from mutation_1 import mutation_1
from mutation_2 import mutation_2
from mutation_3 import mutation_3
from proportion_side import proportion_side
from Bio import SeqIO

import sys

fasta_file = sys.argv[1] 
sequences = list(SeqIO.parse(fasta_file, "fasta"))

seq1 = str(sequences[0].seq)
seq2 = str(sequences[1].seq)
#resultat = analyze_mutations_amb_gaps(seq1, seq2, mutation_1, mutation_2, mutation_3, proportion_side)
resultat = analyze_mutations_amb_gaps_NEW(seq1, seq2, mutation_1, mutation_2, mutation_3, proportion_side)



