from Bio import SeqIO
import sys, os

def count_differences(seq1, seq2):
    nonsyn = 0
    syn = 0
    indels = 0

    codons1 = seq1.split('-')
    codons2 = seq2.split('-')

    for c1, c2 in zip(codons1, codons2):
        if len(c1) != 3 or len(c2) != 3:
            continue
        if c1 == c2:
            continue
        elif '-' in c1 or '-' in c2:
            indels += 1
        elif c1[0:2] == c2[0:2]:
            syn += 1
        else:
            nonsyn += 1

    return nonsyn + indels, syn

def main():
    input_file = "data/sample/example.fasta"
    output_dir = "data/expected_output"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "demo_result.txt")

    records = list(SeqIO.parse(input_file, "fasta"))
    if len(records) != 2:
        print("Error: FASTA must contain exactly two sequences.")
        sys.exit(1)

    seq1 = str(records[0].seq)
    seq2 = str(records[1].seq)

    Nd, Sd = count_differences(seq1, seq2)
    N, S = 100, 100
    dN = Nd / N
    dS = Sd / S if Sd != 0 else 1e-6
    ratio = dN / dS

    with open(output_file, "w") as f:
        f.write(f"Nd: {Nd}\nSd: {Sd}\nN: {N}\nS: {S}\n")
        f.write(f"dN: {dN:.3f}\ndS: {dS:.3f}\ndN/dS: {ratio:.3f}\n")

    print("Demo completed successfully.")
    print(f"Results written to: {output_file}")
    print(f"dN/dS ratio: {ratio:.3f}")

if __name__ == "__main__":
    main()
