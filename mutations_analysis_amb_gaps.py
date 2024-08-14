import math


def analyze_mutations_amb_gaps(seq1, seq2, mutation_1, mutation_2, mutation_3, proportion_side):
   #Agrupar els codons de cada sequencia de de 3 en 3
    def split_sequence(seq):
        return [seq[i:i+3] for i in range(0, len(seq), 3)]
   #Comparar el nombre de diferencies de nucleotids entre el codo de la sequencia de referencia i el codo de la sequencia observada
    def compare_nuc(codon1, codon2):
        changes = 0
        for i in range(len(codon1)):
            if codon1[i] != codon2[i]:
                changes += 1
        return changes
    #Comparar nombre de canvis sinonims i no sinonims entre el codo de la sequencia de referencia i el codo de la sequencia observada en funcio del nombre de diferencies de nucleotids
    def compare_codons(seq1, seq2, mutation_1, mutation_2, mutation_3):
        total_N_1 = 0
        total_S_1 = 0
        total_N_2 = 0
        total_S_2 = 0
        total_N_3 = 0
        total_S_3 = 0
        Nd_gaps = 0  
        codons1 = split_sequence(seq1)
        codons2 = split_sequence(seq2)

        for codon1, codon2 in zip(codons1, codons2):
            changes = compare_nuc(codon1, codon2)
            if codon2 == '---':  
                Nd_gaps += 3 
            elif codon1 == '---':  
                Nd_gaps += 1 
            elif changes == 1:
                for sublist in mutation_1:
                    if codon1 == sublist[0] and codon2 == sublist[2]:
                        total_N_1 += sublist[-2]
                        total_S_1 += sublist[-1]
            elif changes == 2:
                for sublist in mutation_2:
                    if codon1 == sublist[0] and codon2 == sublist[5]:
                        total_N_2 += sublist[-2] / 2 
                        total_S_2 += sublist[-1] / 2   
            elif changes == 3:
                for sublist in mutation_3:
                    if codon1 == sublist[0] and codon2 == sublist[8]:
                        total_N_3 += sublist[-2] / 6  
                        total_S_3 += sublist[-1] / 6  

            # Imprimir valors per a cada parell de codons
            total_Nd = round(total_N_1 + total_N_2 + total_N_3 + Nd_gaps, 2)
            total_Sd = round(total_S_1 + total_S_2 + total_S_3, 2)
            print(f"Codon1: {codon1}, Codon2: {codon2}, Changes: {changes}")
            print(f"total_N_1: {round(total_N_1, 2)}, total_S_1: {round(total_S_1, 2)}")
            print(f"total_N_2: {round(total_N_2, 2)}, total_S_2: {round(total_S_2, 2)}")
            print(f"total_N_3: {round(total_N_3, 2)}, total_S_3: {round(total_S_3, 2)}")
            print(f"Nd_gaps: {round(Nd_gaps, 2)}")
            print(f"Total Nd: {total_Nd}, Total Sd: {total_Sd}")
            print("----------------------------------")

        return (total_N_1, total_S_1), (total_N_2, total_S_2), (total_N_3, total_S_3), Nd_gaps
    
    #Assignar valors de acnvis sinonims i no sinonims esperats per lloc del codo de la sequencia de referencia
    def calculate_reference_values(seq, proportion_side):
        total_N_ref = 0
        total_S_ref = 0
        N_gaps = 0  
        codons = split_sequence(seq)
        for codon in codons:
            if codon == '---':  
                N_gaps += 3
            else:
                for sublist in proportion_side:
                    if codon == sublist[0]:
                        total_N_ref += sublist[-2]
                        total_S_ref += sublist[-1]
        total_N = total_N_ref + N_gaps
        
        return round(total_N, 2), round(N_gaps, 2), round(total_S_ref, 2), round(total_N_ref, 2)
    #Calcul proporcio de canvis sinonims (canvis sinonims obsevats entre canvis sinonims esperats) i proporcio de canvis no sinonims (canvis no sinonims observats entre canvis no sinonims             esperats)
    def pn_ps(Nd, Sd, N, S):
        pn = Nd / N
        ps = Sd / S
        return round(pn, 2), round(ps, 2)

    #Calcul dn amb coreccio Jukes-Cantor i calcul ds amb correccio Jukes-Cantor (s assigna el valor inf quan no es possible el calcul de ds (valor de ds>3/4 i no es pot calcular el log))
    def dn_ds(pn, ps):
        try:
            dn = -3/4 * math.log(1 - ((4 * pn) / 3))
        except ValueError:
            dn = float('inf')
        try:
            ds = -3/4 * math.log(1 - ((4 * ps) / 3))
        except ValueError:
            ds = float('inf')
        return round(dn, 2), round(ds, 2)
    
    #Calcul del ratio dn/ds
    def dnds(dn, ds):
        if ds == 0 or ds == float('inf'):
            return float('inf')
        return round(dn / ds, 2)

    result_mutation_1, result_mutation_2, result_mutation_3, Nd_gaps = compare_codons(seq1, seq2, mutation_1, mutation_2, mutation_3)

    Nd = round(result_mutation_1[0] + result_mutation_2[0] + result_mutation_3[0] + Nd_gaps, 2)
    Sd = round(result_mutation_1[1] + result_mutation_2[1] + result_mutation_3[1], 2)

    N, N_gaps, S, N_ref= calculate_reference_values(seq1, proportion_side)

    pn, ps = pn_ps(Nd, Sd, N, S)

    dn, ds = dn_ds(pn, ps)

    dnds_result = dnds(dn, ds)

    # Print total values
    print(f"Total Nd: {Nd}, Total Sd: {Sd}, Total Nd_gaps: {N_gaps}, Total N: {N}, Total S: {S}")
    print(f"Total N_gaps: {N_gaps}, Total pn: {pn}, Total ps: {ps}")
    print(f"Total dn: {dn}, Total ds: {ds}, Total dnds: {dnds_result}")

    return Nd, Sd, N_gaps, N, S, N_gaps, N_ref,  pn, ps, dn, ds, dnds_result

