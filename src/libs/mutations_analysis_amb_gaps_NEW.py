import math

def analyze_mutations_amb_gaps_NEW(seq1, seq2, mutation_1, mutation_2, mutation_3, proportion_side):
    def split_sequence(seq):
        return [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    def compare_nuc(codon1, codon2):
        changes = 0
        for i in range(len(codon1)):
            if codon1[i] != codon2[i]:
                changes += 1
        return changes
    
    def compare_codons(seq1, seq2, mutation_1, mutation_2, mutation_3):
        total_N_1 = 0
        total_S_1 = 0
        total_N_2 = 0
        total_S_2 = 0
        total_N_3 = 0
        total_S_3 = 0
        Nd_gaps = 0  
        small_value = 1e-7  
        
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
            
        total_Nd = total_N_1 + total_N_2 + total_N_3 + Nd_gaps
        total_Sd = total_S_1 + total_S_2 + total_S_3
        
        if total_Nd == 0:
            total_Nd += small_value
        if total_Sd == 0:
            total_Sd += small_value

        return (total_Nd, total_Sd)
    
    def calculate_reference_values(seq, proportion_side):
        total_N_ref = 0
        total_S_ref = 0
        N_gaps = 0  
        small_value = 1e-7  
        
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
        
        if total_N == 0:
            total_N += small_value
        if total_S_ref == 0:
            total_S_ref += small_value

        return total_N, N_gaps, total_S_ref, total_N_ref
    
    def pn_ps(Nd, Sd, N, S):
        small_value = 1e-6
        Nd += small_value if Nd == 0 else 0
        Sd += small_value if Sd == 0 else 0
        N += small_value if N == 0 else 0
        S += small_value if S == 0 else 0
        pn = Nd / N
        ps = Sd / S
        if ps >= 0.75:
            ps = 0.74999999999  
        if pn >= 0.75:
            pn = 0.74999999999
        print(f"Debug -> Nd: {Nd:.8f}, Sd: {Sd:.8f}, N: {N:.8f}, S: {S:.8f}, pn: {pn:.8f}, ps: {ps:.8f}")
        return pn, ps
    
    def dn_ds(pn, ps):
        try:
            dn = -3/4 * math.log(1 - ((4 * pn) / 3))
        except ValueError:
            dn = float('inf')
        try:
            ds = -3/4 * math.log(1 - ((4 * ps) / 3))
        except ValueError:
            ds = float('inf')
        return dn, ds
    
    def dnds(dn, ds):
        if ds == 0 or ds == float('inf'):
            return float('inf')
        return dn / ds
    
    total_Nd, total_Sd = compare_codons(seq1, seq2, mutation_1, mutation_2, mutation_3)
    total_N, N_gaps, total_S_ref, total_N_ref = calculate_reference_values(seq1, proportion_side)
    pn, ps = pn_ps(total_Nd, total_Sd, total_N, total_S_ref)
    dn, ds = dn_ds(pn, ps)
    dnds_ratio = dnds(dn, ds)
    
    print(f"Final Values -> Total Nd: {total_Nd:.8f}, Total Sd: {total_Sd:.8f}")
    print(f"pn: {pn:.8f}, ps: {ps:.8f}")
    print(f"dn: {dn:.8f}, ds: {ds:.8f}")
    print(f"dN/dS ratio: {dnds_ratio:.8f}")
    
    return pn, ps, dn, ds, dnds_ratio
