# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 19:20:37 2019

@author: manu

A program to generate all possible peptides for a given dna strand

translate_dna() would take the dna strand and directly form the 
peptide chain while translate_rna would translate from rna strand

"""

def read_fasta(fastafile):
    """
    Reads a fasta file and returns a dictionary with sequence
    number as keys and sequence code as values
    """
    sequences = []
    with open(fastafile, "r") as f:
        ls = f.readlines()
        for i in ls:
             sequences.append(i.rstrip("\n"))

    seq_id = []
    for i in sequences:
        if i[0] == ">":
            seq_id.append(i)

    seq_id_index = []
    for i in range(len(seq_id)):
        seq_id_index.append(sequences.index(seq_id[i]))

    seq_dic = {}
    for i in range(len(seq_id_index)):
        if i == (len(seq_id_index) - 1):
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:]
        else:
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:seq_id_index[i+1]]

    seq_dic_2 = {}
    for keys, values in seq_dic.items():
        seq_dic_2[keys] = "".join(values)

    return seq_dic_2
    
    
def complement_dna(dna):
        # Create dictionary of complementing nucleobase pairs
    comp_pairs = {"A" : "T", "T" : "A", "G" : "C", "C" : "G"}
    complementing_strand = ""
        # Generate the complementing strand
    for item in dna:        
        for i in range (len(item)-1, -1, -1):
            complementing_strand += comp_pairs[item[i]]

    return complementing_strand
    
    
def rna_strand(strand):
    rna = ""
        # Generate the RNA string
    for items in strand:
        for i in items:
            # Replace all occurrences of T with U
            if i == "T":
                rna += "U"
            else:
                rna += i

        return rna
    
    
def translate_dna(dnastring):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
    protein = []
    end = len(dnastring) - (len(dnastring) %3) - 1
    for i in range(0,end,3):
        codon = dnastring[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("N")
    return "".join(protein)
    

def translate_rna(rnastring):
    table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
        }
    protein = []
    end = len(rnastring) - (len(rnastring) %3) - 1
    for i in range(0,end,3):
        codon = rnastring[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("N")
    return "".join(protein)
    

# Generates the six possible frames per one sequence

def frame_id(seq):
    '''
    Usage: frame_id(dictionary['key'])
    frame_id(sequences['>Seq1'])
    six_frames = frame_id(sequences['>Seq1'])
    '''
    frames = {'0 skip forward':[],'1 skip forward':[],
              '2 skip forward':[],'0 skip reverse':[],
              '1 skip reverse':[],'2 skip reverse':[]}
    seq_rev = complement_dna(seq)
    for j in range(0,3):
        temp = ''.join([seq[j::]])
        temp_rev = ''.join([seq_rev[j::]])
        seq_trans = translate_dna(temp)
        seq_rev_trans = translate_dna(temp_rev)
        if j==0:
            frames['0 skip forward']=seq_trans
            frames['0 skip reverse']=seq_rev_trans
        if j==1:
            frames['1 skip forward']=seq_trans
            frames['1 skip reverse']=seq_rev_trans
        if j==2:
            frames['2 skip forward']=seq_trans
            frames['2 skip reverse']=seq_rev_trans

    return frames

# Generates all the frames for all the sequences

def gen_frames(dictionary):
    all_dict = {}
    for key, value in dictionary.items():
        all_dict[key] = frame_id(dictionary[key])

    return all_dict

# Find the open frames in the protein sequences

def oframe(amino):
    oframes = []
    for i in range(0,len(amino)):
        if amino[i]=='M':
            temp = ''.join([amino[i::]])
            oframe=temp[0:temp.find('_')+1]
            oframes.append(oframe)
    return oframes

## Finds the longest proteins in each sequence

def find_prots(dictionary):
    prots_dict = {}
    for key, value in dictionary.items():
        poss_protein = []
        for f in value:
            poss_protein += (oframe(value[f]))
            #print key, poss_protein
            c = 0
            result = ""
            for s in poss_protein:
                if len(s) > c:
                    result = s
                    c = len(s)
                else:
                    continue
            prots_dict[key] = result

    return prots_dict
   

def main():    
    
    print("Please provide path to the fasta file: \n")
    
    dna = read_fasta(input())       
    
    #dna = read_fasta("/home/manu/Downloads/sample_dna.txt")
    strands = []    
    for values in dna:
        strands.append(dna[values])
        
    print("Oringinal DNA string:", strands)
    forward_strand = rna_strand(strands)
    complement_dna(strands)
    reverse_strand = rna_strand(complement_dna(strands))
    translate_dna(strands)
    translate_rna(forward_strand)
    translate_dna(complement_dna(strands))
    translate_rna(reverse_strand)

    sequences_frames = gen_frames(dna)
    frames = []
    for items in sequences_frames.values():
        for item in items:
            frames.append(items)
            break
    
    protein = find_prots(sequences_frames)
    this = []
    for items in protein.values():
        for item in items:
            this.append(items)
            break
    
    print("The possible translations are:\n", frames)
    print("\n where the longest protein in given possible proteins is ", this)
    print("\n\nLegend:\n 0/1/2 skip forward: Frame 1/2/3 5' to 3' strand \n 0/1/2 skip reverese: Frame 1/2/3 3' to 5' strand")


if __name__ == '__main__':
    main()    
