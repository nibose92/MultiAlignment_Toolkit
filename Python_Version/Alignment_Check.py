#!/usr/bin/python3


'''
Store residue list depending on type of sequence input.
I avoid nucleotide nomenclatures Y,R,W,S,K,M,D,V,H,N,B for nucleotides because they are ultimately uninformative for
strain/organism level calculations -- If it's not A/C/G/T or an insertion/deletion (shown as a gap '-' in some sequences),
the sequence information is ambiguous.
'''

from Bio import AlignIO
from Bio.Align import AlignInfo


def Alignment_Check(inFile, Limit, Type):

    # Obtain Consensus Sequence from input multiple sequence alignment
    align_all = AlignIO.read(inFile, "fasta")                                                                     # Read input pre-aligned file
    summary_align = AlignInfo.SummaryInfo(align_all)                                                                    # Store information regarding the alignment
    consensus_all = summary_align.gap_consensus(threshold=Limit)                                                        # Store the gapped consensus in a variable
    align_length = align_all.get_alignment_length()

    match Type.lower():
        case "dna" | "rna":
            all_pssm = summary_align.pos_specific_score_matrix(consensus_all,
                                                               chars_to_ignore=["Y", "R", "W", "S", "K", "M", "D", "V",
                                                                                "H", "N", "B"])
            if Type.lower() == "dna":
                residues = ['A', 'C', 'G', 'T']
                Consensus_Header = "Alignment_Position,A%,C%,G%,T%,-%,Consensus_Residue,Consensus_Frequency%\n"
            elif Type.lower() == "rna":
                residues = ['A', 'C', 'G', 'U']
                Consensus_Header = "Alignment_Position,A%,C%,G%,U%,-%,Consensus_Residue,Consensus_Frequency%\n"
        case "prot":
            all_pssm = summary_align.pos_specific_score_matrix(consensus_all)
            residues = ['G', 'A', 'V', 'C',
                        'P', 'L', 'I', 'M',
                        'N', 'W', 'F', 'S',
                        'T', 'Y', 'N', 'Q',
                        'H', 'R', 'K', 'D', 'E']
            Consensus_Header = "Alignment_Position,G%,A%,V%,C%,P%,L%,I%,M%,N%,W%,F%,S%,T%,Y%,N%,Q%,H%,R%,K%,D%,E%,-%,Consensus_Residue,Consensus_Frequency%\n"
        case _:
            raise TypeError()
    return consensus_all, align_length, all_pssm, residues, Consensus_Header