#!/usr/bin/python3

import MultiAlignment_Toolkit
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', '--input', dest='input',
                    help='Enter the .fasta multiple sequence alignment filename which used to generate consensus '
                         'and position specific frequencies')
parser.add_argument('-type', '--type', dest='type',
                    help='Enter whether the sequences are DNA, RNA, or Prot')
parser.add_argument('-limit', '--limit', dest='limit', type=float,
                    help='Enter frequency limit for nucleotide/amino acid to be considered being consensus '
                         '-- default = 0.51', default=0.51)
parser.add_argument('-out', '--output', dest='output',
                    help='Enter output file name without extension.'
                         ' Three files will be generated:'
                         ' (1) a .fasta file of the consensus sequence.'
                         ' (2) a .csv file with position specific frequencies and consensus residue frequency.'
                         ' (3) a .png file of a bar graph with consensus residue frequency at each position.'
                         ' -- default = ./Consensus', default='./Consensus')

args = parser.parse_args()

inFile = args.input
Limit = args.limit
Type = args.type
Consensus = args.output


# Call specific Consensus Function from MultiAlignment_Toolkit. Account for potential errors.
try:
    result = MultiAlignment_Toolkit.Consensus(inFile, Limit, Type, Consensus)
except TypeError:
    print("Error: Entered wrong Sequence Type. Enter 'DNA' or 'RNA', or 'Prot'")
except ZeroDivisionError:
    print("No sequences present in alignment file")
except FileNotFoundError:
    print("The multiple sequence alignment file is either absent or you have entered the filename incorrectly")
except UnboundLocalError:
    print("The sequence type was incompatible in yielding a position-specific frequency")
