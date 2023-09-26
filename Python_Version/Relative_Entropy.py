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
parser.add_argument('-entropy', '--entropy', dest='entropy',
                    help='Enter the output file name for relative entropy.'
                         ' Four files will be generated:'
                         ' (1) a .csv file of the relative entropy for each organism/strain '
                         '     at each position of the multiple alignment.'
                         ' (2) a .png file showing the maximum position-specific relative entropy across organisms/strains (red line - primary y-axis)'
                         '     and the summation of relative entropy at each position across organisms/strains (grey bar - secondary y-axis).'
                         ' (3) a .csv file of theoretical relative entropy values for 0-N strains having identical number of 0-n variations.'
                         ' (4) a .csv file of theoretical cumulative relative entropy values for 0-N strains having 0-n identical variations.'
                         '-- default=./Relative_Entropy', default='./Relative_Entropy')

args = parser.parse_args()

inFile = args.input
Limit = args.limit
Type = args.type
Entropy = args.entropy


# Call specific Relative Entropy Function from MultiAlignment_Toolkit. Account for potential errors.
try:
    result = MultiAlignment_Toolkit.Relative_Entropy(inFile, Limit, Type, Entropy)
except TypeError:
    print("Error: Entered wrong Sequence Type. Enter 'DNA' or 'RNA', or 'Prot'")
except ZeroDivisionError:
    print("No sequences present in alignment file")
except FileNotFoundError:
    print("The multiple sequence alignment file is either absent or you have entered the filename incorrectly")
except UnboundLocalError:
    print("The sequence type was incompatible in yielding a position-specific frequency")
