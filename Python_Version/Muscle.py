#!/usr/bin/python3

import MultiAlignment_Toolkit
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-in', '--input', dest='input',
                    help='Enter the input file path and unaligned fasta file with .fasta extension')
parser.add_argument('-msa_out', '--msa_out', dest='msa',
                    help='Enter the name of the .fasta output file -- default = ./MSA.fasta', default='./MSA.fasta')

args = parser.parse_args()

inFile = args.input
outFile = open(args.msa, 'w')
result = MultiAlignment_Toolkit.Muscle(inFile)
outFile.write(str(result[0]))

print("Completed MUSCLE Alignment!"
      "\nOutput file: " + str(args.msa) +
      "\nAlignment length = " + str(result[1]) +
      "\nNumber of sequences aligned = " + str(result[2]))

outFile.close()