#!/usr/bin/python3

"""

Each function does one specific duty, and each has its own imports:
# 1) Muscle performs multiple sequence alignment according to R. Edgar's algorithm: https://www.drive5.com/muscle/
     The output is a .fasta file with the multiple sequence alignment.
# 2) Consensus_RE takes a multiple sequence alignment (in .fasta format) and gives:
     a) A gapped consensus in .fasta format,
     b) A .csv file which shows the frequency of each residue (nucleotide or amino acid) at each position of the alignment,
        along with the consensus residue and its frequency.
     c) A plot of the consensus residue frequency at each position of the alignment.
     d) A .csv file of the relative entropy per strain per position and a sum of the relative entropy across strains
        at each position
     e) A plot of the Maximum relative entropy at each position and the cumulative sum of relative entropy at each position

"""








# Function to obtain multiple sequence alignment using MUSCLE
def Muscle(inFile):


    # Import necessary functions
    from Bio import SeqIO, AlignIO
    from Bio.Align.Applications import MuscleCommandline
    from io import StringIO


    # Perform multiple sequence alignment
    muscle_cline = MuscleCommandline(input=inFile)  # Perform MUSCLE taking the input file
    stdout, stderr = muscle_cline()  # Generate the standard output and error of the multiple alignment
    align = AlignIO.read(StringIO(stdout), "fasta")  # Read the multiple alignment file
    align_length = align.get_alignment_length()  # Get the length of the alignment
    align_count = 0  # Initialize the count of number of sequences as zero
    for seq_record in SeqIO.parse(inFile, "fasta"):  # Parse through all sequences to get a total number of sequences
        align_count += 1
    return stdout, align_length, align_count  # Return the multiple alignment, alignment length, and count of sequences








# Function to generate a consensus sequence as .fasta and .csv from a multiple sequence alignment in .fasta format
def Consensus(inFile, Limit, Type, Consensus):


    # Import necessary functions
    from Bio import SeqIO
    import Alignment_Check, Graphplot_Functions


    # Open the files that need to be written
    Consensus_Fasta = open((Consensus + '.fasta'), 'w')
    Consensus_CSV = open((Consensus + '.csv'), 'w')


    '''
    Call the Alignment_Check file, taking input file name, the limit for consensus, and type of sequence.
    Relevant values returned will be the consensus sequence, alignment length, position-specific residue count matrix,
    the residues associated with the sequence type, and the header of the consensus file
    '''

    result_Alignment = Alignment_Check.Alignment_Check(inFile, Limit, Type)

    # The returned values from Alignment_Check are assigned to specific variables
    consensus_all = result_Alignment[0]  # variable to assign consensus sequence
    align_length = result_Alignment[1]  # variable to store alignment length
    all_pssm = result_Alignment[2]  # variable to take position-specific residue count matrix
    residues = result_Alignment[3]  # variable to take residues for the sequence-type entered
    Consensus_CSV.write(result_Alignment[4])


    # Get a count of the overall number of sequences and also store the list of headers per organism sample.
    seq_count_all = 0
    for seq_record in SeqIO.parse(inFile, "fasta"):
        seq_count_all += 1


    '''
    Generate the rest of the consensus CSV file.
    Essentially, this for loop goes through each position, records the position number into the file, records the
    frequency of each residue at each position as a percent, including alignment gaps. It then identifies and records
    the maximum frequency at each position as well as the consensus residue.
    '''

    for i in range(align_length):
        Consensus_CSV.write(str(i+1) + ",")
        max_freq = 0.00  # The variable for the max frequency needs to be reset to zero for each position
        ntd_count = 0 # The variable for the sum of residue counts needs to be reset to zero for each position
                      # The sum of frequencies of residues per position is used to calculate the frequencies of gaps
        for n in residues:  # for loop that goes through the frequency of residues at a position
            frequency = round(float((all_pssm[i][n]/seq_count_all)*100), 2)  # Calculate residue frequency as a percent
            Consensus_CSV.write(str(frequency) + ",")  # Write the percent frequency of each residue
            ntd_count = ntd_count + all_pssm[i][n]  # Get the sum of counts of each residue
            if frequency > max_freq:  # Identify the maximum count of all residues
                max_freq = frequency
        gap_count = seq_count_all - ntd_count  # Determine the gap count
                                               # gap count = count of sequences - sum of residue counts
        frequency = round(float((gap_count/seq_count_all)*100), 2)  # Calculate gap frequency as a percent
        Consensus_CSV.write(str(frequency) + ",")  # Write the percent frequency of gaps
        if frequency > max_freq:  # Identify if gaps have the maximum count
            max_freq = frequency
        Consensus_CSV.write(str(consensus_all[i]) + "," + str(max_freq) + "\n")  # Write the consensus residue and frequency
    Consensus_CSV.close()  # Close the consensus csv file, because it is completed!


    # Write the fasta file of the consensus
    Consensus_Fasta.write(">Consensus_length_" + str(align_length) + "\n" + str(consensus_all))
    Consensus_Fasta.close()


    # Call the Consensus_plot function from the Graphplot_Functions file, taking the Consensus output file name
    Graphplot_Functions.Consensus_plot(Consensus)

    print("Consensus Generated!")












def Relative_Entropy(inFile, Limit, Type, Entropy):

    # Import necessary functions
    from Bio import SeqIO, AlignIO
    from Bio.Align import AlignInfo
    import Alignment_Check, Graphplot_Functions
    import pandas as pd
    import numpy as np
    import math, os


    result_Alignment = Alignment_Check.Alignment_Check(inFile, Limit, Type)


    align_length = result_Alignment[1]
    all_pssm = result_Alignment[2]
    residues = result_Alignment[3]


    RE_CSV = open((Entropy + '.csv'), 'w')


    # Get a count of the overall number of sequences and also store the list of headers per organism sample.
    Header_list = []
    seq_count_all = 0
    for seq_record in SeqIO.parse(inFile, "fasta"):
        Header_list.append(seq_record.id)
        seq_count_all += 1


    # Write the column headers of the relative entropy CSV file: "Header ID", and each residue position
    RE_CSV.write("Header_ID")
    for i in range(align_length):
        RE_CSV.write("," + str(i+1))
    RE_CSV.write("\n")


    '''
    Typically, fasta headers of sequences derived from NCBI are in the format "Genbank/RefSeq ID: Sequence-specific info".
    The following bit will take the header list in the multi-fasta file, separate the Genbank/RefSeq ID, and de-duplicate
    the list of headers - essentially, distill down to the organism/strain-specific ID, even if there are multiple copies
    of the gene. This is a critical step for checking relative entropy at each position of multi-copy genes. Also get the
    number of organisms/strains.  
    '''

    df1 = pd.DataFrame(Header_list, columns=['Header_ID'])
    df1[['Header_ID', 'Sep']] = df1['Header_ID'].str.split(':', expand=True)
    df2 = pd.DataFrame(df1.iloc[:,0])
    Org_DF = df2.drop_duplicates(keep='first').reset_index(drop=True)
    Org_count = len(Org_DF)



    '''
    For loop to go through each organism/strain-specific header, calculate the relative entropy at each position.
    We start with the aligned sequences of each organism/strain, which is stored in a temporary file iteratively.
    Then, determine the position-specific frequency at each position for the strain. Determine the Relative entropy of 
    a residue at a position in the gene for a organism/strain relative to all strains/organisms in the population list. 
    '''

    # Start the for loop to go through each organism/strain and derive their position-specific frequencies
    for each in Org_DF['Header_ID']:
        align_count_strain = 0
        temp_fasta = open("./temp.fasta", 'w')  # Create, open, re-create the temporary file of sequences per strain/organism
        for seq_record1 in SeqIO.parse(inFile, "fasta"):  # Parse through all sequences in the population list
            if each in seq_record1.id:  # Identify headers that contain the Genbank/RefSeq ID
                                        # This is specifically required for multi-allelic genes
                temp_fasta.write(">" + seq_record1.id + "\n" + str(seq_record1.seq) + "\n")  # Write their alignment sequences
                align_count_strain += 1  # Count the number of sequences per organism/strain
        temp_fasta.close()
        align_strain = AlignIO.read("./temp.fasta", "fasta")  # Read the organism/strain alignment file
        summary_align_strain = AlignInfo.SummaryInfo(align_strain)  # Store the summary of the alignment
        consensus_strain = summary_align_strain.gap_consensus()  # Obtain the gapped consensus
        if Type.lower() == "dna" or Type.lower() == "rna":  # Get the position-specific frequency depending on whether the sequence type is DNA/RNA/Prot
            strain_pssm = summary_align_strain.pos_specific_score_matrix(consensus_strain,
                                                                         chars_to_ignore=["Y", "R", "W", "S", "K", "M",
                                                                                          "D", "V", "H", "N", "B"])
        elif Type.lower() == "prot":
            strain_pssm = summary_align_strain.pos_specific_score_matrix(consensus_strain)
        RE_CSV.write(str(each))


        # Calculate the relative entropy of the organism/strain relative to all in the population
        for i in range(align_length):
            RE = 0.0  # Relative entropy set and reset to zero
            ntd_count_strain = 0  # Count of number of residues for organism/strain set and reset to zero
            ntd_count_all = 0   # Count of number of residues for all organisms/strains set and reset to zero
            for n in residues:  # Iterate through each residue
                if float(all_pssm[i][n]) != 0.0 and float(strain_pssm[i][n]) != 0.0:  # To avoid Zero division errors, the relative entropy is calculated for only residues present# To avoid Zero division errors, the relative entropy is calculated for only residues present


                    '''
                    The formula for relative entropy = summation(Pi*log(Pi/Qi)).
                    Pi = the frequency relative to total times a residue occurs at a position of the gene in a strain or organism (ranges from 0-1)
                    Qi = the frequency relative to total times a residue occurs at a position of the gene across all copies of strains or organisms (ranges from 0-1)
                    Start with Pi/Qi -> frequency_ratio.
                    Then cumulatively get Pi*log(frequency_ratio).
                    The log is to base 10 - this effects the level of amplification.
                    The final calculation is to four decimal places.
                    '''

                    frequency_ratio = float(((strain_pssm[i][n]) / align_count_strain) /
                                            ((all_pssm[i][n]) / seq_count_all))  # Pi/Qi
                    RE = round((RE + (((strain_pssm[i][n]) / align_count_strain) *
                                      (math.log(frequency_ratio, 10)))), 4)  # summation(Pi*log(Pi/Qi))
                    ntd_count_strain = ntd_count_strain + strain_pssm[i][n]  # Get the sum of counts of each residue for the organism/strain
                    ntd_count_all = ntd_count_all + all_pssm[i][n]  # Get the sum of counts of each residue for all organisms/strains

            if ntd_count_strain < align_count_strain and ntd_count_all < seq_count_all:  # Determine the relative entropy for gaps and add it to the summation
                frequency_ratio = float(((align_count_strain - ntd_count_strain) / align_count_strain) /
                                        ((seq_count_all - ntd_count_all) / seq_count_all))
                RE = round((RE + (((align_count_strain - ntd_count_strain) / align_count_strain) *
                                  (math.log(frequency_ratio, 10)))), 4)
            RE_CSV.write("," + str(RE))  # Write the summed relative entropy for the organism/strain at that position
        RE_CSV.write("\n")               # Write new line for the next organism/strain

    os.remove("./temp.fasta")  # Remove the temporary fasta file for the organism/strain
    RE_CSV.close()


    # Calculate the Max and Cumulative Relative Entropy at each position across all organisms/strains
    df_RE = pd.read_csv((Entropy + '.csv'), sep=',')
    df_RE = df_RE.set_index(df_RE.columns[0])
    Max_RE = []
    Cumulative_RE = []
    Position = []
    for i in range(align_length):
        Position.append(str(i+1))
        Max_RE.append(max(df_RE.iloc[:,i]))
        Cumulative_RE.append(sum(df_RE.iloc[:,i]))

    Graphplot_Functions.RE_plot(Entropy, Position, Max_RE, Cumulative_RE)
    print("Relative Entropy Generated!")


    '''
    Generate CSV files for theoretical Entropy and Cumulative Entropy for the population set.
    The theoretical data will show what the possible relative entropy and cumulative relative entropy values can be at
    a position, based on if 0-N strains had 0-n identical residue changes, where N = the total number of strains 
    in the population, and n = the average number of gene copies per strain. The calculation at a position is:
    Dkl = summation(Pi*log(Pi/Qi))
        = Pi*log(Pi/Qi) + Pi*log(Pi/Qi)
          ^^^^^^^^^^^^    ^^^^^^^^^^^^
          nucleotide 1    nucleotide 2
    cDkl = summation(Dkl at a position)
         = number of strains * Dkl at a position -> since I am only checking theoretical values for multiple strains
                                                    having the same number of allelic variations
    '''

    Avg_Copies = int(round((seq_count_all / Org_count), 0))  # Get the average number of gene copies per organism/strain
    T_Dkl_df = pd.DataFrame(np.empty(((Org_count+1), (Avg_Copies+1))))  # Create an empty dataframe for Theoretical Relative Entropy (Dkl)
    T_cDkl_df = pd.DataFrame(np.empty(((Org_count+1), (Avg_Copies+1))))  # Create an empty dataframe for Theoretical Cumulative Relative Entropy (cDkl)
    for j in range(Org_count + 1):  # Fill out the dataframes with possible values for Relative and cumulative entropy
        for i in range(Avg_Copies + 1):  # i counts through the number of alleles with a mutation, and j counts the number of strains
            if (i*j) == 0:  # i = 0 -> no mutations in the strain (uninformative). j = 0 -> no strains with mutation (uninformative)
                T_Dkl_df.iat[j, i] = 0.0000  # i=0 would give a math error (log 0), and since it is uninformative, Dkl = 0
            elif (Avg_Copies-i) == 0.0:  # To avoid a math error (log 0) in the second nucleotide, that
                T_Dkl_df.iat[j, i] = round(((i / Avg_Copies) * (math.log(((i / Avg_Copies) / ((i * j) / seq_count_all)), 10))), 4)
            else:
                T_Dkl_df.iat[j, i] = round((((i/Avg_Copies) * (math.log(((i/Avg_Copies)/((i*j)/seq_count_all)), 10))) +
                                           (((Avg_Copies-i)/Avg_Copies) * (math.log((((Avg_Copies-i)/Avg_Copies)/
                                                                                     ((seq_count_all-(i*j))/seq_count_all)),
                                                                                    10)))), 4)
            T_cDkl_df.iat[j, i] = j * T_Dkl_df.iloc[j, i]
    Org_count_Index = []  # Initialize an empty list for Dkl and cDkl dataframe rows -> The number of strains
    Allele_count_Headers = []  # Initialize an empty list for Dkl and cDkl dataframe headers -> The allelic variant count
    for i in range(Avg_Copies + 1):
        Allele_count_Headers.append('Allele_count=' + str(i))
    for j in range(Org_count + 1):
        Org_count_Index.append('Strain_count=' + str(j))
    T_Dkl_df.columns = Allele_count_Headers  # Set the column and index names and create the CSV files
    T_cDkl_df.columns = Allele_count_Headers
    T_Dkl_df.set_index(pd.Index(Org_count_Index), inplace=True)
    T_cDkl_df.set_index(pd.Index(Org_count_Index), inplace=True)
    T_Dkl_df.to_csv((str(Entropy) + "_Theoretical_Population.csv"), header=True, index=True)
    T_cDkl_df.to_csv((str(Entropy) + "_Theoretical_Cumulative_Population.csv"), header=True, index=True)

    print("Theoretical Entropy Files Generated!")
