<br>

# Multiple Sequence Alignment Toolkit for Python or R: MUSCLE, Consensus, and Relative Entropy
Tools to perform MUSCLE alignment (in fasta format), generate consensus (in fasta and csv format), and an updated version of Relative Entropy Estimation. These tools work for DNA, RNA, and Protein sequences. Aside from the Relative Entropy tool which is currently only available for Python, all other tools are available in Python and R.\
The tools are not to help Bioinformaticians with multi-alignment analyses, but also to train those delving into Biopython and R for Bioinformatics. The codes have extensive comments to help understand processes.
<br>
<br>
**Author:** Nikhil Bose, Ph.D.\
**Email:** nikhil.bose@ucf.edu

<br>
<br>

# Table of Contents:
1) [General Description](#1-general-description)
2) [Requirements](#2-requirements)
3) [Installation](#3-installation)
4) [Usage (Python Version)](#4-usage-python-version)
5) [Usage (R Version)](#5-usage-r-version) 
6) [Relative Entropy Tool Explanation](#6-relative-entropy-tool-explanation)
7) [Credit Attribution](#7-credit-attribution)
8) [License](#8-license)

<br>
<br>

## 1) General Description
The tools provided perform various functions with multiple sequence alignments.

### a) Muscle alignment (Python or R):
Takes a multiple sequence file in .fasta format, performs multiple sequence alignment using R or Biopython according to R. Edgar's algorithm, and gives the output in .fasta format.
### b) Consensus generation (Python or R)
Takes a multiple sequence alignment in .fasta format and generates a consensus sequence in .fasta format, a .csv file showing the residue and gap frequencies as percents as well as the consensus frequency, and a .png file showing a bar graph of the consensus residue frequencies.
### c) Relative Entropy (RE) estimation (Python ONLY)
Takes a multiple sequence alignment in .fasta format and generates a .csv file showing the relative entropy at each residue position for each organism/strain (see [Relative Entropy Tool Explanation](#6-relative-entropy-tool-explanation)). A graph showing the Maximum RE and Cumulative RE is also generated along with two additional .csv files which show the theoretical values for RE and cumulative RE for the given population.

**NOTE:** If you use the relative entropy code, please cite the paper:\
Bose N, Moore SD. Variable Region Sequences Influence 16S rRNA Performance. Microbiol Spectr. 2023 Jun 15;11(3):e0125223. [doi: 10.1128/spectrum.01252-23](https://journals.asm.org/doi/10.1128/spectrum.01252-23). Epub 2023 May 22. PMID: [37212673](https://pubmed.ncbi.nlm.nih.gov/37212673/); PMCID: [PMC10269663](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10269663/).

<br>
<br>

## 2) Requirements
Depending on the version you want to use, install the requirements for Python or R versions. The requirements and links on how to install them are as follows:

**Python Version:** \
&nbsp; &nbsp; &nbsp; a) [Python (version 3+)](https://www.python.org/downloads/) \
&nbsp; &nbsp; &nbsp; b) [pip (version 3+)](https://pip.pypa.io/en/stable/installation/) \
&nbsp; &nbsp; &nbsp; c) [Anaconda](https://www.anaconda.com/download) \
&nbsp; &nbsp; &nbsp; d) Biopython         (via conda https://anaconda.org/conda-forge/biopython) \
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; (via pip https://biopython.org/wiki/Download) \
&nbsp; &nbsp; &nbsp; e) pandas (comes with anaconda)\
&nbsp; &nbsp; &nbsp; f) numpy (comes with anaconda)\
&nbsp; &nbsp; &nbsp; g) matplotlib (comes with anaconda)\
\
**R Versions:** \
&nbsp; &nbsp; &nbsp; a) [Rstudio (with R base version 3+)](<https://posit.co/download/rstudio-desktop/>)\
&nbsp; &nbsp; &nbsp; b) [Bioconductor package in R](<https://www.bioconductor.org/install/>)\
&nbsp; &nbsp; &nbsp; c) [Biostrings package in R](<https://bioconductor.org/packages/release/bioc/html/Biostrings.html>)\
&nbsp; &nbsp; &nbsp; d) [muscle package in R](<https://bioconductor.org/packages/release/bioc/html/muscle.html>)\
&nbsp; &nbsp; &nbsp; e) [tidyverse (which includes ggplot2)](<https://ggplot2.tidyverse.org/>)

<br>
<br>


## 3) Installation

**Python:** There are six .py files, three that take user inputs and direct them towards the other three computational functions. You need to downloaded all 6 to perform all functions. 

**R:** There are 2 .R files, each perform their respective functions - Muscle\_R.R performs MUSCLE Alignment and Consensus_R.R performs consensus generation. 

**Easy download and install:** Use your mouse/trackpad to download the code, unzip folders ("Python\_Version" or "R\_version" or both) into a folder of your liking. Open terminal in the folder containing python files to run python scripts, and open the .R files in Rstudio to run them.

<br>
<br>

## 4) Usage (Python Version)

### a) Muscle.py:
Takes user input .fasta file of *unaligned* sequences, performs MUSCLE alignment, and generates a .fasta Alignment file.\
Requires **i)** Muscle.py and **ii)** MultiAlignment_Toolkit.py\
\
**Usage in Terminal** (replace /path/to/file with your file path and file name. Both input and output must have .fasta extension): <br>
$ `python3 /path/to/Muscle.py -in /path/to/file.fasta -MSA_out /path/to/output.fasta`

**View help in terminal:** \
$ `python3 file/to/Muscle.py -h`

> usage: Muscle.py [-h] [-in INPUT] [-msa_out MSA]\
> options:\
> -h, --help \
> show this help message and exit\
> -in INPUT, --input INPUT
> Enter the input file path and unaligned fasta file with .fasta extension\
>  -msa_out MSA, --msa_out MSA\
> Enter the name of the .fasta output file -- default = ./MSA.fasta 

<br>

### b) Consensus.py:
Takes user input .fasta file of ALIGNED sequences and generates a consensus.fasta file, a consensus.csv file, and a consensus residue frequency bar plot image as .png.\
Requires **i)** Consensus.py, **ii)** MultiAlignment\_Toolkit.py, **iii)** Alignment\_Check.py, and **iv)** Graphplot_Functions.py\
\
**Usage in Terminal** (replace /path/to/file with your file path and file name. Only input must have .fasta extension):\
$ `python3 path/to/Consensus.py -in /path/to/MSA_file.fasta -type DNA -limit 0.51 -out /path/to/Consensus`

**View help in terminal:** \
$ `python3 file/to/Consensus.py -h`
> usage: Consensus.py [-h] [-in INPUT] [-type TYPE] [-limit LIMIT] [-out OUTPUT]\
> options:\
> -h, --help\
> show this help message and exit\
> -in INPUT, --input INPUT\
> Enter the .fasta multiple alignment filename which used to generate consensus and position specific frequencies\
> -type TYPE, --type TYPE\
> Enter whether the sequences are DNA, RNA, or Prot\
> -limit LIMIT, --limit LIMIT\
Enter frequency limit for nucleotide/amino acid to be considered being consensus -- default = 0.51\
> -out OUTPUT, --output OUTPUT\
> Enter output file name without extension. Three files will be generated: (1) a .fasta file of the consensus sequence. (2) a .csv file with position specific frequencies and consensus residue frequency.
(3) a .png file of a bar graph with consensus residue frequency at each position. -- default = ./Consensus

<br>

### c) Relative_Entropy.py:
Takes user input .fasta file of ALIGNED sequences and generates a consensus.fasta file, a consensus.csv file, and a consensus residue frequency bar plot image as .png.\
Requires **i)** Relative_Entropy.py, **ii)** MultiAlignment\_Toolkit.py, **iii)** Alignment\_Check.py, and **iv)** Graphplot_Functions.py\
\
**NOTE:** In order for this program to work, the fasta headers should be "Organism\_Identifier:Sequence\_details". The **colon (:) is vital** in order for the program to identify multiple sequences that belong to the same organism/strain. This nomenclature for fasta headers is typically followed for fasta sequences downloaded from NCBI. However, some headers may have colons (:) within the organism identifier, in which case those headers need to be relabeled. If you want to cluster sequences of multiple organisms/strains under one, ensure they have the same "Organism_Identifier" before the colon. Details of how the tool works are under [Relative Entropy Tool Explanation](#6-relative-entropy-tool-explanation) section.\
\
**Usage in Terminal** (replace /path/to/file with your file path and file name. Only input must have .fasta extension): \
$ `python3 path/to/Relative_Entropy.py -in /path/to/MSA_file.fasta -type DNA -limit 0.51 -entropy /path/to/Relative_Entropy` 

**View help in terminal:** \
$ `python3 file/to/Relative_Entropy.py -h`
> usage: Relative_Entropy.py [-h] [-in INPUT] [-type TYPE] [-limit LIMIT] [-entropy ENTROPY]\
> options:\
> -h, --help\
> show this help message and exit\
> -in INPUT, --input INPUT\
> Enter the .fasta multiple sequence alignment filename which used to generate consensus and position specific frequencies\
> -type TYPE, --type TYPE\
> Enter whether the sequences are DNA, RNA, or Prot\
> -limit LIMIT, --limit LIMIT\
> Enter frequency limit for nucleotide/amino acid to be considered being consensus -- default = 0.51\
> -entropy ENTROPY, --entropy ENTROPY\
> Enter the output file name for relative entropy. Four files will be generated: (1) a .csv file of the relative entropy for each organism/strain at each position of the multiple alignment. (2) a .png file showing the maximum position-specific relative entropy across organisms/strains (red line - primary y-axis) and the summation of relative entropy at each position across organisms/strains (grey bar - secondary y-axis). (3) a .csv file of theoretical relative entropy values for 0-N strains having identical number of 0-n variations. (4) a .csv file of theoretical cumulative relative entropy values for 0-N strains having 0-n identical variations.-- default=./Relative_Entropy 

<br>

### d) MultiAlignment_Toolkit.py:
Performs computations in functions, depending on which function is called from input python files.

<br>

### e) Alignment_Check.py:
Checks if input alignments are DNA/RNA/Protein and returns values required for computations in MultiAlignment_Toolkit.py. 

<br>

### f) Graphplot_Functions.py:
Takes Computations from the MultiAlignment_Toolkit and creates graphs.

<br>
<br>

## 5) Usage (R Version)

### a) Muscle_R.R:
Takes user input .fasta file of UNALIGNED sequences, performs MUSCLE alignment, and generates a .fasta Alignment file.\
**Step 1:** Make sure the required packages are installed in Rstudio\
**Step 2:** Open the Muscle_R.R file in R studio\
**Step 3:** Edit the inFile, SeqType, and outFile assignments within quotes

        inFile <- "./file/to/path/Unaligned_infile.fasta"       # Change "/file/to/path" and "Unaligned_infile.fasta" to your input file path and file name with .fasta extension
        SeqType <- "DNA"                                        # ONLY POSSIBLE ENTRIES for SeqType are "DNA", "RNA", or "Prot" (not case-specific and within quotes)
        outFile <- "./file/to/path/Aligned_outfile.fasta"       # Change "/file/to/path" and "Aligned_outfile.fasta" to your output file path and file name with .fasta extension

**Step 4:** Run code
> "Completed MUSCLE Alignment!"

<br>

### b) Consensus:
Takes user input .fasta file of ALIGNED sequences and generates a consensus.fasta file, a consensus.csv file, and a consensus residue frequency bar plot image as .png.\
**Step 1:** Make sure the required packages are installed in Rstudio \
**Step 2:** Open the Muscle_R.R file in R studio \
**Step 3:** Edit the infile, SeqType, and outFile assignments within quotes


        inFile <- "./file/to/path/Aligned_infile.fasta"        # Change "/file/to/path" and "Aligned_outfile.fasta" to your output file path and file name with .fasta extension)
        SeqType <-  "DNA"                                      # ONLY POSSIBLE ENTRIES ARE DNA, RNA, or Prot - not case-specific
        Limit <- 0.51                                          # Change value to any decimal between 0 and 1 that you consider as a frequency for a residue being a consensus
        outFile <- "./file/to/path/Consensus_outfile"          # Change "/file/to/path" and "Consensus_outfile" to your output file path and file name WITHOUT EXTENSION. Extensions will be added by the code for the fasta and graph image
**Step 4:** Run code 
> "Consensus Generated!"

<br>
<br>

## 6) Relative Entropy Tool Explanation

This tool is an application of Kullback-Leibler Divergence (*D<sub>KL*) (https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) to provide values at nucleotide positions of your sequences of interest where a rare residue variations in DNA/RNA/peptide sequences reside for an individual relative to that of the population under study. The tool is even more informative if individuals in the population have multiple copies of the gene. Therefore, if the individual has a mutation in a specific position of the sequence among most or all copies of the gene, which no or few other individuals have, then they have evolved signficantly - the same residue changed in genomic copies located at different spots in the genome. This would not be evident from a cursory view of the consensus residue frequency of the population. The reason I wrote this code was to identify rare residue variations that occured across 16 gene copies in a strain or species, relative to a genus (or higher taxonomies). I identified that while 16S gene sequences in strains of the genera *Escherichia* and *Shigella* are over 99% identical, *Shigella dysenteriae*, *Shigella boydii*, and *Escherichia albertii* (pathogenic bacteria) had respective mutations that distinguished them from any other species in the genera (I looked at over 12,000 *Escherichia* and *Shigella* sequences across over 1800 strains). This can be helpful, because typical 16S gene analysis rarely provides species-level identification, especially for Escherichia and Shigella which are typically clumped under a single genus classification. By using the Relative Entropy tool, I compared the nucleotide-to-nucleotide differences in the genes, identified rare strain- and species-specific mutations, and then traced which strains or species they belong to.


The way it works is that it takes a .fasta file containing a multiple alignment of the sequences of interest of a population, separates the sequences belonging to each individual and does a alignment position comparison of the frequency of a residue within the individual relative to the entire population. The formula is:

*D<sub>KL* =  (*Pi * log(Pi/ Qi)*), where

*Pi* = the residue frequency at a position of the gene within the organism or strain (a fraction that gives a decimal value from 0-1; for Escherichia and Shigella 16S genes, they ranged from 1/7 to 7/7) \
*Qi* = the residue frequency at a position of the gene across organisms or strains (a fraction that gives a decimal value from 0-1; for Escherichia and Shigella 16S genes, they ranged from ~1/12000 to ~12000/12000)\
Log is to base 10 \
*summation* implies that *Pi * log(Pi/Qi)* is summed across all observed residues at the position. So if the sequence type is DNA, *Pi * log(Pi/Qi)* at a position is calculated for 'A', 'C', 'G', 'T', and gap ('-') occurrences and summed up.

<br>

**For example**, Say a gene occurs at 10 copies per individual and there are 100 individuals in the population. The total number of gene sequences in the population is therefore 10*100 = 1000

**Case 1:** One individual has a residue change in one gene copy that no other individual in the population has. All other gene copies in the population are identical. \
*Pi* = 1/10 for mutated residue, and 9/10 for the unmutated residue\
*Qi* = 1/1000 for the mutated residue, and 999/1000 for unmutated

*D<sub>KL* at the position = (1/10 * log((1/10)/(1/1000))) + (9/10 * log((9/10)/(999/1000))) \
= 0.2 + (-0.0407) \
= 0.159 &ensp;&ensp;&ensp;&ensp; <- Lower value, hence less informative for organism identification

**Case 2:** One individual has a residue change in all gene copies that no other individual in the population has. All others in the population have identical gene sequences. \
*Pi* = 10/10 for mutated residue, and 0/10 for the unmutated residue \
*Qi* = 10/1000 for the mutated residue, and 990/1000 for unmutated

*D<sub>KL* at the position = (10/10 * log((10/10)/(10/1000))) + (0/10 * log((0/10)/(990/1000)))\
= 2 + (0)**
<br>
= 2 &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;<- Higher value, hence more informative for organism identification

**The second term is set to 0 because logically it is pointless to compute entropy for another residue in the organism when it is known to be absent (https://dsp.stackexchange.com/questions/74057/value-of-0-log0-in-entropy-formula).


Cumulative relative entropy (*cD<sub>KL*) = summation(*D<sub>KL*) at a position across all organisms/strains in the population. So if a residue change is present across multiple individuals, it becomes more informative in the population level (whether the population is all strains in a species or genus). HOWEVER, NOTE that the more inviduals with the position-specific residue variation, the lower the individual *Dkl* (becomes less individual-specific).

The values for *Dkl* and *cDkl* are unbounded (no max and no min, though 0 or negative value would imply uninformative). Also, the values generated cannot generally be normalized because different genes occur in different copy numbers. Therefore, the values completely depend on the number of copy numbers and individuals in each population set (the larger the population, the larger the *D<sub>KL* and *cD<sub>KL* value).

To add importance to values, the Relative Entropy code also generates .csv files showing theoretical *D<sub>KL* and *cD<sub>KL* values for your population set. These files show putative values for a residue mutation in 0-n gene copies for 0-N individuals, where n is the average number of gene copies per individual in your population (rounded to an integer) and N is the total number of individuals in your population. From the observed values for Relative Entropy (either from the .csv or .png file), compare with those of the theoretical values, which will tell you the approximate range of individuals and their gene copies (alleles) with mutations.

<br>
<br>

## 7) Credit Attribution
If you use my code or modify it to your liking, please give credit by citing or linking my github page.
If you use the relative entropy tool, cite the paper: \
Bose N, Moore SD. Variable Region Sequences Influence 16S rRNA Performance. Microbiol Spectr. 2023 Jun 15;11(3):e0125223. [doi: 10.1128/spectrum.01252-23](https://journals.asm.org/doi/10.1128/spectrum.01252-23). Epub 2023 May 22. PMID: [37212673](https://pubmed.ncbi.nlm.nih.gov/37212673/); PMCID: [PMC10269663](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10269663/)

<br>
<br>

## 8) License
The MultiAlignment Toolkit is under an [MIT](https://github.com/nibose92/MultiAlignment_Toolkit/blob/Code-List/LICENSE) License
