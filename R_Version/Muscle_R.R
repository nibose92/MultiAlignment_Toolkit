library(Biostrings)
library(muscle)


# Enter input .fasta file name that is to be aligned, the sequence type, and name of the output multiple sequence alignment file.
# The following four are the only vectors that the user has to change/assign 
InFile <- "Ecoli_Shigella_16Sgenes.fasta"
Seq_Type <- "DNA"
OutFile <- "Ecoli_Shigella_16Sgenes_MSA.fasta"


# Depending on the type of sequence you entered, the sequences need to be read as a DNA/RNA/Amino Acid String Set. also, only consider residues
if (Seq_Type == "DNA") {
  Seq_File <- readDNAStringSet(InFile)
} else if (Seq_Type == "RNA") {
  Seq_File <- readRNAStringSet(InFile)
} else if (Seq_Type == "Peptide") {
  Seq_File <- readAAStringSet(InFile)
} else {
  print("Incorrect sequence type -- Enter DNA/RNA/Peptide")
  return()
}

# Perform multiple sequence alignment using the 'muscle' package, convert the alignment to a .fasta file, and then save the .fasta file
MSA <- muscle(Seq_File)
MSA_fasta <- as(MSA, "DNAStringSet")
writeXStringSet(MSA_fasta, file = OutFile)

print("Completed MUSCLE Alignment")