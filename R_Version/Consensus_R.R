library(Biostrings)
library(ggplot2)


# Enter input multiple alignment file name, types of sequences, the limit to consider a residue as consensus, and name of the output consensus file.
# The following four are the only vectors that the user has to change/assign 
InFile <- "Ecoli_Shigella_16Sgenes_MSA.fasta"
Seq_Type <- "DNA"
Limit <- 0.51
OutFile <- "Consensus"




# Depending on the type of sequence you entered, the sequences need to be read as a DNA/RNA/Amino Acid String Set. also, only consider residues
if (Seq_Type == "DNA") {
  MSA_File <- readDNAStringSet(InFile)
  residues <- matrix(c("A", "C", "G", "T", "-"))
} else if (Seq_Type == "RNA") {
  MSA_File <- readRNAStringSet(InFile)
  residues <- c("A", "C", "G", "U", "-")
} else if (Seq_Type == "Peptide") {
  MSA_File <- readAAStringSet(InFile)
  residues <- c("G", "A", "V", "C",
                "P", "L", "I", "M",
                "N", "W", "F", "S",
                "T", "Y", "N", "Q",
                "H", "R", "K", "D", "E", "-")
  
} else {
  print("Incorrect sequence type -- Enter DNA/RNA/Peptide")
  return()
}

# Create the Consensus sequence based on the input
ConSeq <- consensusString(MSA_File, threshold = Limit, ambiguityMap = 'N')
# Save the consensus as a fasta file. The header has the number of sequences used to generate the consensus and the length of the consensus
cat(paste0(">",OutFile,"_length_",width(ConSeq),"\n",ConSeq),
    file = paste0(OutFile,".fasta"))

# Generate the position-specific frequency matrix, transpose it and remove non-residue columns, store the transposed matrix as a dataframe and convert the frequencies to percents -> consider as dataframe (1)
Res_Mat <- consensusMatrix(MSA_File)
T_Res_Mat <- t(Res_Mat)[, residues]
Res_Freq <- as.data.frame(round(((T_Res_Mat / length(MSA_File)) * 100), 2))
colnames(Res_Freq) <- paste0(colnames(Res_Freq), "%")

# Create a dataframe of the consensus residue frequency -> which is the max frequency at each residue position -> consider as dataframe (2)
Con_Freq <- data.frame("Consensus_Frequency" = apply(Res_Freq, 1, max))

# Create a dataframe of the consensus residue from the sequence -> consider as dataframe (3)
Con_Res <- data.frame("Consensus_Residue" = strsplit(ConSeq, split = "")[[1]])

# Create a dataframe of the residue positions -> consider as dataframe (4)
Res_Pos <- data.frame("Residue_Position" = 1:width(ConSeq))

# Combine dataframes (1), (2), (3), and (4) into one dataframe 'PSSM_DF' and store it under the entered .csv file name. 
PSSM_DF <- cbind(Res_Pos, Res_Freq, Con_Res, Con_Freq)
write.table(PSSM_DF, file = paste0(OutFile,".csv"), sep = ",", row.names = FALSE)

# plot the bar graph with ggplot2
plot <- ggplot(PSSM_DF, aes(x = Residue_Position, y = Consensus_Frequency)) + 
  geom_bar(stat = "identity", fill = "darkblue") + 
  scale_x_continuous(breaks = seq(0, width(ConSeq), 50), expand = c(0.01, 0)) + 
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0.01, 0)) + 
  labs(title = "Frequency of Consensus Residue at each position",
       x = "Residue Position", y = "Consensus Frequency (as %)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

ggsave(paste0(OutFile,".png"), plot = plot, width = 15, height = 8)

print("Consensus Generated!")