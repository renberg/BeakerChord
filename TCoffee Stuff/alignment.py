import Bio
from Bio import AlignIO
from Bio import Align
from Bio.Align import Applications
from Bio.Align.Applications import _TCoffee
from Bio.Align.Applications import TCoffeeCommandline
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import *
from Bio import SeqIO
from Bio.PDB.Polypeptide import three_to_one
#Function definitions
def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()
##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################
def shannon_entropy(list_input):
    import math
    unique_base = set(list_input)                           # Get only the unique bases in a column
    #unique_base = unique_base.discard("-")
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)                        # Number of residues of type i                   
        P_i = n_i/float(M)                                  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    #print sh_entropy
    return sh_entropy
##################################################################
# Function to calculate Shannon's entropy per alignment column for the whole MSA
##################################################################

def shannon_entropy_list_msa(alignment_file):
    shannon_entropy_list = []
    for col_no in xrange(len(list(alignment_file[0]))):
        list_input = list(alignment_file[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))
    return shannon_entropy_list

#START ACTUAL PROGRAM

############################################
#find match
############################################
#read in .pdb sequence
parser = PDBParser()
structure = parser.get_structure('', '1OTH.pdb')
header = parser.get_header()
trailer = parser.get_trailer()
pdbSequence = ''
for residue in structure[0]['A'].get_residues():
    if (residue.get_id()[0]==' '):
        residueName = three_to_one(residue.get_resname())
        pdbSequence += residueName
#check each fasta sequence for match to pdb sequence
bestMatch = ''
bestScore = 0
handle = open("uniprot-ornithine+transcarbamylase-2.fasta", "rU")
for record in SeqIO.parse(handle, "fasta") :
    foundMatch = False
    if foundMatch == False:
        tempFile = open("temp.fasta", "w")
        deleteContent(tempFile)
        tempFile.write(">sp|000000|FAKE HEADER OS=Fakus Faky GN=FAK PE=0 SV=0\n")
        tempFile.write(pdbSequence + "\n")
        tempFile.write(">sp|%s|FAKE HEADER OS=Fakus Faky GN=FAK PE=0 SV=0\n"%(record.name))
        for x in record.seq:
            tempFile.write(x)
        tempFile.write("\n")
        tempFile.close()
        cline = TCoffeeCommandline(infile="temp.fasta",
                                       output="clustalw",
                                       outfile="temp.aln")
        cline()
        tempAlignment = AlignIO.read(open("temp.aln"), "clustal")
        scores = shannon_entropy_list_msa(tempAlignment)
        matchScore = 0
        for num in scores:
            matchScore += num
        if matchScore >= bestScore:
            bestScore = matchScore
            bestMatch = record
        if matchScore == len(pdbSequence):
            foundMatch = True
            print "%s is a match"%(record.name)
handle.close()

#Do alignment
#creates command line wrapper
tcoffee_cline = TCoffeeCommandline(infile="uniprot-ornithine+transcarbamylase-2.fasta",
                                   output="clustalw",
                                   outfile="aligned.aln")
#executes command line wrapper
tcoffee_cline()
#read in the alignment file
alignment = AlignIO.read(open("aligned.aln"), "clustal")
#print results to standard output, this part will be removed in real program
print "Alignment length %i" % alignment.get_alignment_length()
#for record in alignment :
#	print record.seq, record.id

#Trim alignment
bestMatchPos = 0
i = 0
for record in alignment:
    if record.id == bestMatch.id:
        bestMatchPos = i
    i += 1

limit = alignment.get_alignment_length()
i = 0
while i < limit:
    if alignment[bestMatchPos][i] == '-':
        alignment = alignment[:,:i]+alignment[:,i+1:]
        limit -= 1
        print "Edited Alignment length %i" % alignment.get_alignment_length()
    if i < limit:
        if alignment[bestMatchPos][i] != '-':
            i += 1

#calculate Shannon entropy scores (from https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py)
# Add an import path to BioPython
import sys
sys.path.append("~/Desktop/Softwares/biopython-1.65/Bio/")

# import pandas for data frame
import pandas as pd
from numpy.random import randn

###################################################################
## Error message if the lengths of the sequences in MSA are not the same
###################################################################
seq_lengths_list = []
for record in alignment:
    seq_lengths_list.append(len(record))
row_num = len(seq_lengths_list)                                            # Get number of rows in the MSA
seq_lengths = set(seq_lengths_list)                                         # Get unique lengths of the sequences aligned in the MSA
if len(seq_lengths) != 1:
    print "Check you input Alignment!",                                     # Error message if the lengths of the sequences in MSA are not the same


scores = shannon_entropy_list_msa(alignment)

#write scores to .txt file
f = open("entropyScores.txt","w")
f.write("position score\n")
for x in range(len(scores)):
	f.write("%d %f"%((x+1),scores[x]))
	if (x < len(scores)-1):
		f.write("\n")
f.close()



print "done"