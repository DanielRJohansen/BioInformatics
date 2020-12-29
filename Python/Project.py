from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Align import substitution_matrices


BLOSUM62 = substitution_matrices.load("BLOSUM62")

############## TASK 1
def loadSequences(genome):
    sequences = []
    filename = genome + ".fa"
    with open(filename, 'r') as file:
        all_seqs = SeqIO.parse(file, "fasta")
        for seq in all_seqs:
            seq.id = genome     # Makes it possible to distinguish different genomes
            sequences.append(seq)
        print(sequences[0])
    return sequences


def writeFile(filename, sequences):
    records = (SeqRecord(seq=s.seq, id=s.id, name=s.name, description=s.description) for s in sequences)
    SeqIO.write(records, filename +".fasta", "fasta")


def mergeSequences(seq1, seq2):
    return seq1+seq2

def longestSeq(seqs):
    longest = 0
    for seq in seqs:
        if (len(seq.seq) > longest):
            longest = len(seq.seq)
        print(len(seq.seq))
    print("longest ", longest)

def task1():
    human_sequences = loadSequences("human")
    longestSeq(human_sequences)
    #mouse_sequences = loadSequences("mouse")
    #sequences = mergeSequences(mouse_sequences, human_sequences)
    #writeFile("hm_nonpadded_seqs", sequences)


############## TASK 2

def loadSeqDatabase(genome):    # hm_seq, mouse_seq, human_seq
    seqs = []
    filename = "Database\\" + genome + ".fa"
    #filename = "human_mulseq_default.fa"
    with open(filename, 'r') as file:
        all_seqs = SeqIO.parse(file, "fasta")
        for seq in all_seqs:
            seqs.append(seq)
    return seqs

def updateHighestScore(match, seq, genome):
    print("Higher score achieved: ", match.score, " Genome: ", genome, " Sequence id:", seq.id)
    print(match)
    print()

def pairwiseAlignment(genome, database, target_sequence):
    print("Doing pairwise sequence align for: ", genome)
    highest_score = 0
    index = 0
    for seq in database:
        alignments = pairwise2.align.globalxx(seq.seq, target_sequence)
        for match in alignments:
            if match.score > highest_score:
                highest_score = match.score
                updateHighestScore(match, seq, genome)




# Creating sample sequences

def task2():
    #target_sequence = "QVQLQQSGPGPAKPSQTLSLTCAISGDSVSSDSAAWNWIRQSPSRGLEWLGRTYYRSTWYRDYAPSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARDKDSFESNGSLYTAKKMGFDPWGQGTLVTVSETTLTQSPGTLSLSPGERATLSCRASQSVRSSYLAWYQQKPGQAPRLLIYGASSRATGIPERFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSPITFGQGTRLEIKRTVAAPS"
    target_sequence = "QVQLQQSGPGPAKPSQTLSLTCAISGDSVSSDSAAWNWIRQSPSRGLEWLGRTYYRSTWYRDYAPSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARDKDSFESNGSLYTAKKMGFDPWGQGTLVTVS"
    genome = "human"
    database = loadSeqDatabase(genome)
    pairwiseAlignment(genome, database, target_sequence)

    print("\n")

    genome = "mouse"
    database = loadSeqDatabase(genome)
    pairwiseAlignment(genome, database, target_sequence)


def main():
    #print(BLOSUM62)
    #return
    #task1()
    task2()

    # Showing results
    #






















""" Not needed anyway
def padSequences(sequences):
    padded_sequences = []

    longest_length = max(len(s.seq) for s in sequences)
    print("Longest length: ", longest_length)
    for s in sequences:
        while len(s.seq) < longest_length:
            s.seq = s.seq + '-'
        padded_sequences.append(s)
    return padded_sequences
"""




































