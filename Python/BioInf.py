import fastaparser
import numpy as np


class Gene:
    def __init__(self, id, desc, seq):
        self.id = id
        self.desc = desc
        self.seq = seq


def makePWM(genes, categories=('A', 'T', 'C', 'G')):
    seq_length = len(genes[0].seq);
    PWM = np.zeros((seq_length, len(categories)), dtype=int)

    for gene in genes:
        for i in range(seq_length):
            letter_index = categories.index(gene.seq[i])
            PWM[i, letter_index] += 1

    print(PWM)
    highest_prob = np.argmax(PWM, axis=1)
    print(categories[highest_prob[0]])
    for i in range(seq_length):
        print(categories[highest_prob[i]])




genes = []
with open("ex4_sequences.fa") as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    for seq in parser:
        # seq is a FastaSequence object
        gene = Gene(seq.id, seq.description, seq.sequence_as_string())
        genes.append(gene)
makePWM(genes)
