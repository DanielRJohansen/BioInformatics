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
seq_len = len(genes[0].seq)
gene_matrix = np.chararray((len(genes), seq_len),unicode=True)
for i in range(len(genes)):
    gene_matrix[i,:] = genes[i].seq

i = 0
target_len = 19
prob_matrix = np.zeros((target_len,4), dtype=np.float)
while i < seq_len-target_len:
    for j in range(target_len):
        prob_matrix[i][0] = None
        a = np.where(gene_matrix[:,i+j], 'A')
        print(a)
        exit()



print(gene_matrix)
#makePWM(genes)
