
import fastaparser
import numpy as np


class Gene:
    def __init__(self, id, desc, seq):
        self.id = id
        self.desc = desc
        self.seq = seq


def charToPos(char):
    if char == 'A':
        return 0
    if char == 'C':
        return 1
    if char == 'G':
        return 2
    if char == 'T':
        return 3


def startingMotifs(index, n):
    motifs = []
    for i in range(num_genes):
        motif = ""
        for j in range(motif_len):
            motif = motif + gene_matrix[i, index+j]
        motifs.append(motif)
    return motifs


def mostLikelyMotif(prob_matrix):
    likely = ""
    for i in range(19):
        likely = likely + categories[np.argmax(prob_matrix[i])]
    return likely


def modelLikelyhood(gene_sub_seq, prob_matrix):
    P = prob_matrix[0, categories.index(gene_sub_seq[0])]
    i = 1
    while i < len(gene_sub_seq):
        char_likelyhood = prob_matrix[i, categories.index(gene_sub_seq[i])]
        P = P * char_likelyhood
        i += 1
    return P


def modelLikelyhoods(gene_matrix, prob_matrix):
    likelyhood_matrix = np.zeros((num_genes, seq_len-motif_len))
    for i in range(num_genes):
        for j in range(seq_len-motif_len):
            likelyhood_matrix[i, j] = modelLikelyhood(gene_matrix[i, j:(j+motif_len)], prob_matrix)
    return likelyhood_matrix


def motifAtPos(pos):
    motif = ""
    for i in range(motif_len):
        motif = motif + gene_matrix[pos[0], pos[1]+i]
    return motif


def nBestMotifs(likelyhood_matrix, n):
    motifs = []
    for i in range(n):
        best_pos = np.unravel_index(np.argmax(likelyhood_matrix), likelyhood_matrix.shape)
        motif = motifAtPos(best_pos)
        motifs.append(motif)
        #motifs.append(likelyhood_matrix[best_pos[0], best_pos[1]:best_pos[1]+motif_len])
        likelyhood_matrix[best_pos] = 0
    return motifs


def newModel(motifs):
    num_motifs = len(motifs)
    model = np.zeros((seq_len, 4))
    for i in range(num_motifs):
        for j in range(motif_len):
            char = motifs[i][j]
            model[j, charToPos(char)] += 1
    # Normalize model
    model = np.divide(model, num_motifs)
    return model

def emAlgo(gene_matrix):
    model = newModel(startingMotifs(30, 25))    # Initial guess

    most_likely_motif = mostLikelyMotif(model)
    print("Initial guess: ", most_likely_motif)

    for i in range(5):
        model_likelyhoods = modelLikelyhoods(gene_matrix, model)
        motifs = nBestMotifs(model_likelyhoods, 100)
        model = newModel(motifs)

        most_likely_motif = mostLikelyMotif(model)
        print("Initial guess: ", most_likely_motif)


def makeGeneMatrix():
    genes = []
    with open("ex4_sequences.fa") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            gene = Gene(seq.id, seq.description, seq.sequence_as_string())
            genes.append(gene)
    seq_len = len(genes[0].seq)
    gene_matrix = np.chararray((len(genes), seq_len), unicode=True)
    for i in range(len(genes)):
        for j in range(len(genes[i].seq)):
            gene_matrix[i, j] = genes[i].seq[j]
    return gene_matrix


num_genes = 25
seq_len = 60
motif_len = 19
categories = ('A', 'C', 'G', 'T')
gene_matrix = makeGeneMatrix()
emAlgo(gene_matrix)
#makePWM(genes)
